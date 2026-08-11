# ML Warm-Rain Microphysics (TAU Emulator) Integration

This folder documents and provides sources for a neural-network emulator of
the TAU bin-resolved warm-rain (cloud-rain autoconversion/accretion)
microphysics scheme in CESM/CAM, implemented as a self-contained Fortran +
netCDF neural network (no FTorch/PyTorch dependency).

- **Component:** CAM (CESM3), PUMAS microphysics package
- **What's replaced:** the deterministic bin-resolved TAU warm-rain
  calculation (`pumas_stochastic_collect_tau.F90`) with a fast NN surrogate,
  selected via the existing `micro_mg_warm_rain = 'emulated'` namelist option
- **Key idea:** PUMAS already had a namelist-selectable `'emulated'` scheme
  path upstream; this integration swaps in a 9-input neural network (with
  droplet/rain size-distribution shape parameters as inputs and dual
  forward/inverse quantile interpolation) in place of upstream's simpler
  7-input version

## ⚠️ Requires CAM7 physics (`-phys cam7`)

`micro_mg_warm_rain` and the whole `'tau'`/`'emulated'` scheme selection only
exist in the **CAM7** physics driver
(`src/atmos_phys` → `src/physics/cam7/micro_pumas_cam.F90`). They are
**absent** from the CAM6 driver (`src/physics/cam/micro_pumas_cam.F90`).
Using a CAM6-physics compset (e.g. `F2000climo`, which defaults to
`-phys cam6`) will fail at startup with:

```
micro_pumas_cam_readnl:: ERROR reading namelist
```

because the Fortran namelist reader chokes on `micro_mg_warm_rain` as an
unrecognized variable in `&micro_mg_nl`. Use a CAM7-native compset instead,
e.g. `F2000dev` (`2000_CAM70_CLM60%SP_...`).

## Repository Structure

```text
TAU_ML_microphysics/
├── src/
│   ├── module_neural_net.F90       # Dense NN layers, activations, forward/inverse quantile interpolation
│   ├── tau_neural_net_quantile.F90 # TAU emulator driver: builds NN inputs, calls the model, returns tendencies
│   ├── tester.F90                  # Debug utility (write_test_values), used by tau_neural_net_quantile.F90
│   ├── micro_pumas_v1.F90          # Modified: 'emulated' branch call site (NN inputs, in-cloud rescaling)
│   ├── micro_pumas_utils.F90       # Modified: added qsmall_emulator threshold
│   └── ML_fixer_check.F90          # Unmodified from upstream v1.39; included for reference (part of the call chain)
│
├── bld/
│   └── namelist_files/
│       ├── namelist_defaults_cam.xml
│       └── namelist_definition.xml
│
├── docs/
│   └── build_instructions.md
│
└── README.md
```

Note: `pumas_stochastic_collect_tau.F90` (the deterministic bin-resolved TAU
scheme, selected by `micro_mg_warm_rain = 'tau'`) is **not** included here —
it's a separate, mutually-exclusive code path from the ML `'emulated'`
scheme, and was left untouched during this integration (its own
`lcldm`/`precip_frac` handling changed independently upstream between PUMAS
`v1.36` and `v1.39`, unrelated to this ML work).

---

## Prerequisites

- CESM3 (tested with `cesm3_0_alpha09e`, base CAM tag `cam6_4_186`,
  PUMAS `pumas_cam-release_v1.39`)
- CAM7 physics package (`-phys cam7`) — see warning above
- No FTorch/PyTorch dependency — this is pure Fortran + netCDF

## Overview of Integration

```
micro_pumas_v1.F90 (micro_pumas_tend, 'emulated' branch)
         ↓
tau_neural_net_quantile.F90 (tau_emulated_cloud_rain_interactions)
         ↓
module_neural_net.F90 (neural_net_predict, quantile_transform/quantile_inv_transform)
```

### 1. `micro_pumas_v1.F90` — Call Site

**Purpose:** Selects and drives the ML emulator when
`micro_mg_warm_rain = 'emulated'` (mutually exclusive with `'tau'`,
`'sb2001'`, `'kk2000'`).

**Key changes vs. upstream v1.39's `'emulated'` branch:**

- Passes cloud/rain size-distribution shape parameters (`pgam`, `lamc`,
  `lamr`, `n0r` — already computed elsewhere in this subroutine for the
  deterministic scheme) into the NN call, in addition to the in-cloud/
  in-precip state (`qcic`, `ncic`, `qric`, `nric`) upstream already uses:

  ```fortran
  call tau_emulated_cloud_rain_interactions(qcic(1:mgncol,k), ncic(1:mgncol,k), &
                                            qric(1:mgncol,k), nric(1:mgncol,k), &
                                            pgam(1:mgncol,k), lamc(1:mgncol,k), &
                                            lamr(1:mgncol,k), n0r(1:mgncol,k), &
                                            rho(1:mgncol,k), lcldm(1:mgncol,k), &
                                            precip_frac(1:mgncol,k), mgncol, qsmall_emulator, &
                                            proc_rates%qctend_TAU(1:mgncol,k), &
                                            proc_rates%qrtend_TAU(1:mgncol,k), &
                                            proc_rates%nctend_TAU(1:mgncol,k), &
                                            proc_rates%nrtend_TAU(1:mgncol,k))
  ```

- Uses a dedicated `qsmall_emulator` threshold (1e-8) instead of the general
  `qsmall` (1e-18).

- `ML_fixer_calc` is called with in-cloud/in-precip values (`qcic`/`ncic`/
  `qric`/`nric`) rather than grid-mean values.

- After both calls, rescales the resulting tendencies by `lcldm`/
  `precip_frac` to convert from in-cloud/in-precip space back to grid-box
  means:

  ```fortran
  do i=1,mgncol
     prc(i,k)  = -proc_rates%qctend_TAU(i,k) * lcldm(i,k)
     nprc1(i,k)= -proc_rates%nctend_TAU(i,k) * lcldm(i,k)
  end do
  do i=1,mgncol
     if (proc_rates%nrtend_TAU(i,k) > 0._r8) then
        nprc(i,k)  = proc_rates%nrtend_TAU(i,k) * precip_frac(i,k)
     else
        nragg(i,k) = proc_rates%nrtend_TAU(i,k) * precip_frac(i,k)
     end if
  end do
  ```

### 2. `tau_neural_net_quantile.F90` — Emulator Driver

**Purpose:** Builds the 9-element NN input vector, runs it through the
quantile-transformed neural network, and returns tendencies.

**Key difference vs. upstream's 7-input version:**

```fortran
integer, parameter :: num_inputs = 9   ! upstream v1.39: 7

nn_inputs(1,1) = qc(i);   nn_inputs(1,2) = qr(i)
nn_inputs(1,3) = nc(i);   nn_inputs(1,4) = nr(i)
nn_inputs(1,5) = pgam(i); nn_inputs(1,6) = lamc(i)
nn_inputs(1,7) = lamr(i); nn_inputs(1,8) = n0r(i)
nn_inputs(1,9) = rho(i)
```

Upstream's 7-input version uses `precip_frac`/`lcldm` as direct NN inputs
instead of the size-distribution shape parameters; this version instead
applies `lcldm`/`precip_frac` as a rescaling step after the NN call (see
`micro_pumas_v1.F90` above).

### 3. `module_neural_net.F90` — NN Math

**Purpose:** Dense-layer neural network primitives (matrix multiply via
`dgemm`, activation functions) and quantile-based input/output scaling.

**Key difference vs. upstream:** two separate interpolation routines,
`linear_interp_forward` and `linear_interp_inverse`, instead of a single
shared `linear_interp` — used respectively by `quantile_transform` (inputs)
and `quantile_inv_transform` (outputs).

### 4. `tester.F90`

Debug utility (`write_test_values`) — writes NN input/output arrays to a
file for offline validation. `tau_neural_net_quantile.F90` imports this
module but does not currently call it in the main code path.

## Configuration

### Namelist Configuration

```fortran
&micro_mg_nl
  micro_mg_warm_rain = 'emulated'
/

&stochastic_emulated_nl
  stochastic_emulated_filename_quantile     = '/path/to/quantile_neural_net_fortran.nc'
  stochastic_emulated_filename_input_scale  = '/path/to/input_quantile_scaler.nc'
  stochastic_emulated_filename_output_scale = '/path/to/output_quantile_scaler.nc'
/
```

These namelist variables and groups already exist in upstream
`namelist_definition.xml` (`stochastic_emulated_nl`) — no new namelist
entries were needed for this integration, only the underlying Fortran.

## Required Input Files

1. **NN weights file** (`quantile_neural_net_fortran.nc`) — dense-layer
   weights/biases/activations, one dimension pair per layer
   (`dense_NN_in`/`dense_NN_out`). Must have `dense_00_in = 9` to match this
   9-input model.
2. **Input quantile scaler** (`input_quantile_scaler.nc`) — `column`
   dimension must be 9, matching the input feature count.
3. **Output quantile scaler** (`output_quantile_scaler.nc`) — matches the
   3 NN outputs (`qc`/`nc`/`nr` tendencies).

Verify weight-file compatibility before running:
```
ncdump -h quantile_neural_net_fortran.nc | grep dense_00_in
ncdump -h input_quantile_scaler.nc | grep column
```
Both should show `9`.

## Provenance

Ported from the `wkchuang/PUMAS` fork
(https://github.com/wkchuang/PUMAS, branch `1nn_tau_emulator`, commit
`8b76bac`) onto the CCPP-ized `pumas_cam-release_v1.39` baseline (PUMAS
moved from a standalone `src/physics/pumas` submodule to a nested submodule
under `src/atmos_phys/schemes/pumas/pumas` between these two states, and the
fork's `shr_kind_mod` usage was adapted to this codebase's host-agnostic
`pumas_kinds` module).

## References

- Based on the TAU bin-resolved stochastic collection scheme
  (Morrison/Lebo, Gettelman and Chen 2018)
- PUMAS: https://github.com/ESCOMP/PUMAS
- ML fork: https://github.com/wkchuang/PUMAS (branch `1nn_tau_emulator`)
