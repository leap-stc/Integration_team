# Build Instructions (Summary)

Tested against `cesm3_0_alpha09e` (CAM `cam6_4_186`, PUMAS
`pumas_cam-release_v1.39`) on Derecho.

**Requires CAM7 physics.** See the README warning — `micro_mg_warm_rain`
only exists in the CAM7 physics driver. Use a CAM7-native compset (e.g.
`F2000dev`), not a CAM6 one (e.g. `F2000climo`).

1. Clone/checkout CESM3 (with submodules, e.g. via `git-fleximod update`).
   PUMAS in this CESM tag lives at
   `components/cam/src/atmos_phys/schemes/pumas/pumas` (a submodule nested
   inside the `atmos_phys` submodule).
2. Copy the files under `src/` in this folder into that PUMAS directory:
   - New: `tester.F90` (not present upstream)
   - Modified: `module_neural_net.F90`, `tau_neural_net_quantile.F90`,
     `micro_pumas_v1.F90`, `micro_pumas_utils.F90`
   - Unmodified, included for reference: `ML_fixer_check.F90`
3. Copy `bld/namelist_files/{namelist_defaults_cam.xml,namelist_definition.xml}`
   from this folder into `components/cam/bld/namelist_files/` (or confirm
   your tag's files already have the `stochastic_emulated_nl` namelist group
   — these files are unchanged from upstream aside from whatever else
   differs between CAM tags).
4. Create a case with a **CAM7** compset, e.g.:
   ```
   ./create_newcase --case <casedir> --compset F2000dev --res f09_f09_mg17 \
       --project <project> --run-unsupported
   ```
   (`--run-unsupported` needed if this compset/resolution combo isn't yet a
   validated pairing in your CESM tag.)
5. Build as usual — no `USE_FTORCH` or other special `xmlchange` needed,
   this scheme has no FTorch/PyTorch dependency:
   ```
   ./case.setup
   ./case.build
   ```
6. In `user_nl_cam`:
   ```
   micro_mg_warm_rain = 'emulated'
   stochastic_emulated_filename_quantile     = '/path/to/quantile_neural_net_fortran.nc'
   stochastic_emulated_filename_input_scale  = '/path/to/input_quantile_scaler.nc'
   stochastic_emulated_filename_output_scale = '/path/to/output_quantile_scaler.nc'
   ```
7. Before submitting, sanity-check the weight file matches this 9-input
   model (see README "Required Input Files") — a mismatched file will
   either crash or silently produce garbage.

## Symptom → cause reference

| Symptom | Cause |
|---|---|
| `micro_pumas_cam_readnl:: ERROR reading namelist` at startup | Built with `-phys cam6` instead of `-phys cam7`. Rebuild with a CAM7 compset. |
| Crash/NaN in warm-rain tendencies | Check the NN weight file's `dense_00_in` dimension is 9, and the input/output scaler files' `column` dimension matches. |
