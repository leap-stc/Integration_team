# FTorch Deep Convection (YOG) Integration

This repository documents and provides sources for integrating a PyTorch-based
deep-convection scheme (YOG) into CESM/CAM using **FTorch**.

- **Component:** CAM (CESM3)
- **What’s replaced:** ZM/YOG convection tendencies via FTorch TorchScript model
- **Key idea:** Keep CAM physics + vertical remapping intact; swap the NN call with
  a TorchScript forward pass (FTorch), preserving CAM data flow.



## Repository Structure

```text
FTorch_CAM_integration/
├── src/
│   └── cam/                        # Modified CAM physics source files
│       ├── Phys_control.F90
│       ├── physpkg.F90
│       ├── yog_intr.F90
│       ├── nn_interface_CAM.F90
│       ├── nn_convection_flux.F90
│       └── nn_cf_net.F90
│
├── libraries/
│   └── FTorch/
│       └── FTorch_cesm_interface.F90   # Wrapper for FTorch model calls
│
├── docs/
│   ├── build_instructions.md
│   └── troubleshooting.md
│
├── examples/
│   └── user_nl_cam
│
├── MODEL_CARD.md
└── README.md

```

---

## Prerequisites

- CESM 3.0 (tested with `cesm3_0_alpha07f`)
- FTorch (Fortran interface to PyTorch)

## Overview of Integration


The integration follows a modular architecture and the following files are modified during the integration process:

```
CESM Physics Package (physpkg.F90)
         ↓
YOG Interface Layer (yog_intr.F90)
         ↓
CAM Interface (nn_interface_cam.F90)
         ↓
NN Convection Core (nn_convection_flux.F90)
         ↓
FTorch Neural Net (nn_cf_net.F90)
         ↓
FTorch Library (FTorch_cesm_interface)
```


### 1. `phys_control.F90` - Physics Control and Configuration

**Purpose:** Adds YOG scheme control to CESM physics options

**Key Changes:**

- Added `yog_scheme` configuration variable:
  ```fortran
  character(len=16) :: yog_scheme = unset_str
  ```

- Extended `phys_ctl_nl` namelist to include YOG scheme:
  ```fortran
  namelist /phys_ctl_nl/ cam_physpkg, use_simple_phys, cam_chempkg, waccmx_opt, &
                         deep_scheme, shallow_scheme, yog_scheme, ...
  ```

- Added `yog_scheme_out` optional parameter to `phys_getopts`:
  ```fortran
  subroutine phys_getopts(deep_scheme_out, shallow_scheme_out, yog_scheme_out, ...)
      character(len=16), intent(out), optional :: yog_scheme_out
      if (present(yog_scheme_out)) yog_scheme_out = yog_scheme
  end subroutine
  ```

### 2. `physpkg.F90` - Main Physics Package

**Purpose:** Integrates YOG scheme into the physics time-stepping loop

**Key Changes:**

- **Module Import:**
  ```fortran
  use yog_intr, only: yog_tend, yog_init, yog_final
  ```

- **Configuration Retrieval:**
  ```fortran
  call phys_getopts(shallow_scheme_out = shallow_scheme, &
                    yog_scheme_out = yog_scheme, ...)
  ```

- **Initialization** (`phys_init` subroutine):
  ```fortran
  if (yog_scheme == 'on') then
      call yog_init()
  end if
  ```

- **Time-stepping** (`tphysbc` subroutine):
  ```fortran
  if (yog_scheme == 'on') then
      call t_startf('yog_nn')
      call yog_tend(ztodt, state, ptend, pbuf)
      call physics_update(state, ptend, ztodt, tend)
      call t_stopf('yog_nn')
      flx_cnd(:ncol) = prec_dp(:ncol)
      call check_energy_cam_chng(state, tend, "yog_nn", nstep, ztodt, &
                                 zero, flx_cnd, zero, flx_heat)
  end if
  ```

- **Finalization** (`phys_final` subroutine):
  ```fortran
  if (yog_scheme == 'on') then
      call yog_final()
  end if
  ```

### 3. `yog_intr.F90` - YOG Interface Module

**Purpose:** CAM interface layer for the YOG deep convection scheme

**Public Subroutines:**

#### `yog_readnl(nlfile)`

Reads YOG-specific namelist parameters:

```fortran
namelist /yog_params_nl/ yog_nn_weights, yog_nn_scale, SAM_sounding
```

**Parameters:**
- `yog_nn_weights`: Path to PyTorch model weights (.pt file)
- `yog_nn_scale`: Path to NetCDF file with scaling parameters
- `SAM_sounding`: Path to SAM sounding profile data

#### `yog_init()`

- Registers output fields with CAM history buffer
- Initializes the neural network module
- Loads model weights and configuration files

**Output Fields:**
- `YOGDT`: Temperature tendency (K/s)
- `YOGDQ`: Water vapor tendency (kg/kg/s)
- `YOGDICE`: Cloud ice tendency (kg/kg/s)
- `YOGDLIQ`: Cloud liquid tendency (kg/kg/s)
- `PREC_YOG`: Surface precipitation (m/s)
- `YOGDNUMLIQ`: Cloud liquid number concentration tendency (N/s)
- `YOGDNUMICE`: Cloud ice number concentration tendency (N/s)

#### `yog_tend(ztodt, state, ptend, pbuf)`

Calculates physics tendencies using the neural network:
- Calls `nn_convection_flux_CAM` with CAM state variables
- Updates temperature, moisture, and cloud tendencies
- Computes surface precipitation
- Writes output to history buffer

#### `yog_final()`

Cleans up and finalizes the neural network module

### 4. `nn_interface_cam.F90` - CAM-to-NN Interface

**Purpose:** Bridges CAM physics state to neural network input format

#### `nn_convection_flux_CAM_init(nn_filename, metadata_filename, sounding_filename)`

**Inputs:**
- `nn_filename` (char*136): PyTorch .pt model file path
- `metadata_filename` (char*136): NetCDF file with scale/dimension parameters
- `sounding_filename` (char*136): NetCDF file with SAM sounding data

**Actions:**
- Initializes neural network from PyTorch model
- Loads SAM reference sounding and grid data

#### `nn_convection_flux_CAM(...)`

**Inputs:**
- `pres_cam`, `pres_int_cam`, `pres_sfc_cam`: Pressure fields
- `tabs_cam`, `qv_cam`, `qc_cam`, `qi_cam`: Temperature and moisture fields
- `cp_cam`: Specific heat capacity
- `dtn`: Timestep
- `ncol`: Number of columns

**Outputs:**
- `precsfc`: Surface precipitation
- `dqi`, `dqv`, `dqc`: Ice, vapor, and liquid tendencies
- `ds`: Dry static energy tendency

**Process:**
1. Interpolates CAM variables to SAM pressure levels
2. Calls `nn_convection_flux` with SAM-grid variables
3. Returns tendencies on CAM grid

### 5. `nn_convection_flux.F90` - NN Convection Core

**Purpose:** Core neural network convection parameterization

**Key Components:**

#### `nn_convection_flux_init(nn_filename, metadata_filename)`

- Initializes FTorch model interface
- Loads neural network architecture and weights
- Reads scaling parameters from metadata

**Uses:**
```fortran
use FTorch_cesm_interface, only: torch_model
```

#### `nn_convection_flux(tabs_i, q_i, tabs, t, q, rho, adz, dz, dtn, precsfc)`

**Process:**
1. Prepares input features from thermodynamic state
2. Applies input scaling/normalization
3. Calls `net_forward` for neural network inference
4. Applies physical constraints to outputs
5. Updates temperature, moisture, and precipitation fields

### 6. `nn_cf_net.F90` - FTorch Neural Network Interface

**Purpose:** Low-level FTorch interface for neural network operations

**Dependencies:**
```fortran
use FTorch_cesm_interface, only: torch_kCPU, torch_tensor, torch_model, &
                                 torch_tensor_from_array, torch_model_load, &
                                 torch_model_forward, torch_delete
```

**Public Interface:**
- `net_forward`: Execute neural network forward pass
- `nn_cf_net_init`: Initialize neural network from file

#### `nn_cf_net_init(nn_filename, metadata_filename, n_inputs, n_outputs, nrf, model_ftorch, iulog, errstring)`

**Actions:**
- Loads PyTorch model using FTorch:
  ```fortran
  call torch_model_load(model_ftorch, trim(nn_filename), torch_kCPU)
  ```
- Reads input/output dimensions from metadata
- Initializes model on CPU device

#### `net_forward(features, model_ftorch, logits)`

**Process:**
1. Creates FTorch tensor from input features:
   ```fortran
   call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU)
   ```

2. Creates output tensor:
   ```fortran
   call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
   ```

3. Executes forward pass:
   ```fortran
   call torch_model_forward(model_ftorch, in_tensor, out_tensor)
   ```

4. Cleans up tensors:
   ```fortran
   call torch_delete(in_tensor)
   call torch_delete(out_tensor)
   ```

## Configuration

### Namelist Configuration

Add to your CESM namelist file:

```fortran
&phys_ctl_nl
  yog_scheme = 'on'
/

&yog_params_nl
  yog_nn_weights = '/path/to/model_weights.pt'
  yog_nn_scale = '/path/to/scaling_metadata.nc'
  SAM_sounding = '/path/to/sam_sounding.nc'
/
```


## Required Input Files

1. **PyTorch Model File** (`.pt`)
   - Trained neural network weights
   - TorchScript format compatible with FTorch

2. **Metadata File** (`.nc`)
   - Input/output dimensions
   - Feature scaling parameters (mean, standard deviation)
   - Variable normalization factors

3. **SAM Sounding File** (`.nc`)
   - Reference atmospheric profile
   - Pressure level definitions
   - Grid spacing information



## References

- Based on the Zhang-McFarlane deep convection scheme framework
- Adapted for neural network-based parameterization using the Yuval-O'Gorman approach
- FTorch library: https://github.com/Cambridge-ICCS/FTorch
- FTorch CAM integration: https://github.com/addisug/FTorch_CAM


