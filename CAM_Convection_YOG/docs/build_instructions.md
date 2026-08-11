# Build Instructions (Summary)

Tested against `cesm3_0_alpha09e` (CAM `cam6_4_186`) on Derecho. Originally
built against `cesm3_0_alpha07f` — see the "Porting to a New CESM Tag" and
"Build Gotchas" sections in the main README for what changed between tags.

1. Clone/checkout CESM3 (with submodules, e.g. via `git-fleximod update`).
2. Copy the files under `src/cam/` in this folder into
   `components/cam/src/physics/cam/` in your checkout (new files:
   `SAM_consts.F90`, `nn_cf_net.F90`, `nn_convection_flux.F90`,
   `nn_interface_cam.F90`, `yog_intr.F90`; modified: `phys_control.F90`,
   `physpkg.F90` — the latter needs a manual merge if your CAM tag has
   diverged from `cam6_4_186`, see README).
3. Copy `libraries/FTorch/{CMakeLists.txt,FTorch_cesm_interface.F90,buildlib}`
   from this folder into your checkout's `libraries/FTorch/`.
4. Confirm the nested FTorch submodule is actually populated:
   `ls libraries/FTorch/src` should show real files (`CMakeLists.txt`, `src/`,
   etc.), not an empty directory. If empty: `cd libraries/FTorch && git
   submodule update --init src`.
5. Copy `bld/namelist_files/{namelist_defaults_cam.xml,namelist_definition.xml}`
   from this folder into `components/cam/bld/namelist_files/` (or manually
   apply the YOG-related `<entry>` additions if your tag's files have
   otherwise diverged).
6. Create your case as usual, then:
   ```
   ./xmlchange USE_FTORCH=TRUE
   ./case.setup
   ./preview_namelists
   ./case.build --clean-all
   ./case.build
   ```
7. In `user_nl_cam`, set:
   ```
   yog_scheme = 'on'
   yog_nn_weights = '/path/to/model_weights.pt'
   yog_nn_scale   = '/path/to/scaling_metadata.nc'
   SAM_sounding   = '/path/to/sam_sounding.nc'
   ```
   (`yog_lat_max` defaults to 30° — only override this deliberately, see
   README "Build Gotchas").

If `case.build` fails on FTorch with a CMake `ADD_LIBRARY ... SHARED ...
does not support dynamic linking` error, see the Catamount/CMake note in the
README — the fix is already applied in `libraries/FTorch/buildlib` above,
but if you're merging against a different FTorch version, you may need to
reapply it.
