# CLAUDE.md — LEAP / CESM3-MLe Integration repo

## Repo
- Hub for building a hybrid (ML + physics) CESM3 (**CESM3-MLe**): docs, tutorials, and examples for ML parameterizations and calibration. Mostly Markdown plus one Python package.
- Top level: `Start Here.md` (ML-param workflow), `tools/` (functional unit test), `CAM_Convection_YOG/` + `CLM_phenology/` (example integrations), `tuning/` (calibration package — **active development**).
- Git: commits go to `main`. End commit messages with the `Co-Authored-By` trailer. **Stage explicit paths**, not `git add -A` at repo root (stray `.DS_Store`/build artifacts exist; `.gitignore` covers them). Commit only when asked.

## tuning/ — calibration package (pip name `tuning`, src layout)
- Goal: **mix-and-match** CESM calibration. Combine a **Sampler** (how to choose params) × **Component** (which model) × **Runner** (where it runs); swap any one without touching the others.
- One question drives it: *"given the runs so far, which parameter sets next?"* → a Sampler answers.
- Contracts in `src/tuning/core/interfaces.py` (**read first**):
  - Sampler: `ask()` next wave · `tell(params, metrics)` · `is_done()` · `result()`
  - Component: `parameters()` · `apply(case_dir, params)` · `compute_metrics(run_output)`
  - Runner: `run(ensemble, component)`
- A **wave** = one batch of forward-model runs. Flow: ask → run → compute_metrics → tell → repeat → result (`orchestration/loop.py`).
- Plugins self-register by name via decorators (`core/registry.py`), imported in each subpackage `__init__.py`.
- Data: a param set is a plain dict `{name: float}`. `Observations` = targets + uncertainty + `loss()`. `CalibrationResult` = `best_params` + `extras`.

### Two run modes
- `loop.py:run_calibration(sampler, component, runner, max_waves=...)` — synchronous, one process. For runners that finish immediately (`local`); used by tests.
- `campaign.py:Campaign.step()` — one wave per call, state pickled to a file, runner does NOT wait. Re-run once per wave for long HPC jobs (`derecho`).

### Samplers (all ✅ real + tested; pure-Python/sklearn, no CESM)
- `random`, `latin_hypercube` — one-wave designs.
- `eki` — ensemble Kalman inversion (numpy).
- `history_matching` — GP emulator (sklearn) → implausibility → NROY refocus.
- `hmc` — leapfrog Hamiltonian Monte Carlo on a GP emulator, in a normalized [0,1] box; posterior + optional validation wave.
- `nelder_mead` — simplex optimizer as an ask/tell state machine.
- Shared: `samplers/_emulator.py` (`GPEmulator`, lazy sklearn), `samplers/_select.py` (`best()`).

### Components / Runners
- `runners/local.py` ✅ — runs `model_fn(params)` in-process (tests).
- `components/clm.py`, `runners/derecho.py`, `examples/clm_history_matching/` — **DRAFTS** modeled on NCAR/ctsm6_ppe (`gen_paramfiles` + `jobscripts/run_ens.sh`). Untested; need a live CESM/Derecho env.
- `components/cam.py` — stub.

### Conventions
- Docs concise; code beginner-readable, flow obvious without digging.
- Heavy/optional imports (netCDF4, xarray, pandas, sklearn) **lazy inside methods**, so `import tuning` stays light and the core installs with numpy + pyyaml only.
- Mark unfinished work `[STUB]`/`[DRAFT]` and point to the relevant SCOPING phase.
- Status + roadmap: `tuning/SCOPING.md`.

### Commands
- Install: `pip install -e tuning/` (core) · `pip install -e "tuning/[emulator]"` (+ sklearn)
- Test: `cd tuning && pytest` (6 tests, ~1s)

### Next (needs a CESM machine)
- `cam` component; harden the CLM vertical: per-PFT (vector) parameters, real obs-matching diagnostics in `clm.compute_metrics`, qstat wave-readiness check in `derecho`. Gaps listed in `tuning/examples/clm_history_matching/README.md`.

## Response style
Concise and concrete. State assumptions. Verify with the test suite after code changes; report what ran. Don't claim done without running it.
