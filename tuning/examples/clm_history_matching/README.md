# CLM history matching on Derecho (DRAFT)

Calibrate CLM parameters with history matching, running ensembles on Derecho.
This wires `history_matching` + `clm` + `derecho` together. The CLM and Derecho
plugins are **drafts** modeled on [NCAR/ctsm6_ppe](https://github.com/NCAR/ctsm6_ppe)
— they need your machine's paths and a live CESM case, and are untested.

## How it maps to ctsm6_ppe

| ctsm6_ppe | here |
|---|---|
| `gen_paramfiles/*paramranges*.csv` | `CLM.parameters()` reads param names + ranges |
| `gen_paramfiles/ppe_tools` (copy base paramfile, set values) | `CLM.apply()` |
| `jobscripts/run_ens.sh` (clone `--keepexe`, setup, edit `user_nl_clm`, submit) | `DerechoRunner.run()` |
| `ctsm6_calibration` diagnostics | `CLM.compute_metrics()` |
| `analysis_lhc` first wave + emulator refocus | `history_matching` sampler |

## One wave per restart

The runner does **not** wait for jobs to finish. Each run of the script does one
wave and exits; you re-run it after that wave's jobs complete. `Campaign` saves
state in `campaign.pkl` between runs.

```
run 1:  sampler.ask() -> 100 sets (wave 1 = Latin hypercube)
        DerechoRunner clones base case x100, CLM.apply writes each paramfile, submits, exits
        ... you wait for the wave to finish on Derecho ...
run 2:  CLM.compute_metrics reads wave 1's history files -> diagnostic numbers
        sampler.tell() fits the emulator, rules out implausible regions
        sampler.ask() -> wave 2, submitted; exits
        ... repeat until the sampler is done ...
```

## What you must provide

- **Base case** — a CESM case already built (`create_newcase` + `case.build`);
  members are clones with `--keepexe`. Set `derecho.base_case`.
- **Base paramfile** — the default CLM netCDF that `apply()` copies and edits.
  Set `clm.base_paramfile`.
- **Ranges CSV** — from `ctsm6_ppe/gen_paramfiles`. Set `clm.paramranges_csv`.
- **Observations** — targets + uncertainty in `config.yaml`, matching whatever
  `compute_metrics()` returns.
- **Project / paths** — `derecho.project`, `output_root`; CIME scripts on PATH.

## Run — one wave per restart

```bash
pip install -e "tuning/[emulator]"
python tuning/examples/clm_history_matching/run.py    # submits wave 1, exits
# ... wait for the Derecho jobs to finish ...
python tuning/examples/clm_history_matching/run.py    # collects wave 1, submits wave 2
# ... repeat until it prints the best parameters ...
```

Each wave's cases land in `output_root/waveNN/`. State lives in `campaign.pkl`.

## Known gaps (clean up before real use)

- **Scalar params only.** Per-PFT parameters (CSV `min/max` = `pft`/`30percent`)
  are skipped; they need a vector parameter type.
- **`compute_metrics` is a placeholder** — a single global+time mean. Replace
  with the real obs-matching diagnostics.
- **Wave readiness** — the runner never waits, and `Campaign.step()` assumes the
  previous wave is finished when you re-run. If you automate re-runs, add a
  PBS/`qstat` completion check before collecting. `_history_dir` points at
  `run/` — adjust to your archive layout.
