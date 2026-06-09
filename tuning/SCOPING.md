# CESM Tuning Suite — Scoping

**Status:** Draft · **Updated:** 2026-06-08

## 1. Goal

A Python package (`tuning`, in this folder) for calibrating CESM components,
built so people can **mix and match**: combine any *sampler* (how to choose
parameters) with any *component* (which model) and any *runner* (where it runs).

Inspired by
[EnsembleKalmanProcesses.jl](https://clima.github.io/EnsembleKalmanProcesses.jl/dev/)
— the calibration algorithm knows nothing about the model. We do the same in
Python and don't reimplement the math:
[ESEm](https://github.com/duncanwp/ESEm) provides the emulators and samplers
that some methods use.

## 2. The key idea: it's all one question

Every calibration method really answers the same question:

> Given the runs so far, **which parameter sets do we run next?**

That answer is a **sampler**. A *wave* is one batch of forward-model runs. The
sampler proposes a wave, sees the results, and proposes the next:

```
sampler.ask() -> runner.run() -> component.compute_metrics() -> sampler.tell() -> repeat -> result()
```

Same interface, different answers:

| Sampler | Waves | Emulator? | Idea | Backed by |
|---|---|---|---|---|
| `random` | 1 | no | random within bounds | — |
| `latin_hypercube` | 1 | no | spread evenly across ranges | — |
| `eki` | few | no | Kalman nudge toward obs | Python EKI |
| `history_matching` | several | yes | keep "not ruled out yet" region | scikit-learn GP (ESEm optional) |
| `mcmc` / HMC | several | yes | sample the posterior | ESEm |
| `bayesopt` | many | yes | next point most likely to improve | ESEm (GP) |

Emulators are an *internal tool* of some samplers, not a separate concept.
`history_matching` and `eki` use the runs to pick the next wave directly;
`mcmc`/`bayesopt` build an emulator and sample it, after which you can draw N
sets from the posterior and run one more wave — which is just another `ask`.

## 3. The three contracts (`core/interfaces.py`)

- **Sampler** — `ask` (next wave) · `tell` (results) · `is_done` · `result`
- **Component** (`cam`, `clm`) — `parameters` · `apply` (write into a case) · `compute_metrics`
- **Runner** (`local`, `derecho`) — `run`

Swap any one without touching the others. That is the whole design.

## 4. Grounding: `ctsm6_ppe`

[NCAR/ctsm6_ppe](https://github.com/NCAR/ctsm6_ppe) is a CLM-specific version of
this, and seeds the first real implementations:

| `ctsm6_ppe` piece | Maps to |
|---|---|
| `analysis_lhc`, `sparsegrid` | `latin_hypercube` sampler (first wave) |
| `gen_ensembles/gen_paramfiles` | `CLM.apply` (write parameters) |
| `gen_ensembles/jobscripts` | `derecho` runner |
| `ctsm6_calibration/*.ipynb` | `history_matching` sampler (emulator + resample) |
| `diagnostics/`, `postprocessing/` | `CLM.compute_metrics` |

## 5. Layout

```text
tuning/
├── SCOPING.md · README.md · pyproject.toml
├── docs/         architecture · adding_a_sampler · adding_a_component · samplers_overview
├── src/tuning/
│   ├── core/         interfaces · parameters · observations · registry · config
│   ├── samplers/     random ✅ · latin_hypercube ✅ · history_matching ✅ · eki · mcmc · bayesopt
│   ├── components/   cam · clm
│   ├── runners/      local ✅ · derecho
│   └── orchestration/ loop.py
├── examples/clm_history_matching/   config.yaml · run.py
└── tests/test_loop.py
```

## 6. Plan

1. **Phase 1 — scaffold (done).** Core interfaces + loop; working `random` and
   `latin_hypercube` samplers; `local` runner; passing end-to-end test.
2. **Phase 2 — history matching (done).** Real `history_matching` sampler
   (emulate → rule out → refocus) on a swappable GP emulator (scikit-learn
   default; ESEm optional), tested locally. **Remaining:** the real `clm`
   component (`apply`/`compute_metrics`) and `derecho` runner — these need a live
   CESM environment, refactored from `ctsm6_ppe`.
3. **Phase 3 — breadth.** `eki`, `mcmc`, `bayesopt`; `cam`; the `derecho` runner.

## 7. Settled decisions

- Package name: **`tuning`** (matches the folder).
- EKI: **pure Python** (Julia EKP optional later).
- Methods are framed as **samplers** — "which parameter sets next?"
- Components: **CAM, CLM**. Runners: **local, derecho**.
