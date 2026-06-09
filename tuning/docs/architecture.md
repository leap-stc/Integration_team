# Architecture

One question drives everything: **which parameter sets do we run next?**
A *sampler* answers it. The rest is plumbing.

```
           ask()              run()             compute_metrics()       tell()
 Sampler ------->  params -->  Runner -->  output -->  Component -->  numbers -->  Sampler
    ^                                                                                |
    +----------------------- repeat each "wave" until is_done() ---------------------+
                              then result()
```

- **Sampler** — how to choose the next wave: `random`, `latin_hypercube`,
  `eki`, `history_matching`, `mcmc`, `bayesopt`. Some run one wave (random,
  LHC), others many (eki, history matching, mcmc). Some build an emulator
  internally — that's their business, not the loop's.
- **Component** — the model: `cam`, `clm`. Knows its parameters and diagnostics.
- **Runner** — where it runs: `local` (tests) or `derecho` (HPC).

A "wave" is one batch of forward-model runs.

| File | What |
|---|---|
| `core/interfaces.py` | the three contracts (start here) |
| `core/parameters.py` | `Parameter`, `ParameterSet`, `CalibrationResult` |
| `core/observations.py` | targets + uncertainty + `loss()` |
| `core/registry.py` | look up plugins by name |
| `orchestration/loop.py` | the loop above |

## Why these folders

Each folder is the answer to "what must *not* change when I add a plugin?"

- **`core/`** — the shared language. The data types (`Parameter`,
  `Observations`, ...) and the three contracts that every plugin agrees on.
  Mix-and-match only works because a sampler and a component describe a
  "parameter set" the same way. `core/` knows nothing about CESM, EKI, or HPC —
  that's what keeps it stable. (Think wall sockets: one shape, so anything
  plugs into anything.)
- **`orchestration/`** — the conductor. `loop.py` is the *order* the pieces
  fire in (`ask → run → score → tell → repeat`). One place describes the whole
  workflow, and it never changes when you add a sampler or component.
- **`components/`** — the only CESM-aware code. A component translates between
  plain numbers and a real model: `apply()` writes params into a CLM/CAM case,
  `compute_metrics()` turns model output into diagnostics. This translation is
  inherently model-specific, so isolating it lets the *same* sampler calibrate
  CLM today and CAM tomorrow.

Payoff: a new method = one file in `samplers/`; a new model = one file in
`components/`. `core/` and `orchestration/` stay fixed.
