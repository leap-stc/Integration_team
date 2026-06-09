# tuning — mix-and-match calibration tools for CESM

Calibrate CESM by combining a **sampler** (how to pick parameter sets) with a
**component** (which model) and a **runner** (where it runs). Swap any one
without touching the others.

The whole idea in one question: *given the runs so far, which parameter sets do
we run next?* That answer is a sampler.

```
sampler.ask()  ->  runner.run()  ->  component.compute_metrics()  ->  sampler.tell()
                                                                  (repeat each wave, then result)
```

**Scoping & design:** see [SCOPING.md](SCOPING.md).
**The contract (read this first):** `src/tuning/core/interfaces.py`.
**The flow:** `src/tuning/orchestration/loop.py`.

## Status

The `core/` interfaces and loop work; `random`, `latin_hypercube`, and
`history_matching` samplers are real and tested via the `local` runner. The CESM
components and the `derecho` runner are stubs (need a real model environment).

| Slot | Options |
|---|---|
| Samplers | `random` ✅ · `latin_hypercube` ✅ · `history_matching` ✅ · `eki` · `mcmc` · `bayesopt` |
| Components | `cam` · `clm` |
| Runners | `local` ✅ · `derecho` |

## Install

```bash
pip install -e tuning/                 # core
pip install -e "tuning/[emulator]"     # + scikit-learn, for history_matching
```

## Try it

```bash
pytest tuning/tests        # runs the loop end-to-end with the random sampler
```

## Add your own

- New sampler: [docs/adding_a_sampler.md](docs/adding_a_sampler.md)
- New component: [docs/adding_a_component.md](docs/adding_a_component.md)
- Which sampler to use: [docs/samplers_overview.md](docs/samplers_overview.md)
