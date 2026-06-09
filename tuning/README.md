# tuning — mix-and-match calibration tools for CESM

Pick a **sampler** (how to choose parameters), a **component** (which model), and
a **runner** (where it runs), and calibrate. Swap any one without touching the
others.

> One question drives everything: *given the runs so far, which parameter sets
> do we run next?* A **sampler** answers it.

## Install

```bash
pip install -e tuning/                 # core (numpy, pyyaml)
pip install -e "tuning/[emulator]"     # + scikit-learn (history_matching, hmc)
```

## 30-second example (no CESM needed)

Calibrate a toy "model" that echoes its parameter, to hit a target of 3:

```python
import numpy as np
from tuning import Observations, run_calibration
from tuning.core.interfaces import Component
from tuning.core.parameters import Parameter, ParameterSet
from tuning.runners.local import LocalRunner
from tuning.samplers.nelder_mead import NelderMead

class Echo(Component):                          # a stand-in "model"
    def parameters(self):
        return ParameterSet([Parameter("x", low=0, high=10, default=5)])
    def apply(self, case_dir, params): pass     # would write params into a CESM case
    def compute_metrics(self, out): return out  # would read diagnostics from a run

component = Echo()
obs = Observations(targets=np.array([3.0]), uncertainty=np.array([0.5]))
sampler = NelderMead(component.parameters(), obs)
runner = LocalRunner(model_fn=lambda p: np.array([p["x"]]))

result = run_calibration(sampler, component, runner, max_waves=500)
print(result.best_params)                       # -> {'x': ~3.0}
```

Swap `NelderMead` for any other sampler, or `Echo`/`LocalRunner` for a real CESM
component/runner — the rest stays the same.

## How it works

```
sampler.ask()  ->  runner.run()  ->  component.compute_metrics()  ->  sampler.tell()
                                                       (repeat each "wave", then result)
```

- **Sampler** — which parameters to try next.
- **Component** — the model: writes parameters in, reads diagnostics out.
- **Runner** — where it runs: `local` (in-process) or `derecho` (HPC).

Two ways to drive it: `run_calibration` (synchronous, for `local`) and
`Campaign.step` (one wave per restart, for long HPC runs). See
[SCOPING.md](SCOPING.md).

Read next: `src/tuning/core/interfaces.py` (the contracts) and
`src/tuning/orchestration/loop.py` (the flow).

## What's available

| Slot | Options |
|---|---|
| Samplers (✅ tested) | `random` · `latin_hypercube` · `history_matching` · `eki` · `hmc` · `nelder_mead` |
| Components | `clm` (draft) · `cam` (stub) — need a live CESM env |
| Runners | `local` ✅ · `derecho` (draft) |

Which sampler should I use? → [docs/samplers_overview.md](docs/samplers_overview.md).

## Test

```bash
cd tuning && pytest        # 6 tests, ~1s, no CESM needed
```

## Extend

- Add a sampler: [docs/adding_a_sampler.md](docs/adding_a_sampler.md)
- Add a component: [docs/adding_a_component.md](docs/adding_a_component.md)
- Design & roadmap: [SCOPING.md](SCOPING.md)
