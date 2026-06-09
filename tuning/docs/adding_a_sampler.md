# Adding a sampler

A sampler answers: *which parameter sets do we run next?* Subclass `Sampler`,
implement four methods, and register a name.

```python
from tuning.core.interfaces import Sampler
from tuning.core.parameters import CalibrationResult
from tuning.core.registry import register_sampler


@register_sampler("my_sampler")
class MySampler(Sampler):
    def ask(self):
        # return the next wave: a list of param dicts, e.g. [{"x": 1.0}, {"x": 2.0}]
        ...

    def tell(self, params, metrics):
        # learn from this wave's metrics (numpy arrays) to choose the next one
        ...

    def is_done(self):
        return True  # one wave, or stop after N waves / on convergence

    def result(self):
        return CalibrationResult(best_params={"x": 1.0})
```

Put it in `src/tuning/samplers/my_sampler.py` and import it in
`samplers/__init__.py` so it registers. See `samplers/random.py` for a small
working example, or `samplers/history_matching.py` for the emulator pattern.
