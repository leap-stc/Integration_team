"""The calibration loop — the heart of the package.

One wave at a time: the sampler proposes parameter sets, the runner runs them,
the component scores them, the sampler learns. Repeat until the sampler is done.
"""

from ..core.interfaces import Component, Runner, Sampler
from ..core.parameters import CalibrationResult


def run_calibration(
    sampler: Sampler,
    component: Component,
    runner: Runner,
    max_waves: int = 100,
) -> CalibrationResult:
    ensemble = sampler.ask()                            # 1. choose parameter sets
    waves = 0
    while not sampler.is_done() and waves < max_waves:
        outputs = runner.run(ensemble, component)       # 2. run the model
        metrics = [component.compute_metrics(o) for o in outputs]  # 3. score
        sampler.tell(ensemble, metrics)                 # 4. sampler learns
        ensemble = sampler.ask()                        # 5. choose next wave
        waves += 1
    return sampler.result()
