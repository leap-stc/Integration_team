"""History matching should narrow a big parameter space down toward the target.

Toy model: the metrics are just the parameters echoed back, so the target is
known (a=3, b=7) and we can check the search finds it — no CESM needed.
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.loop import run_calibration
from tuning.runners.local import LocalRunner
from tuning.samplers.history_matching import HistoryMatching


class EchoComponent(Component):
    """Toy model whose metrics are just the parameters echoed back."""

    def parameters(self):
        return ParameterSet([Parameter("a", 0, 10, 5), Parameter("b", 0, 10, 5)])

    def apply(self, case_dir, params):
        pass

    def compute_metrics(self, run_output):
        return run_output


def test_history_matching_narrows_to_target():
    component = EchoComponent()
    obs = Observations(targets=np.array([3.0, 7.0]), uncertainty=np.array([0.5, 0.5]))
    sampler = HistoryMatching(component.parameters(), obs, n=30, waves=3, pool=2000, seed=0)
    runner = LocalRunner(model_fn=lambda p: np.array([p["a"], p["b"]]))

    result = run_calibration(sampler, component, runner)

    assert abs(result.best_params["a"] - 3.0) < 1.5
    assert abs(result.best_params["b"] - 7.0) < 1.5
    assert len(result.extras["nroy"]) < 2000   # most of the space was ruled out
