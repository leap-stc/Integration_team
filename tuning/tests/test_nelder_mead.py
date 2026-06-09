"""Nelder-Mead should walk its simplex to the loss minimum — no CESM needed.

Toy model echoes the parameters, so the minimum is at the target (a=3, b=7).
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.loop import run_calibration
from tuning.runners.local import LocalRunner
from tuning.samplers.nelder_mead import NelderMead


class EchoComponent(Component):
    def parameters(self):
        return ParameterSet([Parameter("a", 0, 10, 5), Parameter("b", 0, 10, 5)])

    def apply(self, case_dir, params):
        pass

    def compute_metrics(self, run_output):
        return run_output


def test_nelder_mead_finds_minimum():
    component = EchoComponent()
    obs = Observations(targets=np.array([3.0, 7.0]), uncertainty=np.array([0.5, 0.5]))
    sampler = NelderMead(component.parameters(), obs, max_iterations=80)
    runner = LocalRunner(model_fn=lambda p: np.array([p["a"], p["b"]]))

    result = run_calibration(sampler, component, runner, max_waves=2000)

    assert abs(result.best_params["a"] - 3.0) < 0.2
    assert abs(result.best_params["b"] - 7.0) < 0.2
