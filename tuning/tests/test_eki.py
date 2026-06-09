"""EKI should pull an ensemble toward the observed target — no CESM needed.

Toy model echoes the parameters, so the target (a=3, b=7) is known and we can
check the ensemble mean lands near it.
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.loop import run_calibration
from tuning.runners.local import LocalRunner
from tuning.samplers.eki import EKI


class EchoComponent(Component):
    def parameters(self):
        return ParameterSet([Parameter("a", 0, 10, 5), Parameter("b", 0, 10, 5)])

    def apply(self, case_dir, params):
        pass

    def compute_metrics(self, run_output):
        return run_output


def test_eki_converges_to_target():
    component = EchoComponent()
    obs = Observations(targets=np.array([3.0, 7.0]), uncertainty=np.array([0.5, 0.5]))
    sampler = EKI(component.parameters(), obs, n=50, iterations=10, seed=0)
    runner = LocalRunner(model_fn=lambda p: np.array([p["a"], p["b"]]))

    result = run_calibration(sampler, component, runner)

    mean = result.extras["ensemble_mean"]
    assert abs(mean["a"] - 3.0) < 1.0
    assert abs(mean["b"] - 7.0) < 1.0
