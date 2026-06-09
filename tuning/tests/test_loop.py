"""End-to-end test of the loop with a real sampler — no CESM or HPC needed.

A tiny "model" echoes a parameter x; the random sampler searches for the x that
matches the target. Read this to see how the pieces fit together.
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.loop import run_calibration
from tuning.runners.local import LocalRunner
from tuning.samplers.random import RandomSampler


class ToyComponent(Component):
    """A 'model' with one parameter x; its metric is just x itself."""

    def parameters(self):
        return ParameterSet([Parameter("x", 0.0, 10.0, 5.0)])

    def apply(self, case_dir, params):
        pass  # nothing to write for a toy model

    def compute_metrics(self, run_output):
        return run_output  # the local runner already returns an array


def test_loop_finds_target():
    component = ToyComponent()
    obs = Observations(targets=np.array([3.0]), uncertainty=np.array([1.0]))
    sampler = RandomSampler(component.parameters(), n=200, observations=obs, seed=0)
    runner = LocalRunner(model_fn=lambda p: np.array([p["x"]]))

    result = run_calibration(sampler, component, runner)

    assert abs(result.best_params["x"] - 3.0) < 0.3
