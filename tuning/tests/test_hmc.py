"""HMC's posterior should center on the observed target — no CESM needed.

Toy model echoes the parameters, so the target (a=3, b=7) is known and we can
check the posterior mean lands near it and has the expected number of samples.
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.loop import run_calibration
from tuning.runners.local import LocalRunner
from tuning.samplers.hmc import HamiltonianMonteCarlo


class EchoComponent(Component):
    def parameters(self):
        return ParameterSet([Parameter("a", 0, 10, 5), Parameter("b", 0, 10, 5)])

    def apply(self, case_dir, params):
        pass

    def compute_metrics(self, run_output):
        return run_output


def test_hmc_posterior_centers_on_target():
    component = EchoComponent()
    obs = Observations(targets=np.array([3.0, 7.0]), uncertainty=np.array([0.5, 0.5]))
    sampler = HamiltonianMonteCarlo(
        component.parameters(), obs, n=25, posterior_samples=150, burn_in=50,
        leapfrog_steps=12, step_size=0.06, validation_wave=False, seed=0,
    )
    runner = LocalRunner(model_fn=lambda p: np.array([p["a"], p["b"]]))

    result = run_calibration(sampler, component, runner)

    mean = result.extras["posterior_mean"]
    assert abs(mean["a"] - 3.0) < 1.5
    assert abs(mean["b"] - 7.0) < 1.5
    assert len(result.extras["posterior"]) == 150
