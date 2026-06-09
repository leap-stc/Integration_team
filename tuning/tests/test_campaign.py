"""Campaign should advance one wave per step and survive a restart.

We simulate restarts by building a fresh Campaign for each step that points at
the same state file — exactly what happens when you re-run the script per wave.
"""

import numpy as np

from tuning.core.interfaces import Component
from tuning.core.observations import Observations
from tuning.core.parameters import Parameter, ParameterSet
from tuning.orchestration.campaign import Campaign
from tuning.runners.local import LocalRunner
from tuning.samplers.random import RandomSampler


class ToyComponent(Component):
    def parameters(self):
        return ParameterSet([Parameter("x", 0.0, 10.0, 5.0)])

    def apply(self, case_dir, params):
        pass

    def compute_metrics(self, run_output):
        return run_output


def test_campaign_steps_across_restarts(tmp_path):
    state = str(tmp_path / "campaign.pkl")
    component = ToyComponent()
    obs = Observations(targets=np.array([3.0]), uncertainty=np.array([1.0]))
    runner = LocalRunner(model_fn=lambda p: np.array([p["x"]]))

    def fresh_campaign():
        sampler = RandomSampler(component.parameters(), n=100, observations=obs, seed=0)
        return Campaign(sampler, component, runner, state_file=state)

    # first step submits wave 1 — no result yet
    assert fresh_campaign().step() is None
    # second step (a "restart") collects wave 1 and finalizes
    result = fresh_campaign().step()
    assert result is not None
    assert abs(result.best_params["x"] - 3.0) < 0.5
