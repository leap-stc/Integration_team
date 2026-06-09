"""Run the CLM history-matching experiment — one wave per run.

Each time you run this it does ONE wave: it collects the previous wave's results
(if any), then submits the next wave to Derecho and exits without waiting. Wait
for that wave to finish, then run it again. Repeat until it prints the best
parameters. State is saved in campaign.pkl between runs.

DRAFT: needs a live CESM environment and the paths in config.yaml (see README).
For a runnable demo with no CESM, see tuning/tests/.
"""

from pathlib import Path

import numpy as np

from tuning.core.config import load_config
from tuning.core.observations import Observations
from tuning.core.registry import get_component, get_runner, get_sampler
from tuning.orchestration.campaign import Campaign

here = Path(__file__).parent
cfg = load_config(here / "config.yaml")

component = get_component(cfg["component"])(**cfg["clm"])
runner = get_runner(cfg["runner"])(**cfg["derecho"])
obs = Observations(
    targets=np.array(cfg["observations"]["targets"]),
    uncertainty=np.array(cfg["observations"]["uncertainty"]),
)
sampler = get_sampler(cfg["sampler"])(
    component.parameters(), obs, n=cfg["n"], waves=cfg["waves"]
)

campaign = Campaign(sampler, component, runner, state_file=here / "campaign.pkl")
result = campaign.step()

if result is None:
    print("Wave submitted. Re-run this script after it finishes to start the next wave.")
else:
    print("Done. Best parameters:", result.best_params)
