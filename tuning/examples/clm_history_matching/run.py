"""Wire up an experiment from config.yaml and run it.

The clm / history_matching / derecho plugins are stubs, so this raises
NotImplementedError until Phases 2-3 fill them in. For a runnable end-to-end
demo today, see tuning/tests/test_loop.py.
"""

from pathlib import Path

from tuning.core.config import load_config
from tuning.core.registry import get_component, get_runner, get_sampler
from tuning.orchestration.loop import run_calibration

cfg = load_config(Path(__file__).parent / "config.yaml")

component = get_component(cfg["component"])()
sampler = get_sampler(cfg["sampler"])(component.parameters(), n=cfg["n"])
runner = get_runner(cfg["runner"])()

result = run_calibration(sampler, component, runner)
print(result.best_params)
