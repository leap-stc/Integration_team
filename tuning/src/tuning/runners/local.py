"""Runner that runs the model in this Python process — for testing the loop.

You give it a `model_fn(params) -> run_output`. No HPC needed.
For real CESM runs, use DerechoRunner instead.
"""

from ..core.interfaces import Runner
from ..core.registry import register_runner


@register_runner("local")
class LocalRunner(Runner):
    def __init__(self, model_fn):
        self.model_fn = model_fn

    def run(self, ensemble, component):
        return [self.model_fn(params) for params in ensemble]
