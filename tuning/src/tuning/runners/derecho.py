"""Runner that submits a CESM ensemble to Derecho (PBS).  [STUB — Phase 3]

Will: create a CESM case per parameter set, call component.apply() to write
the params, submit the jobs, wait, then return each run's output path.
"""

from ..core.interfaces import Runner
from ..core.registry import register_runner


@register_runner("derecho")
class DerechoRunner(Runner):
    def run(self, ensemble, component):
        raise NotImplementedError("DerechoRunner is a stub — see SCOPING.md Phase 3.")
