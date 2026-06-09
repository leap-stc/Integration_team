"""CLM/CTSM land component.  [STUB — Phase 2, built from NCAR/ctsm6_ppe]

parameters() shows the shape with two example params. apply() and
compute_metrics() will reuse ctsm6_ppe's paramfile writing and diagnostics.
"""

from ..core.interfaces import Component
from ..core.parameters import Parameter, ParameterSet
from ..core.registry import register_component


@register_component("clm")
class CLM(Component):
    def parameters(self) -> ParameterSet:
        return ParameterSet([
            Parameter("medlynslope", low=1.0, high=6.0, default=4.1),
            Parameter("fff", low=0.02, high=5.0, default=0.5),
        ])

    def apply(self, case_dir, params):
        raise NotImplementedError("CLM.apply is a stub — see SCOPING.md Phase 2.")

    def compute_metrics(self, run_output):
        raise NotImplementedError("CLM.compute_metrics is a stub — see SCOPING.md Phase 2.")
