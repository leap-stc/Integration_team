"""CAM atmosphere component.  [STUB — Phase 3]

parameters() shows the shape with one example param. apply() will write CAM
namelist/xml settings; compute_metrics() will build climate diagnostics.
"""

from ..core.interfaces import Component
from ..core.parameters import Parameter, ParameterSet
from ..core.registry import register_component


@register_component("cam")
class CAM(Component):
    def parameters(self) -> ParameterSet:
        return ParameterSet([
            Parameter("zmconv_c0_lnd", low=0.001, high=0.06, default=0.0059),
        ])

    def apply(self, case_dir, params):
        raise NotImplementedError("CAM.apply is a stub — see SCOPING.md Phase 3.")

    def compute_metrics(self, run_output):
        raise NotImplementedError("CAM.compute_metrics is a stub — see SCOPING.md Phase 3.")
