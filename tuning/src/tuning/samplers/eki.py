"""Ensemble Kalman Inversion (EKI), in Python.  [STUB — Phase 3]

Picks each next wave by nudging the whole ensemble toward the observations with
a Kalman update. Derivative-free and parallel. Pure Python (no Julia).
"""

from ..core.interfaces import Sampler
from ..core.registry import register_sampler

_TODO = "eki is a stub — see SCOPING.md Phase 3."


@register_sampler("eki")
class EKI(Sampler):
    def __init__(self, *args, **kwargs):
        pass

    def ask(self):
        raise NotImplementedError(_TODO)

    def tell(self, params, metrics):
        raise NotImplementedError(_TODO)

    def is_done(self):
        raise NotImplementedError(_TODO)

    def result(self):
        raise NotImplementedError(_TODO)
