"""Bayesian optimization (wraps an ESEm GP emulator).  [STUB — Phase 3]

Fits a GP, then picks the next parameter sets most likely to improve the fit.
Good when each forward-model run is expensive and you have few parameters.
"""

from ..core.interfaces import Sampler
from ..core.registry import register_sampler

_TODO = "bayesopt is a stub — see SCOPING.md Phase 3."


@register_sampler("bayesopt")
class BayesianOptimization(Sampler):
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
