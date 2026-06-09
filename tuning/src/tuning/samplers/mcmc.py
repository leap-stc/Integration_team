"""MCMC / HMC sampling (wraps an ESEm emulator).  [STUB — Phase 3]

Fits an emulator, then samples the posterior over parameters with MCMC (or
Hamiltonian Monte Carlo). result() returns that posterior; you can then draw N
sets from it and run one more forward-model wave.
"""

from ..core.interfaces import Sampler
from ..core.registry import register_sampler

_TODO = "mcmc is a stub — see SCOPING.md Phase 3."


@register_sampler("mcmc")
class MCMC(Sampler):
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
