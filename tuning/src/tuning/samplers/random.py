"""Random sampling: N random parameter sets within bounds. One wave, then done.

The simplest sampler — a good baseline and a working example of the interface.
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.parameters import ParameterSet
from ..core.registry import register_sampler
from ._select import best


@register_sampler("random")
class RandomSampler(Sampler):
    def __init__(self, parameters: ParameterSet, n: int = 20, observations=None, seed: int = 0):
        self.parameters = parameters
        self.n = n
        self.observations = observations
        self.rng = np.random.default_rng(seed)
        self.archive = []        # list of (params, metrics) we have run
        self.waves_done = 0

    def ask(self):
        return [
            {p.name: float(self.rng.uniform(p.low, p.high)) for p in self.parameters.params}
            for _ in range(self.n)
        ]

    def tell(self, params, metrics):
        self.archive = list(zip(params, metrics))
        self.waves_done += 1

    def is_done(self):
        return self.waves_done >= 1   # a single wave

    def result(self):
        return best(self.archive, self.observations)
