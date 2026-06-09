"""Latin hypercube sampling: spread N parameter sets evenly across the ranges.

Better coverage than random for the same N. One wave, then done. Often the
first wave before fitting an emulator (history matching).
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.parameters import ParameterSet
from ..core.registry import register_sampler
from ._select import best


@register_sampler("latin_hypercube")
class LatinHypercubeSampler(Sampler):
    def __init__(self, parameters: ParameterSet, n: int = 20, observations=None, seed: int = 0):
        self.parameters = parameters
        self.n = n
        self.observations = observations
        self.rng = np.random.default_rng(seed)
        self.archive = []
        self.waves_done = 0

    def ask(self):
        unit = _latin_hypercube(self.n, len(self.parameters.params), self.rng)  # in [0, 1]
        sets = []
        for row in unit:
            params = {
                p.name: float(p.low + value * (p.high - p.low))
                for value, p in zip(row, self.parameters.params)
            }
            sets.append(params)
        return sets

    def tell(self, params, metrics):
        self.archive = list(zip(params, metrics))
        self.waves_done += 1

    def is_done(self):
        return self.waves_done >= 1

    def result(self):
        return best(self.archive, self.observations)


def _latin_hypercube(n, d, rng):
    """n points in the d-dimensional unit cube — one sample per equal bin, per axis."""
    points = np.empty((n, d))
    for j in range(d):
        bins = (np.arange(n) + rng.uniform(size=n)) / n  # one point in each 1/n bin
        rng.shuffle(bins)                                # mix axes independently
        points[:, j] = bins
    return points
