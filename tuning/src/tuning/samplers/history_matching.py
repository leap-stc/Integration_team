"""History matching: emulate, rule out, refocus.

How it works, wave by wave:
1. The first wave is a Latin-hypercube design over the full ranges.
2. After each wave we fit an emulator to every run so far.
3. We score a big pool of candidate parameter sets by their *implausibility* —
   how far the emulator says each one lands from observations, allowing for both
   emulator and observation uncertainty.
4. Candidates below the threshold form the "not ruled out yet" (NROY) region.
   The next wave is drawn from there, so the search refocuses each wave.

The emulator is swappable (see _emulator.py); the default uses scikit-learn GPs.
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.parameters import ParameterSet
from ..core.registry import register_sampler
from ._emulator import GPEmulator
from ._select import best
from .latin_hypercube import _latin_hypercube


@register_sampler("history_matching")
class HistoryMatching(Sampler):
    def __init__(self, parameters: ParameterSet, observations, n=50, waves=3,
                 threshold=3.0, pool=5000, emulator=None, seed=0):
        self.parameters = parameters
        self.observations = observations
        self.n = n                       # parameter sets per wave
        self.max_waves = waves
        self.threshold = threshold       # implausibility cutoff (3 is standard)
        self.pool = pool                 # candidates scored by the emulator each wave
        self.emulator = emulator or GPEmulator()
        self.rng = np.random.default_rng(seed)
        self.archive = []                # (params, metrics) for every run so far
        self.waves_done = 0
        self.nroy = None                 # current "not ruled out yet" parameter sets

    def ask(self):
        if self.waves_done == 0:
            return self._to_dicts(self._lhc(self.n))   # first wave: spread out
        if not self.nroy:
            return []                                   # nothing plausible left; loop stops
        pick = self.rng.choice(len(self.nroy), size=self.n)
        return [self.nroy[i] for i in pick]             # later waves: draw from NROY

    def tell(self, params, metrics):
        self.archive += list(zip(params, metrics))
        self.waves_done += 1
        self._refocus()

    def is_done(self):
        return self.waves_done >= self.max_waves or self.nroy == []

    def result(self):
        res = best(self.archive, self.observations)
        res.extras["nroy"] = self.nroy
        return res

    # --- the history-matching step ---

    def _refocus(self):
        """Fit the emulator on all runs, then recompute the NROY region."""
        X = self._to_array([p for p, _ in self.archive])
        Y = np.array([m for _, m in self.archive])
        self.emulator.fit(X, Y)

        candidates = self._uniform(self.pool)
        mean, std = self.emulator.predict(candidates)
        keep = self._implausibility(mean, std) < self.threshold
        self.nroy = self._to_dicts(candidates[keep])

    def _implausibility(self, mean, std):
        """Worst-case standardized distance from observations, over all outputs."""
        obs_var = self.observations.uncertainty ** 2
        distance = np.abs(mean - self.observations.targets) / np.sqrt(std ** 2 + obs_var)
        return distance.max(axis=1)

    # --- converting between param dicts and plain arrays ---

    def _to_array(self, param_dicts):
        names = self.parameters.names()
        return np.array([[p[name] for name in names] for p in param_dicts])

    def _to_dicts(self, array):
        names = self.parameters.names()
        return [{name: float(v) for name, v in zip(names, row)} for row in array]

    def _lows_highs(self):
        lows = np.array([p.low for p in self.parameters.params])
        highs = np.array([p.high for p in self.parameters.params])
        return lows, highs

    def _uniform(self, k):
        lows, highs = self._lows_highs()
        return lows + self.rng.uniform(size=(k, len(lows))) * (highs - lows)

    def _lhc(self, k):
        lows, highs = self._lows_highs()
        unit = _latin_hypercube(k, len(lows), self.rng)
        return lows + unit * (highs - lows)
