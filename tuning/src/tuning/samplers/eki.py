"""Ensemble Kalman Inversion (EKI), in Python.

A derivative-free calibration method: keep an ensemble of parameter sets, run
them all, then nudge the whole ensemble toward the observations using the
ensemble's own covariances (a Kalman update). Repeat for a few iterations.

Each wave = one iteration = running the whole ensemble once.
Based on Iglesias, Law & Stuart (2013).
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.parameters import ParameterSet
from ..core.registry import register_sampler
from ._select import best


@register_sampler("eki")
class EKI(Sampler):
    def __init__(self, parameters: ParameterSet, observations, n=50, iterations=10, seed=0):
        self.parameters = parameters
        self.observations = observations
        self.n = n                        # ensemble size
        self.max_iterations = iterations
        self.rng = np.random.default_rng(seed)
        self.ensemble = None              # current parameter sets, array [n, n_params]
        self.archive = []                 # (params, metrics) for every run
        self.iteration = 0

    def ask(self):
        if self.ensemble is None:
            self.ensemble = self._sample_prior()   # first wave: spread over the ranges
        return self._to_dicts(self.ensemble)

    def tell(self, params, metrics):
        self.archive += list(zip(params, metrics))
        self.ensemble = self._kalman_update(self._to_array(params), np.array(metrics))
        self.iteration += 1

    def is_done(self):
        return self.iteration >= self.max_iterations

    def result(self):
        res = best(self.archive, self.observations)
        res.extras["ensemble_mean"] = self._to_dicts(self.ensemble.mean(axis=0, keepdims=True))[0]
        return res

    # --- the Kalman update: nudge the ensemble toward the observations ---

    def _kalman_update(self, theta, g):
        y = self.observations.targets
        gamma = np.diag(self.observations.uncertainty ** 2)   # observation noise

        d_theta = theta - theta.mean(axis=0)      # spread in parameters
        d_g = g - g.mean(axis=0)                  # spread in outputs
        cov_tg = d_theta.T @ d_g / self.n         # how params co-vary with outputs
        cov_gg = d_g.T @ d_g / self.n             # how outputs co-vary

        gain = cov_tg @ np.linalg.pinv(cov_gg + gamma)        # Kalman gain
        noise = self.rng.multivariate_normal(np.zeros(len(y)), gamma, size=self.n)
        updated = theta + (y + noise - g) @ gain.T            # move each particle toward y

        lows, highs = self._bounds()
        return np.clip(updated, lows, highs)      # keep particles inside the ranges

    # --- helpers ---

    def _bounds(self):
        lows = np.array([p.low for p in self.parameters.params])
        highs = np.array([p.high for p in self.parameters.params])
        return lows, highs

    def _sample_prior(self):
        lows, highs = self._bounds()
        return lows + self.rng.uniform(size=(self.n, len(lows))) * (highs - lows)

    def _to_array(self, param_dicts):
        names = self.parameters.names()
        return np.array([[p[name] for name in names] for p in param_dicts])

    def _to_dicts(self, array):
        names = self.parameters.names()
        return [{name: float(v) for name, v in zip(names, row)} for row in array]
