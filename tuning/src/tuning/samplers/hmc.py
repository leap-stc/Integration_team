"""Hamiltonian Monte Carlo (HMC) on a GP emulator.

Goal: the *posterior* over parameters — the full range of values consistent
with the observations, not just a single best fit.

How it works:
1. The first wave is a Latin-hypercube design over the ranges.
2. We fit a GP emulator (a fast stand-in for the model) to those runs.
3. HMC samples the posterior on the emulator. It treats the parameters as a
   particle: give it a random momentum, simulate frictionless motion over the
   (negative log-)posterior landscape for a few "leapfrog" steps, then accept or
   reject the move so the energy stays honest. This explores efficiently.
4. Optionally, draw N samples from the posterior and run one more forward-model
   wave to check them.

All HMC math runs in a normalized [0, 1] box (each parameter mapped onto its
range) so the step size means the same thing for every parameter.
The emulator is the same GP used by history matching (see _emulator.py).
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.parameters import ParameterSet
from ..core.registry import register_sampler
from ._emulator import GPEmulator
from ._select import best
from .latin_hypercube import _latin_hypercube


@register_sampler("hmc")
class HamiltonianMonteCarlo(Sampler):
    def __init__(self, parameters: ParameterSet, observations, n=30,
                 posterior_samples=500, step_size=0.05, leapfrog_steps=15,
                 burn_in=200, validation_wave=True, emulator=None, seed=0):
        self.parameters = parameters
        self.observations = observations
        self.n = n                              # forward-model ensemble size per wave
        self.posterior_samples = posterior_samples
        self.step_size = step_size
        self.leapfrog_steps = leapfrog_steps
        self.burn_in = burn_in
        self.total_waves = 2 if validation_wave else 1
        self.emulator = emulator or GPEmulator()
        self.rng = np.random.default_rng(seed)
        self.archive = []                        # (params, metrics) for every run
        self.posterior = None                    # array [posterior_samples, n_params]
        self.waves_done = 0

    def ask(self):
        if self.waves_done == 0:
            return self._to_dicts(self._lhc(self.n))           # design wave
        if self.waves_done == 1 and self.total_waves == 2:
            draws = self.rng.choice(len(self.posterior), size=self.n)
            return self._to_dicts(self.posterior[draws])        # validate posterior draws
        return []

    def tell(self, params, metrics):
        self.archive += list(zip(params, metrics))
        self.emulator.fit(self._to_array([p for p, _ in self.archive]),
                          np.array([m for _, m in self.archive]))
        self.posterior = self._sample_posterior()
        self.waves_done += 1

    def is_done(self):
        return self.waves_done >= self.total_waves

    def result(self):
        res = best(self.archive, self.observations)
        res.extras["posterior"] = self._to_dicts(self.posterior)
        res.extras["posterior_mean"] = self._to_dicts(self.posterior.mean(0, keepdims=True))[0]
        return res

    # --- HMC on the emulator (all in the normalized [0, 1] box) ---

    def _sample_posterior(self):
        u = np.full(len(self.parameters.params), 0.5)    # start at the box centre
        samples = []
        for i in range(self.burn_in + self.posterior_samples):
            u = self._hmc_step(u)
            if i >= self.burn_in:
                samples.append(u)
        return self._unit_to_theta(np.array(samples))

    def _hmc_step(self, u):
        momentum = self.rng.normal(size=len(u))
        q, p = u.copy(), momentum.copy()

        # leapfrog: simulate motion over the posterior landscape
        p -= 0.5 * self.step_size * self._grad_potential(q)
        for step in range(self.leapfrog_steps):
            q = q + self.step_size * p
            if step != self.leapfrog_steps - 1:
                p -= self.step_size * self._grad_potential(q)
        p -= 0.5 * self.step_size * self._grad_potential(q)

        # accept or reject so total energy (potential + kinetic) is preserved
        start = self._potential(u) + 0.5 * momentum @ momentum
        end = self._potential(q) + 0.5 * p @ p
        if np.log(self.rng.uniform()) < start - end:
            return q
        return u

    def _potential(self, u):
        """Negative log-posterior. Infinite outside the box (a uniform prior)."""
        if np.any(u < 0) or np.any(u > 1):
            return np.inf
        return -self._log_likelihood(u[None])[0]

    def _grad_potential(self, u):
        """Gradient of the potential, by central differences on the emulator."""
        eps = 1e-3
        points, p = [], len(u)
        for j in range(p):
            up, down = u.copy(), u.copy()
            up[j] += eps
            down[j] -= eps
            points += [up, down]
        loglik = self._log_likelihood(np.array(points))   # one emulator call
        grad = np.empty(p)
        for j in range(p):
            grad[j] = -(loglik[2 * j] - loglik[2 * j + 1]) / (2 * eps)
        return grad

    def _log_likelihood(self, U):
        """Gaussian log-likelihood at unit points U [k, p] using the emulator."""
        mean, std = self.emulator.predict(self._unit_to_theta(U))   # [k, n_outputs]
        var = self.observations.uncertainty ** 2 + std ** 2
        diff = self.observations.targets - mean
        return -0.5 * np.sum(diff ** 2 / var + np.log(2 * np.pi * var), axis=1)

    # --- helpers ---

    def _bounds(self):
        lows = np.array([p.low for p in self.parameters.params])
        highs = np.array([p.high for p in self.parameters.params])
        return lows, highs

    def _unit_to_theta(self, U):
        lows, highs = self._bounds()
        return lows + U * (highs - lows)

    def _lhc(self, k):
        lows, highs = self._bounds()
        return lows + _latin_hypercube(k, len(lows), self.rng) * (highs - lows)

    def _to_array(self, param_dicts):
        names = self.parameters.names()
        return np.array([[d[name] for name in names] for d in param_dicts])

    def _to_dicts(self, array):
        names = self.parameters.names()
        return [{name: float(v) for name, v in zip(names, row)} for row in array]
