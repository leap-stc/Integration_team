"""Nelder-Mead simplex optimization.

A classic derivative-free optimizer. It keeps a "simplex" of n+1 parameter sets
(for n parameters) and each step reflects the worst one through the middle of
the others toward lower loss — then expands, contracts, or shrinks depending on
what it finds. No emulator, no gradients.

It proposes only a few points at a time (usually one), so each wave is small,
and it returns a single best parameter set. The sequential reflect/expand/
contract/shrink logic is run as a small state machine: each tell() decides the
next trial point, each ask() hands it out.
"""

import numpy as np

from ..core.interfaces import Sampler
from ..core.registry import register_sampler
from ._select import best

# standard Nelder-Mead coefficients: reflect, expand, contract, shrink
ALPHA, GAMMA, RHO, SIGMA = 1.0, 2.0, 0.5, 0.5


@register_sampler("nelder_mead")
class NelderMead(Sampler):
    def __init__(self, parameters, observations, max_iterations=100, step=0.1, ftol=1e-6):
        self.parameters = parameters
        self.observations = observations
        self.max_iterations = max_iterations
        self.step = step                  # initial simplex size, as a fraction of each range
        self.ftol = ftol
        self.simplex = None               # list of [point array, loss]
        self.archive = []
        self.iterations = 0
        self.phase = "init"               # what tell() should do with the next results
        self.pending = self._initial_simplex()
        self._scratch = {}                # carries centroid / x_r between sub-steps

    # --- Sampler interface ---

    def ask(self):
        return self._to_dicts(self.pending)

    def tell(self, params, metrics):
        points = self._to_array(params)
        losses = [self.observations.loss(m) for m in metrics]
        self.archive += list(zip(params, metrics))
        getattr(self, f"_after_{self.phase}")(points, losses)

    def is_done(self):
        if self.simplex is None:
            return False
        spread = self.simplex[-1][1] - self.simplex[0][1]   # worst loss - best loss
        return self.iterations >= self.max_iterations or spread < self.ftol

    def result(self):
        return best(self.archive, self.observations)

    # --- one handler per phase: update the simplex, then queue the next point ---

    def _after_init(self, points, losses):
        self.simplex = [[p, loss] for p, loss in zip(points, losses)]
        self._start_iteration()

    def _after_reflect(self, points, losses):
        x_r, f_r = points[0], losses[0]
        f_best = self.simplex[0][1]
        f_second = self.simplex[-2][1]
        f_worst = self.simplex[-1][1]
        centroid = self._scratch["centroid"]

        if f_r < f_best:                              # promising — try going further
            x_e = self._clip(centroid + GAMMA * (x_r - centroid))
            self._scratch.update(x_r=x_r, f_r=f_r)
            self._queue([x_e], "expand")
        elif f_r < f_second:                          # good enough — keep it
            self._accept(x_r, f_r)
        else:                                         # poor — pull back (contract)
            if f_r < f_worst:                         # outside contraction
                x_c = self._clip(centroid + RHO * (x_r - centroid))
                self._scratch["contract_ref"] = f_r
            else:                                     # inside contraction
                worst = self.simplex[-1][0]
                x_c = self._clip(centroid + RHO * (worst - centroid))
                self._scratch["contract_ref"] = f_worst
            self._queue([x_c], "contract")

    def _after_expand(self, points, losses):
        x_e, f_e = points[0], losses[0]
        if f_e < self._scratch["f_r"]:
            self._accept(x_e, f_e)
        else:
            self._accept(self._scratch["x_r"], self._scratch["f_r"])

    def _after_contract(self, points, losses):
        x_c, f_c = points[0], losses[0]
        if f_c < self._scratch["contract_ref"]:
            self._accept(x_c, f_c)
        else:
            self._shrink()                            # contraction failed — shrink everything

    def _after_shrink(self, points, losses):
        self.simplex = [self.simplex[0]] + [[p, loss] for p, loss in zip(points, losses)]
        self._end_iteration()

    # --- simplex moves ---

    def _accept(self, point, loss):
        self.simplex[-1] = [point, loss]              # replace the worst vertex
        self._end_iteration()

    def _shrink(self):
        best_point = self.simplex[0][0]
        shrunk = [self._clip(best_point + SIGMA * (vertex[0] - best_point))
                  for vertex in self.simplex[1:]]
        self._queue(shrunk, "shrink")

    def _start_iteration(self):
        """Sort the simplex and queue a reflection of the worst vertex."""
        self.simplex.sort(key=lambda v: v[1])
        points = np.array([v[0] for v in self.simplex])
        centroid = points[:-1].mean(axis=0)           # middle of all but the worst
        x_r = self._clip(centroid + ALPHA * (centroid - points[-1]))
        self._scratch = {"centroid": centroid}
        self._queue([x_r], "reflect")

    def _end_iteration(self):
        self.iterations += 1
        self._start_iteration()

    def _queue(self, points, phase):
        self.pending = np.array(points)
        self.phase = phase

    # --- helpers ---

    def _initial_simplex(self):
        x0 = np.array([p.default for p in self.parameters.params])
        lows, highs = self._bounds()
        vertices = [x0]
        for j in range(len(x0)):
            vertex = x0.copy()
            vertex[j] += self.step * (highs[j] - lows[j])
            vertices.append(self._clip(vertex))
        return np.array(vertices)

    def _bounds(self):
        lows = np.array([p.low for p in self.parameters.params])
        highs = np.array([p.high for p in self.parameters.params])
        return lows, highs

    def _clip(self, point):
        lows, highs = self._bounds()
        return np.clip(point, lows, highs)

    def _to_array(self, param_dicts):
        names = self.parameters.names()
        return np.array([[d[name] for name in names] for d in param_dicts])

    def _to_dicts(self, array):
        names = self.parameters.names()
        return [{name: float(v) for name, v in zip(names, row)} for row in array]
