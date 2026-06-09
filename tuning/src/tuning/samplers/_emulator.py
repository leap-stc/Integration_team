"""The default emulator for history matching: one Gaussian process per output.

An emulator is a fast stand-in for the model: given parameters, it predicts each
output *and* how unsure it is. History matching uses that uncertainty to decide
which parameter regions can still be ruled in or out.

Any object with the same `fit`/`predict` shape works here, so ESEm (or another
engine) can replace this later.
"""

import numpy as np


class GPEmulator:
    def fit(self, X, Y):
        """X: [n_runs, n_params].  Y: [n_runs, n_outputs] (or 1-D for one output)."""
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel

        Y = np.asarray(Y)
        if Y.ndim == 1:
            Y = Y[:, None]
        self._models = []
        for j in range(Y.shape[1]):
            gp = GaussianProcessRegressor(kernel=RBF() + WhiteKernel(), normalize_y=True)
            gp.fit(X, Y[:, j])
            self._models.append(gp)
        return self

    def predict(self, X):
        """Return (mean, std), each shaped [n_points, n_outputs]."""
        means, stds = [], []
        for gp in self._models:
            mean, std = gp.predict(X, return_std=True)
            means.append(mean)
            stds.append(std)
        return np.column_stack(means), np.column_stack(stds)
