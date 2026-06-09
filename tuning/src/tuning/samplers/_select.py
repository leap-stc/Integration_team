"""Small helper shared by the simple one-wave samplers."""

import numpy as np

from ..core.parameters import CalibrationResult


def best(archive, observations):
    """Pick the lowest-loss parameter set from a list of (params, metrics)."""
    if not archive or observations is None:
        return CalibrationResult(best_params={}, notes=f"{len(archive)} samples",
                                 extras={"archive": archive})
    losses = [observations.loss(metrics) for _, metrics in archive]
    i = int(np.argmin(losses))
    return CalibrationResult(best_params=archive[i][0],
                             notes=f"best loss={losses[i]:.3f}",
                             extras={"archive": archive})
