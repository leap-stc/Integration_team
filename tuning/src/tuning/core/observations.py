"""Observations we calibrate toward, and how far a run is from them."""

from dataclasses import dataclass

import numpy as np


@dataclass
class Observations:
    """Target values and how uncertain each one is."""

    targets: np.ndarray       # values we want the model to match
    uncertainty: np.ndarray   # one std-dev per target

    def loss(self, metrics: np.ndarray) -> float:
        """Lower is better: average squared error, scaled by uncertainty."""
        error = (metrics - self.targets) / self.uncertainty
        return float(np.mean(error ** 2))
