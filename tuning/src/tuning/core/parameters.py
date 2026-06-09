"""Parameters we tune, and the final result of a calibration."""

from dataclasses import dataclass, field


@dataclass
class Parameter:
    """One tunable parameter and the range to search."""

    name: str
    low: float
    high: float
    default: float


@dataclass
class ParameterSet:
    """All tunable parameters for a component."""

    params: list[Parameter]

    def names(self) -> list[str]:
        return [p.name for p in self.params]

    def bounds(self) -> list[tuple[float, float]]:
        return [(p.low, p.high) for p in self.params]


@dataclass
class CalibrationResult:
    """What a sampler returns when it finishes."""

    best_params: dict[str, float]
    notes: str = ""
    extras: dict = field(default_factory=dict)  # e.g. all samples, posterior, emulator
