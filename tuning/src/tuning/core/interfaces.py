"""The three contracts that make calibration "mix and match".

Read this file first — it is the whole design.

Every calibration method really answers one question:

    Given the runs so far, which parameter sets do we run next?

The thing that answers it is a SAMPLER. The calibration loop (see
orchestration/loop.py) only ever talks to three objects, so you can swap any
one without touching the others:

    sampler.ask()             -> the next "wave" of parameter sets to run
    runner.run(...)           -> runs CESM for each set
    component.compute_metrics -> turns model output into numbers
    sampler.tell(...)         -> sampler sees the results, picks the next wave
    ...repeat until sampler.is_done(), then sampler.result()

- Sampler    = how to choose params   (random, latin_hypercube, eki,
                                        history_matching, hmc, nelder_mead)
- Component  = the model               (cam, clm) — its params and diagnostics
- Runner     = where it runs            (local for tests, derecho for real runs)
"""

from abc import ABC, abstractmethod

import numpy as np

from .parameters import CalibrationResult, ParameterSet

# A single parameter set, e.g. {"medlynslope": 4.1, "fff": 0.8}.
ParamVector = dict[str, float]


class Sampler(ABC):
    """Chooses which parameter sets to run next.

    One wave at a time. Some samplers run a single wave (random, latin
    hypercube); others run several (eki, history matching, hmc). Some build an
    emulator internally — that is their business, not the loop's.
    """

    @abstractmethod
    def ask(self) -> list[ParamVector]:
        """Return the next wave of parameter sets to run."""

    @abstractmethod
    def tell(self, params: list[ParamVector], metrics: list[np.ndarray]) -> None:
        """See the results of the wave from `ask` (used to pick the next one)."""

    @abstractmethod
    def is_done(self) -> bool:
        """True when no more waves are needed."""

    @abstractmethod
    def result(self) -> CalibrationResult:
        """Final answer: best params / posterior / selected sets."""


class Component(ABC):
    """A CESM model component (cam, clm)."""

    @abstractmethod
    def parameters(self) -> ParameterSet:
        """Which parameters can be tuned, and their ranges."""

    @abstractmethod
    def apply(self, case_dir: str, params: ParamVector) -> None:
        """Write one parameter set into a CESM case."""

    @abstractmethod
    def compute_metrics(self, run_output) -> np.ndarray:
        """Turn one model run's output into the numbers we compare to obs."""


class Runner(ABC):
    """Where the model runs."""

    @abstractmethod
    def run(self, ensemble: list[ParamVector], component: Component) -> list:
        """Run every parameter set; return one `run_output` per set."""
