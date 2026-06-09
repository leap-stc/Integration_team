"""tuning — mix-and-match calibration tools for CESM.

Start with core/interfaces.py (the design) and orchestration/loop.py (the flow).
"""

from .core.interfaces import Component, Runner, Sampler
from .core.observations import Observations
from .core.parameters import CalibrationResult, Parameter, ParameterSet
from .orchestration.campaign import Campaign
from .orchestration.loop import run_calibration

# Import the plugin packages so they register themselves by name.
from . import components, runners, samplers  # noqa: E402, F401

__all__ = [
    "Sampler", "Component", "Runner",
    "Parameter", "ParameterSet", "CalibrationResult", "Observations",
    "run_calibration", "Campaign",
]
