"""Run a calibration one wave at a time, saving state between waves.

Use this when a wave is a batch of long HPC runs. Call step() to submit a wave
and exit; let it run on its own schedule; then call step() again after it
finishes. State is saved to a file, so each call can be a fresh process — i.e.
you just re-run the script once per wave.

For fast, in-process runs (e.g. the local runner) use run_calibration instead.
"""

import pickle
from pathlib import Path


class Campaign:
    def __init__(self, sampler, component, runner, state_file="campaign.pkl"):
        self.sampler = sampler
        self.component = component
        self.runner = runner
        self.state_file = Path(state_file)

    def step(self):
        """Do one wave.

        Returns None while waves are still running (you should re-run after the
        wave finishes), or the CalibrationResult once the sampler is done.
        """
        sampler, pending = self._resume()

        # 1. if a wave just finished, read its outputs and learn from them
        if pending is not None:
            metrics = [self.component.compute_metrics(loc) for loc in pending["locations"]]
            sampler.tell(pending["params"], metrics)

        # 2. finished?
        if sampler.is_done():
            self._save(sampler, pending=None)
            return sampler.result()

        # 3. submit the next wave (the runner does not wait for it)
        ensemble = sampler.ask()
        locations = self.runner.run(ensemble, self.component)
        self._save(sampler, pending={"params": ensemble, "locations": locations})
        return None

    def _resume(self):
        """Load saved state, or start fresh on the first wave."""
        if self.state_file.exists():
            with open(self.state_file, "rb") as f:
                state = pickle.load(f)
            return state["sampler"], state["pending"]
        return self.sampler, None

    def _save(self, sampler, pending):
        with open(self.state_file, "wb") as f:
            pickle.dump({"sampler": sampler, "pending": pending}, f)
