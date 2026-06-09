"""Runner that builds and submits a CESM ensemble on Derecho (PBS).

Mirrors NCAR/ctsm6_ppe/jobscripts/run_ens.sh: clone a prebuilt base case per
member (reusing the compiled executable), write that member's parameters into
the case, and submit it.

It does NOT wait for jobs to finish. Each wave's cases go in their own
`output_root/waveNN/` directory, and run() returns where each member's output
will appear. You collect a wave's results on the next `Campaign.step()`, after
it has finished on Derecho — so the calibration advances one wave per restart.

DRAFT — assumes a base case is already built (create_newcase + case.build) and
CIME's scripts dir is on PATH. Not runnable here; see the example README.
"""

import glob
import os
import subprocess

from ..core.interfaces import Runner
from ..core.registry import register_runner


@register_runner("derecho")
class DerechoRunner(Runner):
    def __init__(self, base_case=None, output_root=".", project=None):
        self.base_case = base_case        # a prebuilt CESM case to clone
        self.output_root = output_root
        self.project = project

    def run(self, ensemble, component):
        """Submit one wave (no waiting); return each member's history dir."""
        wave_dir = os.path.join(self.output_root, self._next_wave_name())
        locations = []
        for i, params in enumerate(ensemble):
            case_dir = os.path.join(wave_dir, f"member_{i:04d}")
            self._create_case(case_dir)
            component.apply(case_dir, params)   # write this member's paramfile
            self._submit(case_dir)
            locations.append(self._history_dir(case_dir))
        return locations

    def _next_wave_name(self):
        existing = glob.glob(os.path.join(self.output_root, "wave*"))
        return f"wave{len(existing) + 1:02d}"

    def _create_case(self, case_dir):
        # clone the prebuilt base case, keeping its compiled executable.
        # NOTE: `create_clone` is a CIME tool; its scripts dir must be on PATH.
        subprocess.run(["create_clone", "--case", case_dir,
                        "--clone", self.base_case, "--keepexe"], check=True)
        subprocess.run(["./case.setup"], cwd=case_dir, check=True)
        if self.project:
            subprocess.run(["./xmlchange", f"PROJECT={self.project}"], cwd=case_dir, check=True)

    def _submit(self, case_dir):
        subprocess.run(["./case.submit"], cwd=case_dir, check=True)

    def _history_dir(self, case_dir):
        # TODO: point at your run / short-term-archive directory for this case.
        return os.path.join(case_dir, "run")
