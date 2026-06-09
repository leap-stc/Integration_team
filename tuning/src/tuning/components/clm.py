"""CLM/CTSM land component — write parameters into a case, read diagnostics back.

Mirrors the NCAR/ctsm6_ppe workflow:
  - apply():           copy a base parameter file, overwrite the sampled
                       parameters, and point user_nl_clm at the new file
                       (like gen_paramfiles/ppe_tools).
  - compute_metrics(): open the run's history files and reduce them to the
                       numbers we compare to observations (like the
                       ctsm6_calibration diagnostics notebooks).

DRAFT — needs your machine's paths and a real CESM run; see
examples/clm_history_matching/README.md. Not runnable without a live case.
"""

import os
import shutil

from ..core.interfaces import Component
from ..core.parameters import Parameter, ParameterSet
from ..core.registry import register_component


@register_component("clm")
class CLM(Component):
    def __init__(self, paramranges_csv=None, base_paramfile=None,
                 history_var="GPP", paramfile_subdir="paramfiles"):
        self.paramranges_csv = paramranges_csv   # ctsm6_ppe ranges CSV
        self.base_paramfile = base_paramfile      # netCDF to copy and edit
        self.history_var = history_var            # which output to diagnose
        self.paramfile_subdir = paramfile_subdir

    def parameters(self) -> ParameterSet:
        """Read tunable parameters and ranges from the ctsm6_ppe ranges CSV.

        Only scalar (global) parameters are handled here. Per-PFT parameters
        (where min/max are 'pft' or '30percent') need a vector parameter type —
        a TODO for a later version.
        """
        if self.paramranges_csv is None:
            # example defaults so parameters() works without a CSV
            return ParameterSet([
                Parameter("grperc", 0.05, 0.3, 0.11),
                Parameter("wc2wjb0", 0.6527, 1.5, 0.8),
            ])

        import pandas as pd
        table = pd.read_csv(self.paramranges_csv)
        params = []
        for _, row in table.iterrows():
            if row.get("include", 1) != 1:
                continue
            try:
                low, high = float(row["min"]), float(row["max"])
            except (ValueError, TypeError):
                continue  # skip pft / percent rows for now (see docstring)
            params.append(Parameter(row["param"], low, high, (low + high) / 2))
        return ParameterSet(params)

    def apply(self, case_dir, params):
        """Write one parameter set into the case as a fresh CLM paramfile."""
        from netCDF4 import Dataset  # lazy: only needed for real runs

        # 1. copy the base paramfile into the case
        out_dir = os.path.join(case_dir, self.paramfile_subdir)
        os.makedirs(out_dir, exist_ok=True)
        paramfile = os.path.join(out_dir, "params.nc")
        shutil.copy(self.base_paramfile, paramfile)

        # 2. overwrite the sampled parameters
        with Dataset(paramfile, "a") as nc:
            for name, value in params.items():
                nc.variables[name][...] = value

        # 3. point user_nl_clm at the new paramfile
        with open(os.path.join(case_dir, "user_nl_clm"), "a") as f:
            f.write(f"\nparamfile = '{paramfile}'\n")

    def compute_metrics(self, run_output):
        """Reduce one member's history files to the numbers we compare to obs.

        `run_output` is the directory of history (.nc) files for that member.
        Here we take a simple global+time mean of one variable; real cases build
        the obs-matching diagnostics (see ctsm6_calibration notebooks).
        """
        import numpy as np
        import xarray as xr

        ds = xr.open_mfdataset(os.path.join(run_output, "*.clm2.h0.*.nc"))
        field = ds[self.history_var].mean()      # TODO: real spatial/temporal reduction
        return np.atleast_1d(field.values)
