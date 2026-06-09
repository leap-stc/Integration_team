# Which sampler should I use?

Every sampler answers the same question — *which parameter sets next?* — they
just answer it differently.

| Sampler | Waves | Emulator? | How it picks the next wave | Backed by |
|---|---|---|---|---|
| `random` | 1 | no | random within bounds | — |
| `latin_hypercube` | 1 | no | spread evenly across the ranges | — |
| `eki` | few | no | Kalman nudge of the ensemble toward obs | Python EKI |
| `history_matching` | several | yes | keep the "not ruled out yet" region | scikit-learn GP |
| `hmc` | 1–2 | yes | sample the posterior with Hamiltonian Monte Carlo | scikit-learn GP |
| `nelder_mead` | many (small) | no | reflect/expand/contract a simplex downhill | — |

Rule of thumb:
- **Exploring / baseline:** `random` or `latin_hypercube`.
- **Big ensemble you can afford up front:** `history_matching` (the CLM /
  ctsm6_ppe workflow) or `eki`.
- **You want the posterior (uncertainty), not just a best fit:** `hmc`.
- **Few params, just want one best fit, no ensemble:** `nelder_mead`
  (sequential — one run at a time, so slow on HPC unless runs are quick).

`history_matching` and `hmc` build a GP emulator from the runs so far, then
sample it. With `hmc` you can also draw N sets from the posterior and run one
more forward-model wave to check them — that's just another `ask`.
