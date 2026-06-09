# Which sampler should I use?

Every sampler answers the same question — *which parameter sets next?* — they
just answer it differently.

| Sampler | Waves | Emulator? | How it picks the next wave | Backed by |
|---|---|---|---|---|
| `random` | 1 | no | random within bounds | — |
| `latin_hypercube` | 1 | no | spread evenly across the ranges | — |
| `eki` | few | no | Kalman nudge of the ensemble toward obs | Python EKI |
| `history_matching` | several | yes | keep the "not ruled out yet" region | scikit-learn GP (ESEm optional) |
| `mcmc` | several | yes | sample the posterior (incl. HMC) | ESEm |
| `bayesopt` | many | yes | the point most likely to improve the fit | ESEm (GP) |

Rule of thumb:
- **Exploring / baseline:** `random` or `latin_hypercube`.
- **Big ensemble you can afford up front:** `history_matching` (the CLM /
  ctsm6_ppe workflow) or `eki`.
- **Each run is expensive, few params:** `bayesopt`.
- **You want uncertainty, not just a best fit:** `mcmc`.

`history_matching`, `mcmc`, and `bayesopt` build an emulator from the runs so
far, then sample it. For `mcmc`/`bayesopt` you can draw N sets from the
posterior and run one more wave — that's just another `ask`.
