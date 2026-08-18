# Expected metrics

The reference metrics in `metrics.json` were accepted on 2026-08-18 from the
deterministic validation run `20260818-111025`.

The reference environment was:

- FEBio runtime 4.13.0;
- FEBio SDK 4.13.0;
- VFM-engine plugin 1.1.0;
- Python 3.13.14 managed by `uv`;
- validation seed `20260818`.

This baseline covers the deterministic Python `elastic_E` optimization result.
It was reproduced by VFM-engine 1.1.1 in run `20260818-112424` after the legacy
`NLpot_0` path was disabled, and again by VFM-engine 1.1.2 in run
`20260818-150528` after task lifecycle and callback error handling were hardened.
The previous `0xC0000409` plugin crash no longer occurs. Version 1.1.2 preserves
the generated metrics while propagating the failed FEBio solve as a task failure.
The complete case remains a candidate because the FEBio model still reports
failed convergence, negative Jacobians, and error termination.

The tolerances allow small floating-point and optimizer differences while keeping
the optimized modulus tightly constrained. Because the objective value is close
to zero, its comparison is governed primarily by the absolute tolerance.
