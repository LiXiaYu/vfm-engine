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
It does not mark the complete FEBio task as passing: after writing these metrics,
the current plugin continues into the legacy `NLpot_0` path and exits with
`0xC0000409`. The case remains a candidate until that execution failure is fixed.

The tolerances allow small floating-point and optimizer differences while keeping
the optimized modulus tightly constrained. Because the objective value is close
to zero, its comparison is governed primarily by the absolute tolerance.
