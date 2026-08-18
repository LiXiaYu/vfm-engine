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
The Jacobian-semantics refactor reproduced the same metrics in run
`20260818-190423`; its explicit baseline comparison passed.
The complete case remains a candidate because the FEBio nonlinear solve reaches
its maximum retry limit and terminates with failed convergence. One trial state
also reports negative Jacobians, but that report is not the final termination
reason and must not be conflated with VFM-engine's Jacobian rules. In particular,
a virtual test-field determinant may be negative because that mapping is not a
physical deformation and is not used as the integration measure.

The tolerances allow small floating-point and optimizer differences while keeping
the optimized modulus tightly constrained. Because the objective value is close
to zero, its comparison is governed primarily by the absolute tolerance.
