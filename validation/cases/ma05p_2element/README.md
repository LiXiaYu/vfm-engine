# MA05P two-element forward/inverse validation

This case checks the complete numerical path rather than reusing prescribed or
previously cached fields:

1. FEBio solves a two-element model for ten periodic pressure cycles.
2. VFM-engine reads the displacement, stress, and load history captured from
   that same solve.
3. VFM retains the original model configuration as the kinematic reference
   and uses all ten cycles. This is required for a finite-strain constitutive
   comparison: rebasing displacement at a late cycle would discard the
   prestress and material history accumulated before that cycle.
4. VFM compares internal and external virtual work. Parameter identification
   uses the first virtual field and FEBio 4.13's discrete one-term Prony
   recurrence, initialized with zero history at model time zero.
5. The generated metrics are checked against the known material values and,
   independently, against the committed numerical baseline.

## Model

- `Part5`: target viscoelastic neo-Hookean element with known
  `E = 0.3`, `gamma = g1 = 8`, and `tau = t1 = 1`.
- `Part5a`: known isotropic elastic support with `E = 0.03`.
- `Pressure1`: triangular periodic history with a four-second period.
- Both solid elements are `hex8`; the loaded face is `quad4`.

The validation passes only when the FEBio/VFM virtual-work residual is below
`1e-3` in relative L1 norm and each recovered material parameter is within 5%
of its known value. Baseline comparison is an additional regression check; it
does not replace these physical acceptance criteria.

## Reviewed FEBio 4.13 baseline

The reviewed clean run solved 801 states from `t = 0` through `t = 40 s` and
reported a relative virtual-work L1 residual of `4.3851681375e-5`. The recovered
parameters were `E = 0.3098005185`, `gamma = 8.1421802343`, and
`tau = 0.9946489976`; their maximum relative error from the model inputs was
`3.27%`. These numerical values are regression references, while the model
inputs and the 5% recovery limit remain the physical acceptance criteria.

Run from the repository root:

```powershell
.\validation\cases\ma05p_2element\run.ps1
```

Generated solve files, result buffers, debug data, and metrics are written below
`out/validation/ma05p_2element/` and are excluded from Git.
