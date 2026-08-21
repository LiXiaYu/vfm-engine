# MA05P two-element forward/inverse validation

This case validates the periodic viscoelastic constitutive path implemented by
`optim_function_T`:

1. FEBio solves a two-element model for fifty periodic pressure cycles. The
   long forward run is used only to reach a closed periodic response.
2. VFM-engine selects the final complete four-second cycle and rebases its
   displacement to that cycle's first state.
3. `optim_function_T` constructs the periodic `S(t)` response from that final
   cycle; it does not replay the complete fifty-cycle loading history.
4. The three virtual fields from the original MA05P example are used to compare
   internal and external virtual work and identify `E`, `gamma`, and `tau`.
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

The validation runner currently constrains OpenMP to one thread because the
legacy `fun_for_optim_T` implementation accumulates its multi-virtual-field
loss into a shared scalar. This makes the result deterministic until that core
reduction is corrected.

## Reviewed FEBio 4.13 baseline

The reviewed clean run solved 4001 states from `t = 0` through `t = 200 s` and
selected only the 81 states from `t = 196 s` through `t = 200 s`. The final
cycle's endpoint drift in the unit-elastic `S_e_0,xx` response was
`3.47e-13`, compared with a cycle peak of `7.67e-4`. The relative virtual-work
L1 residual was `4.8047241804e-4`. The recovered parameters were
`E = 0.2890757513`, `gamma = 8.2969479304`, and `tau = 1.0072537749`; their
maximum relative error from the model inputs was `3.71%`. These numerical
values are regression references, while the model inputs and the 5% recovery
limit remain the physical acceptance criteria.

Run from the repository root:

```powershell
.\validation\cases\ma05p_2element\run.ps1
```

Generated solve files, result buffers, debug data, and metrics are written below
`out/validation/ma05p_2element/` and are excluded from Git.
