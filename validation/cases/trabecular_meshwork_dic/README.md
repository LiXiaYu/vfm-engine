# Trabecular Meshwork DIC

This case preserves the legacy `xinchen_TrabecularMeshwork_DIC_test` demo in a self-contained validation layout.

## Purpose

The case applies DIC-derived displacements to an FEBio 4 model and runs the VFM elastic modulus optimization callback.

## Inputs

- `input/model/dic_driven_10um.feb`: model selected by the legacy runner.
- `input/model/reference.feb`: related legacy model retained for comparison.
- `input/data/dic_positions.npy`: DIC positions consumed by `scripts/vfm_config.py`.

## Status

This case has an accepted numerical baseline for the deterministic Python
`elastic_E` optimization. It remains a **candidate validation** overall because
the current plugin continues into the legacy `NLpot_0` path after writing the
metrics and exits with `0xC0000409`. The historical run log also ended with
negative Jacobians and failure to converge.

## Run

From the repository root:

```powershell
.\validation\cases\trabecular_meshwork_dic\run.ps1
```

The runner:

1. restores the locked Python environment with `uv`;
2. locates the installed FEBio executable and the current VFM plugin build;
3. creates an isolated directory below `out/validation/`;
4. runs FEBio with this case's VFM configuration;
5. checks generated metrics when a reviewed expected baseline exists.
