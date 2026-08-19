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
`elastic_E` optimization. VFM-engine 1.1.3 runs it as a pure VFM validation:
FEBio still parses and initializes the model, materials, and DataStore records,
but `FEModel::Solve()` is not called. A new result buffer is allocated and the
existing Python displacement-and-pressure setter populates it. Run
`20260819-113318` used FEBio 4.13.0 and plugin 1.1.3, explicitly confirmed that
the solve was skipped, and reproduced the reviewed metrics.

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
5. confirms that the configured FEBio solve was skipped; and
6. checks generated metrics when a reviewed expected baseline exists.
