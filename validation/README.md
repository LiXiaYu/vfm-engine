# Validation

This directory contains small, reusable validation cases for the numerical and physical behavior of VFM-engine. Paper-specific studies and full experimental datasets belong in separate research projects.

## Layout

- `cases/<case>/input/`: versioned FEBio models and compact input data.
- `cases/<case>/scripts/`: VFM configuration and result checking code.
- `cases/<case>/expected/`: reviewed reference metrics and tolerances.
- `cases/<case>/run.ps1`: the case-local Windows entry point.
- `run-all.ps1`: runs every case that provides a `run.ps1` entry point.

Generated files are written below `out/validation/`, which is excluded from Git.

## Cases

- `trabecular_meshwork_dic`: prescribed-field elastic identification.
- `ma05p_2element`: FEBio forward solve to periodic steady state followed by
  final-cycle `optim_function_T` identification and parameter-recovery checks.
- `gauss_integration_solid_elements`: analytical volume and virtual-work
  integration checks for `hex8`, `tet4`, `penta6`, and `pyra5` elements.

## Environment

Validation dependencies are managed with `uv` using `pyproject.toml` and `uv.lock`. The current plugin embeds Python 3.13, so the validation environment is constrained to Python 3.13 as well.

Run one case from the repository root:

```powershell
.\validation\cases\trabecular_meshwork_dic\run.ps1
```

Run all cases:

```powershell
.\validation\run-all.ps1
```

A successful FEBio process is only an execution check. A case becomes a numerical validation only after its reviewed analytical expectations or reference metrics have been committed below `expected/`.
