# VFM-engine v1.1.4

Release date: 2026-08-21

VFM-engine 1.1.4 is a source release that improves virtual-work integration,
constitutive validation, optimization output control, and deterministic
parallel evaluation.

## Highlights

- uses FEBio-provided Gaussian integration weights in solid-element
  virtual-work calculations;
- adds reviewed integration validation for `hex8`, `tet4`, `penta6`, and
  `pyra5` elements;
- adds an MA05P forward-solve and final-cycle periodic viscoelastic inversion
  baseline;
- honors the configured optimization dump path and preserves captured source
  buffers during cached VFM replay;
- makes per-evaluation stress and virtual-work CSV output opt-in through
  `optim_function_output_iteration_stress` and
  `optim_function_output_iteration_virtual_work`;
- adds `optim_function_parallel` and `optim_function_num_threads` controls; and
- removes races from parallel virtual-field loss accumulation while keeping
  FEBio material-point updates that are not safe to share outside the parallel
  region.

No input-file format, output-file format, mapped-buffer layout, or
`timedisplacement`/`timestress` indexing convention is changed by this release.
The legacy in-process `NLpot_0` optimizer remains disabled.

## Validation

The Release plugin was rebuilt against the installed FEBio 4.13 SDK and loaded
by FEBio as version 1.1.4. CTest passed the `jacobian_policy` test.

The Gaussian-integration validation passed for `hex8`, `tet4`, `penta6`, and
`pyra5` (`20260821-113818`). The MA05P case performed a fresh 50-cycle FEBio
forward solve and an 8-thread VFM inversion over the final complete cycle
(`20260821-113835`). It passed the reviewed baseline with:

- virtual-work relative L1 residual: `4.80472418039629e-4`;
- recovered `E`: `0.28907575130477986` (reference `0.3`);
- recovered `gamma`: `8.296947930352085` (reference `8.0`);
- recovered `tau`: `1.0072537749116495` (reference `1.0`); and
- maximum parameter relative error: `0.03711849129401057`.

A separate large private validation case produced the identical objective value
`1.8425449529422005e-13` with one and eight threads, produced no per-evaluation
CSV files under the default output settings, and produced the expected three
files when diagnostic output was explicitly enabled. Private case data and
outputs are not included in this repository.

## License

Project-controlled source in `v1.1.4` is available under either
`AGPL-3.0-or-later` or a separately executed commercial license. The immutable
`v1.0-thesis` source remains available under MIT; this release does not withdraw
or narrow that earlier grant. Third-party components retain their own licenses.
