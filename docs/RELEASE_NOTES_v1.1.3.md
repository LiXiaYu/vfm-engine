# VFM-engine v1.1.3

Release date: 2026-08-19

VFM-engine 1.1.3 is a source release on the post-thesis development line. It
keeps FEBio model and constitutive initialization while allowing VFM workflows
that provide their own data to avoid an unnecessary finite-element solve.

## Highlights

- separates `run_febio_solve` from `reuse_saved_result_buffer`;
- allocates DataStore-derived result storage explicitly after `FEModel::Init()`
  and before `Run()` can invoke `FEModel::Solve()`;
- requires an existing buffer when saved-buffer reuse is requested;
- preserves the old `isRead_FEMresult_fromsavefile` Python property as a
  deprecated compatibility entry point with its legacy coupled behavior;
- splits callback responsibilities into smaller source modules;
- separates physical constitutive Jacobian requirements from virtual-field
  determinant handling; and
- derives the FEBio plugin version from the CMake project version so the reported
  runtime version cannot silently diverge from the source release.

No input-file format, output-file format, mapped-buffer layout, or
`timedisplacement`/`timestress` indexing convention is changed by this release.

## Validation

The Release plugin was rebuilt after cleaning prior build products. CTest passed
the `jacobian_policy` test. The `trabecular_meshwork_dic` case ran in a newly
created directory (`20260819-113318`) with FEBio 4.13.0 and plugin 1.1.3. Its log
confirmed that `FEModel::Solve()` was skipped and that a new result buffer was
created. The reviewed numerical baseline passed with:

- `elastic_E.x = 1.006667430125674e-05`;
- `fun = 1.5168955672099785e-13`;
- `success = true`;
- `nit = 35`; and
- `nfev = 542`.

## License

Project-controlled source in `v1.1.3` is available under either
`AGPL-3.0-or-later` or a separately executed commercial license. The immutable
`v1.0-thesis` source remains available under MIT; this release does not withdraw
or narrow that earlier grant. Third-party components retain their own licenses.
