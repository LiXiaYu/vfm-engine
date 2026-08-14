# VFM-engine v1.1.0

Release date: 2026-08-14

`v1.1.0` is the first post-thesis VFM-engine source release. It promotes the
reviewed output-displacement development work from `develop` to `main` while
preserving the immutable MIT-licensed `v1.0-thesis` baseline.

## Highlights

- reorganizes C++ headers and implementations under
  `FEBio_VFM_Task/include/` and `FEBio_VFM_Task/src/`;
- expands output-displacement processing, callback handling, optimization, and
  FEBio material-point workflows;
- uses contiguous memory-mapped working arrays for time, displacement, stress,
  nodal-force, pressure, and activation data;
- adds Eigen3 and mio to the recorded build dependencies;
- records the provenance of the integrated GitLab development line; and
- fixes the activation-buffer allocation for models with more than one
  nonlinear constraint.

## License

Project-controlled source in `v1.1.0` is available under either
`AGPL-3.0-or-later` or a separately executed commercial license. The historical
`v1.0-thesis` source remains permanently available under MIT. The `v1.1.0` tag
does not withdraw, narrow, or extend the earlier MIT grant.

Third-party components retain their own licenses. See `LICENSING.md` and
`docs/DEPENDENCY_LICENSE_AUDIT.md`.

## Build status

This is a source release. The release-preparation workstation did not contain
the required CMake, vcpkg, compiler, FEBio SDK, or configured Python environment,
so no binary artifact is published and no successful build is claimed. Build
and numerical validation remain required on a configured FEBio development
system.
