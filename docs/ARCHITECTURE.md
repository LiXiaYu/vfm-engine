# VFM-engine architecture

This document describes the current C++ plugin structure. It records implemented
boundaries rather than a promise to remove FEBio as the runtime.

## Runtime flow

1. `VFMTask` loads the Python configuration and registers FEBio callbacks.
2. `VFMTask_result_capture` allocates the memory-mapped result buffers and
   captures time-step displacement, stress, nodal-force, and constraint data.
3. FEBio remains responsible for the model, mesh, material definitions,
   constitutive evaluation, and finite-element solve.
4. `read_solved_information` coordinates VFM post-processing and the Python
   optimization entry point.
5. `VFMTask_laplace` contains the Laplace-domain and related constitutive helper
   calculations used by the post-processing pipeline.
6. `VFMTask_debug_output` contains optional Tecplot diagnostic export.

## Source responsibilities

- `VFMTask.cpp`: task lifecycle and FEBio solve orchestration.
- `VFMTask_result_capture.cpp`: initialization and per-time-step result capture.
- `VFMTask_callback.cpp`: solved-result VFM orchestration. This remains the main
  target for later, behavior-preserving decomposition.
- `VFMTask_laplace.cpp`: Laplace transforms and related internal-work helpers.
- `VFMTask_debug_output.cpp`: diagnostic output that is not part of the solver
  or optimization contract.
- `VFMTask_pybind.cpp`: Python bindings for task configuration.
- `optim.cpp`: reusable optimization functions and objective implementations.
- `common_FEBio.cpp` and `FEBio_refunction.cpp`: shared FEBio-facing numerical
  and compatibility helpers.

The source list is explicit in `FEBio_VFM_Task/CMakeLists.txt`; adding a new
translation unit requires adding it to that list.

## Jacobian semantics

The callback uses several determinants with different meanings. They must not
share one implicit validity rule:

- **Reference mapping Jacobian**: inverted by `domain_init`; it must be finite
  and positive under the current element-orientation convention.
- **Physical deformation determinant** `det(F)`: passed to FEBio constitutive
  calculations and stress conversions; it must be finite and positive.
- **Physical integration Jacobian**: the actual current-configuration Gauss
  integration measure stored in `trueJArray`; both internal and volume virtual
  work use this measure.
- **Virtual-field determinant**: belongs to the test-field mapping used to form
  virtual strain. It is not inverted and is not a physical deformation, so a
  finite negative or zero value is admissible.

`PhysicalDeformationGradient` exposes the physical integration Jacobian.
`VirtualFieldGradient` intentionally does not expose a mapping Jacobian, which
prevents the virtual mapping from being reused accidentally as an integration
measure. The shared policy is covered by the `jacobian_policy` unit test.

## Refactoring constraints

- Structural changes must start from `develop` on a dedicated branch.
- Moving code between modules must not silently change numerical formulas,
  callback order, FEBio material evaluation, or memory layout.
- Each merge back to `develop` requires a clean plugin build and the available
  validation baseline checks.
- Research studies and private case scripts are not part of the plugin modules.
