# VFM-engine

VFM-engine is a Virtual Fields Method (VFM) research software project initiated
and independently developed by Qi Li during his PhD research. The public thesis
baseline integrates VFM workflows with FEBio and includes C++, Python embedding,
and nonlinear-optimization code.

## Thesis baseline

The public `master` commit
`40d99d888409550d334efb5df95139f7dc1f1cda` is frozen as:

- annotated tag: `v1.0-thesis`
- archival branch: `legacy/v1-thesis`

VFM-engine `v1.0-thesis` remains permanently available under the MIT License
contained in that release. Nothing in later governance or licensing documents
withdraws, narrows, or revokes rights already granted for that MIT-licensed
version.

## Current source release

`v1.1.0` is the first post-thesis source release. It promotes the reviewed
output-displacement development work, reorganizes the C++ sources into
`include/` and `src/`, and preserves FEBio as the finite-element and
constitutive runtime.

VFM-engine `v1.1.0` is available under `AGPL-3.0-or-later` or a separately
executed commercial license for project-controlled material. It is not an MIT
continuation of `v1.0-thesis`. No binary artifact is claimed to have been built
or validated as part of this source release.

The current development version is `1.1.2`. It propagates FEBio solve and
post-processing failures, keeps timestep capture state per task instance, and
uses Python `BeforeOptim` as the required optimization entry point. The legacy
in-process `NLpot_0` optimizer remains disabled.

## Licensing

Licensing is version-specific.

- `v1.0-thesis`: MIT License, permanently.
- `v1.1.0`: `AGPL-3.0-or-later` or a separately executed commercial license.
- current post-thesis development line, including `main` and `develop`:
  available under either
  `AGPL-3.0-or-later` or a separately executed commercial license, for material
  that Qi Li or an authorized project licensor has the right to license.

The root [LICENSE](LICENSE) contains the GNU Affero General Public License v3.0
text for the open-source option. Commercial terms arise only from a separate
signed agreement; [`LICENSES/COMMERCIAL.md`](LICENSES/COMMERCIAL.md) is an
explanation, not a license grant.

The license transition on `main` does not withdraw or narrow the MIT rights in
the thesis baseline. In particular, code obtainable from `v1.0-thesis` remains
usable under MIT even when identical code also appears in `main`. See
[LICENSING.md](LICENSING.md), [the transition record](docs/LICENSE_TRANSITION.md),
and [the dependency audit](docs/DEPENDENCY_LICENSE_AUDIT.md).

No v2 release has been published. A third-party dependency, linkage, and source
provenance audit remains a release gate.

## Project governance

- [Provenance and baseline record](PROVENANCE.md)
- [Governance](GOVERNANCE.md)
- [Contributing](CONTRIBUTING.md)
- [Contribution policy](docs/CONTRIBUTION_POLICY.md)
- [Versioning and release policy](docs/VERSIONING.md)
- [Commercial paths](COMMERCIAL.md)

All new contributions require a DCO `Signed-off-by` line. Copyrightable external
contributions to the dual-licensed codebase also require an executed individual
or corporate CLA before merge.

## Development workflow

- `main` is the stable default branch.
- `develop` is the integration branch for reviewed feature work.
- Feature branches should normally merge into `develop` before integrated
  changes are promoted to `main`.

## Source layout

- `FEBio_VFM_Task/include/` contains the public and internal C++ headers.
- `FEBio_VFM_Task/src/` contains the C++ implementation, optimization,
  callback, FEBio-integration, and Python-binding sources.
- `FEBio_VFM_Task/CMakeLists.txt` defines the shared-library target and its
  dependencies.
- `CMakePresets.json`, `vcpkg.json`, and `vcpkg-configuration.json` describe the
  current CMake and vcpkg configuration.

## Current technical dependencies

The current tree references FEBio/FECore, NLopt, pybind11, Python, and OpenMP.
Their inclusion here is not a statement that every future distribution model
has been cleared. Users must comply with each applicable third-party license.

## Commercial collaboration

Commercial use of the already-published MIT v1 does not require an additional
permission from the project. Organizations may nevertheless contract for
official integration, validation, deployment, maintenance, support, or a
separately licensed future version. Commercial licensing does not transfer the
project, repository, roadmap, or underlying copyright.

Draft business documents are in [`commercial/`](commercial/). They are planning
materials, not offers and not ready for signature. Before signing an enterprise
agreement, the documents should be reviewed and adapted by a qualified Chinese
intellectual-property lawyer.

## Citation

See [CITATION.cff](CITATION.cff). Publication and thesis bibliographic details
have intentionally not been guessed and can be added when verified.
