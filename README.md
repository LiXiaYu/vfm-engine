# VFM-engine

VFM-engine is a Virtual Fields Method (VFM) research software project initiated
and independently developed by Qi Li during his PhD work, before his subsequent
postdoctoral appointment. The public thesis baseline integrates VFM workflows
with FEBio and includes C++, Python embedding, and nonlinear-optimization code.

## Thesis baseline

The public `master` commit
`40d99d888409550d334efb5df95139f7dc1f1cda` is frozen as:

- annotated tag: `v1.0-thesis`
- archival branch: `legacy/v1-thesis`

VFM-engine `v1.0-thesis` remains permanently available under the MIT License
contained in that release. Nothing in later governance or licensing documents
withdraws, narrows, or revokes rights already granted for that MIT-licensed
version.

## Licensing direction

Licensing is version-specific.

- `v1.0-thesis`: MIT License, permanently.
- Current source inherited unchanged from that baseline: MIT License unless a
  later release expressly says otherwise.
- `v2+`: planned, but not yet released, under
  `AGPL-3.0-or-later OR Commercial License`.

The v2 dual-licensing plan is not effective until the relevant release is
published with complete license notices. A third-party dependency and source
provenance audit must be completed first. See [LICENSING.md](LICENSING.md) and
[the audit checklist](docs/DEPENDENCY_LICENSE_AUDIT.md).

## Project governance

- [Provenance and baseline record](PROVENANCE.md)
- [Governance](GOVERNANCE.md)
- [Contributing](CONTRIBUTING.md)
- [Contribution policy](docs/CONTRIBUTION_POLICY.md)
- [Versioning and release policy](docs/VERSIONING.md)
- [Commercial paths](COMMERCIAL.md)

All new contributions require a DCO `Signed-off-by` line. Contributions to the
dual-licensed core may also require an individual or corporate CLA.

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
