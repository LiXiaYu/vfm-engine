# Dependency and license audit

**Status: initial inventory only - not legal clearance**

The `main` licensing framework is `AGPL-3.0-or-later OR Commercial License` for
project-controlled material. Release and distribution of a v2 package remain
conditional on a complete source, dependency, linkage, and distribution audit.
This file records the initial inventory visible across the public thesis
baseline and the `v1.1.0` source tree as reviewed on 2026-08-14. It does not
conclude that every commercial or AGPL binary distribution is permitted.

## Initial inventory

| Component or interface | Evidence in tree | Required follow-up |
| --- | --- | --- |
| FEBio / FECore / FEBioMech | Direct header includes and plugin build references | Identify exact FEBio version and governing license/SDK terms; analyze linking, plugin distribution, trademark, and notice requirements |
| NLopt | `vcpkg.json` and `nlopt.hpp` include | Pin exact version; capture license and bundled-algorithm notices; confirm static/dynamic distribution obligations |
| Eigen3 | `vcpkg.json` and Eigen header includes | Pin exact version; verify MPL-2.0 and any bundled third-party notices; confirm that only supported public interfaces are used |
| pybind11 | `vcpkg.json` and header includes | Pin exact version; preserve BSD-style notices; inspect generated/bundled artifacts |
| mio | `vcpkg.json` and `mio/mmap.hpp` include | Pin exact version; preserve the applicable MIT-style copyright and license notice |
| Python | custom vcpkg overlay marked `Python-2.0` and embedding code | Pin exact runtime version; preserve PSF and bundled third-party notices; review redistribution terms |
| OpenMP | `omp.h` and build usage | Identify compiler/runtime implementation and redistribution terms |
| Platform/runtime libraries | C++ standard library, dynamic loader, compiler toolchain | Record each target platform and redistributable runtime license |

## Source audit

For every tracked source and generated file:

- identify author and creation source;
- review Git history and any pre-Git origin;
- locate copied snippets, examples, specifications, generated code, and model
  outputs;
- verify employer, university, funder, and customer rights;
- record copyright and license notices; and
- separate project-owned code from third-party or customer material.

## Distribution audit

Perform the review separately for:

- source release;
- dynamically linked desktop/plugin binaries;
- statically linked binaries;
- Python wheels or embedded Python distributions;
- SDK/source-access packages;
- containers and hosted/SaaS deployments; and
- customer-specific adapters.

The legal result may differ by configuration. Do not assume that permission to
use a dependency internally includes permission to redistribute it or offer it
under a commercial license.

## Release gate

No v2 release should be published until this audit has named versions, verified
source links, stored license texts/notices, resolved all incompatible or unclear
material, and received appropriate technical and legal review. The repository's
license choice for project-owned source does not replace this release gate.

The `v1.1.0` tag is a source release only. It does not represent a cleared or
tested binary distribution.
