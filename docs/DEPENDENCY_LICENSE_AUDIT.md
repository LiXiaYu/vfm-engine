# Dependency and license audit

**Status: initial inventory only — not legal clearance**

The planned `AGPL-3.0-or-later OR Commercial License` model for v2 is conditional
on a complete source, dependency, and distribution audit. This file records what
was visible in the public thesis baseline on 2026-08-13; it does not conclude
that a commercial or AGPL distribution is permitted.

## Initial inventory

| Component or interface | Evidence in tree | Required follow-up |
| --- | --- | --- |
| FEBio / FECore / FEBioMech | Direct header includes and plugin build references | Identify exact FEBio version and governing license/SDK terms; analyze linking, plugin distribution, trademark, and notice requirements |
| NLopt | `vcpkg.json` and `nlopt.hpp` include | Pin exact version; capture license and bundled-algorithm notices; confirm static/dynamic distribution obligations |
| pybind11 | `vcpkg.json` and header includes | Pin exact version; preserve BSD-style notices; inspect generated/bundled artifacts |
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

No v2 dual-license announcement should be treated as effective until this audit
has named versions, verified source links, stored license texts/notices, resolved
all incompatible or unclear material, and received appropriate technical and
legal review.
