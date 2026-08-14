# Licensing policy

## Version-specific licensing

The license for VFM-engine is determined by the exact version or commit received,
not by a general statement about the project name.

### VFM-engine v1.0-thesis

The source at annotated tag `v1.0-thesis`, commit
`40d99d888409550d334efb5df95139f7dc1f1cda`, is licensed under the MIT License
included in that commit.

Those MIT rights are permanent for recipients who comply with the license.
Later project policies, commercial offerings, or differently licensed releases
do not withdraw or restrict the permissions already granted for v1.

### Current main branch

Beginning with the commit titled
`chore: adopt AGPL and commercial dual licensing on main`, VFM-engine material in
`main` that Qi Li or an authorized project licensor has the right to license is
offered under either:

1. **GNU Affero General Public License v3.0 or later**
   (`AGPL-3.0-or-later`), whose full text is the root `LICENSE`; or
2. a **separate commercial license** signed by the applicable licensor and
   licensee.

The commercial option is not a public click-through license and is not granted
by this repository. `LICENSES/COMMERCIAL.md` and the files in `commercial/` are
explanatory material or drafts. Commercial rights exist only under an executed
agreement covering an exact version, package, product, and scope.

### VFM-engine v1.1.0

The source at annotated tag `v1.1.0` follows the current `main` licensing
framework. Project-controlled material is offered under either
`AGPL-3.0-or-later` or a separately executed commercial license.

`v1.1.0` is a post-thesis source release and is not licensed under MIT merely
because its version number begins with `1`. No binary build, bundled dependency,
or redistribution clearance is implied by the source tag.

### Effect of the MIT baseline

The same source can be offered under more than one license by an authorized
copyright holder. Offering the `main` snapshot under AGPL does not cancel the
earlier MIT offer.

Anyone may continue obtaining and using `v1.0-thesis` under MIT. Moreover, a
person relying on code actually obtained under that MIT grant keeps those MIT
rights, including for portions that happen to be identical to content later
present in `main`. The transition does not convert an existing MIT copy into an
AGPL-only copy and does not make an additional commercial license necessary for
ordinary exercise of the v1 MIT rights.

New project-controlled work first added after the transition is not offered
under MIT merely because it shares a repository or history with v1. Its public
license is `AGPL-3.0-or-later` unless a file carries a different explicit notice.

The historical MIT text is preserved at `LICENSES/MIT.txt`. See
`docs/LICENSE_TRANSITION.md` for the recorded boundary.

### v2 release gate

This license transition establishes the intended licensing framework for active
development; it is not an announcement that v2 has been released or that every
possible binary distribution has been cleared. Before a v2 release, the project
must confirm that it has the necessary rights in every included contribution
and that all third-party dependencies, linked components, generated artifacts,
and distribution methods are compatible with both intended paths.

An actual v2 release must include unambiguous license files and notices, identify
the exact covered material, preserve third-party notices, and explain how a
recipient elects a license.

## Contributions and dual licensing

All contributions require DCO sign-off. Copyrightable contributions from anyone
other than the applicable licensor require an executed individual or corporate
CLA before merge into the dual-licensed codebase. The CLA must permit the
project steward to offer the contribution under open-source and commercial
terms. DCO alone does not grant commercial relicensing authority.

Until an appropriate CLA is executed and recorded, an external contribution
must not be merged into the dual-licensed codebase. A pull request or public
discussion does not itself create the required commercial licensing rights.

No contribution will be accepted when its ownership, confidentiality, patent,
or license status is unclear.

## Third-party material

Third-party components retain their own licenses and notices. A VFM-engine
commercial license can cover only material the licensor has authority to
license; it does not replace licenses from FEBio, NLopt, pybind11, Python, or any
other third party.

The initial inventory in `docs/DEPENDENCY_LICENSE_AUDIT.md` is incomplete and is
not legal clearance.

## Commercial licenses and services

A commercial agreement may grant defined rights to `main`-derived material, a
future release, official SDK, source-access package, or support deliverable
without requiring the licensee to use the software under AGPL. It is not needed
merely to exercise MIT rights in v1. Services, warranties, maintenance,
validation, trademark permissions, and future-version rights may still be
separately valuable and separately contracted.

Questions involving a real distribution or enterprise agreement should be
reviewed by qualified counsel.
