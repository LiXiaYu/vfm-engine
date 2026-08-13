# License transition record

## Historical MIT baseline

The public thesis baseline is commit
`40d99d888409550d334efb5df95139f7dc1f1cda`, preserved by annotated tag
`v1.0-thesis`, archival branch `legacy/v1-thesis`, and the historical `master`
branch. It was released under the MIT License contained in that commit. The same
MIT text is preserved at `LICENSES/MIT.txt` on `main` for reference.

Those MIT permissions are not revoked, narrowed, or made conditional on payment.
A recipient may continue using, modifying, distributing, sublicensing, and
selling copies of the MIT baseline subject to the MIT notice condition.

## Active main licensing

Beginning with the commit titled
`chore: adopt AGPL and commercial dual licensing on main`, project-controlled
VFM-engine material in `main` is offered under either:

- `AGPL-3.0-or-later`; or
- a separate commercial license executed by an authorized licensor.

The root `LICENSE` is the unmodified GNU Affero General Public License v3.0 text.
The “or later” election is stated in `LICENSING.md` and uses the SPDX identifier
`AGPL-3.0-or-later`.

The commercial path is contractual. Repository notices and draft agreements do
not grant commercial rights outside AGPL.

## Overlap between v1 and main

The transition creates an additional licensing framework; it does not erase
license history. Source carried forward unchanged from `v1.0-thesis` remains
obtainable from that release under MIT. A recipient who relies on the MIT grant
for such code retains it even if the same text also appears in a later `main`
snapshot offered under AGPL.

New project-controlled work first added after the transition is not automatically
MIT-licensed. Its public license is AGPL unless the file or release contains a
different explicit notice. This distinction becomes commercially meaningful as
post-transition development diverges from the MIT baseline.

## Exclusions and release gate

Only a rightsholder or authorized licensor can grant these choices. Third-party
software, dependencies, data, specifications, names, and marks retain their own
terms. The project grants no rights it does not control.

The transition does not certify dependency compatibility or announce a v2
release. The audit in `docs/DEPENDENCY_LICENSE_AUDIT.md` must be completed before
publishing a v2 source/binary package or making dependency-specific commercial
commitments.

## Contribution boundary

All commits require DCO sign-off. Copyrightable external contributions require
an executed CLA before merge so the project has sufficient copyright and patent
rights for both the AGPL and commercial paths. Until a CLA is executed and
recorded, an external contribution must remain unmerged.
