# Versioning and release policy

## Immutable historical line

`v1.0-thesis` identifies commit
`40d99d888409550d334efb5df95139f7dc1f1cda`. It is the public thesis baseline
and remains MIT-licensed. The tag must never be moved or reused.

`legacy/v1-thesis` is an archival convenience branch pointing to the same
baseline. It is not the active development line and should be protected against
direct commits.

## Current source release

`v1.1.0` is the first post-thesis source release. It is licensed under
`AGPL-3.0-or-later OR Commercial License` for project-controlled material. The
minor version number does not extend the MIT terms of `v1.0-thesis` to new
post-transition work.

The `v1.1.0` tag identifies reviewed source. It does not assert that a binary
package or every possible redistribution configuration has passed the full
dependency and linkage audit.

## Branch roles

`main` is the stable project branch and the GitHub default branch. `develop` is
the integration branch for reviewed feature work before it is promoted to
`main`. Feature branches should use descriptive names such as
`feature/output-displacement` and normally merge into `develop` first.

`main` was created from the thesis baseline without rewriting `master`. The
historical `master` and `legacy/v1-thesis` refs continue to identify the MIT
baseline and must not be force-pushed.

## License transitions

Every release must state the license applying to that exact source. The active
`main` line uses `AGPL-3.0-or-later OR Commercial License` for project-controlled
material beginning with the recorded license-transition commit. A commercial
license exists only when separately signed. This change cannot revoke the prior
MIT release.

A v2 dual-license release must not be cut until:

1. all included code has documented provenance;
2. contributor agreements support both intended license paths;
3. dependency and linkage compatibility is reviewed;
4. third-party notices and source-offer obligations are prepared;
5. release artifacts and source correspond to one immutable commit; and
6. commercial terms have been reviewed by qualified counsel.

## Release checklist

- use a release-specific annotated tag;
- record the full commit SHA and release date;
- build and test from a clean checkout;
- inventory dependencies and preserve notices;
- prepare release notes that identify licensing and compatibility changes;
- create the GitHub release as a draft and attach all intended assets;
- publish only after the draft is complete if release immutability is enabled;
- archive a Git bundle and record its SHA-256 digest; and
- mirror all refs and release metadata to a separately controlled location.

## Suggested naming

- thesis baseline: `v1.0-thesis`;
- first post-thesis source release: `v1.1.0`, with explicit
  `AGPL-3.0-or-later OR Commercial License` notice;
- later MIT-compatible maintenance of the thesis line, if any: use an explicitly
  identified maintenance tag and MIT notice rather than assuming all `v1.x`
  releases are MIT;
- next major public release: `v2.0.0` only after the licensing audit; and
- archived development head, if useful: a neutral tag such as
  `archive/<branch-name>-<YYYYMMDD>` for each independently verified head.

Private evidence refs should not be pushed to a public remote.
