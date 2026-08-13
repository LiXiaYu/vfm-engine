# Versioning and release policy

## Immutable historical line

`v1.0-thesis` identifies commit
`40d99d888409550d334efb5df95139f7dc1f1cda`. It is the public thesis baseline
and remains MIT-licensed. The tag must never be moved or reused.

`legacy/v1-thesis` is an archival convenience branch pointing to the same
baseline. It is not the active development line and should be protected against
direct commits.

## Active line

`main` is the active governance and development branch and the GitHub default
branch. It was created from the thesis baseline without rewriting `master`.
The historical `master` and `legacy/v1-thesis` refs continue to identify the MIT
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
- later MIT-compatible v1 maintenance, if any: `v1.x.y` with explicit MIT notice;
- next major public release: `v2.0.0` only after the licensing audit; and
- private pre-appointment evidence: one annotated
  `pre-postdoc-baseline-YYYY-<branch>` tag per independently verified head.

Private evidence refs should not be pushed to a public remote.
