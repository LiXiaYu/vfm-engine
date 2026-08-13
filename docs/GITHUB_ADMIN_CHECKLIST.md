# GitHub administration checklist

This checklist is intentionally not automated. Complete it only after reviewing
the local refs and governance commit.

## 1. Push the immutable thesis references

Push only the known public baseline refs:

```text
git push origin refs/tags/v1.0-thesis
git push origin refs/heads/legacy/v1-thesis:refs/heads/legacy/v1-thesis
```

Do not use `--force`. Do not push any private or `pre-postdoc-*` ref to the public
remote.

## 2. Push and review main

```text
git push -u origin main
```

Review the branch on GitHub before changing the default. Confirm that
`v1.0-thesis`, `legacy/v1-thesis`, and the old `master` all still resolve to
`40d99d888409550d334efb5df95139f7dc1f1cda`, while `main` contains the later
governance commit.

## 3. Publish the thesis release

In GitHub Releases:

1. enable release immutability in repository settings **before** publishing if
   the feature is available for the repository;
2. create a draft release using existing tag `v1.0-thesis`;
3. title it `VFM-engine v1.0-thesis`;
4. use `docs/RELEASE_NOTES_v1.0-thesis.md` as the starting notes;
5. attach any verified thesis citation, reproducibility inputs, or checksums;
6. verify every asset before publishing; and
7. publish once complete, then confirm the release is shown as immutable.

Release immutability generally applies only to releases published after it is
enabled. Do not plan to move or reuse the tag.

## 4. Change the default branch safely

After `main` is pushed and reviewed, use repository **Settings → Branches** (or
the current default-branch settings page) to change the default from `master` to
`main`. Do not delete `master` as part of this migration. Update local clones
only after the remote default changes.

## 5. Rules and protection

Create branch rules/rulesets for:

- `main`: block force pushes and deletion; require pull requests, review, status
  checks, conversation resolution, and signed commits where practical;
- `master`: block force pushes and deletion; no direct development;
- `legacy/v1-thesis`: block pushes, force pushes, and deletion; and
- release tags such as `v*`: restrict creation/deletion and rely on immutable
  releases for published release tags where available.

Configure a DCO check for pull requests. A cryptographic commit signature and a
DCO `Signed-off-by` line solve different problems; requiring one does not replace
the other.

## 6. Repository metadata

- update the repository description and topics;
- enable issue and pull-request templates only when the project is ready to
  accept outside work;
- require the CLA process before accepting contributions intended for the
  dual-licensed core;
- review collaborator/admin access and enable strong account authentication;
- keep commercial contract drafts visibly marked as drafts; and
- do not label v2 as AGPL/commercially dual-licensed until the release audit is
  complete.

## 7. Independent mirror and evidence archive

- create a separately controlled private mirror with all public refs;
- use a mirror push only to the confirmed empty/dedicated mirror remote;
- create an offline Git bundle and a SHA-256 digest;
- store release notes, tag SHA, bundle digest, and GitHub release evidence in at
  least two independently controlled locations; and
- consider an archival research deposit once bibliographic metadata and release
  contents are final.

Example commands, after setting and verifying a dedicated mirror URL:

```text
git remote add mirror <verified-dedicated-mirror-url>
git remote -v
git push --mirror mirror
git bundle create vfm-engine-public.bundle --all
```

`git push --mirror` makes the destination refs match the source and can delete
destination-only refs. Use it only for a dedicated mirror after inspection; do
not run it against a repository containing independent work.
