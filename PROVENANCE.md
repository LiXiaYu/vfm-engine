# VFM-engine provenance

## Project origin

VFM-engine was initiated and independently developed by Qi Li during his PhD
work. The project and its public Git history predate Qi Li's subsequent
postdoctoral appointment. It is treated as Qi Li's pre-existing background
software and know-how, subject to any rights in identified third-party
components.

This record addresses the VFM-engine platform. It does not determine ownership
of later institution-specific research results, data, patentable inventions, or
commissioned deliverables; those require their own written analysis and
agreements.

## Public-history evidence

The following facts were recorded locally on 2026-08-13 (Asia/Shanghai):

| Item | Recorded value |
| --- | --- |
| Public repository | `https://github.com/LiXiaYu/vfm-engine` |
| Earliest public commit in the recovered history | `2c13f66a45e72a99bbc6f6b81b229bff2a00c365` (2024-02-19) |
| Public `master` thesis baseline | `40d99d888409550d334efb5df95139f7dc1f1cda` (2026-03-16) |
| Baseline license | MIT License, copyright notice for Qi Li |
| Frozen tag | `v1.0-thesis` (annotated) |
| Archival branch | `legacy/v1-thesis` |

The tag and archival branch identify the same commit. The annotation explicitly
records that the MIT grant for this baseline is not withdrawn.

## Local recovery and unpublished refs

At the beginning of the 2026-08-13 governance work, the requested local path
`C:\Users\liqi\programing\mona\vfm-engine` existed but was empty. The public
repository above was verified and cloned into that path before any refs or
documents were created.

The recovered clone exposed only `origin/master`; no remote tags or other remote
branches were advertised. No pre-existing unpublished VFM-engine branch was
found in the recovered clone. This does **not** prove that private branches or
older working copies do not exist on another computer, private remote, backup,
or storage device.

No `pre-postdoc-*` tag has therefore been invented. If a pre-appointment private
branch is later located, record its original commit without rebasing or merging
it first, for example:

```text
pre-postdoc-baseline-2026-<branch-name>
```

Use a separate annotated tag for each distinct head. Do not push private branch
refs to the public remote merely to create evidence. Preserve a full Git bundle,
an offline copy, the tag annotation, the commit SHA, and a SHA-256 digest.

## Provenance categories for future work

Material changes should record one of these origins in the pull request or
equivalent review record:

1. `original`: created by the contributor without restricted third-party input;
2. `third-party`: derived from an identified source under a documented license;
3. `customer-specific`: created under a services agreement and checked against
   its IP allocation;
4. `institution-related`: created during an appointment or funded project and
   cleared in writing before submission;
5. `generated`: produced with automated tools, reviewed by a human contributor,
   and checked for license, confidentiality, and provenance risk.

Repository history is evidence, not a substitute for employment, funding,
institutional, or commercial agreements. Qi Li should obtain a written
background-IP and improvement-rights acknowledgement from any institution whose
rules may affect future VFM-engine development.
