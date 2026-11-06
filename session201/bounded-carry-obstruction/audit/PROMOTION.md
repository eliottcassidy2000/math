# Independent promotion readback

Outcome: **PASS** for the separate bounded-source-error release, including the approved README domain clarification. This supplements, and does not replace or rewrite, REPORT.md, RECEIPTS.json, or any historical input/output.

## Exact final bytes

- PROOF.md: `329151476e7122d36c56f586454fb1e2422e3c2cf5094385b8410af03fb8d4b5`.
- README.md: `47ee66e69cf490c6cbad1d9a226e5964f1c6958bbe9135a068a3811d1f270787`.
- PROMOTION.diff: `b2fc8738e5fdf967df2a8fe5173df392a0d7a93d8e935bdbc69ed6cdcdb15d97`.
- dependencies/CORE_PROOF.md: `49ad1123e80fbf751eac3e72ed0f5ad2ec39e7705f9ee89c77d4826ef5e1066c`, exactly the frozen reviewed central-arithmetic dependency.

The original proof, README, and source manifest in originals/ are byte-identical to the frozen inputs. The promoted proof differs only in audited status, packet-local reproduction commands, the sibling approximation-proof link, and the packet-local product-bound dependency link. The README changes only status and “with R Z-runs” to “with R>=1 Z-runs,” making the already required domain of the displayed division explicit. The full unified diff exactly matches PROMOTION.diff. No mathematics or verifier code was changed. RELEASE.md accurately states the fixed-threshold and both-errors-only limitations.

All historical manifest-listed audit artifacts were compared by both hash and bytes to the original audit. One packaging omission was found and repaired: audit/inputs/session201-bounded-carry-obstruction/MANIFEST.sha256 was absent because the historical audit manifest excluded nested files bearing that basename. The exact unchanged frozen file was restored, hash `5959646a77e11985a59a1611a088a20157a22501c5b0da38d0916e39aa3098f7`. The refreshed audit manifest now excludes only its own root manifest, and includes this nested input. The historical original audit and its manifest were not modified.

## Relocated execution

Executed from the release packet directory:

- `python3 verify.py`
- `python3 -O verify.py`
- `python3 audit/independent_check.py`
- `python3 -O audit/independent_check.py`

All four succeeded. Every output is byte-identical to its corresponding historical exact receipt. New outputs are stored separately as PROMOTION-author-normal.json, PROMOTION-author-optimized.json, PROMOTION-independent-normal.json, and PROMOTION-independent-optimized.json. All INPUTS.json entries resolve and match their frozen hashes. PROMOTION-static-checks.json records the file comparisons, and PROMOTION-receipts.json records final hashes and replay equality.

## Mathematical and packaging scope

For every fixed finite threshold E>=0, only N^{o(1)} actual source indices through N have error at most E. Consequently the twelve fixed forms with BOTH source errors bounded by fixed constants cannot cover an odd integer tail. The R=0 and E=0 boundaries are treated separately in the proof. The result does not exclude a one-sided error bound, a growing threshold, or unrestricted occupancy. **Uniform twelve-template occupancy remains OPEN.**

The sibling link ../hammersley-source-approximation/PROOF.md is intended for the coordinated session201 release layout; this separate packet contains no sibling directory. It is not a runtime dependency of either checker, and the needed bijection/run identities were independently audited here. Final aggregate link resolution and aggregate manifest refresh remain the coordinator's packaging step. The packet-local product-bound dependency resolves and matches exactly.

No source mathematics, verifier code, historical audit report/receipts, repository files, or GitHub state were edited in this promotion review. Only release audit records and the restored frozen audit input were written. The audit manifest inventories all audit files except itself.
