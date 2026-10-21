# Recovery status: exact cross-scale producer, missing independent audit files

2026-10-21. **Current status: original producer proof/program recovered exactly; complete producer receipt regenerated exactly and replayed; original independent audit files are unavailable.** This is not restoration of the former 18-file reviewed publication packet.

Start with the [original mathematical proof](README.md), the [missing-audit notice](review/REVIEW.md), and the [session recovery overview](../session21_recovery_20261021/README.md). The original publication note is preserved byte-for-byte as [ORIGINAL_PUBLICATION.md](ORIGINAL_PUBLICATION.md), solely as pre-incident provenance. Its statements that review files accompany the packet, and its independent-checker replay commands, do not describe this recovered package.

## What is exact and reproducible now

Six files were recovered from an existing unreferenced Git tree and checked against their Git blob identifiers. The seventh, receipt.json, was regenerated using the exact recovered checker and six exact inherited inputs. Its 4,121,366 bytes have the original SHA-256 d2f996a98d6a12b8cc203223cc62f8bc18ec97b105cf9246f21ef612ebde0dd2. The expected seven-file subtree also exists remotely, confirming that the earlier tree-object upload succeeded. No earlier Session21 commit or main-branch update was reached.

The complete 148,959 handoff segments are retained, not replaced by hashes. Normal generation and optimized verification reproduce the original summaries exactly. Both normal and optimized verification were then run from this staged layout against the fresh source checkout. [RECOVERY_VERIFICATION.json](RECOVERY_VERIFICATION.json) records this recovery replay. It is a producer replay, not a new independent mathematical audit.

From this directory, with the mathematical repository root substituted:

    python3 -B check.py --repo-root /path/to/math --verify
    python3 -B -O check.py --repo-root /path/to/math --verify

The checker authenticates its inherited inputs. Verification has one extra saved-receipt equality check relative to generation, so both verification modes match optimized.json. All five scientific file hashes match the original MANIFEST.json.

## Scope and incident

The source-owned fixed-frontier inverse compiler, stated stable scale laws, and finite optimal-capacity certificates are preserved. Uniform native 153N coverage and uniform control of R-prime remain open. The separate native-remainder and tournament claims in the recovery overview are explicitly reconstructed prose with historical finite evidence unless stated otherwise there.

During a long approval wait, the workspace was replaced and the repository and local stages disappeared. After the pending object-creation call returned, saving its response failed because the old output directory no longer existed. Read-only object reconciliation confirmed success of that upload; it was not retried. Only the recovered scientific files and newly labeled recovery material are included here. The former independent checker, handoff checker, review report, audit receipts and audit manifest have not been recovered. [RECOVERY_PROVENANCE.json](RECOVERY_PROVENANCE.json) separates exact recovery from reconstruction.
