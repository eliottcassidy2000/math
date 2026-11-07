# Release verification

Verdict: PASS, independently audited.

The corrected RESULT.md was compared byte-for-byte with REVIEWED_RESULT.md after exactly the approved replacements: both ln(N/(2C)) occurrences become ln(N/(2C^2)); the derivative ln B-1 becomes ln B-2 with its corresponding ln8 comparison; and ln(N/(2C0)) becomes ln(N/(2C0^2)). No other research-text changes occurred before packaging. The frozen original matches its pre-correction SHA-256 receipt.

Packaging additionally promotes only the status sentence to independently audited, and makes audit.py locate check.json relative to its own __file__. The original REVIEW.md and INPUT_SHA256.txt remain historical audit receipts; their original workspace paths describe the audited inputs rather than runtime requirements. CORRECTION.diff is the exact mathematical correction patch. No finite experiment or research conclusion has been expanded.

The complete THM-4887 and THM-4890 dependency files are copied from the pinned Git commit, with original repository filenames, Git blob IDs and SHA-256 hashes in dependencies/PROVENANCE.json. The audit inherits their proved statements rather than claiming to independently recertify every theorem in their dependency graph.

Portable reproduction, from any working directory:

    python3 /path/to/release/check.py > check-rerun.json
    python3 /path/to/release/audit/audit.py > audit-rerun.json

Both commands were run with /tmp as working directory. Their parsed JSON agrees exactly with the recorded check.json and audit/audit.json respectively. The computation remains a finite control; the all-n and all-X proof is analytic.

SHA256SUMS covers every delivered file other than the manifest itself. No publishing was performed.
