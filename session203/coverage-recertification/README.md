# Corrected independent verification of the odd-count prefix

**PASS, independently reviewed in the stated finite scope.** Every odd integer from 1 through 4,294,967,295 except 7 and 21 is certified by the inherited physical one-arc-join families. This recertifies the previously known finite range; it proves no eventual threshold and does not settle the full spectrum conjecture.

An endpoint overrun was found in the inherited independent S78 coverage checker. The intended endpoint contributes no count, but undefined behavior cannot be trusted. The corrected checker was therefore replayed over the full range using F1+F2 alone, with F0 entirely disabled. It again finds exactly 7 and 21 missing. A separate reviewer checked the repair, source bounds, full-run totals and receipts, plus independently reran a smaller prefix and literal graph controls.

- [Concise audit and exact scope](session203-independent-audit/COVERAGE_AUDIT.md)
- [Completed independent review](session203-s78-correction-review/REVIEW.md)
- [Endpoint diagnosis and repair](session203-independent-audit/ENDPOINT_REPAIR.md)
- [Full corrected replay output](session203-independent-audit/CORRECTED_NOF0_32.out)
- [Machine-readable receipt](session203-independent-audit/COVERAGE_REPLAY_RESULT.json)
- [Reviewable patch to the inherited checker](session203-independent-audit/S78_ENDPOINT_FIX.patch)

The replay tested 47,964,383,242 F1 pairs and 134,010,742,236 F2 evaluations in 465.1 seconds using four threads. The 2,732-triple rooted bank was regenerated exactly. The full pass is a corrected replay of the already independent S78 implementation, not a third newly authored full-range sieve. Its separate review does not repeat that full pass.

## Reproduction

Run from this directory, whose two subdirectories retain the audited source layout:

    (cd session203-independent-audit && sha256sum -c COVERAGE_MANIFEST.sha256)
    (cd session203-s78-correction-review && sha256sum -c SHA256SUMS)
    python3 session203-s78-correction-review/check.py
    python3 -O session203-s78-correction-review/check.py

The corrected full-run and bank regeneration commands are preserved in the concise audit. The original files remain unchanged elsewhere in the repository; this packet supplies the isolated corrected source, patch, and certificate. It changes no navigation indexes or theorem status outside this packet.
