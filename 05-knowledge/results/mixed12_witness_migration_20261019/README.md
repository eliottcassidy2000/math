# Frozen mixed-(1,2) witness-migration packet

Date: 2026-10-19. All programs use the Python standard library; Python 3.9 or later is required.

## Exact scope

This packet excludes a **maximum-index-span critical majority root** with one lower slave, two upper slaves, and c≤3 common incomparable neighbors, in a finite poset having no balanced pair. It also covers the order-dual (2,1) case. Exterior size and exterior degrees are unrestricted, subject to retaining this root morphology.

It does **not** settle c=4, all seven-thin posets, or the general 1/3–2/3 conjecture. The computation is a finite rooted-core certificate lifted analytically to arbitrary actual exteriors.

Read `mixed12_witness_migration_20261019.md` for the theorem and proof, and `freeze_review_20261019.md` for the final proof/code audit. The separate independent audit covers core completeness through c=2 and the exceptional six-point inequality, including 200 whole-poset checks; the main replay covers c=3.

## Reproduction

From this directory, check file integrity and verify that normalization made no AST changes beyond replacing assertions:

    python replay_frozen_mixed12_20261019.py --check-only

Replay the full c≤3 certificate and independent exterior audit, both with Python optimization enabled, and require byte-for-byte agreement with the preserved receipts:

    python replay_frozen_mixed12_20261019.py

The full replay can take several minutes. It runs in a temporary directory, leaves frozen files unchanged, and prints a compact JSON report. A digest is not used as a substitute for any generator inequality: the verifier recomputes every selected witness and checks it against every generator before its receipt is compared.

For direct runs, the normalized verifiers retain their original interfaces:

    python -O verify_mixed12_witness_migration_20261019.py --max-common 3
    python -O verify_mixed12_exception_independent_20261019.py

These direct commands write their raw JSON receipts beside the scripts. The main raw receipt is 10,496,542 bytes. Use `--max-common 2` only for a partial replay; it does not certify c=3 and writes a differently named partial receipt.

## Frozen evidence

- `mixed12_witness_migration_audit_20261019.json.gz` is a deterministic gzip of the **exact original receipt bytes**, not a minified or summarized replacement. It decompresses to the original 10,496,542-byte audit.
- `mixed12_exception_independent_audit_20261019.json` preserves the independent receipt exactly.
- `originals/` preserves both original verifier sources and both original notes.
- The two top-level verifier copies replace 11 and 13 assertions, respectively, with explicit `if not ...: raise RuntimeError(...)` checks. All arithmetic, search order, witnesses, and receipt contents are unchanged. The replay driver validates this transformation structurally using Python's AST.
- `freeze_replay_receipt_20261019.json` records the actual optimized replay and source/receipt checksums.
- `manifest.json` records hashes and sizes for every frozen payload file. The manifest excludes itself. The surrounding deterministic archive has a separate SHA-256 checksum.

The scripts at `originals/` are provenance copies; use the top-level normalized scripts for verification. No original research file or Git state was changed by this freeze.
