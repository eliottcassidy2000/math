# Affine fixed-cardinality transfer

Elementary finite-group proofs, exact application identities, and a conditional refinement of the supplied Carmichael counting transfer. Arithmetic conclusions depend on the supplied paper's analytic inputs (S),(A). No independent proof of those inputs, clean Lean build, or historical novelty claim is made.

## Contents

- `FIXED_CARDINALITY_TRANSFER.md`: coefficient extraction, affine pruning, translated Fejer control, quotient cardinality period, and conditional counting exponents.
- `APPLICATIONS_AND_EXACT_CODES.md`: binary control codes, complementary Pascal duality, tournament parity, restricted lifting, and identity-positive reservoirs.
- `EQUATION_CORRESPONDENCE.md`: correspondence to the supplied paper's numbered equations.
- `supporting/PRIME_POWER_BUDGET.md`: exact prime-power envelope and logarithmic bookkeeping.
- `review/AUDIT.md`: independent review, including the resolved retained-coordinate counting issue.
- `review/verify_transfer.py`: author-side finite replay.
- `review/independent_controls.py`: separately implemented exact and numerical controls.
- `SOURCE_METADATA.json`: hashes for the included proof files and pinned provenance of the supplied paper.

## Replay

Requires Python 3.10 or newer; only the standard library is used. From this directory:

    python review/verify_transfer.py
    python -O review/verify_transfer.py
    python review/independent_controls.py
    python -O review/independent_controls.py

The scripts also work when invoked by absolute path from another working directory. The independent script writes `review/controls.json`; each script prints its receipt. Matching normal and optimized receipts are included. All checks use explicit exceptions and remain enabled under `-O`.

No network access or original paper file is needed. The PDF, extracted text, and rendered pages are deliberately absent. Their hashes pin the reviewed source without redistributing it. Included proof files are verified using their actual relative paths.

Integer and rational checks are exact for the tested instances. Floating-point and Decimal controls are numerical evidence. Neither constitutes a proof of the analytic assumptions or a portable formalization audit.
