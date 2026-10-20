# Guarded reset-debt composition

Date: 2026-10-20.

**PROVED:** elementary chronological affine-debt and slack identities,
exact odd-valuation source/image guards, and sufficient arbitrary-finite-word
certificates for the stated first-reset and grouped Collatz banks.

**FINITE-EXACT:** two independent exact validation paths with matching
normal/optimized receipts. **OPEN:** arbitrary-source coverage, the uniform
native signed-row residual estimate, and Collatz convergence.

Start with [the mathematical note](RESET_DEBT_COMPOSITION.md), then read
[the independent audit](PUBLICATION_AUDIT.md). The examples include a growing
reset paid by a later reset, a cylinder descending without a unit-budget
certificate, and a 2-adic infinite address with no infinite positive-integer
source. No general literature-novelty claim is made for affine composition.

## Reproduction

Both checkers require Python 3 with its standard library only. Neither
imports the other checker or needs the inherited research files to run.

    python3 -B verify_reset_debt.py > /tmp/reset-producer.json
    python3 -B -O verify_reset_debt.py > /tmp/reset-producer-optimized.json
    cmp /tmp/reset-producer.json /tmp/reset-producer-optimized.json
    cmp /tmp/reset-producer.json receipt.json
    python3 -B independent_audit.py > /tmp/reset-independent.json
    python3 -B -O independent_audit.py > /tmp/reset-independent-optimized.json
    cmp /tmp/reset-independent.json /tmp/reset-independent-optimized.json
    cmp /tmp/reset-independent.json independent.normal.json
    sha256sum -c SHA256SUMS

The producer records 321,690 explicit checks. The independent audit records
148,619 explicit checks and additionally enumerates 837,930 odd source
representatives. Finite computations validate the formulas; the note's
proofs, not extrapolation, establish the arbitrary finite-length statements.

## Provenance

[SOURCE_INVENTORY.json](SOURCE_INVENTORY.json) identifies the inherited
repository sources and their exact hashes; its paths are repository-root
relative. [VALIDATION.md](VALIDATION.md) describes the unchanged producer
validation universe.

Only two publication sentences were clarified: the sign hypotheses of the
defect formula and the ambient space of the repeated-A intersection. The
[original note](provenance/RESET_DEBT_COMPOSITION.original.md),
[complete diff](publication_overlay.diff), and
[producer hashes](provenance/FROZEN_PRODUCER_HASHES.json) preserve that lineage.
The independent checker is self-contained and the mathematical theorem is
unchanged.
