# Dyadic re-entry packet

- THEOREM.md: exact dyadic layers, full small-array identities, AP threshold, explicit failed lift, and exhaustive divisor counterexamples.
- verify.py: independent top-row enumeration and certificate replay; no imports from the frozen checker.
- certificate.json: full local mask sets and representatives, exact binary compatibility matrices, all 203 candidate results, counts, and witness.
- SOURCE_PINS.json: external theorem/certificate dependencies with hashes.
- SHA256SUMS: this packet's file integrity list.

Reproduce from the workspace root:

    python session55-dyadic-reentry/verify.py --verify
    python -O session55-dyadic-reentry/verify.py --verify

The frozen sibling directory session55-crt-collision-wiring is read-only input. The unbounded completeness claim uses its reduction proof and the pinned prime-field theorem; exhaustive finite tables alone do not establish that reduction. No novelty claim is made.
