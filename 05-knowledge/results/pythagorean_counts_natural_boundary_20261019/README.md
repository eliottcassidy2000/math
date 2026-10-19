# Exact Pythagorean count laws and a natural boundary

Read [the proofs](TERNARY_COUNTS_NATURAL_BOUNDARY.md) and [independent review](review-root/REVIEW.md). The exact count tower recovers its marked ternary wall, has strict discrepancy below 2/3, and its generating function and every fixed-depth progression have the stated natural boundary. These are scoped analytic results using the cited Pólya–Carlson theorem; no transcendence claim for the density constant or global Collatz conclusion follows.

From this directory, run:

    python3 -B verify_ternary_natural_boundary.py
    python3 -B -O verify_ternary_natural_boundary.py
    python3 -B review-root/check.py
    python3 -B -O review-root/check.py

Compare author outputs with the corresponding ternary-natural-boundary receipts and independent outputs with review-root/receipt*.json. The author runs 211,067 explicit checks through digit depth 2,048 and complete tree depth ten. The independent reviewer checks 266 identities, including direct 3-by-3 Berggren matrix enumeration through depth eight. SOURCE_PROVENANCE.json records reviewed source hashes and the scope-only publication edits.
