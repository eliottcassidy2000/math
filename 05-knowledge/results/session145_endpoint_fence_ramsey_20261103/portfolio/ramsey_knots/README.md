# Session 145: sparse Ramsey barriers and honest knot interfaces

Start with [PROOFS.md](PROOFS.md). The tangle scope correction is recorded in [CORRECTION_LINEAGE.md](CORRECTION_LINEAGE.md). The new scoped deductions are:

- exact forest/pseudoforest restricted Ramsey maxima 2,3,4,5 times (k−1), with explicit sharp constructions;
- the support-layer obstruction n≤2r(k−1) when an r-factor attachment-polynomial support is explicitly used as the non-green graph;
- the matching closure polynomial product_{j=0}^{r−1}(z+2j), its one-component coefficient, and contextual matching minimum;
- an explicit same-matching unknot/trefoil hostile, verified with Fox coloring;
- independent recovery of the true Paley-11 roles of 1729 and 17290.

The current source-audited Ramsey ranges are 45–57 and 61–86, with 86 clearly labelled a computer-assisted preprint bound. This session proves neither exact Ramsey value and does not claim an unrestricted improvement.

Run, from this directory:

    python -B verify_structures.py > verification.normal.json
    python -B -O verify_structures.py > verification.optimized.json
    cmp verification.normal.json verification.optimized.json
    python -B independent_checks.py > independent.normal.json
    python -B -O independent_checks.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json
    python -B verify_packet.py

The two scripts use different algorithms and share no implementation. This is independent-algorithm testing, not an external review or formal proof. All correctness gates survive -O. Full historical source copies are pinned in SOURCE_PINS.json; current external claims and their verification boundary are in SOURCES.md.
