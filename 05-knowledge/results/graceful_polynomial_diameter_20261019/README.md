# Polynomial diameter for every graceful upper-tail fiber

Read [the proof](POLYNOMIAL_STAIRCASE_DIAMETER.md) and [independent review](review/INDEPENDENT_REVIEW.md). Tight rank sets must be suffixes. Their chain decomposes the transversal problem into buffered strict-rank blocks, giving general diameter at most s²−s+1 and graceful upper-tail diameter at most r²−3r+3 for r≥3 (at most one for r=1,2).

The strict-rank theorem and its distance bound are attributed to Fon-Der-Flaass. This is a universal path-length theorem; neither optimality nor a general polynomial-time oracle implementation is claimed.

Run from this directory:

    python3 -B verify_polynomial_staircase.py
    python3 -B -O verify_polynomial_staircase.py
    python3 -B review/independent_controls.py
    python3 -B -O review/independent_controls.py

Author receipts record 29,727 checks; independent receipts record 202,884 checks. Each normal/optimized pair agrees exactly. The bounded controls test the implementation; the proof and cited theorem establish the universal claims.
