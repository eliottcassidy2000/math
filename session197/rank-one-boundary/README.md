# Signed root boundary and rank-one endpoint factors

For every rooted tournament, the alternating root-run sum is an explicit signed Hamilton-path sum on the exterior. At crossing matching rank one it is exactly minus a product of two genuine endpoint counts, with no symmetry assumption. Rank-two examples have both signs, so the result is not extended to ternary unit clocks.

- [Proof](PROOF.md) and [independent PASS review](REVIEW.md).
- [Standalone C++ verifier](verify.cpp): compile with C++17 and run; optimized/debug receipts are included.
- [Independent literal Python verifier](independent/check.py): run with Python3, with or without `-O`.

Neither verifier has an external code or data dependency. The source checks 9,774 roots and 6,241 star centers; the independent literal review checks 5,505 roots and 4,975 centers. Frozen source and review bytes are preserved. Original scratch-style paths in the review identify the same-hash files bundled here. No full-spectrum or unproved unit-offset congruence claim is made.
