# Nonlinear residue banks: finite-exact search

**FINITE-EXACT, k=2 through 14 only. No all-k linear-length theorem is claimed.**

Use the one-apex count, modular state, and insulation lemma (15) from `PROOF.md`, section 6. Its lower bound (4) concerns only independent additive disjoint adjacent-switch banks, not this nonlinear search. Its certificate barrier (14) concerns the unchanged S192 magnitude bound.

The verifier exhausts all words and proves the following finite minimal values of `m` (all binary words, not only strong words):

    k:          2  3  4  5  6  7  8  9 10 11 12 13 14
    minimum m:  2  5  6  7  9 10 12 13 14 15 16 17 19
    bank 2k+m:  6 11 14 17 21 24 28 31 34 37 40 43 47

These are **FINITE-EXACT**: the script records all `2^(k-1)` residue witnesses and verifies that the preceding length misses a residue. It tests (15) for every stored witness in three exteriors, and the simultaneous selector for every odd residue at three target sizes for each `k<=14`. At `k=14` the nonlinear insulated bank is 47 bits, versus the sharp independent-switch length 222.

This finite table refutes any attempt to extend (4) to arbitrary banks. It does not prove `m=O(k)` for all `k`, even though that is an attractive next target. Even an all-`k` linear bank would not cross the unchanged certificate barrier (14). The genuine remaining question is a *coupled* decoder controlling the carry or using an exterior responsive to it.

