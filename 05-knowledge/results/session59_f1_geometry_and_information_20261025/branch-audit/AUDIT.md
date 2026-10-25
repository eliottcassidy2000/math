# Session 59 independent audit: branch arrangement (second gate)

**Verdict: PASS**, 2026-10-25. No mathematical correction or scope restriction beyond those already explicit in the packet is needed. No Git operations were performed.

## Frozen targets

- `/tmp/session59-f1-wildcards/branch-addendum.md`: SHA256 `89f858f8871a6608dee911949fc6f91029d4d16f85cef8bc255a4c193028e094`.
- `/tmp/session59-f1-wildcards/check_branch_arrangements.py`: SHA256 `f18abf94bddf6f5b01c04c35cda2cc83a88607986081cebd9b1bdd129c4bf97e`.
- `/tmp/session59-f1-wildcards/branch-certificate.json`: SHA256 `61858eb87888e462c64c91ffa594c21e04aa8cf6c08ad3fb1d720bddeb47ff83`.

## Mathematical checks

1. The universal group magic normal form requires no division by 2 or 3. Summing the two diagonals and the middle column yields `3T = 2T + 3E`, hence `T = 3E`; opposite pairs and the first row/column recover the stated entries. This works for every abelian group, including mixed torsion.
2. Integral Smith transformations act invertibly on powers of every abelian group. A rank-r subsystem with factors d_i therefore has `N^(3-r) product h_(d_i)` solutions. Inclusion-exclusion over the explicitly specified 512 subsets is an exact all-group proof, not empirical interpolation.
3. The signed Smith table and its resulting all-nonidentity formula are correct. A new verifier uses permutation determinants of arbitrary minors, not the author's determinant expression or SymPy. It matches all 512 published invariant-factor records.
4. Independent full center enumeration for 22 finite abelian groups matches the formula. Cases include `C2 x C6`, `C2 x C4 x C4`, `C3 x C9`, `C4 x C4`, and `C5 x C5`. The cyclic and elementary-abelian group structures are not confused.
5. The characteristic-2 and characteristic-3 polynomials follow from the correct additive-group torsion sizes, not cyclic groups of order q. The triple-determinant histogram gives exactly the stated characteristic-dependent changes. The multiplicative residue-class polynomials and the examples q = 3,4,5,7,11 are correct.
6. For the squarefree basis of the nine-root algebra, the trace pairing vanishes off the diagonal and has diagonal entries `512 product z_i^(epsilon_i)`. Its determinant is exactly `2^4608 (product z_i)^256`, including over rings with zero divisors by universal specialization. Independent actual multiplication-matrix trace computations at seven specializations, including zero and negative radicands and rank 32, agree with the analogous general-r formula.
7. After inverting 2, the finite-free cover is etale exactly on the complement of all nine z_i = 0 divisors. In characteristic 2 its relative differential module has nine unconstrained generators at every point, so it is nowhere etale. The transformation `w_i = 1+t z_i` identifies the generic branch loci with w_i = 1 but does not preserve ordinary squarehood.
8. A finite split-torus stratification would give nonnegative coefficients in powers of q-1. The additive count polynomials have a negative quadratic coefficient. For multiplicative counts, even extension degrees fix all prime-to-characteristic 2,3,4 torsion and give the stated three polynomials, each with coefficient -9 on `(q-1)^2`. Infinitely many extension sizes force polynomial equality, yielding the claimed contradiction over any finite base field.

## Cited topology gate

Primary source inspected directly: [Liu, Tran, Yoshinaga, arXiv:1707.04551v3](https://arxiv.org/html/1707.04551v3), formulas (2.1), (2.2), and Corollary 5.6. These give the ordinary and toric Poincare specializations and identify the fully split constituent with the arithmetic characteristic polynomial. For this rank-3 list, substitution produces respectively `1+9s+28s^2+20s^3` and `1+12s+57s^2+105s^3`. Euler characteristics 0 and -59 follow. These are cited applications, not inferred from point-count resemblance. The nonproper family offers no smooth-proper invariance contradiction.

## Reproduction

- Fresh isolated author-checker copies in `replay-normal/` and `replay-optimized/` pass normally and with `python3 -O`. Their certificates agree byte-for-byte with one another and the published certificate.
- `check_independent.py` passes normally and under `-O`; outputs `independent-normal.json` and `independent-optimized.json` agree byte-for-byte.
- All checks use explicit exceptions rather than optimization-sensitive assert statements.

## Scope and open loops

The branch theorem counts labelled magic arrays whose entries are not the group identity, with repeats allowed. These are base-complement counts, not square-root-cover point counts, not injective counts, and not integer square-valued solutions. The split-torus obstruction applies to these open complements only; it is not a general no-F1 theorem. The ordinary-square arithmetic remains separate.

Second-gate blocker: none. Third-gate distinctness addendum was not yet present at this audit's creation and will be reviewed separately; its provisional group formula already passes the independent mixed-group tests.
