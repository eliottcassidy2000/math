# Independent audit: Gaussian orbit compression

## Outcome

PASS for the exact Gaussian orbit, decoder, invariant-ring and joint-frame results, after explicit nondegenerate Berggren boundary guards. This is a scoped proof audit, executable replay and independent finite verification, not a re-audit of every inherited theorem.

## Proof checks

- The disk and shell orbit counts correctly retain the origin and the axis/diagonal stabilizers. The exact fundamental-octant point count is not N/8. The limiting ratio follows because the boundary terms are O(sqrt(M)) while N(M) is bounded below by a constant times M.
- O(2,Z) is exactly the signed-permutation group. This prohibits additional global Euclidean lattice isometries, not stronger arithmetic encodings. The arbitrary-invariant-subset counting lower bound and the radius-only complete-disk specification are different problems and are kept separate.
- Primitive-ray decoding includes the axis and diagonal primitive directions and reconstructs the origin separately. Rational shell rotations use partial lattice maps; they do not falsely enlarge the global integral orthogonal group.
- Q[x,y]^D4=Q[x²+y²,x²y²] follows by sign invariance and the elementary symmetric-polynomial theorem. The two invariants separate real D4 orbits. The integer realizability test checks both squared roots, not only the discriminant.
- The quartic harmonic identity is exact. The first multiple-octant shell is norm 25, with (5,0) and (4,3); this is a finite boundary fact, not a large-scale conjecture.
- Independent orbit compression of summands does not descend addition. The pair (1,1) and pair (1,-1) have identical separate D4 orbit coordinates but unequal sum norms. A simultaneous frame or explicit relative alignment is necessary.
- The modulo-three golden/Gaussian map preserves conjugation and norm as well as multiplication. The modulo-two nilpotent obstruction is valid.

## Normalization corrections

The audit found two minor scope boundary omissions in the original prose and requested their repair without editing the source packet:

1. Ordinary Berggren trees parameterize nondegenerate positive primitive triples. The axis ray (1,0) is not a tree node; it must be handled separately. Thus (1,1) is the smallest missing primitive ray only among rays whose two coordinates are nonzero.
2. The Pell terminal pair (m,n)=(1,1) gives (0,2,2), normalized to the degenerate triple (0,1,1). It may be a descent terminal, but must be excluded from the positive Berggren-ray claim, which uses 0<m<n.

Odd/odd coprime Euclid parameters have content two. That source normalization and the inherited corrected D,D,U chamber word are retained. No arbitrary square-sum island/Berggren branch equivalence is inferred.

## Sources and Bott scope

All four repository source hashes are checked against their pinned commit. The statements about Bott and Clifford periodicity were also checked against the cited primary PDFs:

- [Bott, The Stable Homotopy of the Classical Groups (1959)](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf).
- [Atiyah–Bott–Shapiro, Clifford Modules](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf).

Those papers concern stable classical-group homotopy and Clifford-module/K-theory periodicity, respectively. They do not identify the eight real stabilization degrees with the eight elements of the planar lattice symmetry group. The source's distinction is correct.

## Independent controls

The supplied checker was replayed normally and under -O, with identical output matching its stored receipt. A separate checker, check_independent.py, verifies:

- 251 disk orbit counts by a direct Burnside fixed-point sum;
- both directions of quartic realizability for 84,676 candidate pairs, including unrealizable data;
- the quartic harmonic identity on a signed integer box;
- every modulo-three element's conjugation and norm compatibility;
- 1,101 coprime Euclid parameter pairs' exact primitive content;
- the complete positive Pell locus in the square 1<=m,n<=500;
- all manifest and inherited source pins.

The final independent receipt pins the audited source manifest. Both execution modes match. No Git or other external writes were performed.
