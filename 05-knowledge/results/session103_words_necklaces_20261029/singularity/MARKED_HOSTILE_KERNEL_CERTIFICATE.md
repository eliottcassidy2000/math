# Same marked insertion, different singularity: explicit kernel certificate

**FINITE-EXACT**, independently computed by subset Hamilton-path DP and permutation/path-cut enumeration. This is the parent's separately defined marked witness, **not an identification with the recovered five-letter insertion grammar**.

## Exact objects and observer

The two order-five parent tournament adjacency matrices have rows

    T−: 01000  00100  10010  11001  11100
    T+: 01000  00100  10011  11001  11000.

Thus only the edge between 2 and 4 changes. Apply the same marked insertion to each: add vertex 5 with 5→0 and j→5 for j=1,2,3,4. Both parents and both children have exactly one Hamiltonian cycle. Both parents have sorted score sequence (1,1,2,3,3) and det M2=1. Both children have sorted score sequence (1,1,2,3,4,4), but their native two-run determinants are respectively 0 and 1.

The observer M2 has rows indexed by ordered START pairs and columns by ordered END pairs, in lexicographic order. It counts ordered spanning covers by two nonempty directed paths, singletons allowed. The singular child's dimension is 30, rank is 26, and nullity is 4.

## Full primitive integer kernel

Write q_ij=e_(i,j)−e_(j,i) in the **end-pair column** coordinate space. A basis of the full integral right kernel is

    v1 = q02 − q12,
    v2 = q03 − q13,
    v3 = q02 + q25,
    v4 = q03 + q23 + q35.

Every vector has coefficients ±1 and is primitive. More strongly, these four generate the entire integer kernel lattice: they span the rational kernel, and their coordinate matrix on (2,1),(3,1),(5,2),(5,3) is diagonal with entries 1,1,−1,−1. Any integral kernel vector therefore has integral coefficients in this basis.

All four vectors are antisymmetric under run swap. The plus sector is invertible with determinant 16, while the minus sector has dimension 15 and rank 11. Thus the whole kernel is in the minus sector. This sector calculation alone does **not** establish the native integral Smith form.

## Native torsion: explicit unit pivots

The checker performs integer row/column swaps, row sign changes, and integer row additions, always choosing a pivot ±1. It finds exactly 26 such pivots and leaves a 4×4 zero residual. The pivot block is upper unitriangular up to sign, so integral back-elimination and column elimination reduce the matrix to

    diag(I26, 0_4).

Consequently the native Smith form is **1^26 0^4**, and its cokernel is Z^4: this witness is **singular but torsion-free**. The plus-sector determinant 16 must not be misread as native torsion; the integral sector glue matters. The complete unit-pivot trace is in the receipt.

## Schur residual and exact lifts

For the parent of the singular child T−, H=M1(T−) is unimodular. Order the inserted observer by old ordered pairs first, then (5,u), followed by (u,5), u=0,…,4. Its lower-right block is F=diag(H,H). Put R=C−DF^(-1)E and Z=F^(-1)E.

The 20×20 residual R has nullity four. The following integral residual basis is returned by exact elimination; the table also gives its new-port correction −Zz:

| old-pair vector z | new-port correction | full lifted vector |
|---|---|---|
| −q02 | −q25 | −v3 |
| −q12 | −q25 | v1−v3 |
| q03−q13 | 0 | v2 |
| −q03−q23 | −q35 | −v4 |

Every lift is checked directly against the original 30×30 native observer, not merely against a determinant identity. The universal insertion formula using M1/M2/M3 is separately checked for both parents.

## Reproduction

    python3 -O session103_words_necklaces/singularity/check_marked_hostile.py

Input: `../necklaces/results.json`. Output: `marked_hostile_receipt.json`, containing input SHA-256, exact sparse kernel vectors, all Schur lifts, native Smith form, and unit-pivot trace. All checks use explicit exceptions and remain active under optimization. No publication was performed by this lane.
