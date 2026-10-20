# Independent audit: uniform equivariance obstruction

Date: 2026-10-20 UTC.

## Verdict

**PASS.** The proof in `session12-obstruction/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md` proves both assertions of HYP-9172(a) for every integer n >= 5. No mathematical repair is required. It concerns maps on labelled objects that are equivariant under S_n and induce a bijection of isomorphism types. It does not rule out arbitrary bijections of unlabelled types, and it does not settle HYP-9172(b).

The reviewed file's SHA-256 is `1ffbff76190de3c2603b9512d98aa5192c5454d5d97da358c101fd62ce562db2`.

The audit used independent derivations and fresh code. The producer's computation was not read or reused. All audit files are outside the repository checkout; the repository and producer files were read only. The relevant source definitions in HYP-9172, THM-4531, and THM-4721 were checked.

## 1. Definitions, affine action, and signs

Fix the upward tournament T_0, and encode a tournament T by the binary edge vector t recording reversals from T_0. For g in S_n, let L_g permute unordered edge coordinates and let i_g be the binary vector whose {u,v} coordinate, u < v, is 1 precisely when g(u) > g(v). Then

    g(t) = L_g(t + i_g).

If U is an S_n-invariant binary edge subspace, this affine action descends to the quotient A = V/U. Its size is c = 2^(binomial(n,2) - dim U).

For X in U-perp, the function phi_X([t]) = (-1)^<X,t> is well defined. The c such functions are a basis over R = Z[1/2]: their character matrix multiplied by its transpose is c times the identity. With the contragredient function action,

    g phi_X = (-1)^<i_g,X> phi_(gX).

The exponent is exactly the Royle edge-reversal sign. For a graph stabilizer it is a homomorphism, independent of changing the reference orientation. Changing T_0 only multiplies individual characters by signs and does not change any asserted invariant-basis or determinant conclusion.

For a graph orbit, invariance of the character coefficients permits a nonzero coefficient precisely when all its stabilizer signs are +1. If a stabilizer sign is -1, the coefficient equals its negative, so it vanishes in R. For a sign-even graph, coefficient transport yields a signed orbit sum with coefficients +1 or -1. There is no averaging by |S_n|, no division by odd primes, and no semisimplicity assumption over R.

Thus the signed graph-orbit sums and the source-orbit indicators are two R-bases of the same invariant module. The evaluation matrix M is square, integral, and invertible over R. In particular, det M = +/- 2^e for a nonnegative integer e.

For U = 0, these are ordinary tournaments and sign-even graphs. For U equal to the cut space, A is the usual tournament-switching quotient and U-perp consists exactly of graphs having even degree at every vertex. The target then is precisely the sign-even Eulerian graphs in HYP-9172(a).

## 2. Orthogonality and determinant orientation

Let H_a be a source stabilizer and K_X a graph stabilizer, and use one representative per isomorphism type. Two signed graph-orbit sums have disjoint character supports unless their graph types coincide. Their squared norms over all c labelled quotient elements are c times the number of labelled graph copies. Consequently

    sum_[a] (n! / |H_a|) M[a,X] M[a,Y]
      = c (n! / |K_X|) delta_(X,Y).

After cancelling n!, this is exactly

    M^t D_A^-1 M = c D_B^-1,
    M D_B M^t = c D_A,

with D_A = diag(|H_a|), D_B = diag(|K_X|). Taking determinants gives

    det(M)^2 product_[X] |K_X| = c^N product_[a] |H_a|.

This verifies the direction of every quotient in the proposed argument. It follows that the two stabilizer products have equal valuation at every odd prime. Negative exponents are permitted in their ratio, which is a power of 2.

## 3. Equivariance and simultaneous preservation

For an equivariant map f, every h stabilizing a also stabilizes f(a):

    H_a <= K_f(a).

This is source contained in target, not the converse. The map need not be injective on labelled objects. Equivariance maps a source orbit onto its associated target orbit, and the assumed bijection of orbit sets pairs all source and target types once.

Each odd-prime valuation difference v_p(|K_f(a)|) - v_p(|H_a|) is therefore nonnegative. The product identity makes their sum zero. Each difference is zero. This is valid simultaneously for all odd primes for this one hypothetical map; it does not combine independently chosen prime-local matchings.

The empty graph is sign-even and Eulerian, has stabilizer S_n, and must be hit by a source type. That type's stabilizer H must carry the full 3-part and the full 5-part of n!.

## 4. Source stabilizers and local group argument

A nonidentity involution cannot stabilize a tournament, because any swapped pair reverses its own arc.

For a tournament switching class and fixed vertex v, there is exactly one representative in which v is a source: switch the other vertices pointing into v. If two representatives both have v as source, their difference is a cut containing no edge incident with v, hence the empty cut. Every class automorphism fixing v therefore stabilizes this representative. In particular, a nontrivial involution fixing v is impossible.

For n >= 5, every double transposition has a fixed vertex. Thus both relevant source stabilizers contain no double transposition.

Let H <= S_n have this last property. Two single 3-cycles whose supports meet in two points yield a double transposition after possibly inverting one: (1 2 3)(1 2 4) = (1 3)(2 4), using right-to-left multiplication. If their supports meet in one point, conjugating the first by the second produces a 3-cycle whose support meets the first in exactly two points. Hence distinct supports of single 3-cycles in H are disjoint.

If H has the full ambient 3-part, one of its Sylow 3-subgroups is a Sylow 3-subgroup of S_n. Every such subgroup contains k = floor(n/3) disjoint single 3-cycles: the subgroup generated by any k disjoint 3-cycles extends to an ambient Sylow subgroup, and all ambient Sylow subgroups are conjugate. These supports exhaust the possible single-3-cycle supports in H, since fewer than three vertices remain.

Conjugation by H preserves this collection of k triples and its residual r-set, where r <= 2. Therefore

    H <= (S_3 wr S_k) x S_r,
    v_5(|H|) <= v_5(k!).

For n >= 5, v_5(k!) < v_5(n!): n = 5,6,7 are immediate from the factor 5, and for n >= 8 the interval (k,n] contains at least five consecutive integers. This conflicts with the full 5-part forced by the empty graph. The proof is uniform and requires no classification theorem or unproved additive-prime hypothesis.

At n = 5 specifically, the unique 3-cycle support forces H into S_3 x S_2, whose order is not divisible by 5. This is a complete treatment of the formerly exceptional small case.

## 5. Independent exact controls

Reproduction from the parent workspace:

    python3 -B session12-obstruction-review/independent_controls.py
    python3 -B session12-obstruction-review/group_controls.py

The first program independently enumerates all binary edge subspaces through n = 4, selects the S_n-invariant ones, and tests the exact affine quotient, the signed orbit basis, both integer-scaled Gram identities, and the determinant identity. There are 17 invariant-quotient instances: respectively 1, 1, 2, 4, 9 at n = 0,1,2,3,4. At n = 4 the invariant selection is from all 2,825 subspaces of F_2^6.

For tournament and cut-space quotients through n = 5, it also tests stabilizer containment up to every possible relabelling of the target and computes a maximum bipartite matching. The containment graph has a perfect matching in both cases for every n <= 4, independently verifying the claimed sharpness. At n = 5:

- Tournaments: 12 types; maximum matching 11; source stabilizer orders are seven 1s, four 3s, and one 5.
- Switching classes: 2 types; maximum matching 1; source stabilizer orders 5,3; target stabilizer orders 120,8.

The exact tournament stabilizer list is

    [1, 3, 1, 3, 1, 1, 1, 3, 1, 3, 1, 5],

namely seven 1s, four 3s, and one 5. Its total odd valuations are v_3 = 4 and v_5 = 1, agreeing with the sign-even graph stabilizer product. The quotient determinant is 2^49. For switching classes the total odd valuations are v_3 = v_5 = 1, and the determinant is -8.

The second program checks the explicit double-transposition construction for every ordered pair of overlapping single 3-cycles through n = 7. It also verifies directly on all labelled objects through n = 5 that no tournament has a nontrivial involutory automorphism and no switching class has an involutory automorphism fixing a vertex. All eight labelled switching classes at n = 4 admit a double transposition; none do at n = 5. Thus the n = 4 fixed-point boundary is genuine.

All checks passed. These computations are controls of the proof, not an extension by finite census.

## 6. Artifact checksums

- `independent_controls.py`: `814b304eb19db1927157c50c2eccdf3fdd4824b4d07c81425d3fc33b88474575`
- `independent_controls.json`: `0e96c7ece5c2b4a3098b703867103e69d7b7c189ce23d91a1298c4304797819e`
- `group_controls.py`: `440c3fadb2311eae2d1bbf0ad8af0d7e8d4e916894a5d6e3cba8d2c8571d057e`
- `group_controls.json`: `1c464ca3c37f1201a97b55e28cab0910a3e892f1c1c5faa8757e404b1e7ebd1e`
