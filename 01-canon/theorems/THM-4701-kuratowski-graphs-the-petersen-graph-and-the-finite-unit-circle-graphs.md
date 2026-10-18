---
id: THM-4701
title: "Kuratowski graphs, the Petersen graph and the finite unit-circle graphs G_q = Cay(F_(q^2), mu_(q+1)): any two vertices of G_q have at most two common neighbours, so no G_q contains K_5 or K_(3,3); K_4 occurs iff q is an odd power of 2; G_4 is the Clebsch graph and its ten vertices at distance two from 0, the two-sums of the residue pentagon mu_5, form the Petersen graph as the Kneser graph K(5,2); the golden unit-distance Petersen realization (outer pentagon of circumradius R, inner pentagram of circumradius r turned by a right angle, R^2 + r^2 = 1, R/r = phi) lies in Q(zeta20) and reduces at the prime above 2 exactly onto that distance-two layer; in characteristic 0 the Petersen graph is the unit-sum graph of the pentagon's two-sums; FINITE-EXACT census of Petersen copies in G_q for primes q <= 113"
status: >
  PROVED: statements 1-4 (elementary; exact checks in residue_graph_kuratowski.py, petersen_zeta20.py,
  petersen_mod2_and_unit_sum.py).
  FINITE-EXACT: statement 5, by exhaustive backtracking that is complete up to the edge-transitive symmetry of G_q.
  KNOWN:
  - G_2 = K_4, G_3 = K_3 □ K_3, G_4 = Clebsch, G_8 = H(3,4), and κ(4) = κ(8) = 4 are THM-4624.
  - The two-sum / Kneser picture of the distance-two layer of G_4 is the classical folded-5-cube model of the Clebsch
    graph, whose second subconstituent is the Petersen graph.
  - The golden realization is the 5-fold rotationally symmetric unit-distance drawing of the generalized Petersen graph
    GP(5,2). Žitnik–Horvat–Pisanski (J. Korean Math. Soc. 49, 2012) give rotationally symmetric unit-distance drawings
    of I-graphs, and Kaenthong et al. (2025) study the turn of the inner star. The spoke condition forces a right angle,
    which is R^2 + r^2 = 1.
  - That two unit circles meet in at most two points is the elementary planar fact, used here in its norm-form
    version.
  - Le Anh Vinh (arXiv math/0510092) studies the unit-quadrance graphs on F_q^2, which are G_q for q = 3 (mod 4).
  Not found in the literature checked:
  - the triangle and K_4 laws of statement 2;
  - the 2-adic reduction of the golden realization onto the Clebsch distance-two layer, and the injective reduction at
    q = 19 (mod 20);
  - the unit-sum description in characteristic 0;
  - the census.
  OPEN: the law behind the census (statement 5).
  AUDITED by audit G1 (04-computation/experiments/unit_path_ring_20261018/audit_G1/REPORT.md).
  - Statements 1-5 were reproduced independently, including the census for primes q <= 113 and for q = 9, 25, 27, 49.
  - Its MAJOR item G1-2 (the "pentagon becomes mu_5" remark was false) and minor items G1-4, G1-5, G1-11, N8 were
    applied (MISTAKE-637).
source: mac-mini-2026-10-18-unitpath (owner prompt 2026-10-18: merge the Petersen graph, K_5, K_(3,3) and the ternary tree of primitive Pythagorean triples)
depends_on: [THM-4558 (Lemma R residue colourings; κ(q) = chi(G_q)), THM-4624 (G_q for small q; Clebsch), THM-4700 (unit-path rings, Heron criterion)]
verification:
  - 04-computation/experiments/kuratowski_petersen_20261018/residue_graph_kuratowski.py (+ .out): G_2, G_4, the two-sums, common-neighbour maxima, triangles and K_4 for q in {2,...,32}, Petersen search q <= 31
  - 04-computation/experiments/kuratowski_petersen_20261018/petersen_census.py (+ petersen_census_count47.out, petersen_census_exist113.out)
  - 04-computation/experiments/kuratowski_petersen_20261018/petersen_zeta20.py (+ .out): golden realization, 60 digits, and reductions at q = 19, 59, 79, 139
  - 04-computation/experiments/kuratowski_petersen_20261018/petersen_mod2_and_unit_sum.py (+ .out)
  - 04-computation/experiments/petersen_law_20261018/petersen_count.c (+ petersen_count_400.out): census to 397
  - 04-computation/experiments/unit_path_ring_20261018/audit_G1/ (g1_gq*.py, g1_petersen_zeta20.gp, g1_port_pentagon_mod2.gp): independent reproduction
---

# THM-4701 — Kuratowski graphs, the Petersen graph and the finite unit-circle graphs

## Setting

Let `q` be a prime power. On `F_(q^2)` put `conj(x) = x^q` and `N(x) = x^(q+1)`, so `N` is the norm to `F_q` and is
anisotropic. Let `mu_(q+1) = {N = 1}` and

    G_q = Cay(F_(q^2), mu_(q+1)),

the finite unit-circle graph.

By Lemma R (THM-4558), a conjugation-fixed prime of a CM field whose residue field is `F_(q^2)`, and which is inert over
the real subfield, maps the field's unit-distance graph homomorphically into `G_q`. So `κ(q) = chi(G_q)` bounds that
graph. Throughout, "contains `H`" means `H` is a subgraph, not necessarily induced.

## Statement

**1. Two circles.** Any two distinct vertices of `G_q` have at most two common neighbours. Hence no `G_q` contains
`K_(2,3)`, and no `G_q` contains either Kuratowski graph, `K_5` or `K_(3,3)`. Both of those have two vertices with three
common neighbours.

**2. Triangles and `K_4`.** An edge of `G_q` lies in exactly:
- 2 triangles if `q = 5 (mod 6)`, or if `q` is an odd power of 2;
- 1 triangle if `3 | q`;
- no triangle otherwise.

`G_q` contains `K_4` iff `q` is an odd power of 2. Then every `K_4` is a similarity image of `F_4 ⊂ F_(q^2)`; for
instance `G_2 = K_4`, and `G_8 = H(3,4)` contains `K_4`. In odd characteristic the Euclidean obstruction to `K_4`
persists: the two apexes over a unit edge are at squared distance `N(sqrt-3) = 3`.

**3. The golden residue plane.** `G_4` is the Clebsch graph (THM-4624). Its neighbours of 0 are the residue pentagon
`mu_5 = {1, zeta, ..., zeta^4}`. Its ten vertices at distance two are the two-sums `s_ab = zeta^a + zeta^b`
(`a < b`), and they induce the Petersen graph, with

    s_ab ~ s_cd   iff   {a, b} ∩ {c, d} = ∅,

that is, the Kneser graph `K(5,2)`.

**4. The golden Petersen realization.** In `Q(zeta20)`, put
- `R = (zeta20 + zeta20^-1)/sqrt5`, `r = R/phi`;
- `o_k = R zeta5^k` (outer pentagon) and `i_k = i r zeta5^k` (inner pentagram, turned by a right angle).

Then:
- `R^2 = (5 + sqrt5)/10`, `r^2 = (5 - sqrt5)/10`, and `R^2 + r^2 = 1`. This is a golden Pythagorean triple
  `(R, r, 1)`, with `R/r = phi` and `R r = 1/sqrt5`.
- The fifteen pairs (pentagon sides, spokes `o_k i_k`, pentagram sides `i_k i_(k+2)`) are unit vectors. This is a
  unit-distance realization of the Petersen graph.
- Its rigidity matrix has rank 15, so the real realizations near it form a smooth two-parameter family modulo motions.
- **(a) Reduction at 2.** At the unique prime of `Q(zeta20)` above 2 (residue field `F_16`), the ten points reduce to
  `w·mu_5 ∪ w^2·mu_5`, where `w` is a primitive cube root of unity. These are exactly the ten non-neighbours of 0 in
  the Clebsch graph `G_4`, i.e. the Petersen layer of statement 3.
- **(b) Reduction at `q = 19 (mod 20)`.** At the conjugation-fixed primes of `Q(zeta20)` over every `q = 19 (mod 20)`,
  the reduction is an injective homomorphism into `G_q`. So the Petersen graph is a subgraph of `G_q` for all such `q`.
- **(c) Characteristic 0.** On the ten two-sums `s_ab` of the regular pentagon, `s_ab + s_cd` is a unit vector iff the
  pairs are disjoint, in which case it equals `-zeta5^e`. So the Petersen graph is the unit-sum graph of the two-sums,
  and modulo 2, where sum and difference agree, it becomes the unit-distance layer of statement 3.

**5. Census (FINITE-EXACT).** For primes `q <= 397`, the Petersen graph is a subgraph of `G_q` exactly for
`q ∈ {17, 19} ∪ {41 <= q <= 397, q ≠ 73}`.
- The census to 113 is the Python search, recounted independently by audit G1. From 113 to 397 it is the O(q^3) C
  counter `petersen_law_20261018/petersen_count.c`, which reproduces the Python counts.

The numbers of copies through a fixed directed edge, counted as injective homomorphisms with vertex 0 sent to 0 and
vertex 1 sent to 1, are:

| `q` | Copies |
|---:|---:|
| 17 | 480 |
| 19 | 768 |
| 41 | 240 |
| 43 | 3840 |
| 47 | 480 |

The count is 0 for every other prime `q <= 47`. The prime powers 9, 25 and 27 also admit no copy
(`residue_graph_kuratowski.out` (3)).

## Proofs

**1.** Translate so the two vertices are 0 and `d ≠ 0`. A common neighbour `w` satisfies `N(w) = 1` and `N(w - d) = 1`.
Subtracting gives

    Tr(w conj(d)) = w conj(d) + conj(w) d = N(d),

which is an affine `F_q`-line in the `F_q`-plane `F_(q^2)`. Since `N` is anisotropic, the conic `N(w) = 1` contains no
line, so the line meets it in at most two points.

**2.** Every edge is equivalent to `{0, 1}`, by a translation followed by multiplication by an element of `mu_(q+1)`.
The common neighbours `z` of 0 and 1 satisfy `N(z) = N(z - 1) = 1`, so `z + conj(z) = 1` and `z conj(z) = 1`: they are
the roots of `X^2 - X + 1` lying in `mu_(q+1)`.
- If `p ≠ 2, 3`, the roots are the primitive sixth roots of unity, which lie in `mu_(q+1)` iff `6 | q + 1`.
- If `p = 3`, `X^2 - X + 1 = (X + 1)^2` and there is the single root `-1`.
- If `p = 2`, the roots are the primitive cube roots of unity, which lie in `mu_(q+1)` iff `3 | q + 1`, i.e. iff `q`
  is an odd power of 2.

A `K_4` through `{0, 1}` must use both roots `z`, `z'`, which forces `N(z - z') = 1`. But
`N(z - z') = -(z - z')^2 = 4zz' - (z + z')^2 = 3`, so `3 = 1` and `p = 2`. Conversely, in characteristic 2 the set
`{0, 1, w, w^2}` is a `K_4`, because `w - w^2 = w + w^2 = 1`. Table (2) of `residue_graph_kuratowski.out` confirms
the triangle counts and the `K_4` criterion for every prime power `q <= 32`.

**3.** `Phi_5` is irreducible over `F_2`, since `ord_5(2) = 4`. So `1, zeta, zeta^2, zeta^3` is an `F_2`-basis of
`F_16`, and `zeta^4` is the sum of the other four. The connection set `mu_5` is therefore `{e_1, ..., e_4, e_1+...+e_4}`,
the folded 5-cube.
- *Distance-two layer.* A sum of two distinct elements of `mu_5` is nonzero and not in `mu_5`, because only the sum of
  all five vanishes. The ten such sums are distinct, so they form the whole distance-two layer.
- *Edges.* `s_ab + s_cd` lies in `mu_5` iff the pairs are disjoint, and then it equals the fifth root. If the pairs
  overlap, the sum is a two-sum and hence not in `mu_5`.

**4.**
- `2cos(pi/10) = zeta20 + zeta20^-1`, and `(2cos(pi/10))^2 = 2 + phi`, so `R^2 = (2 + phi)/5 = (5 + sqrt5)/10`. Then
  `r^2 = R^2/phi^2 = (5 - sqrt5)/10`.
- `R^2 + r^2 = R^2 (1 + phi^-2) = (2 + phi)(3 - phi)/5 = 1`.
- Sides: `R^2 |1 - zeta5|^2 = ((5 + sqrt5)/10)((5 - sqrt5)/2) = 1`, and
  `r^2 |1 - zeta5^2|^2 = ((5 - sqrt5)/10)((5 + sqrt5)/2) = 1`.
- Spokes: `|R - i r|^2 = R^2 + r^2 = 1`, since the two vectors are orthogonal.
- Rank 15: the 15 × 20 rigidity matrix at this point has rank 15.
  - Its smallest singular value is 1.19 with rows `d|p_u - p_v|^2` (`petersen_zeta20.py`), or 0.595 with rows
    `p_u - p_v` (audit G1).
  - The realization is induced: exactly 15 of the 45 pairs are at unit distance.

*(a)* At the prime above 2, `i -> 1` because `(1 - i)^2 = -2i`, and `zeta5 -> zeta`, a root of `Phi_5` in `F_16`.
- `zeta20 = i zeta5^4` reduces to `zeta^4`, so `zeta20 + zeta20^-1 -> zeta + zeta^4 = w`. This is a primitive cube root
  of unity, because `(zeta + zeta^4)(zeta^2 + zeta^3) = 1` and `(zeta + zeta^4) + (zeta^2 + zeta^3) = 1`.
- `sqrt5 = 2(zeta5 + zeta5^4) + 1 -> 1`, and `phi -> 1 + w = w^2`.
- Hence `R -> w` and `r -> w^2`. The image of the ten points is `w mu_5 ∪ w^2 mu_5 = F_16^× ∖ mu_5`, which is the set of
  non-neighbours of 0, since `N(x) = x^5 = 1` exactly on `mu_5`.
- Reduction is a ring homomorphism, so edges go to edges. The script checks all of this.

*(b)* For `q = 19 (mod 20)`, `Frob_q = -1` in `(Z/20)^×`, so the prime is `c`-fixed with residue field `F_(q^2)`, and
`zeta20` lands in `mu_20 ⊂ mu_(q+1)`.
- **Injectivity.** The absolute norms `N_(Q(zeta20)/Q)(p_a - p_b)` of all 45 pairwise differences are powers of 5
  (audit G1, `g1_petersen_zeta20.out`). So no difference vanishes modulo a prime not above 5, and the reduction is
  injective at every `q = 19 (mod 20)`.
- Spot checks: `q = 19, 59, 79, 139` (session) and `179, 199, 239, 359` (audit).

*(c)* If `{a, b, c, d, e} = {0, ..., 4}`, then `s_ab + s_cd = -zeta5^e`. Every one of the 45 pairs was checked
numerically: for the overlapping ones `|s_ab + s_cd| ≠ 1`.

**5.** The edge-transitivity used in 2 makes the search over copies with `0 -> 0, 1 -> 1` complete. The script
backtracks along outer cycle, spokes and pentagram. Search vertices 2, 3 and 5 have one placed neighbour each, so they
range over a unit circle. Every later vertex has at least two placed neighbours, so by 1 it has at most 2 candidates.

## Remarks

**The Euclidean picture.**
- The planar unit-distance graph excludes `K_4` (the apexes over a unit edge are `sqrt3` apart) and `K_(2,3)` (two
  circles). So neither Kuratowski graph is a unit-distance graph.
- The Petersen graph, which has both as minors, is one. Its golden realization rests on a right angle, `R^2 + r^2 = 1`,
  a golden Pythagorean triple.
- In `G_q` the circle argument survives in every characteristic, and the `K_4` obstruction survives exactly in odd
  characteristic.

**The census.**
- The real realizations form a 2-parameter family, so `F_q`-points would be expected to number about `q^2` once a
  geometrically irreducible component of the realization variety is defined over `F_q`, with its twisted (anisotropic)
  real structure.
- The data shows gaps (`q <= 37` except 17 and 19, and 73) and irregular counts; for `q > 200`, `count/q^2` lies between
  0.31 and 1.26. Possible causes (HEURISTIC):
  - components permuted by Frobenius in a Chebotarev pattern;
  - Lang–Weil error terms of a high-degree variety;
  - the removal of non-injective points.
- The golden component guarantees copies at every `q = 19 (mod 20)` by (4b). The data suggest copies for every prime
  `q >= 79`.
- Identifying all components and their fields of definition is open.

**Six-colour programme (HYP-9289).**
- The golden pentagon of ports has two-distance graph `K_5` (sides 1, diagonals phi). That is the naive block of 465
  patterns.
- The side vectors of the port pentagon are `mu_5`, and its diagonal vectors are the consecutive two-sums. Modulo 2:
  - the side vectors become the independent neighbourhood `mu_5` of 0;
  - the ten two-sums become the Petersen layer;
  - phi becomes a cube root of unity;
  - the unit-side pentagon itself becomes an induced 5-cycle of `G_4` (audit G1, `g1_port_pentagon_mod2.out`).
- So the `K_5` of the two-distance graph and the Petersen graph of the unit graph are two shadows of the same five
  directions. This is STRUCTURAL; no colouring consequence is claimed. It corrects a pre-audit sentence saying that
  the pentagon itself becomes `mu_5` (G1-2, MISTAKE-637).
