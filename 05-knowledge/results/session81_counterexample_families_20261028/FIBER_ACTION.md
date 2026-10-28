# The full observer fiber is a symmetry action, not a missing invariant

**Status:** PROVED deductions from cited decomposition inputs; FINITE-EXACT
order-12 action census. No literature-novelty claim. This extends
[LIFTING.md](LIFTING.md), not its already established observer equivalence.

## 1. Exact coordinates and actual acting group

Let H be a realizable 3-uniform hypergraph on a fixed finite labelled set V.
Use the hypergraph module convention in Definition 1 of
[Boussaïri–Chergui–Ille–Zaidi](https://arxiv.org/html/1805.04607v1).
Their Theorems 10, 16, 24 and Proposition 51 supply the common strong-module
tree, prime/linear quotient alternatives, prime realization uniqueness up to
converse, and independent quotient-realization coordinates, respectively.
These structural inputs are cited, not newly proved here.

Let P be its prime internal nodes, L its linear internal nodes, and C_v the
set of children of node v. At each prime node choose one reference quotient
Q_v. The complete labelled realization set is naturally in bijection with

    X_H = product_(v in P) {Q_v,Q_v*} × product_(v in L) Orders(C_v).   (1)

Consequently |X_H|=2^|P| product_(v in L) |C_v|!. Empty/singleton sets have
one realization, with an empty product. Coordinates are distinct: the arc
between two leaves is read at their lowest common ancestor.

Put G=Aut(H), the **hypergraph** permutation group. It acts on the tree and
on each transported child set. If g sends v to w, its child bijection sends
Q_v either to Q_w or Q_w*. Write that sign as s(g,v) in F_2. Then the prime
coordinate action is

    (g·epsilon)_(gv) = epsilon_v + s(g,v).

The signs satisfy the cocycle relation
`s(gh,v)=s(g,hv)+s(h,v)`; they are not independent freely chosen group signs.
At a linear node, g simply transports the chosen total order of its children.
Changing reference Q_v conjugates this coordinate action and changes no
counts. Realization and coordinate transport commute. In particular,

    Aut(T_x) = Stab_G(x),           ordinary isotypes = X_H/G.          (2)

Proof: an ordinary tournament automorphism preserves every cyclic triple,
hence lies in G. Inside G, it preserves T exactly when all quotient choices
are preserved, by the injectivity of (1). The same argument applies to
isomorphisms between two realizations on V.

Thus the complete subset-Hamiltonian-path observer determines **a full
space of possible ordinary groups and their distribution**, but not the
chosen group's value. Its loss is the position x in the fiber, not a failure
to determine what fibers are possible. This uses the labelled observer iff
from [LIFTING.md](LIFTING.md); an unaligned bag of subset counts is not enough.

## 2. Fixed-subgroup formula and exact distribution algorithm

For K≤G let K_v={g in K:g(v)=v}. Choose one representative v from each
K-orbit of internal nodes. Let f_v(K_v) be:

* prime v: 2 if every element of K_v has sign s(g,v)=0, and 0 otherwise;
* linear v with k children: k! if every element of K_v fixes every child,
  and 0 otherwise.

Then

    F(K):=|X_H^K| = product_(v in internal nodes/K) f_v(K_v).           (3)

Indeed, a fixed assignment is determined by its choice at one node of each
orbit. Transport is consistent precisely when that choice is fixed by K_v.
A two-point orientation set has both or neither fixed; a finite total order
has no nonidentity order-preserving child permutation. Choices on distinct
node orbits remain independent. Note especially: a permutation cycling
nodes can still have a sign obstruction after completing its cycle.

This gives two finite, exact procedures without enumerating tournament
isotypes:

1. Enumerate coordinate assignments x and compute their stabilizers in the
   known permutation group G using (1). This can still be exponentially
   large; no efficiency or census-free input-generation claim is made.
2. Enumerate subgroups of G, calculate (3), and descend the subgroup lattice:

       E(K)=F(K)−sum_(J>K) E(J).

   Then E(K) counts assignments with **exact** stabilizer K. Equivalently
   use subgroup-poset Möbius inversion, or the table of marks on conjugacy
   classes. The number of realization isotypes with stabilizer conjugate to K
   is `E(K)|K|/|N_G(K)|`. Sum this over conjugacy representatives with a chosen
   group order or abstract group isomorphism type to obtain that distribution.

A complete subgroup action or table of marks can be costly. Aggregate
cycle-index coefficients alone do not supply the required signed node
transports, subgroup incidences or stabilizers. Formula (3) uses full H group
data, not just an ordinary/twisted aggregate prime index. Burnside requires
only F(<g>) to count isotypes, but those element fixed counts alone generally
do not determine a full stabilizer distribution.

A useful mass identity, independent of subgroup enumeration, is

    sum_([T] in X_H/G) 1/|Aut(T)| = |X_H|/|Aut(H)|.                    (4)

This is orbit–stabilizer, not an assertion that all ordinary groups agree.

## 3. The entire old order-12 fiber, not just eight selected tournaments

Let A=singleton→C3 and H_1=C3(C3[A,A,A]). Its canonical tree has one root
prime C3, three linear nodes with distinguishable singleton/C3 children,
and three leaf prime C3 nodes. Therefore

    |X_H1|=2^(1+3+3)=128,
    Aut(H_1) ≅ (S3)^3 ⋊ S3,       |Aut(H_1)|=6^4=1296.

Here the outer S3 permutes the three four-vertex blocks; each inner S3
permutes the vertices of its leaf triangle and fixes its isolate. The linear
children cannot be exchanged, since their cardinalities differ.

Write coordinates `(r; l_0,l_1,l_2; c_0,c_1,c_2)`: r orients the root,
l_i chooses source/sink placement, c_i orients the leaf triangle. If an H
automorphism has block permutation pi and local triangle permutations rho_i,

    r' = r + parity(pi),
    l'_(pi(i)) = l_i,
    c'_(pi(i)) = c_i + parity(rho_i).

An element fixes some configurations iff pi is even and the sum of local
triangle signs around every pi-cycle is zero. In that case its fixed count
is `2·4^(number of pi-cycles)`. The exact fixed-count histogram is

| Number of group elements | Number of fixed configurations |
|---:|---:|
| 1053 | 0 |
| 216 | 8 |
| 27 | 128 |

Burnside gives `(216·8+27·128)/1296=4` isotypes. The action kernel consists
of independent rotations of the three leaf triangles, of order 27; the
fiber action is not faithful despite the faithful vertex action of G.

| Number of isotypes | Ordinary group | Labelled orbit size each |
|---:|---|---:|
| 2 | (C3)^3 semidirect C3, with cyclic factor permutation | 16 |
| 2 | (C3)^3 | 48 |

The first two types have all source placements or all sink placements; the
other two have one or two sink placements. Thus 32 labelled realizations
have automorphism order 81 and 96 have order 27, and

    2/81+2/27 = 128/1296 = 8/81.

`check_fiber.py` independently computes G by global hyperedge-preserving
permutation backtracking, then checks all 128 literal arc matrices, all
1296 element fixed counts, all 128 ordinary automorphism orders by a separate
arc-preserving backtracking, and orbit closure and mass. It also checks the all-height mass specialization
against the polynomial recurrence for h=0,...,4. Exhaustiveness of
the 128 choices follows from (1), not from closure alone. All guards remain
active under Python -O.

## 4. All-height upgrade: the earlier restricted family meets every isotype

Let H_0=C3(A) and recursively let H_(h+1) be the hypergraph obtained by
inflating the single-triple root with three copies of H_h. Let F_h be the
fixed-root-orientation family in [LIFTING.md](LIFTING.md). Although its
labelled representatives are only a subfamily, it meets **every ordinary
isotype in the full realization fiber of H_h**.

Proof by induction. At height zero, each realization is a singleton placed
before or after a triangle; either triangle orientation is isomorphic to the
fixed one. Thus the two isotypes are A_0 and A_1. At a higher root, a
realization chooses either orientation of a C3 quotient and arbitrary child
realizations. If its orientation is opposite the chosen convention, exchange
two whole blocks, possible because all child hypergraphs are copies of H_h.
This relabelling normalizes the quotient; it does not reverse child arcs.
Normalize each child recursively. The result lies in F_(h+1). Between two
normalized representatives the root quotient is the same directed C3, so
only cyclic permutations are allowed, not arbitrary S3 permutations. Hence
the cyclic-triple classification in the earlier recurrence stays exact.

Put L_h=3^h and P_h^nodes=(3^(h+1)−1)/2; these count linear seed nodes and
all prime nodes. The complete labelled fiber and hypergraph group orders are

    |X_Hh|=2^(L_h+P_h^nodes),         |Aut(H_h)|=6^(P_h^nodes).         (5)

The group recursion is G_0=S3 and G_(h+1)=G_h^3 semidirect S3. The independent
seed linear orders do not add hypergraph symmetries: their children differ
in size. This proves (5), not merely numerically matches it.

Avoiding notation collision, write p_h(z) for the earlier exponent
polynomial. Its recurrence now enumerates **the full fiber's isotypes**:

    p_0(z)=2,
    p_(h+1)(z)=(p_h(z)^3−p_h(z^3))/3 + z p_h(z^3).

For each class, |Aut(T)|=3^(L_h+e(T)), where e counts upper nodes with three
identical full child tournament types. In particular the full-fiber isotype
counts at h=0,1,2,3 are 2,4,24,4624. The complete possible order spectrum is
`3^L_h,...,3^(L_h+(3^h−1)/2)` with no gaps, by the earlier support proof.
The mass specialization is an independent exact identity:

    3^(−L_h) p_h(1/3) = 2^L_h / 3^(P_h^nodes).                       (6)

There is also a direct algebraic explanation of (6): at z=1/3 the two
p_h(z^3) terms cancel, leaving
`p_(h+1)(1/3)=p_h(1/3)^3/3`. Starting from 2 yields
`p_h(1/3)=2^(3^h)/3^((3^h−1)/2)`. Thus the stabilizer correction and orbit
mass agree by an exact cancellation, not an empirical numerical fit.

More precisely, if c_(h,e)=[z^e]p_h(z), then the number of **labelled**
realizations with ordinary automorphism order 3^(L_h+e) is

    c_(h,e) · 6^(P_h^nodes) / 3^(L_h+e).                              (7)

For a uniformly random labelled realization in this fixed fiber, its excess
symmetry exponent therefore has probability

    Pr(e)=c_(h,e)·3^(−e) / p_h(1/3).                                  (8)

Uniform isotype sampling instead gives c_(h,e)/p_h(1). These are distinct
sampling models: the labelled distribution weights a class inversely by its
automorphism order. The entire distribution is determined by the same
polynomial, not merely its attainable order interval.

The normalized family has only 2^L_h labelled constructions; (5) is usually
much larger. Full isotype coverage must not be confused with full labelled
coverage. Nor does this special C3-tree normalization assert that arbitrary
prime quotient orientations can always be normalized: non-self-converse
prime quotients are a direct obstruction.

## Reproduction

    python3 check_fiber.py > check_fiber.out
    python3 -O check_fiber.py > check_fiber.optimized.out
    cmp check_fiber.out check_fiber.optimized.out

Inheritance: [Session80 C3 families](../session80_tournament_symmetry_20261027/C3_FAMILIES.md), [group law](../session80_tournament_symmetry_20261027/GROUP_LAW.md), and
[LIFTING.md](LIFTING.md). Concept board: observer quotient / signed
transport / subgroup stabilizers / normal forms / orbit mass. The new bridge
is from fixed-H realization coordinates to exact ordinary symmetry
**distribution**, not a recovery of one unrecorded realization.

## Independent audit

A separate independent audit reports PASS for the theorem statements,
all-height normalization and exact-stabilizer conjugacy count. Its separate
`check.py` tests the actual C3 equations on all 8,192 parity-reduced arc
assignments at order 12: exactly the 128 product choices survive, independently
confirming exhaustiveness. It also checks all 1,296 affine coordinate actions,
80 two-generator subgroup fixed-count cases, and every literal fiber through
order five. Its normal/optimized receipts agree. This supplements rather than
replaces the cited all-order structural inputs.
