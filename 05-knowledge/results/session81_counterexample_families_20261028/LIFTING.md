# A procedural collision factory: full subset path data versus symmetry

**Status:** PROVED deductions with cited structural inputs; FINITE-EXACT checks.
Ordinary directed tournaments on a fixed finite labelled set; isomorphisms
preserve arcs, not converse. No literature-novelty claim. The order-12 witness
is inherited, not new. The contribution is its exact all-height family law,
gapless symmetry-order range, and identification of the observer quotient.

## Inheritance / concept board

Anchor: generate counterexamples, not merely repeat a seed. Niche: identify the
precise path observer's information. Wildcard: obtain every intermediate group
order inside one observer fibre. Board: module reversal / labelled subset
paths / canonical child types / cyclic color stabilizers / endpoint markings.

Closest mechanism: Session80 `GROUP_LAW.md`, decorated canonical stabilizers.
Hostile seed: its `C3[A,A,A]` versus `C3[A*,A,A]`, orders 81 and 27.
Corrected near miss: equal C3, child group, or cycle index never guarantees
same parent symmetry. Least-used sidecar: THM-3121,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`,
retains full path covers when scalar H does not compose. Here the entire
*label-addressed* path-cover profile still forgets tournament symmetries.

## 1. Run-reversal lemma (elementary)

A module M of T is a vertex set for which every outside vertex beats all of M
or loses to all of M. Let T^M reverse precisely the arcs internal to M.
For each directed path, reverse the vertex order of each maximal consecutive
run lying in M. Internal run arcs are now valid in T^M; all junctions remain
valid by the module property. Runs and their membership pattern do not change,
so this is an involutive bijection.

It preserves the exact vertex set of each path, its length, the ordered word
of inside/outside memberships, every outside vertex's position, and the
number and vertex sets of path-cover components. Thus for every S⊆V and d,

    H(T[S]) = H(T^M[S]),       pc(T[S],d) = pc(T^M[S],d).

Here H counts directed Hamiltonian vertex sequences, H(empty)=1; pc counts
unordered spanning collections of d nonempty vertex-disjoint directed paths.
Both endpoints are preserved if outside M; an endpoint inside M generally
changes. Inclusion/deletion marks are safe, orientation-sensitive endpoint
marks are not. This is not a bijection preserving the full individual path.

A nonmodule has no such guarantee: reversing one arc of C3 produces a
transitive tournament, changing H from 3 to 1. At n=2, a reversed arc already
changes its labelled endpoint path matrix despite preserving H and C3.

## 2. Exact observer equivalence (all finite orders)

For tournaments T,U on the **same** finite labelled V, these are equivalent:

1. C3(T)=C3(U), as actual sets of vertex triples.
2. U is obtained from T by finitely many reversals of sets that are modules
   **at the step at which they are reversed**.
3. H(T[S])=H(U[S]) for every S⊆V.
4. pc(T[S],d)=pc(U[S],d) for every S⊆V and d.

The structural input is [Boussaïri–Chergui–Ille–Zaidi (2018)](https://arxiv.org/html/1805.04607v1),
Theorems 10, 16, 24 and Proposition 51: common strong-module tree; prime or
linear tournament quotients; prime realization unique up to converse; and
independent realization choices at nodes. This uses their Definition 1 of
hypergraph module, not the alternative in Remark 26. Checked in the primary
text 2026-10-28; inherited route: Session80 `C3_FAMILIES.md`.

**Proof of 1⇒2.** Work on the common strong-module tree. At a prime node,
reverse its entire vertex set, then reverse each child set. Internal child
arcs are restored while all quotient arcs are reversed. At a linear node,
two adjacent child blocks have module union. Reverse that union and reverse
each of the two children; this swaps their order without changing either
child internally. Adjacent swaps realize any linear ordering. All named sets
remain modules throughout each move. Applying these operations at the nodes
matches every target quotient and therefore all target arcs. This proves the
claim also for degenerate n=0,1,2, with vacuous or linear operations.

**Other implications.** The lemma gives 2⇒4, and d=1 gives 4⇒3.
A triple has H=3 exactly when cyclic, otherwise H=1, so 3⇒1.
Alternatively 3⇒4 follows without decomposition: for nonempty S,

    pc(T[S],d) = sum_{π partition of S into d blocks} ∏_{B∈π} H(T[B]).

Thus the exponentially large labelled subset path observer contains *exactly*
the cyclic-triple hypergraph's information, no more. Equal total number of
cyclic triples is not sufficient. Nor is equality of an unaligned multiset of
subset H values the assertion proved here. After relabelling, the same iff
holds only with one common relabelling aligning every subset.

## 3. Ternary collision factory and exact symmetry polynomial

Put A0=1→C3 on labels 0,1,2,3, with 0 the source and
1→2→3→1; put A1=A0* (literal converse on the same labels).
Both groups are C3, but their score multisets differ. Define F_0={A0,A1} and

    F_(h+1) = { C3[B0,B1,B2] : B0,B1,B2∈F_h }, modulo ordinary isomorphism.

The underlying labelled construction has 3^h addressed leaves, each a
four-vertex block with an independent bit, so 2^(3^h) labelled choices at
order n_h=4·3^h. All upper quotients have the fixed orientation C3. Any pair
of choices differs by reversing leaf modules. Consequently all choices have
identical labelled C3 and every labelled subset path/path-cover count.

These are genuine canonical prime-root substitutions at every upper node.
Thus isomorphism classes correspond exactly to cyclic triples of child
isotypes: reflection is forbidden. Let e(T) be the number of upper nodes
whose three children have the *same full tournament isotype*. Then

    |Aut(T)| = 3^(3^h + e(T)).

This is a node-count formula, not a claim that the group is elementary
abelian. Its full group is recursively `(G0×G1×G2)⋊H`, H=C3 for three equal
isotypes and H=1 otherwise, with actual cyclic permutation action.

Define P_h(z)=Σ_[T∈F_h] z^e(T). The exact procedure is

    P_0(z)=2,
    P_(h+1)(z) = (P_h(z)^3-P_h(z^3))/3 + z P_h(z^3).          (★)

Explanation: remove constant triples from ordered triples, divide the
remaining free rotation orbits by three, and put each constant triple back
with one extra symmetry node. Actual child isotype equality, not equality of
e or group order, determines which triples are constant.

In particular a_h=|F_h| satisfies a_0=2 and
`a_(h+1)=(a_h^3+2a_h)/3`; a_0,a_1,a_2,a_3=2,4,24,4624.
The construction does not list every **labelled** realization of its C3
hypergraph. The subsequent [full-fiber normalization theorem](FIBER_ACTION.md)
proves that it nevertheless meets **every ordinary isotype in this particular
fiber**. These are not counts of all tournaments or all prime tournaments.

### Every intermediate symmetry order occurs

For every h≥0 the support of P_h is exactly

    {0,1,...,M_h},       M_h=(3^h−1)/2.

Hence one labelled observer fibre contains tournaments with **every** order

    3^(3^h), 3^(3^h+1), ..., 3^((3^(h+1)−1)/2).

Proof: exchanging all leaf bits preserves e and pairs distinct isotypes,
since the odd number of leaves prevents equal counts of A0 and A1.
Thus each nonzero coefficient is at least two. If the support at h is
0..M, any sum in 0..3M can be realized by a nonconstant triple: choose child
exponents summing to it, and if all chosen types would coincide, choose the
other type of the same exponent in one slot. Such a triple has no upper
rotation. A constant maximum triple realizes 3M+1. No greater exponent is
possible because there are only M_(h+1)=3M+1 upper nodes.

The uniform construction realizes the maximum. A concrete minimum pair is
D_0=A0,E_0=A1 and

    D_(h+1)=C3[D_h,D_h,E_h],
    E_(h+1)=C3[E_h,E_h,D_h].

Their numbers of A1 leaves differ (indeed by one), so they are distinct
isotypes at every level. All upper stabilizers are trivial, giving
Aut(D_h)≅Aut(E_h)≅C3^(3^h). The maximum/minimum order ratio is 3^M_h.
For h=1 this recovers the old 81-vs-27 witness; h≥2 gives the all-height
amplification and the whole intervening order spectrum.

## 4. Loss ledger and decisive tests

Source: labelled tournaments. Map: full S↦H(T[S]) or full (S,d)↦pc(T[S],d).
Preserved: exactly C3 realization class, including its common strong-module
vertex sets. Destroyed: quotient orientation choices, SCC block order,
ordinary child isotypes and consequent stabilizers. Needed sidecar for the
exact group: quotient arcs, child isotypes, and coherent transports as in
Session80 `GROUP_LAW.md`; H alone does not repair these.

For example the four-vertex seeds have H=3, but end-vertex count vectors
(0,1,1,1) and (3,0,0,0). No theorem here preserves prescribed inside endpoints,
individual directed paths, ordinary isomorphism type or automorphism group.

## 5. Reproduction / exact universe

    python3 check_lifting.py > check_lifting.out
    python3 -O check_lifting.py > check_lifting.optimized.out
    cmp check_lifting.out check_lifting.optimized.out

Standard library only; no assert-based validity controls. The checker:

* exhausts every labelled tournament of orders 1..6 (33,867 total), grouping
  by literal C3 sets and independently calculating all subset H by endpoint DP;
* verifies the direct run map on every nonempty module of every labelled
  tournament through order four and every directed nonempty partial path;
* independently computes full automorphism orders by global arc-preserving
  backtracking for all eight depth-one leaf colorings (order 12), and compares
  all 4,096 subset H values and actual cyclic-triple sets;
* enumerates every binary leaf word through h=2 (2+8+512 choices), canonicalizes
  by cyclic child rotation, and checks the class exponent polynomial (★);
* calculates (★) through h=3, checks gapless support, and rejects endpoint and
  nonmodule overclaims with explicit hostile controls.

The all-order theorem is proved above, not inferred from these finite ranges.
No full permutation group is enumerated at order 36 or beyond. No assertion
of minimal order for unequal symmetries is made.
