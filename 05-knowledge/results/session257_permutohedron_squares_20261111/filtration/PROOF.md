# Stabilizer marks, fixed permutohedra, and the valuation wall

Status: PROVED elementary constructions below; CITED inherited Monster orders,
centralizer and Paley-normalizer inputs; FINITE-EXACT bounded controls in
`check.py`. No new sporadic-group classification or moonshine explanation is
claimed. This is a common functorial carrier, not an identification of Monster
with a twelve-vertex square-sum graph.

## 0. Inheritance and concept board

Sources: `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`;
Session32 `SYLOW_FIXED_POINT_PARTITION.md`, `VALUATION_CYCLE_AND_QUOTIENT.md`,
`FISCHER_SUPPORT_SHELLS.md`, `paley/PALEY_MONSTER_BRIDGE.md`; Session256
`SYNTHESIS.md`, especially its finite-response interfaces and the inherited
endpoint-symmetry obstruction. Their large-group inputs are inherited, not
reaudited. Session32 assembled copies supply the latter four files locally.

Closest mechanism: subgroup containment detects fixed points in a coset action.
Hostile: equal p-adic orders need not give equal fixed-point marks. Corrected
near miss: B is not the actual stabilizer; 2.B is. Least-used sidecar: the
entire orbit partition, rather than merely prime support or orbit count.
Concept board: stabilizer depth; marked singleton ports; ordered partitions;
valuation weights; signed reversal. Anchor: merge the three prime shells with
permutohedral geometry. Niche: exact normalization wall. Wildcard: port
responses usable by square-sum attachment and signed Paley presentations.

## 1. Fixed-face transport theorem

Let X be a finite set of n elements, let a finite group G act on X, and let
Pi(X) be the convex hull of all bijections X -> {1,...,n}, in R^X. For a
subgroup P <= G, let its orbits be O_1,...,O_k, of sizes m_1,...,m_k.

**Theorem 1.**

1. Faces of Pi(X) correspond to ordered partitions B_1|...|B_r of X. Such a
   face is P-invariant iff every B_i is a union of P-orbits.
2. The inclusion-minimal P-invariant faces have the individual O_i as their
   blocks in any order. There are k! of them, each combinatorially
   Pi(O_1) x ... x Pi(O_k), of dimension n-k.
3. The fixed section Pi(X)^P = Pi(X) intersect (R^X)^P is a weighted
   realization of a (k-1)-dimensional permutohedron. Its face poset is that
   of Pi({O_1,...,O_k}); it has k! vertices. For an orbit ordering
   O_(i1)|...|O_(ik), its vertex has constant coordinate

       (2s + m_(ij) + 1)/2 on O_(ij),   s=sum_(h<j) m_(ih).

4. The singleton P-orbits are exactly X^P. Thus an ordered list of r distinct
   fixed vertex-ports has (|X^P|)_r possibilities, the falling factorial.
   The singleton flags on the orbit blocks recover this port response from
   the fixed-face construction. With an invariant allowed port set E,
   replace X^P by E^P.

**Proof.** Maximizing a linear functional partitions coordinates by equal
coefficients. Rearrangement assigns consecutive ranks to each coefficient
block, proving the usual ordered-partition face description and dimension
sum(|B_i|-1). An action preserves this ordered partition exactly when it
preserves each block, equivalently each block is a union of orbits. Refining
an invariant partition until each block is one orbit proves (2).

Average any point over P. This Reynolds operator maps Pi(X) onto Pi(X)^P.
A functional on the fixed subspace has one coefficient per orbit. If these
coefficients are distinct, rearrangement puts the entire orbits into
successive rank intervals. Averaging within each orbit gives the displayed
vertex. Different orbit orders give distinct vertices. Ties merge precisely
the corresponding consecutive orbit blocks. Hence the normal fan is the
braid fan on the k orbit coefficients modulo constants, proving (3),
including the entire face lattice, not just a vertex count. Item (4) is the
definition of a fixed vertex and an ordered injection. QED.

**Important distinction.** A P-invariant face need not be fixed pointwise.
For transitive nontrivial P, the whole permutohedron is its unique minimal
invariant face, while its fixed section is a single point. The ambient
permutohedron has a P-fixed vertex iff P acts trivially on all X, not merely
iff X^P is nonempty. The marked singleton orbit is the correct port test.

**Naturality.** Equivariant bijections transport all these objects. Subgroup
inclusion P<=Q coarsens the orbit partition and gives the literal inclusion
Pi(X)^Q subset Pi(X)^P. A marked vertex x survives precisely for subgroups
P contained in G_x. This is the usable map from a subgroup filtration to
fixed polytopes plus ports. Abstract combinatorial type discards orbit sizes
and fixed vertices; restore the orbit partition and singleton marks.

### Exact identification with the square-sum zonotope

The weighted realization in Theorem1 is exactly the companion chart's
zonotope, not merely a polytope of the same combinatorial type. Write r_i
for the common coordinate on orbit O_i and m_i=|O_i|. Define the invertible
affine coordinate change on the orbit-constant subspace

    y_i = m_i (r_i-(m_i+1)/2).

At the fixed-section vertex indexed by an orbit order, r_i is the average
of the consecutive ranks L_i+1,...,L_i+m_i. Therefore y_i=m_i L_i. Taking
convex hulls and applying the weighted carrier theorem in
`../permutohedron/PROOF.md`, Section2, proves the literal affine equivalence

    Pi(X)^P  ≅  Z(m) = sum_(i<j) m_i m_j [e_i,e_j].

The inverse is r_i=y_i/m_i+(m_i+1)/2. Orbit weights are positive, so no
coordinate is lost by this transformation. This is the same coordinate
map as block-position sums: m_i r_i equals the sum of ranks in orbit i.
It is independent of any graph structure. A graph contributes admissible
ports, orientations and seams as additional data over this exact carrier.

**Averaging is the companion refinement projection.** Suppose a group
P's orbits refine those of a larger coordinate-permutation group Y. Group
fine orbits a into coarse blocks I_i, with weights u_a and
w_i=sum_(a in I_i)u_a. The Reynolds average R_Y sends Pi(X)^P onto Pi(X)^Y:
it maps the large permutohedron into itself and onto its Y-fixed section,
and every Y-fixed point was already P-fixed. In zonotope coordinates this
map is precisely

    (Lz)_i = sum_(a in I_i) z_a - sum_(a<b in I_i) u_a u_b.

Indeed the coarse constant rank is sum_(a in I_i)u_a r_a/w_i; substitute
the two affine coordinate formulas and use
w_i(w_i+1)/2-sum_a u_a(u_a+1)/2=sum_(a<b)u_a u_b.
Thus the companion cut-refinement projection and subgroup averaging are
the same affine map on their shared polytopes. One interpretation uses path
pieces; the other uses permutation orbits. This equality is exact even
when the underlying path seams are infeasible.

Do not reverse its direction: subgroup inclusion gives the inclusion
Pi(X)^Y subset Pi(X)^P, while averaging gives a surjection the other way.
These maps are a retraction pair in rank coordinates, not inverse
bijections. Fine vertices may project to interior coarse points.

## 2. Monster application: prime shells become port-depth regimes

Take X=M/C_M(t), t in class 2A, with H=C_M(t)=2.B and
n=97,239,461,142,009,186,000. These inputs are CITED from the inherited notes.
For any prime p let a_p=v_p(|M|), b_p=v_p(|H|). Then

    max {log_p |P| : P a p-subgroup, X^P nonempty} = b_p.

Indeed X^P nonempty iff P lies in some conjugate of H; a Sylow p-subgroup
of H attains the bound. Subgroups of each intermediate p-power order exist.
This gives a new geometric realization of the inherited three regimes:

- p in {17,19,23,31,47}: a=b=1. Every Sylow subgroup's orbit-labelled fixed
  permutohedron has at least one marked singleton block.
- p in {2,3,5,7,11,13}: 0<b<a. A p-subgroup at every order up to p^b can
  have a singleton port; none above that order can. In particular the full
  Sylow fixed section has no singleton mark, even though some order-p
  fixed section has one.
- p in {29,41,59,71}: b=0, a=1. Every nontrivial p-subgroup has no singleton
  marks; its orbits all have size p. Its fixed section has dimension n/p-1,
  and its minimal invariant faces are products of n/p copies of Pi_p.

Thus the Fi22 prime support identifies the intermediate port-depth regime;
the two outer shells are full-depth and zero-depth. This uses actual
subgroups of M and the actual coset action, not a fictional chain Fi22<B<M.
It does not turn the order coincidence with Fi22 into a new embedding.

For any order-p subgroup P, write f_P=|X^P|. The exact fixed-section rank is

    dim Pi(X)^P = f_P + (n-f_P)/p - 1.

The unmarked dimension plus known n,p recovers f_P; but the depth vector
alone does not determine this dimension. The missing datum is the subgroup
mark f_P, computable in principle by

    f_P = |{g in M : g^(-1) P g <= H}| / |H|.

For P=<u> cyclic, this is |C_M(u)| |u^M intersect H| / |H|.
These formulas specify an exact next computation rather than guessing
fixed-point counts from valuations.

**Hostile control.** In S4, let H1=< (12) > and H2=< (12)(34) >. Both coset
sets have size12 and the same order valuations a_2=3,b_2=1. For P=<(12)>,
the marks are respectively2 and0, since |C_S4((12))|=4. Thus the fixed
permutohedra have dimensions6 and5. No valuation-only rule can recover the
full fixed-face response, even in this tiny faithful transitive setting.

## 3. Q12: three colors, orbit faces, and endpoint ports

Here Q12 means the graph on 1,...,12 with an edge iff the sum is a square,
not a twelve-dimensional cube. Its path components, with explicit traversal,
are

    A=(2,7,9), B=(11,5,4,12), C=(8,1,3,6,10).

For this paragraph only, label the component classes A,B,C. These are not
the inherited three birth-role edge colors used in the companion note.
The full graph automorphism group is C2^3,
generated by independent path reversals; sizes are distinct so components
cannot be permuted by a graph automorphism. Its orbit sizes on the 12
vertices are 1,1,2,2,2,2,2. Therefore its fixed permutohedron section has
dimension6, with 7!=5040 vertices, while its minimal invariant ambient
faces have dimension5 and type Pi_2^5 (a five-cube). This is not a cube
claim about Q12 itself.

Now enlarge the coordinate-permutation group to the Young group
Y=S_A x S_B x S_C. This is a subgroup of S12 containing the graph
automorphism group P=C2^3; it is **not** a graph automorphism group, since
arbitrary permutations within a path destroy its edges. Its three orbits
are exactly the components of weights3,4,5. Hence

    Pi_12^Y ≅ Z(3,4,5),

with exactly the six vertices listed in the companion note. Meanwhile
Pi_12^P ≅ Z(1,1,2,2,2,2,2), with weights assigned to the actual named orbits.
In rank space the hexagon is literally a section inside the six-dimensional
P-fixed polytope; averaging within each component retracts that polytope
onto the hexagon. In weighted coordinates it is the explicit refinement
projection just proved, grouping the reversal-orbits by their components.
The grouping is essential; the sorted list of seven sizes alone does not
specify the map.

This produces an actual common polytope for the group-action and square-sum
charts. It simultaneously shows exactly what the larger symmetry loses:
the two singleton midpoint marks disappear since all Y-orbits have size
at least3, the endpoint set is not Y-invariant, and path orientation and
internal adjacency are not Y-invariant either. Reynolds averaging retains
component position sums, not endpoint identities or seam feasibility. Keep
the six named endpoints and three reference traversals outside the quotient.

Two fixed singleton vertices are the odd-path midpoints 7 and3. Neither is
an endpoint. For the allowed port set E of the six path endpoints, E^G is
empty. Thus an unmarked fixed vertex says too little for endpoint gluing:
port incidence is essential. Choosing one endpoint in each component gives
8 choices; C2^3 acts freely and transitively on them. The kernel is trivial,
so orienting all three paths completely removes these graph symmetries.

For each ordering of the three color blocks, the ambient Pi_12 has a face
Pi_3 x Pi_4 x Pi_5. Its vertices are words with the blocks contiguous. Of
these, exactly 2^3 words traverse every component as a path. Across the six
block orders there are 48 such words, carrying a free C2^3 action and six
ordinary orbits. This is an exact bridge between a parabolic permutohedral
face bank and the three-color path-cover interface. It does not declare
all 3!4!5! vertices of one face to be valid square-sum words.

The incidence data needed to glue the colored components are their endpoint
pairs and the actual available cross edges. A complete operation that adds
one allowed edge between successive chosen path endpoints turns one of the
48 words into a Hamilton path exactly when both required edges exist.
Preserving only color sizes, or only the number8 of endpoint selections,
loses this predicate. The endpoint-marked ordered partitions preserve it.
This is the same fixed-port interface as in Section2, with a small graph's
actual allowed port set replacing the Monster coset set. It supplies a
shared exact construction, not a canonical bijection between their vertices.

## 4. The valuation cycle has a genuine permutohedral wall nearby

On the ordered prime coordinates (2,3,5,7,11,13), the quotient |M|/|B|
has small-prime valuation vector

    v=(5,7,3,4,1,2).

Replacing B by the actual stabilizer 2.B gives

    w=v-e_2=(4,7,3,4,1,2).

There is an exact geometric meaning to this normalization. Form the orbit
polytope P(v)=conv{sigma v : sigma in S6}; similarly P(w). P(v) has 720
vertices and the braid normal fan. P(w) has 360 vertices because equal
entries at the original prime coordinates2 and7 give a stabilizer S2.
The continuous path v(t)=v-t e_2, 0<=t<=1, stays in one strict weight
chamber until it hits exactly the wall v_2=v_7 at t=1.

**Theorem 2 (normalization contraction).** Exactly360 disjoint edges of
P(v(t)) collapse at t=1: the pairs obtained by exchanging the entries
originally attached to primes2 and7. They are edges because their values
5-t and4 are consecutive in the sorted list
1<2<3<4<5-t<7. Their length is sqrt(2)(1-t). No two other distinct vertices
coalesce. Equivalently, the vertex specialization is the quotient of S6 by
the order2 stabilizer of w, with fibers of size2.

**Proof.** Rearrangement gives the braid normal fan whenever all entries are
distinct. Exchanging neighboring values is an edge. At t=1 the sole repeated
value is4, so a vertex equality occurs iff permutations differ by exchanging
these two entries. The displacement has two nonzero coordinates of opposite
size1-t. QED.

This is a concrete operation on a permutohedron forced by the central
factor2, and it identifies exactly which prime pair collides: {2,7}.
The polytope alone forgets which coordinates are primes and how the
normalization arose. Keep the prime labels, the valuation vector, and the
chosen central-factor direction as sidecars. It is an S6 orbit polytope,
not a claimed Weyl group of the Monster.

The user's chain 13->2->5->3->7->4 becomes a five-cycle after the specified
set identification4=13. Its permutation matrix on five states has
characteristic polynomial z^5-1. Under the same identification, normalization
changes only2->5 into2->13: there is now a two-cycle13<->2 with tail
5->3->7->13. Its functional-graph matrix has characteristic polynomial
z^3(z^2-1). The same one-coordinate normalization therefore causes both a
permutohedral wall contraction and loss of invertibility of the auxiliary
cycle operator. These are two exact consequences of the same input change;
they are not equivalent invariants. The noninjective second operator cannot
be an automorphism of any finite permutohedron vertex set carrying these
five states faithfully. One must allow a noninvertible state map or retain
a larger history state.

This preserves the legitimate chain interpretation without inventing a
mod9 map. The additional arrow11->1 obstructs such a full congruence map,
and even the five-state restriction cannot be represented by an integer
polynomial modulo9, as Session32 already proved. The wall construction does
not use that invalid arithmetic quotient.

## 5. Orientation/converse makes the other inherited bridge compatible

The CITED normalizer quotient N_M(<7A>)/He is AGL(1,7), whose square-slope
subgroup is Aut(Paley7)=7:3. A nonsquare affine map is an anti-automorphism,
so on directed Hamilton words it acts only together with position reversal.

On Pi_7 vertices regarded as rank assignments r:X->{1,...,7}, position
reversal is the affine involution r -> 8-r. Hence a signed permutation
(g,epsilon) acts by coordinate permutation followed, when epsilon=-1, by
this central reflection. It acts on the entire permutohedron and preserves
the Hamilton-word subset exactly with the sign correction. Ordinary
coordinate permutation alone would send a directed word to a reverse-directed
word. Thus the actual Monster-local signed quotient embeds into a genuine
permutohedron symmetry action preserving the marked Paley Hamilton subset.
The square-sum component reversals instead preserve undirected traversal
without this sign correction. The sign is data, not a free identification.

## 6. Consequences and strongest next questions

1. **Compute marks, not new prime coincidences.** For representatives of
   order-p subgroups in the small-prime shell, compute f_P from class fusion
   into2.B. This would determine actual fixed-permutohedron dimensions and
   all ordered fixed-port counts. Orders alone provably cannot do it.
2. **Port-restricted marks.** Given a distinguished attachment relation on a
   G-set, replace X^P by fixed allowed ports and then by fixed incident port
   pairs. Which finite collection of such marks determines every bounded
   attachment response? Session256's finite-response viewpoint supplies the
   right target, but no universal finite interface independent of arity is
   claimed here.
3. **Filtration beyond support.** Can the level1/2/6 modular correspondences
   recover any of these subgroup marks, rather than only prime support?
   This is stronger and more testable than asking whether two prime lists
   match. No such map has been constructed in this note.
4. **Normalization family.** For other centralizer covers, subtract the
   valuation vector of the central factor and track the walls hit. Do the
   colliding prime coordinates correlate with subgroup-mark changes? This
   is an open comparison, not implied by Theorem2.
5. **Signed port responses.** Can the three fixed and three paired ordinary
   Paley Hamilton presentations be recovered by a marked fixed-face
   inventory for the signed action, retaining endpoint orientation? This
   would connect the inherited local quotient to actual response data,
   rather than merely repeating its group order.

Reproduce: `python3 -B check.py > exact.json` in this directory. The finite
universe is S4 cosets, all Q12 components and their48 traversal words, all
720 valuation permutations, explicit small orbit-section vertices, and all
42 signed affine maps on the189 Paley7 Hamilton words. An unsigned nonsquare
map has zero overlap with those directed words, a decisive sign control.
The large-group statements are applications of proved generic formulas to
inherited input, not giant permutation enumerations. No LRC or Chan–Pak
coverage consequence is asserted.
