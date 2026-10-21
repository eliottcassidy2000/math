# Which symmetry classifications meet, and which operations do not commute

Session24, 2026-10-21. Proofs below use the independently audited signed-exception classification and current repaired metagraph canon. No literature-priority claim or complete-wiring theorem is made.

## 1. Exact intersections of three different symmetry conditions

Let E denote binary signed exceptionality: no vertex subset has trivial stabilizer in Aut±. Let R denote Hamilton-path rigidity h=H/|Aut|=1. Let P denote pure-blue status, equivalently the converse involution fixes every ordinary Hamilton-path orbit. These predicates act on different carriers.

Current inheritance, not historical THM2450 prose:
- THM-2477-rigid-strong-classification-through-sixty-two.md, current all-order addendum, and merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md classify strong path-rigid atoms as {1,C3}.
- THM-2454-pure-blue-alphabet-and-center-law.md classifies strong pure-blue atoms as {1,C3,R5}, where R5 is the regular cyclic five-vertex tournament.
- The new signed-exception proof shows every E tournament is strong, and classifies all of them as X[C3^[k]], X=C3,R5,P7.

Therefore, up to isomorphism and among tournaments of order at least two:

    E intersect R = {C3},
    E intersect P = {C3,R5}.

Proof: exceptionality forces strong connectivity; apply the two strong-atom classifications and check the listed seeds. The singleton is not exceptional under the faithful permutation action. This proves that the two symmetry notions are incomparable rather than universally aligned. P7 is exceptional but not pure-blue. Every transitive tournament of order at least two is pure-blue but not exceptional. Higher exceptional towers are never pure-blue.

## 2. A genuine cross-order map preserving converse-path structure

Let A be any nonempty tournament and T=A[C3]. The bottom copies are exactly the cyclic three-vertex modules: a triple using two vertices of one copy and one of another is not cyclic, while a cyclic triple using three distinct copies is distinguished by an unused vertex from one of those copies. Hence every automorphism or anti-automorphism preserves the bottom blocks. Thus Aut(T)=C3^|A| semidirect Aut(A), and T is self-converse exactly when A is. Write Hcal(A) for ordinary Hamilton-path orbits.

There is a canonical injection

    L: Hcal(A) -> Hcal(A[C3]).

Replace each vertex of a representative path by a Hamilton path through its triangle, keeping all vertices of each triangle consecutive. All three local triangle paths are equivalent under the local C3 action, so the resulting orbit does not depend on these choices. An automorphism taking two such lifted paths to one another induces an automorphism taking their base paths to one another; hence the map is injective.

For self-converse A, the injection commutes with converse followed by path reversal. This follows by projecting block order and observing that the local triangle has one ordinary Hamilton-path orbit, fixed by converse. Thus it embeds the complete involution action, not just its cardinality. In particular, writing h for ordinary path orbits and b for converse-fixed orbits,

    h(A[C3]) >= h(A),  b(A[C3]) >= b(A),
    (h(A[C3])-b(A[C3]))/2 >= (h(A)-b(A))/2.

The final inequality follows by injecting the two-element involution orbits, not by subtracting the first two inequalities.

For a non-self-converse merged pair, the fundamental black degree is h, so the ordinary injection gives degree monotonicity there too. Consequently triangle substitution injects ALL converse-merged node classes across orders n to3n, preserves self-converse versus paired status, and never decreases fundamental black degree. This gives cross-order substitution links between path-orbit carriers. It does NOT assert a homomorphism of the fixed-path-complement metagraph.

## 3. No signed exception is a black leaf

Use the fundamental paired-black multigraph convention of current `THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`, Statement 3: at a self-converse node degree=(h-b)/2, with loops contributing two. The following exact seeds were independently checked in the signed and path-cover packets:

    C3: degree0; R5: degree0; P7: degree3;
    C3[C3]: degree14; R5[C3]: degree80646.

The injection of Section2 makes these degrees nondecreasing under further triangle substitution. Together with the complete exceptional classification, it follows that every signed exception is either C3 or R5, isolated in the BLACK subgraph, or has black degree at least3. No signed exception is a black-degree-one leaf.

This is not a bridge classification. A bridge can be incident to vertices of degree greater than one, and neither this argument nor the degree data excludes that possibility. Blue edges are not included in the black-degree statement.

## 4. Exact obstruction to commuting with fixed-path complementation

Let kappa(T,P) reverse every arc not in the specified Hamilton path P, retaining its path arcs. Triangle substitution with a block-contiguous lifted path does not commute with kappa, even up to unlabelled isomorphism.

Take A=B=C3 and P=(0,1,2). Lift to A[B] with vertex order0,1,...,8, grouped into three consecutive triangles. Exact Hamilton counts are

    H(A[B]) =3159,
    H(kappa(A[B], lifted P)) =2983,
    H(kappa(A,P)[B]) =27,
    H(kappa(A,P)[kappa(B,P)]) =1.

The unequal counts certify nonisomorphism. These counts are independently checked by subset DP and all9! vertex permutations in check.py.

The mechanism is a boundary-port defect, not merely a numerical mismatch. Between consecutive lifted blocks, only ONE cross arc belongs to the lifted Hamilton path: the last vertex of the first block to the first vertex of the next. kappa preserves that arc but reverses every other cross arc. The cross-block relation is therefore no longer uniform whenever at least one of the two blocks is nonsingleton (equivalently their size product exceeds one). An external endpoint distinguishes the path endpoint inside any adjacent nonsingleton block from its other vertices. Hence original nonsingleton module blocks cease to be modules after this operation.

In contrast, applying kappa to the quotient and then substituting retains uniform cross-block relations. Even reversing the internal gadgets cannot repair the cross-boundary defect. Any future wiring transport must retain the actual boundary ports; an unmarked module quotient discards essential information.

## 5. Connection contract

Source: canonical triangle-substitution towers and ordinary Hamilton-path orbits.
Target: larger path-orbit carriers with their converse involutions.
Preserved: ordinary orbit identity, fixed versus paired converse type, and injective two-cycle count.
Lost: fixed-path-complement adjacency and unmarked block uniformity after a flip.
Required sidecar for adjacency work: actual entry/exit vertices at every block boundary.
Cheapest decisive hostile: C3[C3] above, with2983 versus27 (or1 if internal gadgets are flipped too).

The intersection and no-leaf corollaries are all-order consequences of the cited classifications. The port obstruction is elementary and the numerical witness is finite-exact. Complete metagraph adjacency and bridges not incident to leaves remain open here.
