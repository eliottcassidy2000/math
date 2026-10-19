# Saturated core vertices annihilate boundary-tree increments

Date: 2026-10-19. Status: proved structural lemma. The subsequent finite c=4 exploration is separate and proves no c=4 exclusion theorem.

## 1. Setup and forced exterior directions

Let x∥y in a finite poset P, and let Q=I(x)∪I(y), including x and y themselves. Every v∈Q is incomparable with at least one root. Every exterior element is comparable with both roots and hence lies either in

    A={a outside Q: a<x and a<y},
    B={b outside Q: x<b and y<b}.

The mixed alternatives would compare x and y. Every a∈A lies below every b∈B.

Suppose every vertex of P has incomparability degree at most d. Write d_Q(v) for v's incomparability degree in the induced core Q, and let

    S={v∈Q: d_Q(v)=d}.

**Saturation lemma.** For each s∈S and every a∈A,b∈B,

    a<s<b.

**Proof.** The d incomparable core neighbors of s exhaust its entire permitted incomparability degree, so s is comparable with every exterior element. If s<a, then s<a<x and s<a<y, contradicting that s is incomparable with at least one root. Hence a<s. If b<s, then x<b<s and y<b<s, giving the same contradiction. Hence s<b. This argument includes either root when it belongs to S. ∎

More generally, a core element cannot lie below any element of A or above any element of B, whether or not it is saturated. Its allowed comparison directions with the exterior are A<core<B.

## 2. Prefix increments: every later increment also vanishes

For an initial core word u that can precede the first root in a core extension, let L(u) be the number of interleavings of A with the ordered word u respecting P, divided by e(A). Let L(empty)=1. Define

    U_u=L(u)−L(u with its last letter removed),
    U_empty=1.

Suppose u contains a saturated letter s. Write u=q s r, where q is the part strictly before that occurrence. All A elements precede s. Therefore every allowed interleaving consists of an allowed interleaving of A with q, followed by the fixed word s r. Deletion and appending of s r are inverse bijections. Appending is valid because no core letter can lie below an A element; the core word already respects all internal relations. Consequently

    L(q s r)=L(q).

It follows that U_u=0 whenever u contains a saturated letter. In particular, along a fixed full pre-root word, the increment at the first saturated letter and **all subsequent prefix increments** vanish, not merely the increment at that letter.

## 3. Suffix increments: the dual statement

For a terminal core word v following the last root, define the normalized interleaving count R(v) with B, and

    V_v=R(v)−R(v with its first letter removed),
    V_empty=1.

If v=q s r contains a saturated letter, every B element follows s. Every valid interleaving is the fixed prefix q s followed by an allowed interleaving of r with B. Deletion and prepending are inverse bijections, giving

    R(q s r)=R(r).

Thus V_v=0 whenever v contains a saturated letter. Along suffix-tree growth by prepending letters, the increment when the last saturated letter is first included and every later growth increment vanish. Equivalently, in the original word, increments beginning at or before its last saturated letter are zero.

## 4. Consequence for the exact generator cone

The ambient extension-count vector has the boundary-tree expansion

    sum_(u,v) U_u V_v G_(u,v),

where G_(u,v) counts core extensions beginning with u and ending with v. The saturation lemma gives the exact vanishing rule

    U_u V_v=0 whenever (letters(u) union letters(v)) intersects S.

Therefore a proof for degree≤d posets need check only generators whose initial and terminal words avoid S. This removal is sound for every actual exterior; no independence or realizability assumption is made about the remaining nonnegative coefficients.

The retained cone can still be larger than the realizable set. A nonsaturated core vertex with degree d_Q(v)<d has only d−d_Q(v) available exterior incomparable neighbors. Those additional finite slack budgets are not enforced merely by deleting saturated-generator columns. A counterexample to this pruned cone would therefore be a limitation of that relaxation, not a counterexample poset.

## 5. Application to mixed (1,2), c=4

Here |Q|=9. The roots have core degrees 6 and 7. The lower slave already has its comparability with y, and each upper slave with x, so each already has degree≤7. The only way to violate the seven-thin core budget is a common neighbor with no comparability at all inside Q. Thus every common neighbor must have at least one core comparable neighbor.

A common vertex with exactly one core comparable neighbor is saturated and must be avoided in every nonzero prefix/suffix increment word. The lower slave is also saturated if y is its only core comparable neighbor. An upper slave is saturated if x is its only core comparable neighbor. Root saturation causes no extra boundary restriction because roots never occur in the pre-root or post-root words.

If every nonroot boundary-capable letter is saturated, all nonempty boundary increments vanish, leaving only the empty/empty generator. Then every ambient comparison probability within Q equals the bare-core probability. For example, this occurs when the common poset is two disjoint two-element chains, there are no cross-relations to slaves, and the two upper slaves are incomparable. In that example A lies below all of Q and B above all of Q.

## 6. Separate census observation

The canonical naturally labeled c=4 family has 96,428 cores; 79,121 survive all-vertex core-degree≤7, eliminating 17,307 (approximately 17.95%). This count does not include the generator pruning above and is not itself a proof of a c=4 majority exclusion.

The excluded count also has an inclusion-exclusion check using the already verified smaller family counts F(0)=7,F(1)=40,F(2)=357,F(3)=4824. Requiring a specified common vertex to be isolated leaves exactly a smaller family. Hence

    4F(3)−6F(2)+4F(1)−F(0)=17,307.

The numbers with exactly 1,2,3,4 over-budget common vertices are 15,464; 1,704; 132; 7 respectively. The proof of the zero-increment lemma does not depend on these computations.
