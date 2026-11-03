# A finite upper certificate for the order-ten response dimension

This is an additive certificate to the Session 146 order-ten lower-bound packet. It uses the same rational linear response language: all labelled old tournament laws and all fixed source-oblivious attachments. It does not assert a bound for source-dependent attachments, conditioning, shared-source products, or a restricted S/B-only word bank.

**Theorem.** The order-ten pooled-parity carrier of dimension 347626 has at least 2520 independent invisible directions. Together with the separately certified lower bound 345106, the exact response dimension is **345106** (or 345105 varying coordinates for normalized laws). The upper certificate has been independently audited by a second diagram enumerator; its independent receipt and controls accompany this packet.

## Reduction to a single 735-column block

Fix the full-support old signatures with degree-one set E={0,...,7} and degree-two set {8,9}. They have four paths, six old forest edges, and 735 pooled-parity signatures. For an exterior coefficient context G, if G meets either old degree-two vertex, every coefficient in this block vanishes. The same holds if any vertex in E has G-degree at least two. Otherwise each vertex in E is untouched by G or has G-degree one.

Suppose v in E is untouched. Add a new path v−x−y with fresh labels x<y larger than all previously used labels. On every old forest in this block, the augmented union adds two vertices and two edges to an existing path. Linear-forest validity, all path parities, the component count and nontrivial-component count are unchanged. The new two-edge tail has traversal sign + in either direction; on an even path reversing canonical orientation also has sign +. Thus the entire restricted coefficient row is unchanged. Repeating gives a saturated context in which each vertex of E has G-degree exactly one. This applies also to zero rows.

## Finite saturated diagrams

A saturated G with a nonzero coefficient row is a disjoint union of paths. Since every old vertex in E has G-degree one, each component of G that meets E is either (i) an exterior chain joining two distinct vertices of E, or (ii) an exterior chain from one vertex of E to a new leaf. Thus the incidence data are a partial matching of E, with a stub at every unmatched endpoint. Record the parity of the number of edges on each joining chain and each stub. Give a joining chain canonical length 2 for even parity or 3 for odd; give a stub canonical length 2 for even or 1 for odd. All internal vertices in these chains are new and distinct. There are exactly 32400 such matching/stub/parity diagrams.

Additional G-only cycles or odd-length paths make the entire row zero. G-only even paths and isolated new vertices contribute only a nonzero rational factor, independent of the old signature, and may be omitted for row-span purposes. Replacing any chain by the canonical chain of the same parity preserves the contracted incidence graph, its cycles, all enlarged path parities and the sign changes upon reversing a chain. Its intrinsic traversal sign changes by a fixed factor independent of the old forest, so the restricted row is rescaled by a rational scalar. The number of resulting components and nontrivial components may change due to discarded G-only components but, for a fixed G, they are independent of the old forest whenever the union is admissible. Therefore every saturated row is a scalar multiple of one of these 32400 diagram rows, or zero.

More explicitly, every old representative has six edges, so for a fixed exterior G, any admissible union has c=N−|G|−6 components. All old vertices are non-isolated; the number of isolated new vertices is fixed by G. Hence the number of nontrivial components is also independent of the representative. The Hamilton coefficient law gives the same positive magnitude in all nonzero entries of a fixed diagram row, with signs ±1. Consequently a sign-normalized 0/±1 row is an exact rational row representative, not a modular approximation. For these finite canonical diagrams all N<101, and the common magnitude is nonzero modulo 101; normalizing modulo 101 therefore recovers the exact sign row.

## Exact integer kernel certificate

Concatenate `kernel.part00`, then `kernel.part01`, to obtain a 735-by-56 integer matrix U, with SHA-256 `91c68094d7c5d74467a9de29e949b6bbabe844edd0e2c0193bca9d5b241522fb`. Rows are the 735 canonical parity-signature representatives in the same order as the order-ten r=4 lower certificate. `verify_upper.cpp` reconstructs these representatives, checks rank(U)=56 modulo 101, generates all 32400 canonical diagrams, and checks that every exact 0/±1 sign row times U is the integer zero vector. It finds 17042 all-zero diagram rows. Thus the rational rank of the entire block is at most 735−56=679. The independent lower certificate already gives 679 independent rows, so its exact block rank is 679.

Relabel E in all C(10,8)=45 ways. The 56-dimensional kernels are supported on disjoint signature sets and hence form 2520 independent kernel vectors of the full response matrix. A signed relabelling only permutes/scales coordinates. Thus the full response dimension is at most 347626−2520=345106. The previous lower bound gives equality.

This is a linear kernel statement, not a claim that every kernel vector is a probability law; sufficiently small signed perturbations around the uniform law can be normalized to strictly positive indistinguishable source laws. The proof uses all fixed attachment responses through exact exterior-variable interpolation, not a finite list of actual S/B words.

## Blindness persists at every higher order

For every n ≥ 10,

    d(n) ≤ c(n) − 2520 C(n,10),

where c(n) is the pooled-parity carrier count. To see this, select any ten-element subset A of the n old labels and embed the 2520 order-ten kernel vectors on signatures whose non-isolated support is exactly A, with all other old labels isolated. For any n-order exterior context, restrict its coefficient row to these embedded columns. Regard the n−10 old vertices outside A as additional exterior vertices, while leaving every old-old edge incident with them unselected. The restricted row is an allowed order-ten exterior coefficient row, up to the fixed sign gauges arising from relabelling. It therefore annihilates every embedded kernel vector. Different ten-subsets have disjoint non-isolated signature supports, so all 2520 C(n,10) embedded vectors are linearly independent. This proves the upper bound and shows the pooled-parity carrier is not minimal at any order n ≥ 10.

Combining with the order-ten lower-bound packet gives

    L9(n) + 45405 C(n,10) ≤ d(n) ≤ c(n) − 2520 C(n,10),

where L9 is defined there. For example 1654181 ≤ d(11) ≤ 1874511 and 6796879 ≤ d(12) ≤ 10611580. These include the constant coordinate; subtract one for normalized laws. Equality at n ≥ 11 remains open.
