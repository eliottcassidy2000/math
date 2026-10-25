# Classification of all bare unique-sum frames for an affine square plane

## Scope and exact results

Let F=F_(3^n), q=|F|, and let A=C+U⊆F be a fixed nine-point affine F3-plane. A frame is an unordered pair {B,D} of three-element subsets of F such that the nine sums b+d are distinct and constitute A. The frame is disjoint when B∩D is empty. If A consists of square values, every disjoint frame supplies a square-sum K3,3 cross skeleton. The graph need not be induced, and vertex labels need not themselves be squares. Square values include zero unless A is explicitly required to avoid it.

There are exactly:

- 30q frames before requiring disjointness: 6q AP/AP and 24q AP/non-AP.
- 30(q−9) disjoint frames: 6(q−9) AP/AP and 24(q−9) AP/non-AP.
- 30 free gauge orbits, of which six are AP/AP and 24 mixed; every orbit has q frames, nine bad and q−9 disjoint.
- 1080(q−9) distinct unoriented simple six-vertex paths belonging to these cross skeletons.
- 1080(q−9) distinct unordered pairs of unoriented six-vertex Hamilton paths whose union is one of these skeletons.

Exactly one fifth of the pre-disjointness frames have both parts AP; the remainder have exactly one AP part. The same fractions hold among disjoint frames when q>9; at q=9 all disjoint counts are zero. These assertions count only frames with the prescribed nine distinct cross sums A, not arbitrary square-sum bicliques. Every frame has the same magic-admitting entry set A, including the mixed frames. The distinction is whether a particular additive chart has both AP guards, not whether its entry set admits some magic arrangement.

## 1. Every frame has an AP part

Fix d∈D. Because B+d⊆A, every difference of two elements of B lies in U. Likewise D−D⊆U. Uniqueness of the nine cross sums is equivalent to

    (B−B)∩(D−D)={0}.

Indeed, a nonzero shared difference b−b'=d−d' gives the repeated sum b+d'=b'+d, and conversely a repeated cross sum gives such a difference.

A three-element set whose differences lie in U is either an affine F3-line, or, after translation, {0,u,v} with u,v independent. In the second case its nonzero differences are precisely ±u, ±v, ±(u−v), six distinct elements. The two remaining nonzero elements of U form the single direction ±(u+v). Hence if B is non-AP, all differences of D lie in that remaining one-dimensional subspace. D has three elements, so it is an entire affine line. This proves that at least one part is AP and rules out a non-AP/non-AP frame.

## 2. Count with a distinguished AP part

There are four directions L⊂U. For each, B may be any of the q/3 affine cosets b+L in F. The possible D are exactly the transversals of the three L-cosets in the affine plane C−b+U: one point from each coset. Necessity follows because sums must stay in A and two elements of D in the same L-coset produce overlapping translates B+d. Sufficiency follows because the three translates then partition A. There are 3³=27 transversals.

Thus 4(q/3)27=36q ordered frames with a distinguished first AP part. Among the 27 transversals, exactly nine are affine lines: choose one of the other three directions M⊂U and one of its three affine cosets in C−b+U. Hence 12q of these ordered frames are AP/AP, and 24q are AP/non-AP.

When forgetting the order of the parts, AP/AP frames have been counted twice, while mixed frames have been counted once because their AP part is unique. The unordered total is therefore 12q/2+24q=30q, split as 6q+24q.

The choice of representative b for B does not affect C−b+U, because changes of b lie in L⊂U. Likewise replacing C by any other representative of A changes nothing. This count contains no hidden base-point multiplicity.

## 3. Gauge orbits and disjointness

First fix a gauge-invariant orientation: put the AP part first for mixed frames; for AP/AP frames put their distinct directions in a fixed order. With this convention the additive group of F acts by

    {B,D} → {B+g,D−g}.

Without this convention the displayed operation would depend on the chosen ordering, although its full orbit set would not: exchanging the initial ordering replaces g by −g.

For ordered parts, a stabilizer g requires g to be a period of both B and D. Any nonzero such period belongs to both difference sets, which uniqueness forbids. For unordered parts there is one additional possibility: the action might exchange B and D. That would make D a translate of B, so their difference sets would be equal and nontrivial, again impossible. Therefore this action is free even on unordered frames. AP status is preserved, so there are 30 orbits, split into six AP/AP and 24 mixed orbits.

A cross-part collision after gauging is equivalent to

    b+g=d−g, or 2g=d−b.

The nine differences d−b are distinct. Equality d−b=d'−b' would give the equal cross sums d+b'=d'+b, forcing b=b' and d=d' by the unique-sum hypothesis. Multiplication by 2 is bijective in characteristic three. Thus there are exactly nine bad gauges in each free orbit, giving 30(q−9) disjoint frames and the asserted split.

## 4. Distinct paths, not merely incidence counts

A K3,3 on six fixed labeled vertices has 3!·3!=36 Hamilton paths modulo reversal: start on a designated part and permute each part along the alternating positions. Every unoriented six-vertex path determines its unordered alternating bipartition, by taking its odd and even positions. Consequently two distinct frames cannot contribute the same path. Additional edges in the ambient square-sum graph do not affect this uniqueness.

It follows that the frame/path incidence count 36·30(q−9)=1080(q−9) is also a count of distinct ambient simple six-vertex paths in this precisely described subclass. It is not necessarily the total number of six-vertex paths in the entire square-sum graph. The AP/AP subclass similarly has 216(q−9) distinct paths.

Exactly 36 unordered pairs of the 36 unoriented Hamilton paths cover all nine K3,3 edges; this finite claim is independently enumerated in the verifier. A covering pair determines its union graph, and the connected K3,3 graph determines its unordered bipartition. Thus covering pairs from distinct frames also cannot coincide. This yields the same 1080(q−9) exact count for covering pairs. None of these six-vertex paths is being called a Hamilton path of the entire q-vertex field graph.

## 5. Corrections and relationship to the frozen AP-only packet

The preceding AP-only packet's frame counts and 72-fold labeling quotient are unchanged. Its warning against treating the path count as a deduplicated ambient count was unnecessarily weak, and any suggestion that distinct frames could repeat the same six-vertex path is incorrect: the alternating bipartition recovers the frame. The corrected statement is exact distinctness within the counted subclass, without exhaustiveness in the ambient graph.

Both the pinned F81 and F729 planes contain zero, so their square certificates use the zero-allowed convention. The F2187 plane avoids zero and therefore works for both square conventions. This supplement does not alter the frozen earlier packet.

## 6. Reproduction and finite witnesses

The standard-library verifier uses independent base-three polynomial field arithmetic and the pinned witness planes. It enumerates all AP parts and all 27 transversals per AP part, deduplicates unordered frames, and checks the exact type and disjointness counts for q=81,729,2187.

For F81 it independently examines all binomial(81,3)=85,320 possible first triples, without imposing an AP or direction restriction. It forms every possible second triple whose cross sums lie in A and separately tests nine distinct sums. This finds exactly the same complete frame set. It then explicitly enumerates the 30 gauge orbits and verifies 81 members and nine bad members in each, and enumerates every resulting six-vertex path while checking no frame collisions.

The three disjoint frame counts are 2,160, 21,600 and 65,340; their path counts and unordered two-path-cover counts are 77,760, 777,600 and 2,352,240. The q=81 AP/AP subcount is 432 frames and 15,552 distinct paths. Extension-field integers encode polynomial coefficients; arithmetic modulo the integer q would be wrong. No integer lift or new existence theorem is asserted.
