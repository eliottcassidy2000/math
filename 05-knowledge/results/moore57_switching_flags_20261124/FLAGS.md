# Exact local surgery on the asymptotic partial-flag construction

**CITED** for the underlying family and **PROVED** for the exact degree, generic inverse and equivariant thinning below. The q=2 nonequivariant result is **FINITE-EXACT**, with replay code. These have 1085 vertices at q=2, not the 3250 required by Moore(57,2).

[Cames van Batenburg and Korsky, §3](https://arxiv.org/html/2608.03965v1#S3) build their diameter-k graphs by retaining even-dimensional flags and passing through the odd-dimensional flags in a bipartite compatibility graph. Odd/even sorting gives paths, not a unique Moore address. Here is the degree-diameter 2 case with the witness data retained.

## The finite-field chart

The vertices are (U,W) with U a 2-space inside a hyperplane W in a 5-space over F_q. Distinct flags (U,W),(V,Z) are adjacent exactly when U,V share a line and U+V ⊆ W∩Z (this includes U=V). This is equivalent to compatibility with a common odd flag (line in 3-space): choose a shared line and any intermediate 3-space between U+V and W∩Z.

Write [i]_q=1+q+...+q^(i−1). The order is [5]_q(q²+1)[3]_q. Count neighbors of a fixed flag in three disjoint types:

| Type | Number |
|---|---:|
| same U, different W | q(q+1) |
| different U, same W | q(q+1)² |
| both coordinates different (D) | q²(q+1)² |

For same U, choose a different hyperplane among the [3]_q containing U. For same W, choose the common line in U and a different plane through it in W. For type D, choose a line of U (q+1), a 3-space through U in W (q+1), a different hyperplane through the 3-space (q), then a plane in that 3-space through the selected line but not equal to U (q). Thus the unthinned degree is q(q+1)((q+1)²+1).

For disjoint U,V let A be any line of U∩Z and B any line of V∩W; both exist. Then (A+B,U+V) is a common neighbor of the endpoints. On the cell where dim(U∩Z)=dim(V∩W)=1 this is the unique midpoint: the candidate plane must meet each of U,V in those respective lines, and its hyperplane must contain U+V. Those targets are nonadjacent. For a fixed (U,W) there are q⁴+q³ choices of Z not containing U. For each, let A=U∩Z and H=W∩Z. The admissible V are the planes in Z neither containing A nor contained in H. Count them by inclusion-exclusion: Gaussian(4,2)−2[3]_q+[2]_q=q⁴+q³. So exactly (q⁴+q³)² targets per root belong to this explicit unique-midpoint cell, a proportion 1−O(q^(-2)) as q grows. This is a computable inverse, not just a density bound. It does not constrain the other pairs to have one midpoint.

## Every D edge is indispensable, and all same-U edges can go

All D edges form one GL(5,q)-orbit: adapt a basis to the common line, the two plane extensions, and the two distinct hyperplane extensions. One representative occurs on the unique generic path

    (12,1234)  —  (13,1235)  —  (35,1345),

where digits index basis vectors, and the first and last flags are disjoint-plane endpoints with both cross intersections lines. Therefore every D edge must remain in **every spanning subgraph of diameter two**. The subgraph that keeps just D edges fails. For disjoint U,V in the *same* hyperplane W, every midpoint hyperplane must be W=U+V, so both links of every route are same-W edges.

Keeping D and all same-W edges, and dropping all same-U edges, suffices. The disjoint case uses the midpoint above, whose plane differs from both U and V. For U,V sharing a line L, choose a plane X through L in W∩Z different from U,V (there are at least q+1≥3 candidates), then a hyperplane Y containing U+V+X (dimension at most four). Both legs change plane. Thus the minimum-degree GL(5,q)-invariant spanning diameter-2 subgraph has degree q(q+1)³; no edge deletion from its two required orbits is GL-equivariant. At q=2 that degree is 54 (from 60).

For comparison, at the thin set-flag limit, replace subspaces by subsets of [5]. The 30 even flags are ordered (2,2,1)-partitions. Two are adjacent iff choose permutation representatives agreeing at odd prefixes (1,3), i.e. one block of parallel swaps at positions 2 and 4 of the five-symbol permutation. The alternating sorting network gives diameter two after passing to these flags. Degree is 10 by the three counts at q=1, **not** a single tournament-arc flip. This bridge preserves the compatible shared representative, not spectrum or path counts.

## Sharp q=2 thinning when same-U edges are forbidden

Fix a hyperplane W. The 35 planes U in W are the vertices of its repair graph (only same-W edges). Each line A in W gives a seven-point Fano pencil: its points are planes through A, and its lines are the three planes through A contained in a particular 3-space H in W. Each candidate same-W edge is in exactly one pencil, for A equal to the intersection of its two planes.

**Necessary two hits per Fano line.** Take a pencil point U and a Fano line H not containing it, so U∩H=A. Choose W' with W∩W'=H. For each of the three pencil points V on H, consider the pair (U,W),(V,W'). It is not adjacent. A common neighbor (X,Y) must have X in H meeting U, so A⊂X; X must also differ from V since same-U edges are forbidden. If X≠V then (U+V)+X=W, as (U+V)∩H=V. Hence Y=W, and the first leg is a retained same-W edge U--X. To succeed for every V, the retained neighbors of U must meet the three-point line H at least twice.

So the omitted graph in each pencil has max degree two, and whenever a vertex has two omitted neighbors, those three points are collinear in the Fano plane. Indeed two omitted neighbors not collinear with the vertex would form a deficient line avoiding it, and three distinct omitted neighbors cannot all be pairwise collinear with the vertex. Its components are isolated points, edges, two-edge paths, or Fano triangles: a longer path would force two lines sharing two points to coincide. There cannot be two triangle components, as two Fano lines intersect. Thus at most five of the 21 pencil edges are omitted. At least 16 per pencil are needed, and each W needs at least 240 repairs. Its 35 flags require max repair degree at least ceil(480/35)=14. Adding the 36 forced D neighbors gives max degree at least 50 anywhere in this model.

**Attainment.** The deterministic [certificate](flag-check.py) omits a Fano triangle and a matching on its complement in each pencil, choosing the 15 patterns jointly to cover all 280 skew plane pairs in W and to place every plane in at least one omitted triangle. The retained degree distribution per W is 14 on 26 planes, 13 on eight, and 12 on one. Transport the same template to each W by a linear isomorphism. The script verifies every one of the 1085 binary flag vertices reaches every other within two steps, by integer bitsets, without relying on a generic-case assumption. It gives 19,530 forced D edges plus 7,440 repairs, and max degree 50, both optimal within the specified no-same-U model. The construction does not claim an overall extremal degree-diameter record.

The logical distinction matters for transport to the Moore search: generic pairs force most edge classes here, but the nongeneric same-W defect calls for pencil-specific, overlapping local repairs. Independent choices that hit one pencil do not guarantee that two skew planes can share a midpoint. In a Moore graph there is no spare witness to pay for either type of duplication.
