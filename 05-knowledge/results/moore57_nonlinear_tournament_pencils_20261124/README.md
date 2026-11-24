# Degree 57: from parallel tournament palettes to nonlinear pencils

## The results

1. **Hoffman–Singleton's uniform rule fails in the projective degree-57 model.** On the known graph's projective chart, all 35 line-centers have parallel local partitions. In PG(3,7), not only can two adjacent centers not both be parallel: they cannot both use the broader *plane-slice* pattern. The finite test rejects all 46,144 slice pairs up to the specified normalization. Consequently, **at least 2,450 of the 2,850 line-centers must use genuinely nonlinear line pencils** in any degree-57 completion using only regular spreads in this model.
2. **Single-center feasibility is not the obstruction.** There are 6,031,872 exactly two-direction partitions at a fixed center, but two centers with that form cannot be adjacent. In the same model, all 721 relative regular-spread anchors pass the first global Hall test, each matching 1,920 distinct cross-pair demands to 1,920 witness lines. The matching is not a solution: it omits simultaneous compatibility of the witness spreads.
3. **A tournament target independent of this model.** Two short partial blue words exactly detect the girth defects in any 57-regular graph on 3,250 vertices. Their weighted count hits 6,021,051,175,776,000 if and only if the graph is Moore. The same objective works for the smaller degrees, but a verified degree-seven example already stalls under strict single-switch descent. See [blue words and switch experiments](alternating-blue.md).

**Status:** local normal forms are proved; plane-slice rejection, Hall checks and known-case replays are finite exact; all degree-57 conclusions are conditional; the compatible pencil conjecture is open. No degree-57 graph is claimed. The local pencil results require the PG(3,7) point/line model and regular (Desarguesian) pointed spreads. They do not rule out irregular spreads, other Steiner designs, or a graph without a tight coclique.

## Inherited geometry and what changes

The [graph-spine session](../moore57_graph_spine_tiles_20261124/) gave a cheap hostile case for simply replacing the field in the 5-cohort affine formula: mixed midpoint counts force field size five. The [seam session](../moore57_metagraph_seams_20261124/) identified the complementary failure of choosing branch contacts independently. This session changes the carrier from a fixed affine cohort partition to pointed projective spreads, and replaces the single tournament palette by a constrained family of local direction words. Forgetting the spread labels leaves the ordinary tournament; preserve the labels as a sidecar to evaluate the two-endpoint compatibility rules. The alternating-blue test is the separate order-independent check.

## Pointed spreads and tournaments

Start with the known projective-space representation of Hoffman–Singleton and its previously proposed [PG(3,7) analogue](https://doi.org/10.37236/1224). The size-400 coclique in this analogue is a modeling assumption. Point-type vertices form the coclique; line-type vertices have their projective incidence neighbors, together with an auxiliary line graph H. Every line B selects a *pointed spread* P_B: 50 mutually skew lines partitioning the 400 points, one of which is B. Its spokes P_B minus B are the H-neighbors of B. The term regular below refers to the Desarguesian geometry of the spread, not merely to vertex degrees.

For a fixed order of the line labels, record each proposed H-contact by reversing the corresponding initially transitive tournament arrow. A same-label source/observer copy can protect the reference order if a fixed directed Hamilton spine is wanted. The pointed spreads couple those arrow choices at both endpoints; they are not independent tiles. The precise inverse is a pairwise constraint on which pointed-spread labels can sit at the tournament's vertices:
- P_B and P_C must select each other or neither.
- If they select each other, their intersection is exactly {B,C}; otherwise they intersect in at most one line.

These conditions are sufficient as well as necessary. In the assembled graph, point-line incidence and the spread partition eliminate every triangle or square involving a point. The two overlap rules eliminate triangles and squares entirely within H. The regular graph has 3,250 vertices, exactly the size of a collision-free radius-two neighborhood at degree 57. An alternative cover-only formulation, spectral checks and resolution-bank rules are in [pointed spreads](pointed-spreads.md).

### The chart at one center

Use a general finite field with q elements, taking q=2 or 7 here. Fix B as the vertical line in a decomposition of the four-dimensional vector space into two 2-spaces. The q⁴ lines skew to B are graphs of linear maps M in the affine space M₂(F_q). Their skewness test is det(M−M′)≠0. In this chart every regular spread containing B is the vertical line plus an affine 2-flat M₀+U, where U=A F, A is invertible and F is a copy of F_(q²) in M₂(F_q). To see the formula, use the standard field spread and the affine chart transformations M↦P M Q+R preserving the vertical line. The factor Q can be absorbed by conjugating F and replacing A. Conversely these operations exhibit each such flat as a regular spread.

Let the finite portion S of P_B be the central 2-flat. For each of its q² lines Y, P_Y shares B and Y with P_B. Two different Y have nonadjacent centers (otherwise they form a triangle with B), so their spreads share only B. **The finite portions of P_Y therefore partition all q⁴ matrices into q² affine 2-flats**, each meeting S in exactly its center. Choose S through zero by a translation.

The directions of any two disjoint 2-flats meet nontrivially, since complementary directions would force a unique intersection. Thus the directions either share a common 1-space L or all lie in a common 3-space W. For the proof, take two distinct directions U,V, put L=U∩V and W=U+V. Any direction not containing L meets U and V in distinct lines so lies inside W; that direction excludes every direction through L that leaves W. The parallel case is included in both descriptions, and the two descriptions can overlap.

- **Plane slices:** If the directions lie in W, each of its q affine cosets is partitioned by q parallel 2-flats (two different 2-directions inside a 3-space cannot be disjoint). The central S meets each slice in a line, and the chosen direction in that slice must be transverse to S. This yields at most q local direction states.
- **Line pencil:** Quotient by L. The image of S is a plane, and the q² flats become affine lines with q² different basepoints y on that plane. Set the slice coordinate to t and write the other coordinates as y+t f(y). The partition condition is **y↦y+t f(y) is a permutation for every t∈F_q**. This is also sufficient. Allowed slopes come from the regular directions containing L and transverse to S. In tournament language every contact between two basepoints y,z must forbid exactly the nonzero proportional slope differences f(y)−f(z)∈F_q^×(y−z); fixed orientation of that contact is just a way to store this symmetric constraint.

After scaling L to the span of the identity at q=7, regular directions through it are the 21 trace-zero projective directions (a,b,c) for which a²+bc is a nonsquare. The [standalone checker](spread-check.py) enumerates the 57 projective lines and finds that they contain 0, 3, or 4 such directions. Hence no nonconstant affine slope function can work: rank two would require 49 directions; rank one would use seven collinear directions in the affine slope chart. A constant function is parallel and is already plane-slice. This is a structural distinction, not just a different field size.

## Why plane slices cannot touch at q=7

For adjacent B,C, use B=vertical, C=horizontal. Normalize the finite chart S of P_B to a copy of F_49. The finite chart T of P_C meets S only at zero. The candidates are the 2-spaces A F: 21 embedded fields (the conjugacy orbit under GL₂, whose normalizer has size 2·48), each with 42 cosets of its multiplicative group, hence 882 regular spreads through the pair. There are 721 options for T. Indeed each of the 48 nonzero vectors of S is in 21 chart spaces; after removing S, each intersecting space accounts for six vectors, giving (48·21−48)/6=160, hence 882−1−160=721. The choice of S as any listed field-coset direction loses no possibilities: changes of coordinates fixing the two axes act transitively on these directions.

If the chart at B has plane slices, the slice through C has direction T since it contains P_C, and forces P_D to have finite portion D+T for the six nonzero D on a scalar line of S. In the chart at C the change is inversion of nonsingular M; the finite portion of P_B in that chart is again S. The slice through B forces P_E to have finite portion E⁻¹+S in that chart for six E on a scalar line of T (inversion takes a scalar line to a scalar line). There are eight possibilities on each side. For each pinned D,E, neither is an H-neighbor of the other, since an edge would close a square through B,C. Their selected spreads must intersect once if D,E are skew, and zero times otherwise. To see necessity of the zero/one distinction, an intersecting pair already has a point witness; a skew nonedge needs one and only one line witness. Thus the exact test is
$$
\#\{X\in D+T:\det X\ne0,\ X^{-1}-E^{-1}\in S\}
 = [\det(D-E)\ne0].
$$
The vertical and horizontal lines cannot contribute, since S and T are transverse and D,E nonzero.

[anchor-hall.cpp](anchor-hall.cpp) enumerates all 721×8×8 choices. The **best** has 12 mismatches among the 36 pinned pairs. Even all-parallel neighborhoods alone fail at least 1,752 of 2,304 cross pairs over every anchor.

Therefore the plane-slice centers form an independent set of H. Since H is 49-regular and its smallest eigenvalue is at least −8 by interlacing with the full Moore adjacency spectrum, the ratio bound is 2,850·8/(49+8)=400. For completeness, test the Rayleigh quotient of the centered indicator of an independent set against −8 to obtain the bound. The other **at least 2,450** centers fall in the line-pencil case but not in the plane-slice case, so have a nonaffine slope function. At q=2 the same exact CSP constructs the Hoffman–Singleton graph from 15 points and 35 lines, and all 35 centers give parallel local partitions. This is the precise failure of the simplest extrapolation to 57.

## A live constructive test: adjacent charts and the 1,920 witnesses

The two adjacent anchored stars have 48 noncentral spokes apiece. Every spoke D on one side meets exactly eight lines of the other spread (one disjoint spread line per point on D), all among the 48 opposing spokes. Thus 48·40=1,920 cross pairs D,E are skew. These are H-nonedges, so each needs one common line-neighbor X. A single X cannot witness two pairs: it would then meet two spokes at the same center, creating a four-cycle. The 1,920 witnesses are exactly the invertible graph matrices outside the two anchored spreads. Indeed they must be skew to both centers (otherwise that center and X have a point witness as well as the common spoke); exactly 2,016−96=1,920 such lines are available. The residual 832 lines divide into 384 left-only witnesses, 384 right-only witnesses, and the 64 lines meeting both B and C, naturally indexed by an 8×8 grid of points on the two centers. For the one-sided count, each of the 48 spokes has 48 neighbors besides its center, 40 of which are already used; a common neighbor on that side cannot be reused. The q=2 replay has 0 two-sided witnesses and an edgeless 3×3 residual grid for all 70 line edges; it does not assert that the q=7 residual grid is independent.

For each of the 1,920 slots (D,E), connect it in a bipartite graph to X only if some regular spread containing D,E,X meets the anchored spreads precisely at D and E. Fix D,E as coordinate axes and use the matrix transform (M−D)(M−E)⁻¹; enumerate the 882 charts and discard those containing forbidden anchored lines. The original vertical line maps to I. An exact maximum matching passes for **all 721 anchors**. Each has 1,080,192 to 1,134,624 admissible slot-witness incidences. This is a local test only: each matched X still has to choose a spread agreeing with the other witnesses.

**Working conjecture to test next.** Some regular anchor supports a choice of two genuinely nonlinear line-pencil charts that agree on all their shared incidences and admit a *compatible* 1,920-witness assignment. Here compatible means witness spreads at least satisfy the pair constraints within each common-spoke row or column, not just each slot separately. This remains much weaker than completing all 2,850 centers. An implementable next test is an exact-cover or SAT layer indexed by (slot, witness, regular chart), reserving each X once and forbidding two charts sharing a spoke and any further line. In parallel try the other tier: solve the two adjacent pencil constraints first with slope maps on F₇² and use the Hall layer to prune. The main heuristic change is to search nonlinear permutation tables, not to extend one small tournament's uniform rule.

## Reproduction

```sh
python spread-check.py
c++ -O3 -std=c++17 anchor-hall.cpp -o anchor-hall
./anchor-hall plane-slices
./anchor-hall parallel
./anchor-hall 0 721   # exhaustive independent Hall layer; about two minutes
python moore-alternating-blue-check.py
python moore-objective-cubic.py
python moore-objective-d7-fixture.py
python moore-skew-square-check.py
```

[Signed skew squares](moore-skew-square.md) gives another independent tournament interface: complement signs determine the underlying edge orientations up to total converse when local star products and global parity are both satisfied. The 2026 [Smith–Montemanni paper](https://www.mdpi.com/2227-7390/14/23/4318) excludes a different approach using a group of derangements; neither that result nor the tests here settle unrestricted degree 57.
