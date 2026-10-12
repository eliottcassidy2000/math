# Exact eight-colour obstruction for whole Stechkin tiles, and a local split escape

2026-10-12. **PROVED, scoped:** the fixed-scale trihexagonal partition below requires exactly eight colours if every tile interior is monochromatic, with no periodicity assumption. Its triangle-only region graph requires exactly four. **FINITE-EXACT:** a source-bound fourteen-region obstruction, a seven-colour finite split witness, and an eight-colour obstruction for one specified periodic split. **OPEN here:** a pointwise six-colouring of the unit-distance plane. These are graphs of region conflicts, not claimed eight-chromatic graphs of unit-distance points.

## 1. Inheritance and the changed obligation

The [Stechkin coverage packet](stechkin_whole_plane_obstructions_20261011g.md) recovers the normalized mixed-distance type $(1,1,1,1,1/2,1/2)$. It proves that the two triangle classes each contain unit equilateral triangles, and that separately rescaling those classes leaves a region of density $3/16$ uncovered. The [square-cell packet](square_cell_conflicts_20261011g.md) proves an eight-colour bound for its specified forty-cell, fixed-period square model. The [diagonal holonomy packet](diagonal_palette_holonomy_20261011g.md) excludes its particular commuting-reflection palette mechanism.

Those results do not settle larger colour periods on the Stechkin tiles. This packet examines that obligation directly. The closest positive mechanism is exact convex-region conflict compilation with a complete context graph. The canonical hostile is mistaking a quotient or region clique for simultaneous unit distances between chosen representatives. The corrected near miss is confusing six sets avoiding different distances with six sets avoiding one distance. The useful retained coordinate is **position inside a tile**, which an indivisible-region colour erases.

The active objects are the original tile, its open interior, its lattice address, its exact distance interval to another tile, and its colour constraint. The map from geometry to the conflict graph preserves the possibility of a colouring constant on those interiors. It discards the pair of points witnessing each edge. Such witnesses may differ from edge to edge.

The incoming [eleventh-root layers](eleven_port_layers_20261011.md) and [rhombus overlap barriers](rhombus_overlap_elliptic_20261011.md) concern different point-graph and field interfaces. No result below circumvents them or uses a field escape as a chromatic lower bound.

## 2. Exact tile and conflict compiler

Use triangular coordinates in the basis

$$u=(1,0),\qquad v=(1/2,\sqrt3/2),\qquad q(a,b)=a^2+ab+b^2.$$

For $(i,j)\in\mathbb Z^2$, the three open regions are translates by $iu+jv$ of the interiors of these convex polygons:

| Kind | Vertices in triangular coordinates |
|---|---|
| $H$ | $(1/2,0),(1/2,-1/2),(0,-1/2),(-1/2,0),(-1/2,1/2),(0,1/2)$ |
| $T_-$ | $(1/2,0),(1/2,1/2),(0,1/2)$ |
| $T_+$ | $(1/2,1/2),(1,1/2),(1/2,1)$ |

Write $H(i,j)$ and $T_\pm(i,j)$ for the translated interiors. These are the side-$1/2$ regular hexagons and the gaps of the inherited construction. The hexagon has diameter one. Its open interior contains no unit pair; neither does either open triangle.

Two region vertices are adjacent iff **some** pair of points, one in each region, has distance one. For full-dimensional convex polygons $P,Q$, form the closed difference hull $D=\overline Q-\overline P$ from all vertex differences. Let

$$m=\min_{z\in D}q(z),\qquad M=\max_{z\in D}q(z).$$

Then their open interiors conflict exactly when

$$m<1<M. \tag{1}$$

Indeed $Q-P$ is the interior of $D$. The squared norm on this connected open set has neither a positive minimum nor a maximum; continuity supplies every value strictly between its infimum and supremum. Its maximum over $D$ is attained at a vertex. Its minimum is zero if $0\in D$, otherwise it is found by exact projection to each edge using the Gram form $q$. This gives a rational, complete compiler, not a numerical distance tolerance.

All base coordinate differences lie in $[-3/2,3/2]$. If one coordinate of an integer address difference has absolute value at least three, the corresponding displacement coordinate has absolute value at least $3/2$. Since

$$q(a,b)\ge\tfrac34a^2,\qquad q(a,b)\ge\tfrac34b^2,$$

its squared distance is at least $27/16>1$. Thus all possible conflicts occur among the twenty-five offsets in $[-2,2]^2$. The complete directed counts, in order $H,T_-,T_+$, are

$$\begin{pmatrix}12&12&12\\12&6&3\\12&3&6\end{pmatrix}. \tag{2}$$

The six same-orientation triangle offsets are

$$(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),$$

and the three $T_-\to T_+$ offsets are

$$(-1,-1),(-1,1),(1,-1).$$

The twelve $H\to H$ offsets are the six nearest triangular-lattice vectors together with

$$(-2,1),(-1,-1),(-1,2),(1,-2),(1,1),(2,-1).$$

The script keeps all mixed offsets and computes every window directly. Two useful endpoint controls are $H(0,0)$ versus itself, with window $[0,1]$, and $H(0,0)$ versus $H(2,0)$, with window $[1,9]$. Neither is an **open-interior** conflict. Boundaries require a separate ownership rule.

## 3. Triangle-only regions require exactly four colours

Give both orientations the colour $(i\bmod2,j\bmod2)$. Every same-orientation conflict changes this parity pair, while each cross-orientation conflict changes both bits. This is a four-colouring of the entire triangle-region graph.

For the lower bound, take the six $T_+$ addresses

$$(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),$$

and add $T_-(0,0)$. The six first vertices form the triangular lattice patch of side two, with nine edges. Its edge-connected triangles make any proper three-colouring unique up to palette permutation. Its three corner vertices $(-1,-1),(-1,1),(1,-1)$ have distinct colours. The added region conflicts with exactly those three corners. The resulting seven-region, twelve-edge graph is not three-colourable.

The companion independently generates and checks a complete three-colour branching refutation, with 103 proof nodes, as well as the explicit four-colouring. This is stronger than the previous unit-triangle obstruction when the remaining regions must stay indivisible. It does not prove that the actual point graph on the union of all triangle interiors has chromatic number four: a pointwise colouring may split a region.

## 4. Fourteen regions forbid seven colours at every period

The following ordered source regions have exactly 65 conflict edges:

| Index | Region | Index | Region |
|---:|---|---:|---|
| 0 | $H(0,0)$ | 7 | $T_+(-1,0)$ |
| 1 | $H(0,1)$ | 8 | $T_+(-1,1)$ |
| 2 | $H(1,0)$ | 9 | $T_+(0,-1)$ |
| 3 | $H(1,1)$ | 10 | $T_+(0,0)$ |
| 4 | $T_-(0,0)$ | 11 | $T_+(0,1)$ |
| 5 | $T_-(0,1)$ | 12 | $T_+(1,-1)$ |
| 6 | $T_-(1,0)$ | 13 | $T_+(1,0)$ |

Here is a short certificate using only incidences authenticated by (1).

1. Vertices $1,2,3$ form a clique and are adjacent to every other vertex. Their colours cannot occur on the remaining eleven vertices $W$.
2. On $W$, the six vertices in the cyclic order $(7,8,11,13,12,9)$ induce $C_6$. Vertex 10 is adjacent to the entire rim. Vertices $4,5,6$ form a clique, and each has a neighbour in each of the two alternating rim triples.
3. Consequently every independent set in $W$ has size at most three. If it avoids 0 and 10, it takes at most one of $4,5,6$ and at most three rim vertices; a rim triple must be alternating and cannot accompany any of $4,5,6$.
4. The two special vertices 0 and 10 are adjacent. The non-neighbours of each within $W$ form a clique, so an independent set containing either special vertex has size at most two.

If $W$ had four colours, the two special vertices would belong to different colour classes. The total capacity of those four independent classes would be at most

$$2+2+3+3=10<11.$$

Thus $W$ needs at least five colours and the full graph needs at least eight. This is an elementary finite certificate, not an unrecorded SAT conclusion. The script rechecks every geometric premise and independently enumerates the independent-set counts in $W$: $11,26,14,0$ at sizes one through four.

**The lower bound has no period hypothesis.** Any assignment of at most seven colours constant on all original tile interiors fails on this one finite patch. It applies even if its colour pattern is nonperiodic or its boundary ownership is altered arbitrarily.

**Metric hostile.** For example, the centres of $H(0,0)$ and $H(1,1)$ have squared distance three, whereas their region window is $[3/4,7]$ and they conflict. Replacing a region by its centre destroys edges. Even a region clique supplies different witnesses for different pairs, not one simultaneous set of unit-separated representatives.

## 5. Eight colours suffice, including every boundary

Use disjoint palettes of four colours each, with

$$c(H(i,j))=2(i\bmod2)+(j\bmod2),$$
$$c(T_\pm(i,j))=4+2(i\bmod2)+(j\bmod2). \tag{3}$$

The complete offset table proves (3) proper on open interiors. The two-period quotient has twelve vertices and 54 edges, and carries the same positive certificate.

To colour the full plane, keep triangle interiors open. Every other point belongs to a closed hexagon. A nonvertex hexagon boundary point belongs to just one hexagon. At a vertex, two hexagons meet; assign that vertex to the lexicographically smaller lattice centre. This assigns every point exactly once.

The only pairs at distance one within a closed regular hexagon are opposite vertices. Lexicographic ownership assigns exactly one vertex from each opposite pair to that hexagon. For different hexagons of the same parity colour, centre differences lie in $2\Lambda$. Their shortest nonzero length is two. At length two, the two closed hexagons can reach distance one only at their facing vertices, and the lexicographic rule assigns one of those vertices to the intervening hexagon, which has a different parity colour. All remaining same-colour centre differences have length at least $2\sqrt2>2$, since their squared triangular norm is an integer larger than one before doubling. Their closed hexagons are therefore separated by more than one. Different hexagon and triangle palettes cannot cause a same-colour cross conflict.

This proves an actual eight-colouring of all points of the plane constant on each tile interior. Together with Section4, the optimum for this fixed partition is exactly eight. The rational `locate` API implements the ownership convention; the proof applies to real points as well. Exact rational unit-pair tests exercise both interiors and facing-vertex boundary cases.

## 6. The omitted coordinate repairs one finite patch, but not yet the plane

Bisect only region 0 of Section4 along the symmetry axis $a+b=0$. In the new ordered list, first put the half $a+b>0$, then $a+b<0$, then old vertices $1,\ldots,13$. The complete fifteen-piece conflict graph has 70 edges and the following proper seven-colouring:

$$ (0,5,6,4,5,0,3,2,1,3,3,2,1,1,0). \tag{4}$$

All cross-piece windows are recomputed, including constraints between pieces which formerly had different colours. No edge is discarded on the basis of a preferred witness. The internal bisecting line can be assigned to the first half, retaining the original open hexagon boundary. This does not introduce a new same-colour conflict: every other region remains open, so the difference with that other region is open with the same closure and the same strict test (1). The two halves have different colours; each remains inside the original unit-free open hexagon. Thus (4) colours the union of these fourteen original interiors pointwise, including the internal cut.

This is a constructive use of the previously lost within-tile coordinate. It is only a finite patch. Splitting finitely many original tiles cannot repair the whole infinite model: translate the fourteen-region obstruction away from that finite modified set.

**One periodic escape test, with complete context retained.** Split every $H(i,j)$ with $i,j$ both even along the same local axis, and require a period-two colour pattern. This hits every translated old witness, whose four hexagon addresses form a $2\times2$ block. The quotient now has thirteen pieces, 63 edges and no loops, but still needs eight colours. Its explicit eight-clique consists of the negative half of $H(0,0)$, the three unsplit hexagon residues, and all four $T_-$ residues. Every clique edge carries a concrete translation lift and an exact distance window. This is a quotient clique, not a set of eight simultaneously unit-separated points. The old eight-colouring remains a positive witness by giving both halves the old hexagon colour.

The periodic quotient is complete: a translation by $2(i,j)$ with $\max(|i|,|j|)\ge2$ puts parent-centre address differences outside $[-2,2]^2$, so the inherited exterior bound excludes all omitted interactions. The tested periodic split therefore has optimum eight. This last lower bound is restricted to that split pattern and period; the unsplit lower bound of Section4 was unrestricted in colour period.

The next open construction needs a different infinite refinement or a finer colour rule inside tiles, and must retain all translated cross constraints. Merely removing one familiar obstruction, or introducing a new within-tile bit, is insufficient evidence for a global improvement.

## 7. Exact certificate boundaries and reproduction

The production interfaces retain exact rational polygons, tile type and original lattice address. `verify_eight_obstruction` authenticates the marked capacity proof against freshly compiled geometric conflicts. `split_packet` retains every pair in the finite split context; `periodic_split_packet` also retains a translation lift for every quotient edge. All public coordinate and tile inputs reject boolean/float aliases; no floating tolerance or external SAT solver is used by the saved verifier.

| Source | Map / retained predicate | Lost data and needed correction |
|---|---|---|
| Open tile interiors | Exact region graph; whole-region colourability | Individual edge witness points differ |
| Whole region | Centre or one representative | Unit-conflict edges can disappear |
| Infinite unsplit partition | Finite fourteen-region restriction | Enough for a lower bound at every period |
| Split finite patch | Complete fifteen-piece context | Does not cover all translates |
| Periodic split | Quotient with all edge lifts | Palette periodicity is an explicit extra hypothesis |
| Open interior colouring | Full-plane colouring | Boundary ownership must be proved separately |

[Script](../../04-computation/experiments/stechkin_refined_conflict_graph_20261012a.py), [saved output](stechkin_refined_conflict_graph_20261012a.out):

```text
python -B -X utf8 04-computation/experiments/stechkin_refined_conflict_graph_20261012a.py
python -B -O -X utf8 04-computation/experiments/stechkin_refined_conflict_graph_20261012a.py
```

The declared finite universe is all nine ordered type pairs and address offsets $[-4,4]^2$ (with the proved bound making $[-2,2]^2$ complete), every pair of the fourteen-region witness, every pair in the fifteen-piece split and complete thirteen-piece periodic quotient, all independent subsets of $W$ through size four, twenty-five translated marked witnesses, and exact rational point/boundary controls. Finite type and wrong-source hostiles are included. The arbitrary-period and full-plane claims rest on the proofs above, not extrapolation from these controls.

Normal and optimized runs agree with the saved output at **7,008 always-active exact checks**; LF SHA256 `07c5814733cd85743ba4476c36de8a14cd2246de7a7af87c961bfd4866dbcffe`.
