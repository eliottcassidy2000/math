# What Fano triangle moves and surface geometry retain about flows

2026-10-04. **PROVED** for the elementary flow quotient, exact fibres,
face-boundary formulas, and conditional Euler inequality below.
**FINITE-EXACT** for the declared flow census. **CITED** for the named
external results, with primary sources checked on this date. Tutte's
5-flow existence conjecture remains **OPEN**. No historical novelty is
claimed for triangle/star operations or surface flow/color duality.

## 1. Inheritance and two selected connections

[THM-4529, Petersen and Heawood families in Paley coordinates](../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md)
and its [full result, sections 2--3](procgen_petersen_20261001_petersen_heawood_paley.md)
already construct the Heawood graph from the seven edge-disjoint Fano
triangles of the Paley tournament on seven vertices. They also recover
Petersen by deleting a Heawood point and suppressing three degree-two
vertices. The old Hamiltonian-path parity does not survive these graph
moves. The Fano colouring dictionary and snark obstruction already appear
in [the Paley/snark note, section 3](collatz_label_map_paley_snarks_20261001.md).
They are inherited, not new numerological connections.

The present geometric input is the companion rotation construction:
the same labelled Heawood graph has the two rotation systems specified
in section 3 below, with seven hexagonal faces on a torus or three
14-gonal faces on a genus-three surface. The script reconstructs their
dart permutations independently. No graph automorphism census is used.

The closest conserved quantity is signed boundary flow. The hostile is a
nonzero circulation around a triangle whose boundary current is zero.
The corrected near miss is planar flow/color duality used without its
homology hypothesis. The least-used sidecars are the seven triangle
circulations and the surface homology class. The board is Paley arcs,
Fano incidence, group flows, face colors, homology, and curvature.

Only two connections are pursued:

| Source -> target | Explicit map and preserved property | Lost information / repair |
|---|---|---|
| Paley-oriented `K7` flows -> Heawood flows | Replace each triangle by its three outgoing boundary currents; Kirchhoff conservation survives | One circulation per triangle; retain seven group elements |
| Surface face colors -> graph flows | Signed differences across the two incident faces; proper coloring gives nowhere-zero flow | Global color translation, and most graph flows are not boundaries; retain the embedding and homology class |

## 2. Seven triangle circulations are exactly the lost coordinates

Let `A` be a finite abelian group, written additively, of order `q`.
A flow assigns group values to oriented edges with zero total outgoing
minus incoming value at every vertex. Reversing an edge negates its value.

Orient a triangle `x->y->z->x`, with respective values `a,b,c`.
Its outgoing currents at the old vertices are

```text
beta_x=a-c, beta_y=b-a, beta_z=c-b; their sum is zero.       (1)
```

Replace the triangle by a new vertex `t` and three edges directed from
`x,y,z` to `t`, carrying these currents. This preserves conservation at
the old vertices, and the new vertex is conserved by their zero sum.
Conversely, every zero-sum triple has exactly `q` triangle lifts:

```text
(a,b,c)=(h,h+beta_y,h-beta_x), h in A.                     (2)
```

If all three star currents are nonzero, the three forbidden values
`0,-beta_y,beta_x` are distinct. Hence a nowhere-zero star has exactly
`q-3` nowhere-zero triangle lifts. This includes the sharp failure at
`q=3`: a nonzero star exists but none of its triangle lifts stays nonzero.
For a general boundary triple, the number is `q` minus the number of
distinct forbidden values. Zero-current patterns must not be discarded.

Now take points modulo seven and `D={1,2,4}`. The Fano line
`L_t={t+1,t+2,t+4}` is the directed triangle
`t+1 -> t+2 -> t+4 -> t+1` of the Paley orientation. These seven triangles
partition all 21 edges of `K7`. Performing (1) on every triangle produces
the Heawood incidence graph, oriented from point `p` to line `L_t` when
`p in L_t`. The result is an exact sequence of abelian groups

```text
0 -> A^7 -> Flow_A(K7) -> Flow_A(Heawood) -> 0.             (3)
```

Surjectivity follows by choosing one `h` in (2) independently for every
line. The kernel consists exactly of constant circulations on the seven
triangles. The connected-graph cycle ranks agree:
`21-7+1=15`, `21-14+1=8`, and `15-8=7`.

For `q>=3`, write `NZ_A(H)` for the number of nowhere-zero Heawood flows.
Then the number of nowhere-zero `K7` flows **whose projection also stays
nowhere-zero** is exactly

```text
(q-3)^7 NZ_A(H).                                         (4)
```

Equation (4) counts a guarded subset of `K7` flows, not all of them.
For any nonzero `g in A`, giving every directed Paley arc the value `g`
is a nowhere-zero flow: each of the seven triangle circulations has
zero boundary. Its Heawood projection is identically zero. This is the
explicit hostile to transferring a certificate merely because the graph
move is defined. For `q>=4`, the reverse direction does work: a
nowhere-zero Heawood flow always has nowhere-zero triangle lifts.

## 3. One graph, two surfaces, different boundary subspaces

Use points `P_p` and lines `L_t`, indexed modulo seven. At a point, order
the three neighbours as

```text
(L_(p-1), L_(p-2), L_(p-4));
```

at a line, order them as `(P_(t+1),P_(t+2),P_(t+4))`.
Let `alpha` reverse a dart and `sigma` advance around these orders.
Faces are cycles of `sigma alpha`. Keeping both displayed orders gives
the first map below; reversing every line order gives the second.

| Point/line signs | Face boundaries | Euler characteristic / genus | Embedded dual |
|---|---|---|---|
| `+,+` | seven 6-cycles | `0`, genus 1 | `K7` |
| `+,-` | three 14-cycles | `-4`, genus 3 | three vertices with seven parallel edges between each pair |

The dual statements are literal: no dual loop occurs, and the displayed
edge multiplicities are exact. The script checks all 42 darts and 21
dual incidences. The two occurrences of `K7` in this note have different
roles: the triangle-move source in section 2 and the torus dual here.
Their equal abstract graphs do not make those two maps equal.

For a connected cellular embedding in a closed orientable surface, assign
a color `h_f in A` to each face. The signed face-boundary sum assigns an
edge the difference of its two incident face colors. Around each vertex
the differences telescope, so it is a flow. It is nowhere-zero precisely
when adjacent faces have different colors. The dual is connected, so two
colorings give the same flow exactly when they differ by a common additive
constant. Thus face-coloring counts must be divided by `q`.

It follows, for every finite abelian group of order `q`, that the numbers
of nowhere-zero **boundary** flows in the two Heawood embeddings are

```text
torus:   (q-1)(q-2)(q-3)(q-4)(q-5)(q-6),
genus 3: (q-1)(q-2).                                     (5)
```

The first formula is zero for `q<7`, and first becomes positive at seven.
The second is already positive at three. Consequently the abstract graph
can carry small nowhere-zero flows even when one chosen embedding has no
corresponding face coloring. This is the exact obstruction to treating
surface flows as if all of them came from planar dual colors.

For a field of `q` elements the flow space has dimension eight. The face
boundary subspace has dimension `F-1`: its only face relation is the sum
of all oriented boundaries, as follows from connectedness of the dual.
The quotient therefore has dimensions `8-6=2` and `8-2=6`, respectively;
these are the surface homology dimensions `2g`. The same free-group
description gives quotient `A^(2g)` for arbitrary finite abelian `A`.
One can obtain this description by collapsing a primal spanning tree
and a disjoint dual spanning tree; the remaining `2g` surface generators
have no abelian relation on an orientable closed surface.
In this note the binary vector groups use only their additive structure.

The distinction is visible at four states. A nowhere-zero `F_2^2` flow
on a cubic graph has all three different nonzero values at every vertex,
so it is exactly a proper three-edge-coloring. Heawood has 48 such flows.
None is a boundary in the torus embedding; six are boundaries in the
genus-three embedding. The Four-Colour Theorem is a planar theorem, so
there is no contradiction; its primary published statement is
[Robertson--Sanders--Seymour--Thomas (1997)](https://doi.org/10.1006/jctb.1997.1750).

## 4. Exact flow census and the role of 5, 6, and 7

The experiment solves Kirchhoff equations along a spanning tree. Its eight
cotree edges are independent cycle coordinates; requiring them to be
nonzero and rejecting a zero tree edge gives every nowhere-zero flow
exactly once. Face-boundary membership is checked independently by row
reduction of the oriented face cycles.

| Additive group | All nowhere-zero Heawood flows | Torus boundaries | Genus-three boundaries |
|---|---:|---:|---:|
| `Z/3` | 2 | 0 | 2 |
| `F_2^2` | 48 | 0 | 6 |
| `Z/5` | 1692 | 0 | 12 |
| `F_2^3` | 765072 | 5040 | 42 |

The counts in the last two columns also follow from (5), giving a second
verification route. A separate permutation enumeration finds 24 perfect
matchings, each leaving a single 14-cycle. Choosing the first edge color
as a matching and alternating the other two around its complement gives
`24*2=48` Tait colorings, independently checking the four-state total.
At four states the torus census occupies six homology
classes with eight flows each; in genus three the zero class has six flows
and 42 nonzero classes have one flow each. More genus does not mean these
different boundary spaces are nested.

There is a separate, precise classical connection to the requested
`5--6--7` boundary. [Seymour's 1981 theorem](https://doi.org/10.1016/0095-8956(81)90058-7)
proves that every bridgeless graph has a nowhere-zero integer 6-flow.
Tutte asks for five instead. The current primary
[Esperet--Hendrey--Lagoutte--Marseloo--Norin--Steiner preprint, v4, 3 July 2026](https://arxiv.org/html/2512.17342v4),
section 1, still states this existence problem as a conjecture. Its
counterexamples concern **reconfiguration** of existing 5-flows, not their
existence. That is a different predicate.

The published abstract of
[Möller--Carstens--Brinkmann (1988)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190120208)
reports that a 5-flow counterexample minimal within a fixed-surface-genus
class has no face boundary shorter than seven, and gives the orientable
order bound `28(g-1)`. This short-face reduction is **CITED**, not proved
by the Euler calculation below or asserted for arbitrary graphs. The
abstract also reports the 5-flow result through orientable genus two.

Here is the elementary geometric implication with its own hypotheses.
For a connected cubic cellular map on a closed orientable genus-`g`
surface, if every face boundary has length at least seven, counting with
multiplicity gives `2E=3V` and `7F<=2E`. Therefore

```text
2-2g = V-E+F <= V-3V/2+3V/7 = -V/14,
so V <= 28(g-1).                                       (6)
```

For equal face size `ell`, the per-vertex combinatorial curvature is
`1-3/2+3/ell`. It is positive at `ell=5`, zero at six, and negative at
seven. The three roles must remain typed: five is the sought flow bound,
six is the proved general bound and separately the flat cubic face size,
and seven is the cited minimal-counterexample face threshold. The numerical
overlap does not identify flow magnitudes with polygon lengths.

The Heawood example is a decisive hostile to a stronger inference. It has
an embedding with 14-gonal faces and negative curvature, yet it already
admits a nowhere-zero 3-flow: orient every edge from points to lines and
put value `-2` on the perfect matching `p-t=1`, and `+1` on the other
two matchings `p-t=2,4`. At every vertex the signed sum is zero over the
integers. Its other embedding is toroidal. A
chosen negative-curvature embedding does not obstruct small graph flows.
The script also checks constant value two as a `Z/6` flow and constructs
a `Z/7` boundary flow from seven distinct torus face colors.

## 5. Reproduction, finite universe, and boundaries

Run from the repository root:

```text
python -X utf8 -B 04-computation/experiments/geometry_flow_obstructions_20261004.py
python -O -X utf8 -B 04-computation/experiments/geometry_flow_obstructions_20261004.py
```

The [script](../../04-computation/experiments/geometry_flow_obstructions_20261004.py)
and [output](geometry_flow_obstructions_20261004.out) use standard-library
integer arithmetic and checks that remain active under optimization.
Normal and optimized output agree. The explicit universe is:

* all 283 zero-sum boundary triples over `Z/2,...,Z/8`, `F_2^2`, and
  `F_2^3`, including zero-current hostiles;
* both displayed Heawood rotation systems, all face darts and dual edges;
* all `7!` point-to-line permutations for the independent perfect-matching
  and three-edge-color count;
* all `2^8`, `3^8`, and `4^8` nonzero cotree assignments for groups of
  orders three, four, and five;
* all `7^6` normalized assignments for `F_2^3`: the first two incident
  edge values are fixed to 1 and 2. Every nowhere-zero cubic star has an
  ordered independent pair, and `GL(3,2)` acts transitively on its 42
  possibilities, so multiplying the normalized count by 42 is exact;
* every nowhere-zero triangle lift of one Heawood witness per group:
  respectively `0,1,128,78125` lifts, independently reprojected to the
  original witness and checked for conservation and nonzero edges.

Neither the quotient theorem nor the surface counts prove any new flow
conjecture. They isolate two concrete failures of transport: nonzero
edges can acquire zero boundary, and a conserved flow can carry nonzero
surface homology. The retained circulation and homology coordinates are
what a lossless graph-to-geometry carrier needs.
