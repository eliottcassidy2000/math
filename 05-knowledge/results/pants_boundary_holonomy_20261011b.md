# Three cuffs, triangular seam counts, and the holonomy needed for gluing

**Status: PROVED elementary statements below; FINITE-EXACT companion controls.**
This package uses classical surface and graph constructions. It claims no new
classification of surfaces or laminations, no decision of plane 6/7-colourability,
and no additional Collatz ROOT coverage. Its contribution here is an explicit
marked interface and exact tests of the information lost by boundary-only data.

## 1. Inheritance and the three objects

The closest prior mechanism is the checkerboard gluing character in
[glued chessboard rings](glued_chessboard_rings_scaffolds_20261006.md), especially
its audited sections 5–7: a pair of pants occurs as a surface piece, but its
topological type alone does not determine how two colour classes glue.
[Geometry and flow obstructions](geometry_flow_obstructions_20261004.md) already
separates local face-boundary flows from global surface homology. The current
[boundary palette interface](boundary_palette_compression_20261011.md) retains
realizable equality partitions and their actual colouring lifts.

The type warning in
[THM-3354, inequivalent H1 carriers](../../01-canon/theorems/THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan.md)
applies: an integer arc count, a binary surface character, and a graph-colouring
constraint are different objects. We give the maps between specified versions
of them, rather than identifying objects because they have three ports.

A pair of pants is a sphere with the interiors of three disjoint disks removed.
We label its boundary circles, or cuffs, $A,B,C$. A genus-$g$ closed orientable
surface, $g\ge2$, decomposes into $2g-2$ pants along $3g-3$ cuffs. These standard
counts and conventions are recalled in [Hatcher, *Pants Decompositions of
Surfaces*, p. 1](https://pi.math.cornell.edu/~hatcher/Papers/pantsdecomp.pdf).
The standard presentation is
$\pi_1(P)=\langle a,b,c\mid abc=1\rangle$; see
[Goldman, *Trace Coordinates on Fricke Spaces of Some Simple Hyperbolic
Surfaces*, §3.2](https://arxiv.org/abs/0901.1404).
The arguments below use this elementary presentation, not a character-variety
classification.

## 2. A fixed three-band seam model gives triangular numbers

Fix three mutually disjoint standard bands of arcs, joining $AB$, $AC$, and
$BC$. Allow $x,y,z\ge0$ parallel copies in those bands, respectively. This is a
restricted family: no same-cuff returning arcs, boundary twists, or unmarked
changes of seam system are included.

The cuff endpoint counts are

$$
(a,b,c)=(x+y,x+z,y+z).
$$

The integer matrix of this map has determinant $-2$. Its inverse is

$$
(x,y,z)=\left(\frac{a+b-c}{2},\frac{a+c-b}{2},
                   \frac{b+c-a}{2}\right).
$$

Consequently a nonnegative integer endpoint triple is realizable in the fixed
model **if and only if** it satisfies the three triangle inequalities and
$a+b+c$ is even. This is both a decoder and a rejection certificate: the parity
condition alone is insufficient, and the inequalities alone are insufficient.
For example, $(1,1,1)$ satisfies the inequalities but has odd total. The triple
$(2,0,0)$ is outside the three-band image although an essential arc returning
to cuff $A$ can have those endpoints on a pair of pants. Thus this theorem
does not classify all properly embedded arc systems.

For exactly $N=x+y+z$ arcs there are

$$
\binom{N+2}{2}
$$

distributions. A lossless two-cut encoding is

$$
(x,y,z)\longmapsto (x+1,x+y+2)
       \in\{(i,j):1\le i<j\le N+2\}.
$$

Its inverse is $(i-1,j-i-1,N+2-j)$. This is the precise triangular-number
connection: lattice points of a two-dimensional simplex are in bijection with
two cut positions. The bijection supplies no orientation on pairs and hence no
tournament without an additional binary observable.

Reducing the endpoint equation modulo two gives $a+b+c=0$. One geometric
interpretation is to take two copies of the pants cut along the arcs and
cross-glue the sheets: crossing an arc toggles the sheet, and a cuff-parallel
loop detects the parity of endpoints on that cuff. This gives the local
binary-character interface used next; it forgets the full integer counts.

## 3. Boundary characters omit half the global binary data

Fix a connected closed orientable surface $\Sigma_g$, $g\ge2$, and a marked
pants decomposition. The surface gluing is fixed and uses orientation-reversing
cuff identifications. Let $\Gamma$ be its connected cubic dual multigraph:
vertices are pants and indexed edges are paired cuffs. Parallel edges and
loops are allowed; a loop has two incidences at its vertex.

We work only with rank-one sign local systems, equivalently characters
$\pi_1(\Sigma_g)\to\{\pm1\}$, written additively in $\mathbb F_2$.
On one pair of pants, $abc=1$ says that the three cuff bits have even sum.
Those three bits determine its local character. Compatibility across a glued
cuff makes its two bits equal, since inverse signs are equal in $\mathbb F_2$.
Therefore the global cuff vector $\epsilon\in\mathbb F_2^{E(\Gamma)}$ satisfies

$$
\partial\epsilon=0,
$$

the incidence-parity equation. Its allowed space is the graph cycle space
$Z_1(\Gamma;\mathbb F_2)$, of dimension

$$
|E|-|V|+1=(3g-3)-(2g-2)+1=g.
$$

These cuff bits are not a complete global character. Choose based paths and
local frames, and fix a reference identification on each seam for each
compatible local system. Relative seam identifications are bits $\tau_e$.
Changing the frame on pants $v$ by $s_v$ changes them by

$$
\tau_{uv}\longmapsto\tau_{uv}+s_u+s_v.
$$

Thus the independent seam classes form
$C^1(\Gamma;\mathbb F_2)/\delta C^0(\Gamma;\mathbb F_2)
=H^1(\Gamma;\mathbb F_2)$, also of dimension $g$.
For each fixed allowed cuff vector there are exactly $2^g$ gluing classes.
Without the reference choices, the fibre is a torsor for this group rather than
a preferred copy with a distinguished zero.

Equivalently, there is an exact sequence

$$
0\longrightarrow H^1(\Gamma;\mathbb F_2)
\longrightarrow H^1(\Sigma_g;\mathbb F_2)
\longrightarrow Z_1(\Gamma;\mathbb F_2)
\longrightarrow0. \tag{1}
$$

Here the last map records restrictions to the marked cuffs. To see surjectivity
and the kernel directly, apply van Kampen to the pants pieces, choose a spanning
tree of $\Gamma$, and retain one stable letter for every non-tree edge. In the
abelianization over $\mathbb F_2$, the pants relations are exactly the incidence
equations on the cuff variables. Edge conjugacies identify the two cuff
variables and impose no condition on the stable letters. The $g$ non-tree
stable-letter bits are free. This proves (1), its dimension $2g$, and the
description of its kernel. A splitting requires the marking choices.

These seam bits are discrete sign transport. They are not continuous metric
twist lengths, and varying them does not alter the previously fixed surface
gluing. Similarly, the fact that traces forget the affine carry in
[THM-4604, reducible Collatz characters](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md)
is a separate marked-representation loss, not an identification of its carry
with the group in (1).

### The smallest closed-pants hostile and a positive repair

Glue two pants along three cuffs. The dual graph is the theta multigraph,
with two vertices and three parallel edges, and the result has genus two.
The cuff vectors are $000,110,101,011$. For each there are four seam classes,
giving sixteen global characters.

In particular, cuff vector $000$ with seam vector $000$ and cuff vector $000$
with seam vector $001$ have identical cuff restrictions but different global
holonomy. The first character is trivial. The second detects a dual cycle.
The associated double cover is disconnected in the first case and connected
in the second, because any nontrivial character into $\mathbb F_2$ is
surjective. This is the smallest hostile among closed orientable pants
decompositions; no minimality claim is made for punctured or bordered surfaces.

For a lossless codec, retain the marked graph and a spanning tree $T$:

1. Record the cuff bits on the $g$ non-tree edges. Recover the tree bits
   uniquely by eliminating the incidence residual from leaves towards the root.
2. Normalize the seam vector to zero on $T$ by fixing the root frame and
   propagating its frame changes. Record the remaining $g$ seam bits.

The resulting two length-$g$ bitstrings recover the full gauge class. Recording
only the first string loses exactly $g$ bits. This does not discard every local
frame while still claiming a based path decoder: the graph, tree, edge indices,
and reference paths remain part of the interface.

## 4. A directly testable colouring version

Let $F_1,\ldots,F_m$ be connected bipartite graph fragments on marked vertex
labels. Supply an actual proper binary colouring $c_i$ of each fragment.
The graph to be coloured has exactly the union of the declared edges.
Assume their overlap graph is connected.

For overlapping fragments $i,j$, all shared vertices must give the same value

$$
\tau_{ij}=c_i(v)+c_j(v)\pmod2. \tag{2}
$$

If two shared ports disagree in (2), their local equality partitions are
incompatible. Otherwise place this bit on the overlap edge. A global binary
colouring exists **if and only if** every overlap-cycle has zero sum of
$\tau$. A spanning-tree normalization both decides this and exports a colouring.

Indeed, a connected bipartite graph has exactly two proper binary colourings,
differing by a constant flip $s_i$. Local colourings glue precisely when
$\tau_{ij}=s_i+s_j$. The tree equations determine the flips up to one common
flip; the non-tree equations are exactly the cycle tests. Agreement on shared
vertices then defines a global map, and every declared edge is proper because
it belongs to a fragment. Conversely any global colouring supplies such flips.
Thus the test is an equivalence, not merely a sufficient choice of local frames.

Three single-edge fragments forming $K_3$ are the minimal simple-graph hostile:
every piece is bipartite, but their overlap cycle has nonzero holonomy.
Four single-edge fragments forming $C_4$ pass, for all choices of their local
frames. Two connected pieces can also fail before the cycle test: the edge
$0$--$1$ and the path $0$--$2$--$1$ disagree on their two common ports.

This is a binary specialization of the prior boundary-palette mechanism.
It can check the two-colourable parts of an already supplied ternary partition.
It does not construct a partition of every plane unit-distance graph into
three bipartite parts. For a geometric union, actual coordinate identifications
and every additional cross-edge must be retained and tested; an abstract
overlap graph cannot certify that they are absent.

## 5. Exact implementation and scope

Companion:
[pants_boundary_holonomy_20261011b.py](../../04-computation/experiments/pants_boundary_holonomy_20261011b.py).
The public graph representation supports indexed parallel seams and loops;
loops contribute twice to incidence. All bit, vertex, and multiplicity inputs
require exact integers, rejecting boolean and float aliases.

The finite universe is:

- theta and separating-cuff genus-two duals, the tetrahedral genus-three dual,
  $K_{3,3}$ at genus four, and the cube dual at genus five;
- every cuff vector and seam vector on these five graphs, and every resulting
  two-string codec; explicit gauge-orbit comparison for both genus-two graphs
  and the tetrahedral graph;
- all seam layers $N=0,\ldots,16$, checking both inverse maps;
- all connected simple graphs on one to four labelled vertices with every
  binary edge labelling: 647 instances, compared with exhaustive vertex flips;
- actual triangle, square and multiple-port fragment controls, plus malformed
  typed inputs.

The surface theorem is proved above; enumeration checks its finite
presentation and codec. No geometric realization, scalar trace, colour count,
or local success is silently substituted for the actual global lift.

Reproduction:

~~~text
python -B -X utf8 04-computation/experiments/pants_boundary_holonomy_20261011b.py
python -B -O -X utf8 04-computation/experiments/pants_boundary_holonomy_20261011b.py
~~~

Both runs perform 14,857 exact checks; saved output is
[pants_boundary_holonomy_20261011b.out](pants_boundary_holonomy_20261011b.out).
