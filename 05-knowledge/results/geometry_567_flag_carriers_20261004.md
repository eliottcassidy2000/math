# The 5–6–7 curvature boundary needs a rotation carrier

2026-10-04. **PROVED** permutation and rotation-bit mechanisms below;
**FINITE-EXACT** for the exhaustive Heawood rotation/automorphism census and
the two matrix-group controls. The classical Klein-map identification is
**CITED** separately. No new classification or historical priority is claimed.
No Collatz convergence statement is transferred.

Artifacts: [script](../../04-computation/experiments/geometry_567_flag_carriers_20261004.py)
and [output](geometry_567_flag_carriers_20261004.out).

## 1. Recovered mechanisms and the overlooked coordinate

The prior [30/42 triangle-group discussion](ogg_triangular_chowla_20260927.md#3-thirty-and-forty-two-one-formula-proved-classical)
already identifies the spherical/Euclidean/hyperbolic transition through
the sign of (1/2+1/3+1/q-1), and identifies the Paley group (7:3)
as an index-eight subgroup of the Fano/Klein group of order168.
[THM-2996 — prime-modular affine-defect trichotomy](../../01-canon/theorems/THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness.md)
provides another exact triangle-group quotient and warns that a finite
quotient is not the spherical, Euclidean or hyperbolic group itself.

[THM-4529 — Petersen and Heawood families in Paley coordinates](../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md)
already proves the Fano-line constructions: replacing all seven Paley
triangles in (K_7) by stars gives the Heawood incidence graph; deleting one
point and suppressing its three resulting degree-two line vertices gives
the Petersen graph. It also refutes preservation of Hamiltonian-path parity
through these moves. The seven Petersen-family members are not an action
of the seven Paley translations.

The closest recent carrier result is
[noble subdivision with retained flags](noble_subdivision_flags_20261004.md):
incidence ranks matter, extra dualities appear after forgetting them, and
cell transitivity does not imply flag transitivity. MISTAKE-559 in
`01-canon/MISTAKES.md` already rejects the unsupported identification of
small tournament-spectrum numbers with polyhedral symmetry data.

Our board is now **graph / tournament / rotation / face / curvature / flow**.
The canonical hostile is the same Heawood graph carrying different surfaces.
The corrected near miss is “the same group of order21 or168 determines the
surface.” The least-used coordinate is the local cyclic order around each
vertex. Recovering it produces the following exact test.

## 2. The faithful map carrier is two permutations

For a finite graph, a dart is an ordered edge ((v,w)). Edge reversal is the
fixed-point-free involution $\alpha(v,w)=(w,v)$. A rotation system supplies
a permutation $\sigma$, cyclically ordering the darts at each vertex.
With our composition convention,

\[
 \psi=\sigma\circ\alpha
\]

is the face permutation. The cycles of $\sigma,\alpha,\psi$ count vertices,
edges and faces. Gluing oriented disks along the resulting face cycles gives
the corresponding closed orientable surface; connectedness of the graph
gives connectedness of this surface. This construction retains the surface
embedding combinatorially, without requiring coordinates in three-space.

For a cubic graph whose face cycles all have length (q), let (d) be the
number of darts. Then

\[
 V=d/3,\qquad E=d/2,\qquad F=d/q,\qquad
 \chi=d\left(\frac12+\frac13+\frac1q-1\right).
\]

There are $2d=4E$ flags $(v,e,f)$. A dart records an oriented edge;
a flag additionally selects an incident face-side. Thus the transition
at (q=5,6,7) is exact:

| face size (q), cubic valence | triangle excess | curvature sign |
|---:|---:|---|
|5| (1/30) | positive |
|6| (0) | zero |
|7| (-1/42) | negative |

Regular (q)-gons with angle $2\pi/3$ can be glued along the map edges:
spherical for (q=5), Euclidean for (q=6), hyperbolic for (q>6).
Three angles meet to give $2\pi$, so no vertex curvature defect remains.
This is an intrinsic surface metric, not an assertion of a faithful
Euclidean polyhedron with all abstract symmetries.

The existence of a tournament on the graph's labels supplies neither
$\sigma$ nor $\psi$. These must be specified and checked.

## 3. One relative rotation bit changes the Heawood surface

Work modulo seven. Put (D=(1,2,4)), use points (P_p), lines (L_t),
and declare incidence when $p-t\in D$. These are exactly the inherited
Paley/Fano coordinates. The graph has 14 vertices, 21 edges and 42 darts.
The Paley group consists of

\[
 P_p\mapsto P_{ap+b},\qquad L_t\mapsto L_{at+b},
 \qquad a\in\{1,2,4\},\quad b\in\mathbb Z/7.
\]

It has order21 and preserves the tournament relation on the point labels.
Choose the reference cyclic neighbor lists

\[
 P_p:\ (L_{p-1},L_{p-2},L_{p-4}),\qquad
 L_t:\ (P_{t+1},P_{t+2},P_{t+4}).
\]

At points advance this list by $\varepsilon_P\in\{+1,-1\}$; at lines
advance by $\varepsilon_L\in\{+1,-1\}$. All four resulting rotations
are Paley-equivariant. They are the only such rotations: the Paley group is
transitive on each of the two vertex parts, preserves each reference local
cycle, and each cubic star has only two cyclic orders.

**PROVED rotation-bit theorem.** Equal signs give seven hexagonal faces and
genus one. Opposite signs give three 14-gonal faces and genus three.

Here is the mechanism, rather than just the counts. Encode a point-origin
dart $P_p\to L_{p-d}$ by ((p,d)), with $d\in D$. Set

\[
 a=2^{\varepsilon_L}\pmod7,\qquad
 b=2^{\varepsilon_P}\pmod7.
\]

Two face steps are

\[
 \boxed{\psi^2(p,d)=(p+(a-1)d,\ ab\,d).}
\]

For equal signs, $ab\in\{2,4\}$ has order three and
$1+ab+(ab)^2=0\pmod7$. Thus every $\psi^2$-orbit has length three
and every face has length six. For opposite signs, (ab=1), so the
second coordinate is fixed and the first is translated by nonzero
((a-1)d); every orbit has length seven and every face has length14.
The even face length is forced by the graph's bipartition. Euler now gives

\[
 14-21+7=0,\qquad14-21+3=-4.
\]

The displayed face lists in the output also verify that every face boundary
is a simple cycle. The three long faces are Hamiltonian cycles. Their face
size is14, not seven: this genus-three Heawood map is not the Klein
({7,3}) map.

Global sign reversal reverses the chosen surface orientation and preserves
genus. Thus, within these four marked systems, one relative bit determines
the surface type; the remaining sign retains orientation. Without the
equivariance restriction, a complete rotation system needs all14 local bits.

There is also a positive **canonical tournament-to-map construction**.
Keep the point and line copies matched by their original tournament labels,
and at every star follow the directed three-cycle induced by the Paley
tournament on its neighbor labels. An out-neighborhood at $L_t$ follows
the displayed reference order, giving $\varepsilon_L=+1$. An
in-neighborhood at $P_p$ follows its reverse, since negation reverses
Paley edges, giving $\varepsilon_P=-1$. Thus this explicit rule selects
the genus-three map. The script checks all42 local successors against the
tournament relation. Reversing just the in-star convention instead selects
the torus. This rule uses the tournament and the matched labels; the bare
incidence graph does not retain that rule or its output rotation.

## 4. Same graph, different dual and symmetry carrier

The complete finite audit enumerates every permutation of the seven points
preserving the seven line sets. There are168. Every graph automorphism
either preserves both bipartite parts or swaps them. The explicit duality

\[
 P_p\leftrightarrow L_{-p}
\]

supplies all swaps, giving336 graph automorphisms. This is an exhaustive
enumeration, not an inference from a desired group order.

On darts, a graph automorphism (g) is orientation-preserving for a marked
map exactly when $g\sigma=\sigma g$, and orientation-reversing exactly
when $g\sigma=\sigma^{-1}g$. All336 candidates are tested:

| signs | faces | genus | preserving | reversing | dual graph |
|---|---|---:|---:|---:|---|
|((+,+)), ((-,-))| seven hexagons |1|42|0| (K_7) |
|((+,-)), ((-,+))| three 14-gons |3|21|21| seven parallel edges on each pair of (K_3) |

Both full map groups have order42 and are transitive on vertices, edges,
and faces. Neither is flag-transitive: there are84 flags and the flag
action is free, hence two flag orbits. The torus examples are orientably
regular and have no reversing map automorphism. The genus-three examples
have reversals, but their preserving subgroup has two dart orbits. These
claims concern intrinsic/combinatorial maps, not Hill's Euclidean geometric
definition of noble polyhedra.

The dual difference is exact and important. In the torus map each pair of
faces meets in one edge. In the genus-three map each pair meets in seven
edges. There are no dual loops in either case. Keeping only the graph or
the common Paley group loses this information completely.

For an independent finite scope check, every labelled rotation system on
this cubic graph was enumerated: exactly $2^{14}=16384$. Its genus census is

| genus |1|2|3|4|
|---|---:|---:|---:|---:|
| labelled systems |16|1008|10880|4480|

These are labelled, oriented rotation systems, not isomorphism classes.
Exactly four preserve the fixed Paley subgroup. The census's genus bounds
also have elementary controls: girth six bounds (F\le7), hence genus at
least one; (F\ge1) bounds genus at most four.

## 5. The 5–6–7 comparison with its complete flags

For (p=5,7), the script independently constructs $G=\mathrm{PSL}(2,p)$
as determinant-one matrices modulo ({I,-I}). Put

\[
 s=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
 t=\begin{pmatrix}0&-1\\1&1\end{pmatrix}.
\]

Their projective orders are two and three; (st) is projectively an upper
unipotent of order (p). On the dart set (G), use right multiplication
by (s) for (alpha), and by (t) for $\sigma$. The face permutation
is right multiplication by (st). The two generators generate all of
(G), as checked by exhaustive closure; alternatively their upper and
lower unipotents generate the determinant-one matrices. Left multiplication
acts transitively and freely on these darts and commutes with the two
structure permutations. The graph and face boundaries are separately
checked to be simple.

| carrier | (V,E,F) | darts / flags | genus | face size |
|---|---|---|---:|---:|
| $\mathrm{PSL}(2,5)$ | (20,30,12) | (60/120) |0|5|
| equal-sign Heawood | (14,21,7) | (42/84) |1|6|
| $\mathrm{PSL}(2,7)$ | (56,84,24) | (168/336) |3|7|

The first is the spherical dodecahedral carrier; the last is the standard
Klein-map carrier. **CITED context:** Scholl–Schürmann–Wills describe the
Klein quartic's dual maps with24 heptagons/56 vertices and56 triangles/24
vertices, and distinguish their full combinatorial symmetries from the
smaller symmetry groups of three-dimensional realizations; see
[their primary geometric account, “Basic facts”](https://math.ucr.edu/home/baez/klein_quartic_scholl.pdf).
The present verification uses the explicit finite permutations, not a
numerical identification with the quartic equation.

For comparison, the opposite-sign Heawood map has (d=42,q=14), so its
excess is (-2/21) and (chi=42(-2/21)=-4). Equal genus does not identify
its graph, map, triangle signature, symmetry action, or complex structure
with the Klein carrier. The distinction is visible before any sophisticated
invariant: its graph has14 vertices, not56.

## 6. Useful transfer and stopping boundary

| source to target | map and preserved data | loss and restoration |
|---|---|---|
| Paley tournament to a marked Heawood map | incidence stars plus induced neighbor three-cycles select $(-1,+1)$, genus3 | keep the matched point/line labels and local-cycle rule; the bare graph forgets them |
| marked rotation map to bare graph | forget $\sigma$, retain darts and reversal | face cycles, genus, dual adjacency and admissible map automorphisms disappear |
| map to flag carrier | retain vertex/edge/face incidence and orientation conventions | no combinatorial loss; geometric coordinates remain separate |
| group quotient to surface map | retain the marked generators (s,t), acting on darts | abstract group order alone loses the cycle partitions |

The graph's cycle/flow space is unchanged by the rotation bit, whereas the
cellular boundary space changes with the faces. For any field, the connected
Heawood graph has cycle dimension (E-V+1=8); the face-boundary dimensions
are (F-1=6) and2. The remaining first-homology dimensions are therefore2
and6, respectively. This supplies an exact interface to the
[companion flow audit](geometry_flow_obstructions_20261004.md): transporting
a planar face-coloring/flow argument requires the
rotation and homology sector, not only the common graph. It also explains
why a group-symmetry match does not select a curvature regime.

There is no map here from a positive integer or a Collatz certificate to a
chosen rotation system preserving an arithmetic return predicate. The
four-system example is a controlled obstruction to forcing such a transfer,
and a concrete storage prescription when surface data really is required.

## Reproduction and scope

```text
python -X utf8 -B 04-computation/experiments/geometry_567_flag_carriers_20261004.py
python -O -X utf8 -B 04-computation/experiments/geometry_567_flag_carriers_20261004.py
```

The script checks all336 graph automorphisms, all16384 local rotations,
all four Paley-invariant systems, their actual face lists and dual
multiplicities, vertex/edge/face orbits, and both complete matrix-group
dart carriers. It uses exact integers and rational triangle excesses.
All checks survive optimization. The positive controls are the sphere,
torus and genus-three regular carriers; the hostile is the same Heawood
graph with a changed relative rotation bit. No wider census of surfaces,
tournaments, or Collatz routes is asserted.
