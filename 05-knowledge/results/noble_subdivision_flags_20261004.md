# Noble polyhedra, retained flags, and two different triangle refinements

**2026-10-04. Status: PROVED elementary mechanisms; FINITE-EXACT examples;
CITED geometric definitions.** No new classification or tournament claim.

The useful connection is an equivariant change of representation, with an
explicit warning: preserving every old symmetry does not preserve transitivity
on the newly enlarged set. Barycentric refinement remembers incidence flags;
edgewise refinement remembers a simplex lattice and its scale. They are
different operations even when their numbers of triangles agree.

Companions: [script](../../04-computation/experiments/noble_subdivision_flags_20261004.py),
[output](noble_subdivision_flags_20261004.out),
[mixed refinement and faithful overlay](mixed_triangle_refinement_20261004.md).

## 1. Inheritance and scope

The closest recovered mechanisms are
[the corrected Hill model study](noble_polyhedra_hill_fissary_20261002.md),
[typed medial reconstruction](medial_pair_tiles_20261004.md), and the
[affine-versus-metric dissection boundary](prime_shells_20260921_dissections.md).
MISTAKE-559 in `01-canon/MISTAKES.md` retracts stronger claims about the
Hamiltonian-path spectrum and unsupported number matches. Nothing here uses
those matches. Four fissary duals must not be counted as two of Hill's 146
isolated nobles, and coincident geometric points must not replace abstract
vertices. The least-used useful sidecar here is the **rank of an incidence**.

Hill's noble predicate is geometric vertex and face transitivity. His
nondegeneracy conditions require faithful realization, planar faces, and no
coplanar faces sharing an edge. Geometric symmetry is induced by Euclidean
isometries, a possible proper subgroup of abstract map automorphisms.
These definitions are [CITED: Hill, Definitions 2.2--2.4](https://arxiv.org/html/2607.28711v1).
The script recomputes the two selected model combinatorics; it does not rerun
the geometric classification or infer geometric symmetry from graph symmetry.

The concept board is: flags / orbit types / simplex lattice / metric /
duality / retained parent labels. Each operation below states which coordinates
it keeps and which it forgets. No intrinsic tournament relation is needed.

## 2. Flag refinement preserves the group and destroys vertex transitivity

**PROVED.** Let M be a finite connected closed polyhedral map. Assume its graph
has no loops or multiple edges, each face boundary is a simple polygon of at
least three sides, every edge has two distinct incident faces, and each vertex
link is a cycle of length at least three.
Use abstract cell identities even if a geometric realization self-intersects.
Write S(M) for the order complex of its proper cells: vertices are old vertices,
old edges, and old faces, with ranks 0, 1, and 2. Triangles are flags v < e < f.

If M has counts (V,E,F), then

\[
 (V_S,E_S,F_S)=(V+E+F,6E,4E).
\]

There are 2E comparable pairs of each rank pair 01, 12, and 02, and 4E flags.
Hence the Euler characteristic is unchanged. The new degrees are

\[
 \deg_S(v)=2\deg_M(v),\qquad \deg_S(e)=4,\qquad
 \deg_S(f)=2|f|.
\]

Consequently **S(M) is never vertex-transitive, even for its full uncolored
abstract automorphism group**: edge vertices have degree 4; both other ranks
have degree at least 6. This is an obstruction to transitivity, not a failure
of symmetry transport.

Every automorphism of M acts faithfully on S(M), and every rank-preserving
automorphism of S(M) reconstructs an automorphism of M. Thus

\[
 \operatorname{Aut}_{\rm rank}(S(M))\cong\operatorname{Aut}(M).
\]

There is a useful stronger reconstruction. The degree-4 vertices intrinsically
identify rank 1. Delete them. The remaining graph is the radial incidence graph
between old vertices and faces; it is connected and bipartite. Connectedness
follows by replacing any old edge-path with vertex--incident-face--vertex
paths and then adjoining every face. Its bipartition is unique up to exchange.
An uncolored automorphism therefore either preserves ranks 0 and 2 globally
or exchanges them globally. The latter is exactly a duality of M. Conversely,
every duality induces such an automorphism. The order complex is a flag
simplicial complex, so there is no distinction here between its graph
automorphisms and simplicial automorphisms. Therefore

\[
 [\operatorname{Aut}(S(M)):\operatorname{Aut}(M)]
 =\begin{cases}2& M\text{ is self-dual},\\1&\text{otherwise}.\end{cases}
\]

For any specified inherited group G, triangular-face orbits of S(M) are
exactly flag orbits of M. Its vertex-orbit count is the sum of the three old
cell-orbit counts. In particular, if G is vertex- and face-transitive with k
edge orbits, S(M) has k+2 vertex orbits. Noble does not imply flag-transitive.
For a connected map, an automorphism fixing a flag fixes its three adjacent
flags and then every flag, so the flag action is free. Hence its number of
flag orbits is 4E/|G|.

**FINITE-EXACT controls.** Counts and map groups are independently reconstructed
from literal (vertex, edge, face) incidences; small refined graphs are also
enumerated by a generic graph-isomorphism implementation.

| Original map | S(M) counts | Refined degree multiplicities | Rank-preserving / full abstract group |
|---|---:|---|---:|
| tetrahedron | (14,36,24) | degree 4: 6; degree 6: 8 | 24 / 48 |
| cube | (26,72,48) | degree 4: 12; degree 6: 8; degree 8: 6 | 48 / 48 |
| octahedron | (26,72,48) | degree 4: 12; degree 6: 8; degree 8: 6 | 48 / 48 |
| Hill D-4 | (200,720,480) | degree 4: 120; degree 8: 60; degree 24: 20 | 240 / 240 |
| Hill D-5 | (170,540,360) | degree 4: 90; degree 6: 60; degree 18: 20 | 120 / 120 |

The full D-4/D-5 entries follow from the proved reconstruction and the computed
old groups, not from enumerating all their refined-graph permutations. Both
have V=20 and F=60, precluding a duality. Their old edge-orbit sizes under the
full abstract group are (60,60) and (30,60), respectively. Their refined
triangles have **exactly 2 and 3 orbits under the full abstract group**.
They therefore supply noble inputs whose flag refinements fail both vertex
and face transitivity. For the regular tetrahedron placed at
(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1), an exact rational squared-distance
test retains only 24 of the 48 abstract refined automorphisms: the extra
dualities exchange vertices and face centers at different radii.

Literal barycenters of a planar face place adjacent small triangles in the
same plane, whenever these triangles are nondegenerate. Thus even the regular
tetrahedron's flat refinement fails Hill's adjacent-face condition. A projection
or bending procedure is an additional geometric operation; none is implicit
in these abstract group statements.

The typed-medial statement in the earlier note is the parallel operation:
old edges become vertices, while vertex-star and face-boundary cells retain
their two types. Forgetting their types can add dualities. The present flag
proof makes the same information loss explicit in the triangle subdivision.
For primary context see
[Hubard et al., *Medial symmetry type graphs*](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v20i3p29/pdf/).

## 3. Lattice refinement: congruence is weaker than transitivity

**PROVED.** The degree-r edgewise triangle has vertices

\[
 T_r=\{(i,j,k)\in\mathbb Z_{\ge0}^3:i+j+k=r\},
\]

with root edges differing by e_i-e_j. Its upward triangles are
{b+e_1,b+e_2,b+e_3} for sum b=r-1; downward triangles are
{b+e_1+e_2,b+e_1+e_3,b+e_2+e_3} for sum b=r-2. Empty negative-sum sets contribute
nothing. Counting gives

\[
 (V_r,E_r,F_r)=\left(\binom{r+2}{2},\frac{3r(r+1)}2,r^2\right),
 \qquad |\partial E_r|=3r.
\]

Coordinate permutations S_3 preserve cell kind and act on b by permutation.
Thus the number of triangular-face orbits under the full symmetry group of an
equilateral parent triangle is

\[
 p_{\le3}(r-1)+p_{\le3}(r-2),
\]

where p counts partitions into at most three parts and p(negative)=0.
For r=1,...,12 the exact counts are
1,2,3,5,7,9,12,15,18,22,26,30. Already r=2 has two orbits of sizes 3 and 1:
the corner cells and the central cell. All four triangles are congruent; no
symmetry preserving the whole parent triangle carries a corner to its center.

For a closed triangular map the compatible facewise refinement has

\[
 V_r=V+(r-1)E+\frac{(r-1)(r-2)}2F,\quad E_r=r^2E,\quad F_r=r^2F.
\]

An old vertex retains its old valence; every new vertex has valence 6. For
tetrahedral, octahedral, and icosahedral maps, whose old valences are 3, 4, and
5, this prevents vertex transitivity for every r>=2, independently of any
metric. The script checks the first two map families for r=1,...,12. The
icosahedral statement follows from the general local-degree proof.

Subdivision at degrees r and then s equals subdivision at degree rs in this
triangle lattice. Within each upward or downward cell, scaled barycentric
coordinates produce exactly the finer parallel-line grid at denominator rs;
their root directions agree, so both the vertices and the little triangles
agree. This is an exact combinatorial identity, with all 36 pairs r,s=1,...,6
also checked as sets of cells. The broader simplex construction is classical:
[Edelsbrunner--Grayson, *Edgewise Subdivision of a Simplex*](https://www.graysonfamily.org/dan/Papers/p24-edelsbrunner.pdf).
No new subdivision construction is claimed.

The cumulative chart C(i,j,k)=(i,i+j,r) sends this grid to 0<=u<=v<=r. Root
directions become horizontal, vertical, and simultaneous diagonal steps.
It transports incidence and the S_3 action by conjugation. Euclidean symmetry
requires the transported metric; coordinate integrality alone supplies none.

## 4. Cheap hostiles and the information that repairs them

Barycentric refinement of one triangle has 7 vertices and 6 triangles;
degree-2 edgewise refinement has 6 vertices and 4 triangles. More generally,
S^k has 6^k triangles and 2^k segments on each marked original side. Equality
with a degree-r edgewise refinement would force r=2^k from the boundary and
r^2=6^k from the faces, impossible for k>0. Even the coincident face count 36
does not help: S^2 has (V,E,F)=(25,60,36) and boundary length 12, whereas
degree 6 has (28,63,36) and boundary length 18.

Putting the barycentric vertices in a common lattice does not repair adjacency.
At denominator 6 an edge midpoint is (3,3,0), the centroid is (2,2,2), and
their joining median has direction (-1,-1,2). No triangular-lattice edge is
parallel to it: a root direction always has one zero coordinate. Arbitrarily
fine edgewise subdivision therefore cannot turn that median into a straight
path in its one-skeleton. This does not prohibit a further common refinement.

The [mixed-refinement companion](mixed_triangle_refinement_20261004.md)
addresses the positive repair: retain both parents in an intersection-complex
overlay and triangulate it. Its one-level relation and its hostile to automatic
higher-level iteration must not be upgraded to commutation of S and E_r.
The combinatorial cure is a common refinement with parent labels; the metric
and global group action remain additional data.

For an information-bearing object, useful storage is therefore
**(cell, rank/type, parent cell, scale, group action)**. The first four fields
support reconstruction and refinements; the last specifies which transitivity
claim is being tested. Quotienting to shape, triangle count, or an unlabeled
orbit can erase exactly the coordinate needed for the next operation. The
cheapest check is to count cell orbits or vertex degrees before discarding it.

## 5. Reproduction, independent paths, and limits

Run from the repository root:

```text
python 04-computation/experiments/noble_subdivision_flags_20261004.py
python -O 04-computation/experiments/noble_subdivision_flags_20261004.py
```

The output is identical with optimization; checks use explicit exceptions.
An independent agent audited the degree obstruction and the reconstruction
up to duality, including its closed-map assumptions; no correction was needed.
Dependencies are Python 3 and NetworkX. The finite universe consists of three
regular maps, two pinned noble models, local degrees 1--12, tetrahedral and
octahedral surface degrees 1--12, and 36 composition pairs. Positive controls
are exact group reconstruction and multiplicative edgewise composition;
hostiles are orbit splitting, degree splitting, abstract/geometric duality,
equal face counts, and a non-root median direction.

For small maps, color-preserving flag propagation is compared with generic
graph automorphisms. The triangular graph is independently regenerated from
all root-difference pairs. The geometric tetrahedron test uses exact rational
distances. Hill's data are not vendored: the two OFF files are read from a
temporary cache or fetched at pinned commit
`a801da7582fa927ba07e0af74d7c9c39445b7264`; Git blob hashes are checked against
`eac135ab499faf4437ca05e7f9d05be2d954e5dc` (D-4) and
`5978dbe15987115b8264c871d4203cfa3e5c8cd1` (D-5). The model repository is
[Plasmath/noble-tools-revised](https://github.com/Plasmath/noble-tools-revised/tree/a801da7582fa927ba07e0af74d7c9c39445b7264),
licensed GPL-3.0. This script derives incidences from those files and does not
import the older experiment's implementation.

The proved universal statements use the stated closed-map assumptions.
Boundary maps, repeated incidences, degree-2 links, collapsed geometric cells,
and untyped naked medial graphs require separate treatment. Nothing here
establishes a Collatz coverage theorem, a Hamiltonian-path spectrum completion,
or an identification of tournament classes with a triangle lattice.
