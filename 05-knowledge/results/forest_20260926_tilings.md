# Excursion forests: flat stars, boundary closure, and unbounded tiling addresses

**Status:** PROVED elementary statements and explicit constructions below;
FINITE-EXACT independent censuses; CITED production-system framework; OPEN
transport to Collatz. These are research results, not a proof of Collatz and
not a claim of a new classification of Euclidean tilings.

![Finite crops of two valid regular-polygon tilings: periodic sparse hexagons and a Sierpinski selection of sparse hexagons.](forest_20260926_tilings.png)

[Standalone vector figure](forest_20260926_tilings.svg). Both panels show
actual unit regular-polygon edges. The right panel is a crop of the infinite
selection, not a claim of a finite-k fractal tiling.

## Inheritance and the active board

Closest proved mechanism: regular-polygon angle deficit, recovered through
[HYP-2943, regular-solid/tiling recursion carrier](../hypotheses/HYP-2943-lrc14-regular-solid-tiling-recursion-carrier.md).
Its Gaussian/Eisenstein norm labels are separate from ordinary patch counts.
The corrected near miss is [HYP-3061, geometry-regime audit](../hypotheses/HYP-3061-lrc14-geometry-regime-archive-audit.md):
spherical, Euclidean, and hyperbolic curvature are different regimes; a small
count shared with tournaments is not a map between objects. The hostile
example here is the perfectly angle-balanced star `3.8.24`, which cannot
occur in an edge-to-edge tiling of the whole plane by regular polygons.
The least-used sidecar is the *word of neighbors around a face*, not the
angle sum around one vertex.

The board is: (1) local curvature, (2) cyclic boundary words, (3) refinement
rules, (4) global symmetry orbits, (5) unbounded hierarchical addresses,
(6) exact Collatz affine carry. Each new construction below is checked
against these six coordinates.

## 1. The counts 1,3,7,10 have an exact source

For a flat star with polygon side numbers `p_1,...,p_d`,

```
sum_i (p_i-2)/(2p_i) = 1,
sum_i 1/p_i = (d-2)/2.                                  (T1)
```

Every polygon angle is at least one sixth of a turn, so `3<=d<=6`.
Sorting the side numbers and solving (T1) gives:

| Valency d | Unordered multisets | Cyclic words up to rotation/reflection |
|---|---:|---:|
| 6 | 1 | 1 |
| 5 | 2 | 3 |
| 4 | 4 | 7 |
| 3 | 10 | 10 |
| Total | 17 | 21 |

The complete lists are in the [retained exact output](forest_20260926_tilings.out).
The computation has no guessed denominator cutoff: if `r` denominators
remain with reciprocal sum `s`, the next denominator is at most `r/s`.
This recursively exhausts the search. The resulting maximum denominator
is 42. A separate integer-angle search through denominators 3..42 then
checks the same census. Quotienting all permutations by the dihedral group
produces the 21 cyclic types.

This classical count is also recorded in
[Grünbaum--Shephard, *Tilings by regular polygons*](https://faculty.washington.edu/moishe/branko/BG108.Tilings.by.regular.polygons.pdf).
Here the arithmetic census is independently reproduced. The 21 are local
stars, not 21 globally realizable tilings and not 21 symmetry-orbit types.

For a regular map `{p,q}`, flatness reduces to `(p-2)(q-2)=4` and gives
`{3,6},{4,4},{6,3}`. Strictly positive curvature gives
`{3,3},{3,4},{4,3},{3,5},{5,3}`, the five Platonic cases.
These are exact classifications of the corresponding integer equations;
they do not identify the plane with a torus or identify polygon duality
with tournament complementation.

### The F=U+S comparison retains an important mismatch

The recovered [divisor-balance result, DB1](arithmetic_braids_20260917_divisors.md)
defines `F` as proper nontrivial divisors, `S` as those divisors that are
squarefree, and `U` as proper prime divisors. At `N=p^2 q r`, with distinct
primes, their sizes are `F=10,S=7,U=3`; the unit supplies a further `1`.
The mechanism is the exponent box `{0,1,2} x {0,1} x {0,1}`.

But `U subset S subset F`, whereas the four vertex-valency classes are
disjoint. A bijection carrying these membership predicates to the four
valency classes is therefore impossible. A tagged union of the divisor
sets has size 21, but the tags destroy the original inclusions. The shared
numbers merit investigation; they do not yet supply a structure-preserving
map. This is a concrete failed implication, not just a request for caution.

## 2. Six local stars fail a boundary-closure test

Every edge-to-edge tiling by regular polygons has a common side length
across adjacent faces. At a fixed face, record in cyclic order the polygon
on the other side of each edge. The complete local census forces:

| Locally flat star | Odd central face | Forced alternating neighbors |
|---|---:|---|
| 3.7.42 | 7 | 3,42 |
| 3.8.24 | 3 | 8,24 |
| 3.9.18 | 9 | 3,18 |
| 3.10.15 | 15 | 3,10 |
| 4.5.20 | 5 | 4,20 |
| 5.5.10 | 5 | 5,10 |

For example, any corner incident with both a triangle and an octagon must
have a 24-gon as the remaining face; every corner incident with a triangle
and a 24-gon similarly forces an octagon. Around the triangle, neighbors
would have to alternate `8,24,8,24,...`, which cannot close after three
edges. The other five rows have the same proof. The script checks that
*every* census entry containing each indicated face pair gives the forced
third face, so this is not assuming all vertices have the original type.

Thus angle zero does not imply global extension. The missing invariant is
boundary closure. This makes the analogy to an arithmetic source condition
precise at the level of *which information a local quotient forgets*.
It does not identify the geometric boundary word with a Collatz carry.

## 3. A small exact refinement graph

Splitting one angle into two polygon angles requires

```
1/q + 1/r - 1/p = 1/2.
```

With `3<=q<=r` and finite `p>=3`, the only solutions are

```
6 -> (3,3),   12 -> (3,4),   30 -> (3,5).
```

Proof: `q>=4` gives at most `1/2` for the first two reciprocals, so
`q=3`; then positivity forces `r=3,4,5`. A 30-gon never occurs in the
flat-star census. Allowing both orders of the `(3,4)` replacement on the
21 cyclic stars yields an exact graph with **14 directed edges**. Valency
increases by one, so it is acyclic. The graph preserves the local angle
equation and cyclic order. Its depth is at most three; it cannot by itself
encode arbitrarily deep Collatz excursions.

The `6 -> (3,3)` rule has a global geometric realization: dissect a regular
hexagon into six unit equilateral triangles. At each old hexagon corner,
the old 120-degree angle becomes two 60-degree angles. However, choosing
arbitrary local replacements independently is not a proof of a globally
compatible refinement. The boundary edge/neighbor word remains necessary.

## 4. Two local stars, unbounded k-uniformity: an explicit construction

Let `Lambda=Z e_1+Z e_2` be the unit triangular lattice, with the angle
between `e_1,e_2` equal to 60 degrees. Six unit triangles incident with a
lattice vertex form a regular hexagon. For any subset `A` of `Lambda` whose
distinct points are at least distance three apart (in particular any subset
of `3 Lambda`), replace each cluster centered at `a in A` by that hexagon. Distinct
clusters are separated, since their centers are at least distance three
and their circumradius is one. The result is an edge-to-edge tiling of
the whole plane by unit triangles and regular hexagons.

There are only two possible vertex stars:

```
3^6,              3^4.6.                               (T2)
```

The six vertices bordering a replacement have the second type; all other
remaining vertices have the first. Chosen centers themselves are removed.
This supplies a geometrically valid carrier for an arbitrary selection
address, with no appeal to an unproved extension of locally legal stars.

Now choose `A=N Lambda`, for an integer `N>=3`. The symmetry group is
exactly `N Lambda semidirect D_6`: every symmetry must preserve the hexagon
centers, and all these lattice symmetries preserve the construction.
The quotient of vertices by translations is `(Z/NZ)^2 \ {(0,0)}`.
Therefore the number of vertex orbits is

```
k_N = (N^2+6N-10+2 gcd(N,3)+gcd(N,2)^2)/12.             (T3)
```

This is a PROVED formula, independently checked by explicit orbit
enumeration for `3<=N<=60`. Burnside's lemma gives it directly: identity
fixes `N^2` sites; the rotations through 60 and 300 degrees fix one each;
120 and 240 degrees fix `gcd(N,3)` each; 180 degrees fixes `gcd(N,2)^2`;
each of the six reflections fixes `N`. Subtract the orbit of the removed
center. In coordinates, generators are

```
R(a,b)=(-b,a+b),     S(a,b)=(a+b,-b).
```

The first values for `N=3,...,10` are `2,3,4,6,7,9,11,13`.
In particular `k_N >= (N^2-1)/12` tends to infinity although the list of
local vertex figures remains exactly (T2). The parameter k counts symmetry
orbits, not distinct angle words. This distinction is essential when
"fractalizing k-uniform tilings."

## 5. Fractal selection with genuine unbounded address depth

Use the same valid hexagon replacement, now at centers

```
A = {3(i e_1+j e_2): i,j>=0,  i bitwise-AND j = 0}.
```

Let `Q` be the corresponding index set. It obeys the exact disjoint recursion

```
Q = 2Q + {(0,0),(1,0),(0,1)}.                            (T4)
```

Consequently `|Q intersect [0,2^r)^2|=3^r`; the normalized finite address
patterns give the familiar Sierpinski hierarchy. This is a hierarchy of
*hexagon-center selections*. The unit-sided tiling itself is not claimed
to be invariant under dilation, and its tile boundaries are not fractals.
Every finite stage and the entire infinite selection still define valid
plane tilings by the geometric argument in section 4.

This example has infinitely many vertex orbits: in the negative quadrant,
the distance to the nearest hexagon is unbounded, and that distance is
preserved by every tiling symmetry. More strongly, for every prescribed
finite observation radius, vertices sufficiently far into that quadrant
have identical triangular neighborhoods while this distance is unbounded.
Thus any finite family of bounded-radius local features can be fixed while
the needed global address information remains unbounded.

This answers the geometric construction question affirmatively, but it
does **not** retain finite k. Finite k and aperiodic hierarchical selection
are distinct restrictions. Periodic finite patterns give finite-k
approximations; their limiting hierarchy need not have finite k.

## 6. What transfers to the excursion forest

For a finite polygonal disk patch, put

```
K(v)=1-deg(v)/2+sum_{faces f incident with v} 1/|f|.
```

Summing counts vertices, edges, and faces exactly, so `sum_v K(v)=1`.
All flat interior vertices contribute zero. The entire value is carried
by the boundary; for a patch with h holes the sum is `1-h`. This is an
elementary discrete Gauss--Bonnet identity. It survives compatible patch
refinement, but it does not check the forced-neighbor closure of section 2.

The useful transfer contract is:

| Coordinate | Tiling carrier | Collatz forest requirement |
|---|---|---|
| Local equation | angle sum | affine block identity |
| Composition | glue ordered boundaries | compose affine matrices in actual orbit order |
| Missing data in a scalar | cyclic neighbor word | source residue and affine carry |
| Hierarchy | chosen sparse-center addresses | unbounded-depth excursion ancestry |
| Realization test | boundary closure and a valid plane patch | every intermediate state is the required positive integer |
| Global target | compatible whole-plane tiling | a well-founded descent or bounded total reset debt |

These rows preserve the *logic of composition with a boundary condition*.
No letter map identifying Euclidean curvature with Collatz drift has been
constructed. Root's affine curvature `(A-1)/B`, for `x -> A x+B`, is a
different invariant: its composition is a weighted average and descent
at source n requires `(A-1)/B < -1/n`. A zero polygon angle defect cannot
pay that source-dependent arithmetic threshold.

A relevant primary framework is
[Goodman-Strauss, *Regular Production Systems and Triangle Tilings*, sections 1--2](https://strauss.hosted.uark.edu/papers/TARP.pdf).
It retains oriented edge pairings and compatible strip-production relations;
not every locally allowed production extends to an orbit. The cited role
here is the boundary-language framework, not the paper's historical claims
about then-open decidability problems and not a theorem about Collatz.

**Strong next question (OPEN):** can an exact Collatz excursion certificate
have a finite set of local composition rules but an unbounded boundary
word, with a computable boundary obstruction that must eventually fire
for each fixed positive source? The six odd-face obstruction certificates
show what such a theorem would look like. The sparse-hexagon construction
shows why finitely many local states alone cannot supply it.

## Reproduction and scope

```
python -B 04-computation/experiments/forest_20260926_tilings.py
python -O -B 04-computation/experiments/forest_20260926_tilings.py
python -B 04-computation/experiments/forest_20260926_tilings.py --figure 05-knowledge/results/forest_20260926_tilings.svg
```

All checks use integer or rational arithmetic and active `require` calls.
Controls: two independent complete star censuses; all forced-pair boundary
extensions; all 14 refinement edges; the five/three regular-map equations;
58 direct D6 orbit enumerations against (T3); independent bitwise and
substitution constructions through address depth eight. The proofs above
give the unbounded statements; finite computations are corroboration.
The PNG preview was rendered from the SVG in headless Microsoft Edge at
1400 by 720 pixels and visually inspected. SVG generation is dependency-free.
