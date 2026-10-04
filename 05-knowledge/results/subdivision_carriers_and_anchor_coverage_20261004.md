# Retained subdivision addresses and Fourier-separated route families

2026-10-04 (America/Denver).

**PROVED:** the scoped geometric reconstruction/refinement mechanisms and
the all-period repeated-chart separation and density-existence statements
linked below. **FINITE-EXACT:** the declared model, chart and route censuses.
**OPEN:** universal positive Collatz convergence, a complete integer-to-golden
passport construction, and LRC(14). The work supplies useful representations
and certified families, not those global conclusions.

## Inheritance, portfolio, and live board

The anchor was the recommended next step from
[golden phase trees](difference_diagonals_golden_phase_recursion_20261004.md):
measure what the 48 rational charts actually add, retaining endpoint proofs.
The niche was the owner's cumulative-coordinate triangle rule; the wildcard
was the connection to noble polyhedra and the six-pair octahedron. The niche
produced a precise symmetry repair, while the anchor produced a general
Fourier separation theorem beyond its initial finite atlas.

Closest proved mechanisms: the source-preserving rational return compiler,
the exact golden phase lift, and typed medial/flag reconstruction. Canonical
hostiles: equal denominator ideals with different arithmetic realizations;
equal mesh counts with different incidences; and a legal route whose endpoint
has not yet repaid the original source. Corrected near misses: independent
absolute signs in a cumulative lattice, blindly iterating barycentric
containment, and confusing summed cylinder measure with natural density.
Least-used sidecars: the exact first-descent clock, local simplex order,
and a pair of retained parent cells.

The six live concepts were: source cells / first-descent clocks / completed
suffixes / cumulative addresses / incidence flags / coordinate symmetry.
The decisive comparisons and their limits are recorded in the companions:

1. [Added rational-chart coverage and completed suffixes](anchor_coverage_refinement_20261004.md).
2. [All-period chart separation by a Fourier character](periodic_chart_separation_20261004.md).
3. [Cumulative edgewise coordinates and exact decoding](edgewise_cumulative_coordinates_20261004.md).
4. [Mixed triangle refinements and retained-parent overlays](mixed_triangle_refinement_20261004.md).
5. [Noble polyhedra and reconstructible flag subdivision](noble_subdivision_flags_20261004.md).
6. [Tetrahedral pair charts and their symmetric completion](tetrahedral_pair_refinement_20261004.md).

## 1. The supplied coordinate rule has a meaningful sign boundary

On the degree-r triangle, x+y+z=r and x,y,z are nonnegative integers.
The transform C(x,y,z)=(x,x+y,x+y+z) is an invertible integer shear.
The standard edge rule, excluding equal points, is

    C(v)-C(u) in {0,1}^3 OR C(u)-C(v) in {0,1}^3.

The sign is common to the entire vector. This is the condition in
[Athanasiadis, section 4](https://arxiv.org/pdf/1310.0521), independently
derived here for the triangle. It yields root steps e_i-e_j and counts

    (V,E,F)=((r+1)(r+2)/2, 3r(r+1)/2, r^2).

Independent absolute signs instead admit (0,2,0)--(1,0,1) at r=2; its
cumulative difference is (1,-1,0). That chord crosses a legitimate edge.
This is the smallest failure. Only four of the eight displayed absolute
binary triples are available on the fixed plane, because the final
cumulative difference is always zero.

The corrected graph retains more than its picture suggests. Its corner
distances recover every coordinate: x_i=r-dist(v,r e_i). Thus the graph
plus an ordering of its three corners reconstructs all integer addresses.
The ambiguity without corner labels is exactly S3. This is a concrete
reconstruction theorem, not a claim that all unlabelled graphs remember
their original arithmetic.

## 2. Scale and recursive flags remain different coordinates

Edgewise degrees multiply: E_s E_r=E_(rs), with compatible child ordering
in higher dimensions. Barycentric depth adds: S^j S^k=S^(j+k). On a
triangle, E_r gives r^2 faces and S^k gives 6^k; their boundary lengths are
3r and 3*2^k. Both can be recursive.

The two mixed orders E_r S and S E_r even share their entire f-vector,
yet differ at r=2. Each has (19,42,24), but S E_2 has three degree-seven
vertices and E_2 S has none. Counting cells therefore cannot recover
the order of operations.

A positive theorem survives: S E_r refines both E_r and S on a triangle.
The proof uses the affine coordinate-permutation group. Barycentric
subdivision turns every fixed median into a subcomplex, dividing the
small triangles among the original six flag chambers. This does not
iterate automatically; an exact S^2 E_2 versus S^2 hostile is saved.

For finite geometric triangulations of the same simplex, intersect their
cells and subdivide the
intersection complex, retaining both parent labels and coordinate maps.
Overlaying every coordinate-permuted version first extends the symmetry
repair to every finite dimension. This construction makes two different
representations jointly available; it does not declare them identical.

## 3. The six-pair octahedron reappears inside the tetrahedron

Degree-two refinement of a tetrahedron has four original vertices and
six edge midpoints. Those six points are the six unordered pairs of four
labels, exactly the octahedral pair graph from the previous tile work.
Its three opposite pairs are the three perfect matchings of four labels.
An edgewise chart chooses one corresponding internal diagonal, reducing
the coordinate-permutation group from 24 elements to eight.

Insert the tetrahedron center. The central octahedron now has eight star
tetrahedra, with four corner tetrahedra outside it: 12 cells in total,
with f=(11,30,32,12). This is the minimal common refinement of the three
specified degree-two charts. It restores all 24 endpoint permutations.
The price is unequal volumes: four cells of volume fraction 1/8 and
eight of 1/16. The exact boundary assumptions and minimality proof are
in the companion note.

Use the three matching cuts

    u=x1+x2-x3-x4,
    v=x1+x3-x2-x4,
    w=x1+x4-x2-x3,

alongside s=sum(x_i). This is a Hadamard transform with determinant -16
and inverse divided by four. It sends the six midpoints to the three
positive/negative coordinate axes. The eight central star cells now
have the eight sign triples as addresses. Each old diagonal chart
forgets one sign; any two distinct charts recover all three. Three
charts have only eight consistent joint states out of 4^3 formal tuples.

This is also an exact Fourier representation: after labelling the four
endpoints by the Klein four group, the Hadamard rows are its four real
characters. The three pairing cuts are precisely its three nonconstant
Walsh coefficients. This spatial Fourier transform and the cyclic clock
character in section 5 act on different objects; both have explicit
decoders or annihilation predicates, rather than a shared numerical motif.

Here the user's binary cube has a precise second role: it labels eight
central cells. It is different from the incorrect independent-sign edge
test on the original two-dimensional triangle. The coordinate transform
also requires its mod-four lattice conditions; determinant 16 does not
identify this construction with an unrelated parameter -29/16.

## 4. Noble symmetry is transported but transitivity can fail

For the stated connected closed polyhedral maps with simple face boundaries
and vertex degree at least three, barycentric subdivision creates vertices
of three incidence ranks. Old edge centers have degree four; old vertices
and face centers have degree at least six. Hence the refined graph is never
vertex-transitive, even under its full abstract automorphism group.

Nevertheless it retains the original map. Degree four recognizes old
edges; removing them leaves the connected bipartite vertex-face incidence
graph, determining the other two ranks up to one global exchange. Thus
the only additional untyped symmetries are original-map dualities.
For the tetrahedron this distinguishes 24 typed automorphisms from 48
untyped ones; the regular tetrahedron's literal barycenter realization
has 24 geometric symmetries.

The exact Hill D-4 and D-5 controls have respectively two and three
orbits of refined triangular faces under their full abstract groups.
Their refinements fail both transitivity properties. Literal planar
refinement also makes adjacent triangles coplanar, which fails Hill's
geometric polyhedron convention. An additional bending or projection
operation would need its own proof.

Flag ranks, three matching cuts, and the earlier marked Zeckendorf colors
therefore play different roles. Each is useful because its labels retain
specific information. Reusing a three-color palette is not an identification
between their state spaces; the extra-unit provenance in the Zeckendorf
construction remains a separate retained coordinate.

## 5. First-descent clocks turn a finite atlas into a general theorem

Let w be a marked primitive positive valuation word of length p, with
every prefix multiplier 3^j/2^(a1+...+aj)>1. Its negative rational fixed
point is -h/d. Put Q=2^(a1+...+ap). The inherited compiler supplies a source cylinder for
every repetition m>=2; h<Q^2 proves that all such repetitions lie in its
domain. Every source in that cylinder has exact first descent at pm,
and all valuations before its final step equal the repeated nominal word.

**PROVED universal separation.** Distinct marked primitive words give
disjoint source cylinders at all repetitions m,n>=2. For a common source,
the first-descent times must agree, T=pm=qn. The two length-T repeated
words then agree except possibly at the last entry. Let zeta be a primitive
T-th root of unity. Each repeated word has Fourier coefficient zero at
zeta, since its repetition factor is a nontrivial geometric sum. Their
difference is supported in one position; its coefficient cannot vanish
unless that last difference is also zero. The two repeated words are
therefore equal and have the same marked primitive root.

This is a precise Fourier-to-route transfer. It needs only the vanishing
test for one complex coefficient; even its magnitude suffices for this
test. It does not contradict the earlier guard examples where retaining
only magnitudes loses the legal operation. A fixed cosine component alone
has a quarter-turn blind spot. The repetition condition is essential:
word (1,2) at m=1 and word (1) at m=2 both compile to 3 mod16.

The complete all-period atlas has a natural density, justified by a
uniform no-descent bound on the omitted periods, rather than by assuming
countable additivity. Exact enumeration through p=12 and m=30 gives

    0.11862865 < density < 0.11862913.

These are fractions of all positive integers. This is the mass of the
whole structured atlas, including inherited families, not new coverage
relative to the old bank. Its sources have certified first descent;
this statistic does not say that every source already has a saved route
to 1. Increasing the period cutoff cannot cover the source 7: its exact
first-descent valuation prefix (1,1,2,3) has no twice-repeated nominal
word agreeing through its third entry.

## 6. The initial 48-chart task is also completed, with a separate metric

Against the frozen 65-cylinder bank and five named infinite families,
the 48 negative rational charts contribute 47 disjoint all-height additions.
Their added natural density is 0.00033038280540161839..., with a strict
tail smaller than 2.087e-59 after m=30. This subtraction metric is different
from the universal-atlas mass above.

The strongest two in that 48-chart comparison are 1113 and 1122, centered
at -65/17 and -73/17. Their first added cylinders are 719 mod8192 and
6983 mod8192. A shared parameter -121/64 in the previous quadratic-norm
work did not identify their cycles; their distinct source sets now make
the lost word coordinate operationally visible.

For endpoint storage, 432 explicitly chosen sources all reach 1. Their
suffix graph uses 20,774 distinct edges instead of 32,027 separately
stored suffix edges, with every edge arithmetically checked. General hub
sealing yields 192 symbolic infinite families with already certified
suffixes. For chart1122 at repetition two, targeting 17 and reusing
17->13->5->1 lowers the least terminal exponent from 3102 to 153.
These sealed families have density zero for each fixed chart, repetition
and hub; their benefit is completed routes, not additional cylinder mass.

## Recommended next steps and stopping boundaries

1. Use nonperiodic heads to address structural misses such as 7. Repeating
   a larger period cannot fix the exact time-four obstruction. Keep a
   source-relative descent proof and compare added mass with the complete
   named baseline before ranking a proposed head.
2. Optimize certified hub choice jointly with terminal exponent and suffix
   cost. Preserve the actual integer at every shared DAG node; coincident
   residues or equal word lengths do not authorize merging routes.
3. Use the three tetrahedral charts as a small exact model of joint-view
   reconstruction. Test further quotients against its eight compatible
   sign states and retain lattice congruences. For higher degrees, the
   finite symmetry-orbit overlay is available but its size needs measuring.
4. Keep geometric, combinatorial, and metric symmetry predicates separate
   when projecting or bending a subdivided noble model. Current results
   settle the flat and abstract operations only.

These are new precise obligations, not promises that a representation
alone resolves the global conjecture. The scripts and independently
audited proofs in the six companions preserve the successful constructions
and the cheapest hostile witnesses for the routes that stopped.
