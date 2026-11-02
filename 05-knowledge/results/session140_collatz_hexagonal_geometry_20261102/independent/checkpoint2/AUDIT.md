# Session140 independent audit — checkpoint 2

**Outcome: PASS within the scopes below.** Checkpoint1 remains byte-frozen.
This separate audit accepts the elementary synthesis, proper Euclidean
refinement, and explicit hyperbolic control; it does not re-certify every
inherited theorem named in the synthesis. Exact local bytes are pinned in
`SOURCE_PINS.json`. A notation ambiguity was repaired before freeze: in the
site-conjugacy section C now explicitly denotes the standard map T defined
at the beginning of the synthesis. No mathematical correction was required.

## 1. First returns and exact-clock site addresses

Theorem S is correct, with its load-bearing **first** positive return
hypothesis. Positive integral return times have unbounded cumulative sums,
so every visit to the marked set occurs at a stated cumulative time. This
proves reachability and first-hit correspondence, least periods as sums of
route clocks, and eventual-periodicity equivalence. The extra hypothesis
that every ambient state eventually enters the marked set rules out hidden
orbits remaining forever outside it. Endpoint equations without first
returns fail on the stated two-state swap hostile. The canonical suspension
has the required first-return property; clocks are functions `tau:X->N+`,
not necessarily uniformly bounded constants.

The ring enumeration gives a genuine bijection from positive integers to
triangular-lattice sites: ring r has 6r sites, and the first r rings plus the
origin have `1+3r(r+1)` sites. The exact inverse radius is
`ceil((sqrt(12n-3)-3)/6)`. Conjugating standard Collatz by this bijection
retains original clock, indegrees, components, hits and periods, but arrows
are nonlocal jumps. Radial distance alone forces jumps at even n of length
at least `r(n)-r(n/2)`, asymptotic to
`(1-1/sqrt(2))*sqrt(n/3)`. Occupying every site is not a single orbit visiting
every site, not local mesh adjacency, and not intrinsic uniqueness.

The explicit bitlength coloring was independently proved and controlled in
checkpoint1; its use here remains valid. The distinction between this source
coloring and a coloring of every plane unit-distance pair is retained.

## 2. Other elementary synthesis principles

The quotient/fibre test, legality-saturated operation descent, and
operation-closed continuation congruence are valid as stated. The linear
boundary-rank criterion applies to all response vectors in its declared
vector space; legal geometric contexts must actually realize the response
vectors used to claim sharpness. Joint boundary extension counts compose
by summing products over the same assignment, not by multiplying marginals.
The first-hit target-fibre condition and nonnegative ranking equivalence
are sound and do not construct a Collatz ranking algorithm.

The integral repair example is correct: the triangle Laplacian repairs
`(1,-1,0)` over the rationals but an integral repair would require
`3(z0-z1)=2`. The arithmetic-jet induction is also correct. From
`A>=L`, `A+C±2D>=L` and `|x|<=1`, the stated decomposition bounds the new
energy below `(1+x²)L`. The lower-dimensional u and u±v coefficient vectors
remain in the required congruence class. This is not a proof that arbitrary
Pfaffian coefficients satisfy stronger-than-parity congruences.

The historical inventory is a selected connection inventory, not an
exhaustive novelty census. Inherited rigidity, tournament, Selmer, modularity,
and order-42 results keep their existing status. The audit does not upgrade
them merely because their summaries occur beside accepted elementary lemmas.

## 3. Proper induced unit-distance refinement

The refinement is valid. A generic rotation gives a cycle seed distinct
x-coordinates, allowing disjoint root intervals. Component strips with gaps
greater than one prevent all cross-component unit edges; choosing left
endpoints tending to infinity ensures bounded regions meet finitely many
strips, even when widths are unbounded.

Each finite set of children receives disjoint subintervals whose closures
lie in the parent's interval and within horizontal displacement 1/4.
Their available upward unit-circle arcs are open and nonempty. Only finitely
many old points and old-circle intersections are excluded at each stage,
so a legal point exists. Its vertical increment exceeds `sqrt(15)/4`.
Finite rooted levels and this uniform positive increment imply properness
of the vertex set. Unit edge length and local finiteness then imply that
compact regions meet only finitely many drawn edges as well. Crossings are
allowed: this is not a crossing-free straight-line graph embedding.

The exponential Collatz inverse family still excludes uniform positive
separation. This does not force finite accumulation: in a proper discrete
set, each finite subset has positive distance from all other vertices.
Choosing closer and closer pairs outside successively larger finite subsets
produces nonadjacent pairs escaping every compact set. The chosen upward
construction height is not necessarily decreasing along Collatz and is not
a convergence argument.

## 4. Hyperbolic curvature control and induced-distance corollary

The standard metric inputs were independently reopened in
[Anna Felikson's university-hosted notes](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf),
Section 6.2 and Theorem 6.22, printed pages 95–96: the half-plane distance
formula and boundary-orthogonal circle/vertical-line geodesics. All remaining
construction claims follow by elementary calculations, not an appeal to a
hyperbolic tessellation theorem.

For `z_(k,j)=((j+1/2)2^(-k),2^(-k))`, every parent–child distance has
`cosh d=21/16`, so `L=arcosh(21/16)`. Different levels give the height-ratio
bound. Adjacent-level horizontal offsets are odd multiples of h/4, with the
minimum exactly on parent–child pairs; same-level and nonadjacent-level
bounds are strictly greater. Thus all distinct vertex distances are at
least L, with equality **if and only if** the vertices are tree neighbors.

**Immediate accepted corollary:** this is an induced fixed-distance-L
realization of the abstract binary tree. Scaling all metric distances by
1/L gives an induced unit-distance realization, with curvature `-L²` rather
than `-1`. Properness, noncrossing and uniform separation survive scaling.

The endpoint-derived circle center is `a-s*11h/8` and radius squared is
`185h²/64`. Both endpoints lie on the same monotone side of the circle, so
the actual arc stays in its endpoint horizontal interval and height band
`[h/2,h]`. Separate parent intervals and separate open height bands prevent
nonincident crossings. This checks the noncrossing claim without assuming
that arbitrary hyperbolic arcs stay below endpoint heights. Compact subsets
of the open half-plane have a positive minimum height, hence meet finitely
many levels and edges. Boundary accumulation on `[0,1]×{0}` is outside the
hyperbolic plane.

The construction is **not quasi-isometric**: the two depth-k vertices
straddling the root split have tree distance 2k but constant ambient distance
`arcosh(3/2)`. It models only the abstract binary tree, not every standard
Collatz state or edge, and constructs no regular hyperbolic hexagonal mesh.
The positive contrast is metric capacity, not a solution of Collatz or HN.

## 5. Reproduction, portability, and stopping boundary

`verify.py` is a separate standard-library exact checker with explicit
exceptions. Normal and Python `-O` outputs agree. Controls comprise 8,001
hyperbolic vertex pairs and 7,875 edge pairs, sharp-distance equality iff
adjacency, generic endpoint-derived circle equations, 64 root-split hostile
levels, a separately constructed 31-vertex **rational** induced upward tree,
30,301 lattice addresses, 10,000 exact-clock site steps, and 2,625
three-variable arithmetic-jet checks. Finite evidence does not replace the
all-level arguments above or the infinite finite-exclusion construction.

The producer principles and curvature verifiers were additionally replayed
in normal and optimized modes; saved JSON pairs are identical. Portability
is tested in a temporary relocated packet tree with the same relative
paths. Preserve this directory as `independent/checkpoint2` alongside
`principles`, `curvature_control`, and `hadwiger_nelson_refinement`.
The portable selector includes no PNG or large artifact.

No Git or remote change was made by this audit. No assertion of Collatz
convergence, geometric uniqueness, all-site local filling, external HN
six-bound verification, or full inherited proof recertification is made.
