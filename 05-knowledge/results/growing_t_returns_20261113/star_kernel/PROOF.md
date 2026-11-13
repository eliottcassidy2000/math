# A reusable full-count compiler for selective-star repairs

**PROVED, scoped:** the exact newborn response on every unique-start star
face, its positive network representation, and the chain of minimizing
repairs within each face. **FINITE-EXACT:** the881 repair and its immediate
one-star continuation obstruction. A selector producing a unit-ready
source at every405-ray step is **REFUTED** by the incoming all-order k17
exclusion; see [the current synthesis](../README.md#exact-inverse-progress).
Constructors with a changed interface and general count coverage remain
**OPEN**. The kernel itself and all finite repairs below are unchanged.

The one-vertex quadratic and positive cut geometry are inherited from
[session124](../../session124_tournament_cut_response_20261031/README.md)
and [session136](../../session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md).
They are not new here. The additional ingredient is the ending-cover field
needed for a *second* birth and the complete restriction to unique-start
faces. This keeps the full Hamiltonian target and readiness simultaneously.

## 1. Supply the deletion tournament, then choose its star

Let U be a nonempty tournament. Add a vertex v, with v->j exactly when
x_j=1. Write P_x for this tournament. All arcs inside U stay fixed;
selectively changing x is a genuine star repair. Define on U:

* h=H(U), E_i and S_i the Hamiltonian endpoint and starting counts;
* F_i the number of ordered spanning two-path covers with first nonempty
  path ending at i, allowing the second path to be empty;
* Q_ij the number of ordered spanning pairs of nonempty paths with first
  path ending at i and second path starting at j; set Q_ii=0.

These are nonnegative integers. Cutting a Hamiltonian path after i proves
F_i>=h. The inherited pair-deletion inequality gives
Q_ij+Q_ji>=2H(U-{i,j})>=2 for distinct i,j.

Deleting v from a Hamiltonian path gives three disjoint cases: v was first,
last, or internal. Therefore the inherited one-ear formula is

    H(P_x)=sum_j x_j S_j + sum_i(1-x_i)E_i
            +sum_ij(1-x_i)x_j Q_ij.                         (1)

The other coordinates have especially simple forms:

    S_v(P_x)=sum_j x_j S_j,
    H(P_x-v)=h,
    F_v(P_x)=h+sum_i(1-x_i)F_i.                            (2)

For the last formula, either the first path is the singleton v, leaving
a Hamiltonian path on U, or delete v from that first path. What remains
is exactly a cover counted by F_i, with i->v. No second-path sector is
discarded.

Now add y with v->y and y beating every other vertex of P_x, marking y.
The inherited birth identity gives

    H_new=H(P_x)+F_v(P_x)-S_v(P_x),
    F_y-H_new=S_v(P_x).                                   (3)

Combining(1)--(3), the exact full target is

    N(x)=h+sum_i(1-x_i)(E_i+F_i)
              +sum_ij(1-x_i)x_j Q_ij.                    (4)

There is no bounded-defect assumption, and h may vary with the supplied U.
The empty-U singleton boundary is deliberately separate: its starting
count is1 rather than the empty sum in(2).

## 2. Every unit-ready star belongs to an explicit cube face

By(3), the final newborn is unit-ready iff sum_j x_j S_j=1. Thus choose
exactly one vertex w with S_w=1, force x_w=1, force x_j=0 at every other
vertex with S_j>0, and freely choose x_j at all zero-start vertices.
These faces are disjoint and exhaust every possible star of every width.
If no vertex has S_j=1, there is no such repair at this deletion graph.

This gives an exact target recognizer:

    choose a unique-start face; solve N(x)=prescribed target on that face;
    build P_x and its canonical newborn.

It is necessary and sufficient for this specified two-stage construction.
It is not a recognition theorem for all tournament constructors. Strongness
and a desired physical root interface can be imposed and checked separately;
the formulas and the finite tests retain weak candidates as well.

## 3. Positive cuts, minima, and the limits of optimization

Equation(4) is a directed cut with an external field. Put a source s and
sink t around the vertices of U. The source side means x_i=1. Add arcs

    s->i of capacity E_i+F_i,
    j->i of capacity Q_ij.

Then N(x)-h is precisely the cut capacity. Force x_w=1 and each other
positive-start vertex to0 with capacities greater than the sum of all
finite capacities. Hence, once the path kernel is supplied, the minimum
ready newborn count is obtained by one ordinary directed minimum-cut
problem per possible w. Computing the kernel is not claimed polynomial
in the tournament order.

For two masks A,B on the same face,

    N(A)+N(B)-N(A intersection B)-N(A union B)
      =sum_{i in A\B, j in B\A}(Q_ij+Q_ji).              (5)

The right side is strictly positive for incomparable masks. Intersection
and union stay on the face. Therefore all minimizing masks on one face
form a chain, of size at most z+1 when z is the number of zero-start
vertices. The linear F-field changes the minimizers and target levels,
but not this strict curvature. This is the scoped transfer from the older
ear-cut geometry; no assertion of interval coverage follows.

The construction also gives a cheap lower bound. On a unique-start face,
the inherited weighted-cut identity yields

    N(x)=2h-1+sum_{i:x_i=0}F_i+2c_b(x),
    N(x)>=(2+|{i:x_i=0}|)h-1,                            (6)

where c_b is the nonnegative weighted outgoing cut in U. This bound can
reject an oversized deletion graph before a target search. It does not
force h or the number of zero-start vertices to remain bounded in an
infinite construction.

Exact level membership is stronger than minimum-cut optimization. There
can be holes above the minimum, and a successful repair need not lie at
the minimum. Also, a minimum in one face need not compare with the minimum
in another face. The receipt tests these statements only in their declared
scopes.

In fact the inherited unique-start component grammar makes the unrestricted
minimum on each face still simpler. If a face exists, U=C->TT_a, where C
is its first strong component (possibly a singleton), and S_C(w)=1.
Indeed S_U(w)=S_C(w)H(U-C)=1 forces the suffix to have count1. Every vertex
of C has positive start count; all zero-start vertices are precisely that
transitive suffix. The unique minimizing mask sends v to *every* suffix
vertex. The completed newborn graph is then the fixed graph C+v+y followed
by TT_a, with unchanged count. Any other mask has some suffix vertex z->v,
and contains the triangle v->w->z->v. The odd-cycle packing formula retains
all packings on C+v+y and gains at least this positive triangle contribution;
the count strictly increases. Thus a generic minimum-cut solver is optional
for this unrestricted face problem. The positive network remains useful
for additional fixed-port constraints and for exposing the exact target
equation. This refinement prevents mistaking a cheap minimum computation
for the difficult missing target selector.

There is also a quantitative target filter. If r suffix vertices point to
v rather than away from it, then

    N(x)-N_min >= 2r (H(C)+H(C-w)).                       (7)

For each changed suffix vertex i, the completed graph contains both
triangles v->w->i->v and v->y->i->v. The complement of the first has count
H(C-w), since y is then a source and all remaining suffix vertices form
a dominated transitive tail. The complement of the second has count
H(C). Each triangle contributes twice that complement count to the
odd-cycle packing formula. These families are disjoint (their cycles
through v are different), and none occurs in the baseline graph. All
baseline packs persist. This proves(7), and bounds how many suffix ports
may point inward when a particular near-minimum target is prescribed.
The odd-cycle formula and its authorship are inherited via the positive
cut note above; this is an application, not a new proof of that formula.

## 4. The881 repair and the failure of its immediate repetition

The supplied A2 core has rows

    (518,388,456,1011,487,455,387,257,1,502), H=881.

Choose physical vertex0 and reverse its comparisons with2,7,9. The parent
becomes

    (130,388,457,1011,487,455,387,256,1,503),
    (H,S_0,H(-0),F_0)=(153,1,27,1539).

The actual newborn has rows

    (1154,388,457,1011,487,455,387,256,1,503,1022),
    (H,S_10,H(-10),F_10)=(1691,152,153,1692).

The full target is153+1539-1=1691 and the readiness is1692-1691=1.
The input881 count is intentionally lowered by the star repair; this is
how it escapes the pure-birth obstruction. The desired final count is
preserved as the explicit target, not at every intermediate operation.

Across all ten possible ports of A2 there are exactly12 unique-start star
choices. Exactly one has newborn target1691: the three-arc repair above.
Across all eleven ports of its literal1691 child there are exactly6 such
choices, and none has newborn target3311. These are exhaustive statements
at every star width for the two specified inputs, including weak parents.
They do not exclude two successive star repairs, changed deletion graphs,
extra vertices, or a different1691 representative.

The separate [repair family](../repairs/PROOF.md) varies an actual module
and gives the run-kernel response. Its existence is not a repeatability
claim; the exact target equation still has to be solved at each step.

## 5. Connection contract and checks

Source: a deletion tournament U with its actual endpoint and path-cover
kernel. Target: integer cuts in a positive capacity network, restricted
to explicit unique-start cube faces. Map: equations(1)--(4). Preserved:
every repaired parent count, newborn count, and newborn readiness.
Discarded: most graph adjacency and future responses after another core
change. Sidecar: U itself, chosen physical port, and the complete cut mask.
Cheapest hostile: the1691 child has no one-star continuation to3311 despite
the successful881 repair. This is a concrete geometry of feasible integer
choices, not an identification with geometric number theory or a claim
about rational/transcendental dichotomies.

Run `python -B .../star_kernel/check.py`, then with `-O --output receipt.optimized.json`.
The explicit universe is all75 labelled nonempty deletion tournaments
through order4, all1098 stars and their newborns, and all135 minimizing
faces. Separate endpoint counts of the literal parent and child verify
the formulas. The larger named repairs are checked independently at every
feasible star. No strongness filter is inherited silently.
