# Conditional indexed Hamiltonian observers on native hierarchical cells

23 October 2026. **PROVED scoped**: a source-specific compiler for partial
indexed path observations in the mixed-depth ternary family. **FINITE-EXACT**:
the independent actual-tournament controls in `check.py`. No unobserved
kernel lift is selected, no arbitrary Hamiltonian path is assumed to lie in
the tested orbit, and no electrical invariant is replaced by a group count.

## Inheritance and the new object

The incoming [mixed-depth/shared-flags theorem](../../session42_shared_lattices_square_colors_20261023/hierarchical-kernels/MIXED_DEPTH_AND_SHARED_FLAGS.md)
at `8b71b746f7` proves the group order, actual-owner base, and complete-path
arc-flag cost for named hierarchical cells. Our earlier
[kernel factorization](../../attachment_lattices_square_sum_20261023/KERNELS.md)
and [wide-move/root-state bridge](../../attachment_lattices_square_sum_20261023/MOVES.md)
retain independent kernel freedoms. These are inherited mechanisms.

The canonical hostile is the incoming actual regular order-63 source: its
15 required owner slots need nine indexed arc flags, not eight. The corrected
near miss is a quotient that retains only the number of lifts: it can count
them but cannot identify an actual path owner. The least-used sidecar here
is the precise set of *input-prefix coordinates exposed by an observation*.
The live concepts are shared coordinates, partial pins, conditional observer
fibres, actual arc availability, and equivariant selection.

The new target is conditional and source-specific: after receiving arbitrary
consistent indexed owner/arc answers, determine exactly which further path
observations remain possible, how many lifts they retain, and the minimum
additional flags needed for a specified subset of owners. Whole-frame
existence is not re-proved as a new result.

## 1. The authenticated native carrier

There are `m>=1` named cells, each the depth-`k>=1` lexicographic directed-triangle
tournament on words `x in {0,1,2}^k`. An intercell edge of visibility `ell`
uses equality of length-`ell` prefixes, or its complementary transpose
according to its retained outer orientation. Depth zero is uniform.
All actual blocks and source vertex identities are supplied.

For `r=0,...,k-1`, form the graph `G_r` of cell edges with `ell>=r+1`.
Let `C_r(i)` be the named connected component containing cell `i`. The full
cell-preserving automorphism group `H` has independent *parameter choices*

    a_(r,C,w) in F3,
    C a component of G_r, w an INPUT prefix of length r.

They act on a source owner `(i,x)` by

    g_a(i,x) = (i,y),
    y_r = x_r + a_(r,C_r(i),x[:r]) mod3.             (1)

There are `D=sum_r 3^r c(G_r)` coordinates and `|H|=3^D`. Their group law
is semidirect, not coordinatewise addition. Independent parameter choices
and independent group multiplication are different assertions.

Retain a supplied actual directed Hamiltonian path `P` of the source `T`.
Every `g_a P` is an actual Hamiltonian path of that SAME tournament.
The orbit has exactly `3^D` paths, because an ordered spanning path has
trivial pointwise stabilizer. It need not be the whole set of Hamiltonian
paths of `T`. The pairwise observable is the actual tournament arc relation;
the fixed named cells and source path are retained data, not cosmetic ties.

## 2. Exact compiler for indexed owner answers

An answer at position `j` is an actual target owner `(i',y)`. Let
`P[j]=(i,x)`. If `i'!=i`, reject: the allowed charts preserve cell names.
Otherwise (1) gives the `k` equations

    a_(r,C_r(i),x[:r]) = y_r-x_r mod3.               (2)

An indexed arc answer at clock `j` contributes the two owner answers at
`j,j+1`; both endpoints and the integer clock are kept. Let `F` be the map
of forced coordinates collected from all answers. Reject iff one shared
coordinate is assigned two different residues (or a requested owner/index
is outside the declared carrier).

**Theorem.** If the answers are consistent, the exact number of remaining
actual Hamiltonian paths in `H P` is

    3^(D-|F|).                                     (3)

Every assignment to the unforced coordinates extends to one such path,
and distinct assignments give distinct full paths. Thus (2) is both
necessary and sufficient, including when different paths' ports share
cells and prefix components.

Proof: the native group parameterization (1) has no additional conditions.
The answers fix exactly the specified coordinate values. Every completion
is an actual allowed automorphism, and freeness of its action on `P` proves
injectivity. The transporter perspective remains valid: these paths form
a coset of the stabilizer of the observed actual owners, although the
coordinate box is not being identified with an abelian subgroup. QED.

This detects injectivity failures as well: asking two different source
owners to have the same image produces incompatible equations at their
first differing input digit, once the relevant cell condition is retained.

## 3. Conditional observation table without listing every frame

For a requested tuple `Q` of indexed owner positions, let `C(Q)` be the
union of the coordinates used by those source addresses in (1). Put

    q = |C(Q) minus keys(F)|.

**Theorem.** Exactly `3^q` ordered answer tuples to `Q` remain possible.
Each such tuple has precisely `3^(D-|keys(F) union C(Q)|)` full-path lifts.
An exact representation of the answer set is (1) restricted to the `q`
new shared coordinates. It must not duplicate a coordinate separately for
each occurrence. Distinct assignments to these coordinates give different
answer tuples because each coordinate is read by at least one requested
digit and can be recovered by subtracting that source digit.

In particular the answer to `Q` is already determined iff
`C(Q) subset keys(F)`. The conditional exposure function is a coverage
function: overlap of the newly exposed coordinate sets is exactly the
dependence that multiplying separate answer counts would discard. Counts
are finite cardinalities; a probabilistic interpretation requires choosing
the uniform distribution on the remaining frames explicitly.

For a single requested owner, forced coordinates along its input address
form an initial segment. A deeper forced coordinate comes from an observed
owner sharing that prefix and lying in the same deeper component; nested
components and shared shorter prefixes then force all its ancestors.
Its possible images are therefore a single ternary-prefix cylinder. This
is exact address information, rather than a three-color summary.

## 4. Minimum additional flags for a TARGETED observer

A bottom slot is `(component of G_(k-1), input prefix of length k-1)`.
One owner answer in a slot fixes its last coordinate and all ancestors
needed for that owner. No owner in another slot can fix this last coordinate.
Let `S_F` be the slots already touched by the received owner answers, and
let `S_Q` be the slots of the requested owner tuple. Set `R=S_Q minus S_F`.

Use the ACTUAL consecutive arcs of `P` to form the slot-incidence graph
`J`: endpoints in different slots give an edge; endpoints in one slot give
a singleton flag. These are positions on one actual path, not hypothetical
edges of an abstract cover. Assume each required slot has an available flag;
this always holds for the supplied spanning path on these nontrivial cells.

**Theorem.** The minimum number of additional indexed arc flags needed to
determine the requested tuple, after the current answers, is

    |R| - nu(J[R]),                                 (4)

where `nu` is maximum matching size in the induced graph on the still
required slots. A flag joining a required slot to a nonrequired slot can
serve as a singleton cover. When every slot has an internal path arc, one
can use those actual arcs for all remaining singletons.

Proof: each uncovered required slot has an independent bottom rotation
`k` fixing the observed SOURCE owners. For any current solution `g0`,
the solution `g0 k` preserves every received answer while changing a
requested owner. This right-action convention also covers nonzero answers.
Thus every required slot must be touched. Touching every required slot
is sufficient because it fixes every coordinate needed by the requested
owners. Pairing two required slots saves one flag exactly along an actual
edge of `J[R]`; the usual edge-cover/matching proof gives (4), with flags
to outside slots counted as singleton covers. The cover is realized by
retained actual path clocks and endpoints. QED.

The same optimum holds for adaptive querying within this pure orbit model.
The coordinate/slot set exposed by a chosen clock does not depend on the
answer values. Every consistent branch must cover `R`, while a nonadaptive
minimum cover already achieves (4). Adaptivity may matter after additional
non-invariant filters; no such extension is asserted here.

Equation (4) specializes to the incoming full-path formula when `Q` is
the entire path and there are no earlier pins. It adds arbitrary consistent
answers and a specified target observation. It does not claim that nine
flags distinguish the order-63 path from every Hamiltonian path outside
its `H` orbit.

## 5. Compatible pins leave a native Hamming move graph

Let `r=D-|F|`. For consistent pins, join two candidate frames when their
parameter lists differ in exactly one free coordinate. This graph is
exactly the Hamming graph `H(r,3)`. Every vertex is an actual native
automorphism and an actual Hamiltonian observer `g_a P`. Every move keeps
all indexed pins, including the endpoints and clocks of arc answers.
Its diameter is `r`; for `r=0` it is a single determined frame.

The assertion has a geometric realization, not only a coordinate count.
Suppose coordinate `(s,C,w)` changes by `delta`. The earlier coordinates
send the input prefix `w` to a common target prefix `v` in every cell of
`C`: all these cells lie in the same earlier threshold components. The
actual map `g_(a') g_a^(-1)` adds `delta` to digit `s` precisely on owners
in cells `C` whose TARGET prefix is `v`; every other digit and owner is
unchanged. This is a subtree rotation on the genuine tournament. Blocks
of depth at most `s` cannot see it; blocks of greater depth have both
cells in `C`, so see the same rotation on both sides. Internal arcs are
preserved by the cyclic digit rotation. The source-prefix coordinate
description remains nonabelian under group multiplication.

These moves change actual owner lifts and port alignments, but keep every
ordered arc bit `T(P[i],P[j])` unchanged. In particular they keep the same
backward-edge free word along the path. Any equivariant construction of a
target graph from `(T,P)` is transported by `g` and retains its isomorphism
class. Thus this is kernel-lift connectivity inside one such native owner
class; it does not cross free-word or target-owner class boundaries.

For the explicitly chosen random-scan heat-bath chain, choose one of the
`r` free coordinates uniformly and replace it uniformly by a residue in
`F3`, including its current value. When `r>0`, the exact spectral gap is
`1/r`. Indeed the one-coordinate averaging operators commute, and their
tensor decomposition has eigenvalues `1-j/r`, with multiplicity
`binomial(r,j) 2^j`. The uniform measure here is a chosen search law, not
an inherited distribution on all tournament Hamiltonian paths.

Arbitrary extra guards may disconnect the surviving move graph or change
its diameter and gap. Rechecking every realized move makes the search
sound but does not retain connectivity. A constant guard retains all the
statements above. These are native automorphism-orbit moves; without an
additional map they are not asserted to be elementary tile flips, path
reversals, or edges of a source metagraph.

**Authenticated observer-footprint corollary.** Suppose a supplied guard
has been proved to factor through the actual indexed observer `Q`:
`Gamma(g_a P)=gamma(Obs_Q(g_a P))`. Let `q=|C(Q) minus keys(F)|` as above.
The conditional table identifies its possible observations with `F3^q`;
let `A` be the assignments accepted by `gamma`. The exact guarded frame
graph is

    H(q,3)[A] Cartesian-product H(r-q,3).           (5)

Consequently its path count is `|A| 3^(r-q)`. Its connected components
are precisely those of the visible induced graph `H(q,3)[A]`, each with
all `3^(r-q)` invisible lifts retained. It is connected iff `A` is nonempty
and that visible graph is connected. Inside any component, its diameter
is the visible component's diameter plus `r-q`. For empty `A` there are
no candidates and no components; a zero-dimensional Hamming graph is a
singleton. Proof: acceptance depends only on the `q` exposed coordinates;
a one-coordinate edge changes either one exposed coordinate within `A`
or one arbitrary invisible coordinate. This is exactly the Cartesian
product adjacency. The standard product shortest-path decomposition gives
the component and diameter claims.

This reduces an authenticated guard calculation to its actual observer
footprint, without recycling or selecting invisible lifts. It does not
assert that a native electrical or mod-eight guard has such a footprint.
That factorization is a separate proof obligation. Multiple guards can
use the union of their authenticated footprints; correlations between
their visible coordinates remain in `A`.

## 6. Why an equivariant selector cannot discard the remaining lifts

For this paragraph the input is the **unpointed candidate fibre plus its
received observations**. The reference `P` used to calculate coordinates
is only an auxiliary chart, not a distinguished selector input. If one
retains `P`, a parameter origin, or an external label order as extra data,
it can break the symmetry and the obstruction below does not prohibit a
chart-dependent completion algorithm.

Suppose these unpointed data leave more than one path in the orbit. The
stabilizer of the observed TARGET owners acts freely and transitively on
those full ordered paths. It is conjugate by any current solution to the
source-pin stabilizer. A deterministic rule
that selects a single candidate and is equivariant under every remaining
source automorphism would have to select a path fixed by that stabilizer.
There is no such path. Thus unbroken kernel lifts cannot be removed by a
covariant choice depending only on the received data. An external label
order can choose a representative, but supplies an extra gauge; it does
not prove that the original observer determined that path.

More generally, an equivariant section of `H -> H/K` for the left action
would send the `K`-fixed coset `K` to a `K`-fixed point of the free action
on `H`, so it exists only when `K` is trivial. This elementary obstruction
explains why compatible frame existence is weaker than an actual port
selector.

## 7. Converse and exact native guards

The same automorphisms preserve `T^op`. Reverse the actual path and replace
each clock `j` by `|P|-1-j` for owner observations (and `|P|-2-j` for arc
flags, reversing its endpoints). The coordinate equations refer to the
same actual owners, so all counts and compatible completions are unchanged.
The checker verifies this using literal arcs of `T` and `T^op`.

If an additional native guard `Gamma(T,P_g,actual ports,integer clocks)`
is required, retain it as a predicate on the realized candidate. The
compiler's counts (3)-(4) apply after that guard ONLY if it is certified
constant on the relevant remaining orbit, or if its extra restriction has
been explicitly incorporated. A source- or path-label dependent guard
need not be constant. A genuine isomorphism invariant of the fully marked
pair is constant under the corresponding transported isomorphism, but
that fact must be established for the actual guard being used.

In particular no native mod-eight electrical fingerprint is replaced here
by the order of a ternary kernel, a parity, or the clock index. This new
family has not been identified with the old order-eight metagraph carrier.
An exact guard adapter may enumerate the remaining finite completions or
consume an authenticated observer table; it must not accept an untested
branch merely because the port equations are consistent.

## Reproduction

    python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/compiler/check.py
    python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/compiler/check.py

`receipt.json` records the finite universes and actual controls. The
unbounded claims above follow from the supplied native coordinate action;
finite tests do not assert coverage of arbitrary tournament metagraphs.

The frozen receipt is reproduced by 117,691 explicit checks in both normal
and optimized Python. Controls include all 9! permutations of the R9
carrier, all 81^3 internal-chart triples of the mixed order-27 carrier,
all 32,768 required-slot subsets of the regular order-63 hostile, every
30,618 directed coordinate move in the mixed order-27 orbit, and all 512
indexed owner-clock subsets in R9. Every one of the 512 guards on a single
indexed R9 owner is compared with its exact visible/invisible product;
connected and disconnected cases retain literal component diameters.
The heat-bath check uses a complete
rational tensor eigenbasis through dimension four. It checks the stated
finite universes, not all arbitrary source tournaments or guard predicates.
