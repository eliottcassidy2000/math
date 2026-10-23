# Exact boundary compatibility through triangular recursion

23 October 2026. **PROVED:** the connected boundary-state classification,
minimality, triangular composition, arbitrary-label extension, and the
triangle-to-circle section-existence equivalence below. **FINITE-EXACT:**
the 64-state, 262,144-composition controls in `check.py` and `receipt.json`.
**CITED / RECOVERED:** the historical routes and Hopf-bundle background.

The new increment is an exact compositional carrier for parallel-owner
compatibility. It does not classify all tournaments or authenticate all
Hamiltonian-port partners. A valid tournament may have no parallel owner
section. Existing affine-holonomy theory and the incoming circle embedding
are credited rather than rediscovered as new results.

## 1. Inheritance and the live concepts

The anchor is preserving actual owner compatibility while recursively
contracting connected pieces. The niche is the distinction between additive
closure and square-sum birth graphs. The wildcard is the triangle/circle
extension and its precise separation from Hopf-bundle topology.

The five live concepts are actual boundary names, parallel transport,
simultaneous holonomy, overlap gluing, and the observer being preserved.

Closest proved mechanism:

- [Session39 affine connections](../session39_signed_phases_and_prime_cylinders_20261023/marked-cells/AFFINE_CONNECTIONS.md)
  decodes actual equitable singleton/triple tournament partitions, retaining
  block polarity and the jointly based affine holonomies. It inherits
  [THM-2622, affine torsor holonomy](../../../01-canon/theorems/THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary.md).
- The circle representation already appears in
  [Session39 SYNTHESIS, section 1](../session39_signed_phases_and_prime_cylinders_20261023/SYNTHESIS.md).

Canonical hostile: the incoming
[girth-phase tournament family](../session39_signed_phases_and_prime_cylinders_20261023/marked-cells/GIRTH_PHASE_HOSTILES.md)
has valid nonzero attachment phases that local incidence observations lose.
Two individually solvable reflection constraints can also have no common
fixed point. Independent cycle summaries cannot replace a common root.

Corrected near miss: equal face counts do not make barycentric and edgewise
subdivisions commute, and a one-level common-refinement theorem does not
iterate automatically. The exact correction and the repair by retained
parent pairs are in
[mixed triangle refinement](../mixed_triangle_refinement_20261004.md).

Least-used relevant sidecar: the old
[summand-fibre audit](../arithmetic_braids_20260917_summand.md)
retains both parents and the integer clock. It distinguishes the actual
arithmetic rule from its much less informative directed shadow.

## 2. The square-sum and summand objects are different

[THM-2422, operation fibres and summand closure](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md)
proves that distinct-positive-parent additive closure from seeds `{2,3}` is

    Cl^<({2,3}) = N \ {1,4,6}.

This is an AND-of-ORs dependency hypergraph. Its missing states are not the
connected components of an ordinary directed graph. Allowing equal parents
restores 4 and 6. The existential directed additive shadow has every ascending
arc except `x -> 2x`; its complement is the doubling forest.

The separate square-sum graph has vertices `1,...,n`, with an undirected edge
`x--y` exactly when `x+y` is a square and `x!=y`. The corrected
[square-filter note](../collatz_mod6_20260922_w6_summand_closure_square_filter.md)
proves:

- its component founders are **1,2,4**;
- it has exactly three components only for `4<=n<=12`;
- vertices 13 and 14 perform two successive merges;
- it is connected for every `n>=14`.

At 12 the components are `{1,3,6,8,10}`, `{2,7,9}`, and `{4,5,11,12}`.
In particular 6 belongs to the component founded at 1. The set `{1,2,4,6}`
therefore touches two different constructions; it is not evidence for four
founders or a repeated three-component decomposition. For a general target
filter, the exact birth condition is that it misses `[n+1,2n-1]`. This is an
interval test with its scale retained, not a phase or connectivity invariant.

## 3. The connected affine boundary-state theorem

Let a finite connected support graph carry a copy of `F3` at every vertex.
Each edge has a specified bijection

    pi_ij(x) = epsilon_ij x + delta_ij,
    epsilon_ij in {+1,-1}, delta_ij in F3,

and reverse transport is its inverse. This is exactly the normalized
matching carrier of the incoming triple-cell theorem. Its three labels are
actual cyclic coordinates, not three newly interchangeable colors.

Retain `b>=1` distinct, named boundary vertices. A section assigns one label
at every support vertex and satisfies every edge constraint. Its boundary
relation R is the set of the resulting ordered b-tuples.

**Theorem.** There are exactly

    N(3,b) = 1 + 3^b + 6^(b-1)                              (1)

possible boundary relations. Every one occurs on a connected simple support.
Every admissible boundary tuple has exactly one interior completion.

**Proof.** Root a spanning tree at the first boundary vertex. Tree transport
expresses every label as an invertible affine function of its root value x.
The admissible root values are the common fixed points of the based chord
holonomies. A nonidentity translation fixes nothing, a reflection fixes
one point, and the identity fixes all three. Intersecting these fixed sets
therefore leaves 0, 1, or 3 values.

The empty case gives one relation. A singleton can be any of the `3^b`
named tuples. In the three-point case, the relation is

    {(x,g_2(x),...,g_b(x)): x in F3},

where each g_j is one of the six affine permutations. Different choices
give different relations, yielding `6^(b-1)` possibilities. Once the first
boundary label is fixed, connected tree transport uniquely fixes every
interior label. Thus no multiplicity is lost by passing to this relation.

To realize these states, join the root to the other b-1 boundary vertices by
the desired transports, and attach a triangle `root--u--v--root`. Use identity
transport on its first two edges and holonomy h on the last. Choose h to be
the identity, a reflection with the required fixed root value, or a nonzero
translation. Singleton coordinates can be supplied by translations on the
boundary star. All supports are simple and connected. QED.

For three boundaries this gives

    1 empty + 27 singleton + 36 three-point states = 64.     (2)

An enumerated state needs exactly six bits. The sequence for b=1,...,5 is
`5,16,64,298,1540`; the apparent powers-of-four pattern at b=2,3 already
fails at b=4. A complete four-boundary state needs nine bits, not eight.

**Minimality for the stated observer.** If R and S differ, choose a tuple
p in their symmetric difference. A connected piece exposing only `{p}`
exists by the realization above. Glue its b named boundaries to the original
boundaries. To retain a simple union support, first subdivide every edge of
the pin context joining two boundary vertices through its own fresh interior
vertex, using identity transport followed by the original transport. All
other context interior vertices are fresh as well. This preserves the
relation and unique completion, and prevents overlapping boundary edges
from creating parallel or conflicting edge prescriptions. The result has a
section for exactly one of R,S. Therefore any
summary required to decide section existence under every such boundary
context must distinguish all the states in (1). This is an information bound
for this compatibility observer, not for whole-tournament reconstruction.

## 4. Three-piece gluing is an exact recursion

Arrange three connected pieces with named boundary triples

    (a,x,z),  (x,b,y),  (z,y,c),

where x,y,z are shared corners and a,b,c remain external. Define

    Glue(R,S,T) = {(a,b,c): there exist x,y,z with
                   (a,x,z) in R, (x,b,y) in S, (z,y,c) in T}.  (3)

The union support is connected: every pair of pieces shares a named corner.
It therefore satisfies the theorem and exposes one of the same 64 states.
Equation (3) is exact elimination of the three internal corner variables.
No admissible assignment is invented and none is discarded. Connectedness
also implies unique completion for each surviving external tuple.

Consequently arbitrary-depth recursion with this three-copy incidence
pattern can use a fixed 64-state transition table. Compositions remain exact
at every finite depth. This is the triangular overlap pattern of the
Sierpinski gasket; it is not an identification with a three-adic tournament,
a geometric fractal limit, or a dyadic Pascal operator.

The fixed table has `64^3=262144` entries and can be stored with one byte per
entry. Building one summary per node of a supplied recursion tree takes
linear time in that tree. General affine-support section existence already
has a direct linear propagation algorithm; the point of the table is exact
caching and local recomposition, not an improvement over that whole-graph
algorithm. It does not construct the recursion tree or authenticate missing
port identities.

### The cheapest hostile: counts alone give three different answers

Take R=S to be the diagonal relation `{(x,x,x):x in F3}`. All three choices
below for the third child also have cardinality three:

    T_0 = {(x, x,   x): x in F3},
    T_1 = {(x,-x,   x): x in F3},
    T_2 = {(x, x+1, x): x in F3}.

The parent relations have respectively 3, 1, and 0 elements: the closed
constraint is `x=x`, `x=-x`, or `x=x+1`. Thus three separately nonempty
children with identical counts need not glue to a nonempty parent. T_0 and
T_2 even have the same positive transport signs; the phase is essential.

This gives a precise version of “blocks must succeed together.” The shared
owner equations, not the individual success counts, determine compatibility.

### Quotienting the boundary charts too early also fails

Independent cyclic relabelings of the three boundary coordinates reduce the
64 states to six orbits, of sizes `1,27,9,9,9,9`. T_0 and T_2 are in the same
orbit: shift their second coordinate. However, that coordinate is a shared
corner in (3). Changing its chart must also change the neighboring piece.
Independently choosing a representative for each child loses precisely the
transition data needed for gluing. Keep actual boundary coordinates, or keep
every overlap chart transformation alongside the quotient state.

### Scope boundary

Do not discard a disconnected component that has no retained boundary. An
unseen component can have no sections, or can multiply the number of sections
by three. Our recursion never creates such a component: every input support
is connected, the pieces share corners, and an external boundary remains.
Arbitrary contractions need additional connectivity and multiplicity data.

## 5. Extension toward seven-label cells

The same argument works for an arbitrary r-element label set when every
support edge is a specified permutation. It does not require affine maps.

For `r>=3`, the exact number of connected b-boundary relations is

    N(r,b) = sum_{k=0}^{r-2} binom(r,k) ((r)_k)^(b-1)
             + (r!)^(b-1),                                 (4)

where `(r)_k=r!/(r-k)!`. A subgroup of S_r cannot fix exactly r-1 points:
every permutation fixing those also fixes the last. Every other fixed set D
is realizable, using one cycle on its complement (or the identity if D is
the whole set). For a fixed root domain D of size k, each other boundary
coordinate is any injection from D into the label set, hence `(r)_k`
choices. These injections extend to permutations. Rooting at the first
boundary makes the count unique. The star-plus-triangle construction realizes
every state, and the same pin-context argument proves minimality.
Use fresh subdivisions on boundary-to-boundary context edges as in section 3;
the general-r minimality claim also concerns simple union supports.

For r=7 the one-, two-, and three-boundary counts are

    121, 95642, 185036888.

Thus a universally complete three-boundary compatibility code for arbitrary
seven-label permutation connections needs at least **28 bits**. This count
is a theorem from (4); the checker evaluates the formula, not an exhaustive
185-million-state test. A small independent enumeration gives 113 states
for r=4,b=2.

This identifies two separate advantages of triples. Every permutation of
three labels is affine over F3, and every regular nonuniform 3-by-3 binary
block is a matching or its complement. Neither statement extends to seven:
regular 7-by-7 blocks of row sum 2 through 5 are not selected matchings.
The seven-label extension therefore applies only after actual edge
bijections have been supplied. The internal seven-vertex tournaments,
nonmatching blocks, polarities, cell identities, and marks remain sidecars.

## 6. Triangle and circle: an exact existence equivalence

Use the inherited Session39 embedding, with `omega=exp(2*pi*i/3)`:

    x in F3 -> omega^x in S1,
    (epsilon,delta): z -> omega^delta z^epsilon.             (5)

For epsilon=-1, z^epsilon is complex conjugation on S1. Formula (5) is a
faithful copy of the dihedral group D3 in O(2), and composition is exactly
the affine transport law. The new observation is a precise preserved
predicate beyond merely embedding the group.

**Theorem.** On any connected support, the discrete triangular-label
connection has a parallel section if and only if its circle extension (5)
has one.

**Proof.** Both section sets are fixed sets of the same root holonomy
subgroup H<=D3. If H is trivial, there are three discrete sections and a
whole circle of continuous ones. If H is generated by one reflection,
there is one discrete section and two antipodal circle sections. Otherwise
H contains a nonidentity rotation: two distinct reflections multiply to
one, and the remaining nontrivial subgroup is the rotation subgroup.
A nonidentity rotation fixes no point on the circle or the triangle.
These are all six subgroups of D3. QED.

Existence is therefore preserved; multiplicity and the actual section space
are not. A full frame, a root mark, and a choice of one antipodal fixed point
are stronger data. In particular this theorem does not make a circle and a
triangle homeomorphic.

### Actual order-31 control: all 420 complete matching systems are empty

The concurrent [endpoint overlap packet](../endpoint_overlap_ports_20261023/README.md)
grounds 420 retained ordered-mark cases whose actual partitions have seven
singleton cells and eight triples. This is a count of cases, not of isotypes.
Their matching support is K8. With the retained first cell as root and its
star as spanning tree, the 21 chord holonomies have exactly these profiles:

| Actual systems | Fixed-point-free translations | Reflections | Identity loops |
|---:|---:|---:|---:|
| 168 | 9 | 12 | 0 |
| 252 | 21 | 0 | 0 |

Thus every one of the 420 systems has **zero common sections**, and a single
literal triangle loop already witnesses that fact. The separate
`actual_control.py` reconstructs all 28 selected matchings directly from
the incoming adjacency rows and retained cell lists, compares all three root
assignments with all 21 loop tests, and saves one obstructing triangle per
case. It independently obtains 6,804 fixed-point-free loops and 2,016
reflections. This is a finite statement about this fixed authenticated host
and these actual partitions, not every order-31 tournament.

This is a decisive limit of the 64-state observer: contracting the entire
connected K8 to three external cells gives the same empty state in all 420
cases. Such a summary cannot reconstruct these sources or serve as a rule
rejecting empty native systems. K8 also comes with no asserted width-three
Sierpinski decomposition. The compatibility carrier can still summarize
declared smaller pieces, but the complete connection must retain more.

The 168 reflection cases show what that extra information does. The fixed
owner of one specified negative loop is unique even though it satisfies no
global parallel section. It supplies the intrinsic mark used by the incoming
decoder. The circle extension doubles that loop's fixed set to two antipodal
points, so it loses the unique-owner property. Retaining the actual embedded
triangle selects the unique fixed point in its three-element fibre. Exact
section-existence transfer and exact owner recovery are different goals.

### The same holonomy has two different outputs: fixed points and centralizer

The inherited joint decoder also identifies the residual rotational
automorphisms. Let H<=S3 be the subgroup generated by the jointly based
holonomies, and let C3 denote the rotations of the actual directed root
triple. Then

    parallel owner sections = Fix_F3(H),
    cell-preserving rotational automorphisms = C_C3(H).     (6)

These are different constructions from the same retained transport data.
The first asks for a root value fixed by every loop. The second asks for a
root rotation commuting with every loop. Tree transport propagates such a
rotation uniquely to each other cell; conjugating a rotation by an affine
map is still a rotation. The chord equations hold precisely when the root
rotation centralizes H. The joint decoder then preserves every original
block, its polarity, and all retained singleton attachments.

In the actual 420 cases the recorded loop profiles determine H completely:

| Ordered-mark cases | H | Fix_F3(H) | C_C3(H) | Additional actual owners after triple partition |
|---:|---|---|---|---:|
| 252 | C3 | empty | C3 | 1 |
| 168 | S3 | empty | trivial | 0 |

In the first row every chord is a nonidentity translation, so the generated
group is C3 and all root rotations commute with it. A nontrivial rotation
moves every vertex of each triple; an actual owner is necessary to remove
that residual symmetry, and one suffices by connected matching propagation.
In the second row there are both nonidentity translations and reflections,
which generate S3. A translation `x->x+a` commutes with a reflection only if
`2a=0` in F3, hence only if a=0. The incoming negative-loop decoder obtains
an intrinsic fixed owner and resolves these 168 cases without an additional
externally supplied owner. Rigidity alone would not prove that last
algorithmic statement; the existing decoder does.

Thus the same empty fixed-section summary hides two different symmetry
regimes and two different owner budgets. Equation (6) is an application of
the inherited affine decoder, not a new general holonomy theory. It explains
why the full transport representation remains useful after its common
fixed-point set has become empty.

This diagnostic does **not** supply a closed 65-state recursion by splitting
the empty state into two symmetry types. An empty piece with H=C3 can still
have a reflection in a tree transport between its boundary cells. That
relative chart is invisible to both its empty fixed set and its root
centralizer. Gluing three such pieces can create either an orientation-even
new loop or a reflection loop; in the latter case the combined holonomy
becomes S3. The boundary transports must therefore accompany any proposed
symmetry-sensitive compression. The 64-state theorem remains exact for its
declared compatibility predicate; reconstruction retains the full transport
tuple and block data.

## 7. Hopf, Bott, and the operation-bearing sidecar

The recovered old octonionic route is
[octonions at the center, section 4](../../../07-reflections/octonions-at-the-center-s6-bott-lyapunov-and-class-rank-frontiers-root-20260824.md).
It records `S7 -> S15 -> S8` as a homogeneous Hopf fibration, and explicitly
separates common dimensions from identical structures. The separate
[Hopf/S6 boundary audit](../../../07-reflections/hopf-s6-triangle-family-conductor-clutch-and-crossfrontier-root-20260824.md)
rejects transferring an A2-specific boundary-forest theorem to an unrelated
quotient merely because both have a normalization and boundary branches.
The repository's Hopf-algebra work is a different construction again.

For primary background, Baez's
[The Octonions, projective lines](https://math.ucr.edu/home/baez/octonions/node9.html)
describes canonical line bundles through their two hemisphere charts and
equatorial transition functions, then derives the Hopf sphere bundles.
The useful lesson is to keep the chart gluing map. Our boundary-chart hostile
in section 4 is an exact finite counterpart of that bookkeeping principle.

There is no asserted map identifying our support graph with a Hopf base or
our finite owner labels with Hopf fibres. In particular, a U(1) bundle over a
graph has no degree-two Chern obstruction, while its connection can still
have nontrivial loop holonomy and no parallel section. Indeed one can choose
a section at graph vertices and extend it along each edge using connectedness
of the circle fibre; being parallel imposes additional transport equations.
The Hopf circle bundle
over S2 has a topological section obstruction; graph-loop transport is a
different predicate. The finite D3 action (5) supplies a real map; invoking
Hopf topology without supplying a base, fibre, projection, and connection
would not.

The earlier native modulo-eight character guard likewise lives on a different
carrier: a fixed labelled Laplacian potential. Its three-bit residue rejects
all three nonnative members of one certified feasible switch family. The
64-state carrier here retains a different six-bit compatibility object.
These guards can be applied jointly after both maps have been constructed;
the shared small powers of two do not identify their invariants or establish
Bott periodicity. See the inherited
[Hall/Bott packet](../metagraph_hall_bott_20261023/README.md) and incoming
[converse obstruction](../session39_signed_phases_and_prime_cylinders_20261023/linking/CONVERSE_IS_NOT_TIME_REVERSAL.md).

## 8. Map and loss ledger

| Source | Map / target | Preserved predicate | Lost information / required sidecar | Cheapest test |
|---|---|---|---|---|
| Connected named affine matching piece | Restrict sections to boundary relation | Compatibility with every external pin context; unique completion | Full tournament, individual holonomies, internal geometry; retain block and transport tables for reconstruction | Three 3-state children give 3,1,0 outcomes |
| Three compatible connected pieces | Existentially eliminate shared corner labels | Exact section relation at every finite recursion depth | If pieces are disconnected, boundary-free components need a separate factor | All 64^3 joins against literal support propagation |
| Actual boundary relation | Quotient independent boundary rotations | Its isolated gauge orbit | Shared chart alignment; retain overlap transformations | T_0 and T_2 share an orbit but glue differently |
| Discrete D3 connection | Embed triangle vertices in the unit circle | Existence of a parallel section | Multiplicity 3/1/0 becomes infinite/2/0 | Enumerate the six holonomy subgroups |
| Connected r-label permutation carrier | Root domain plus boundary injections | Exact compatibility | General regular adjacency blocks are not bijections | Seven-label formula; four-label two-boundary enumeration |
| Hopf chart presentation | Retain transition maps | Its own bundle gluing | No direct graph-to-Hopf identification is supplied | Distinguish topological section from parallel section |

The next useful actual-port query retains a specified loop's action and its
owner, rather than replacing the whole connection by the common section set.
Boundary-state comparisons can be used on supplied connected patches, with
their actual overlap coordinates; no width-three decomposition is inferred
from the eight triple cells. The test must compare to the native state and
must not reject every empty state. For exact reconstruction continue to
retain the full jointly based holonomy tuple and original block data.
This note does not supply a missing uniform
Collatz payment bound or an all-size metagraph classification.

## 9. Reproduction and proof/experiment boundary

From the repository root:

    python 05-knowledge/results/hopf_sierpinski_port_recovery_20261023/check.py
    python -O 05-knowledge/results/hopf_sierpinski_port_recovery_20261023/check.py

`--write` regenerates the receipt. Ordinary execution compares the generated
and frozen receipts after JSON normalization; optimization does not disable checks.

The consumer checks all 64 literal five-cell realizations by both brute
assignment and propagation, all 2016 pairs of distinct states against a pin
context on a literal simple nine-cell union, and all 262144 three-piece joins against independent propagation on
their literal 12-cell connected supports. The exact composition counts are:

    empty: 186571; singleton: 67797; three-point: 7776.

The full byte transition table has SHA-256
`9e3db55b8142aa6cc2406a005b3afaeb23875262ff42a429e230c21b8867f018`.
Additional checks cover the boundary gauge orbits, all six D3 subgroups,
the circle action law, the 3/1/0 hostile, and the stated small state counts.
The general b and r formulas and arbitrary finite-depth recursion are proved
above; finite enumeration is a control rather than a substitute for those
quantifiers.

The separate actual-host consumer reads the frozen endpoint packet and
reconstructs its matching maps:

    python 05-knowledge/results/hopf_sierpinski_port_recovery_20261023/actual_control.py
    python -O 05-knowledge/results/hopf_sierpinski_port_recovery_20261023/actual_control.py

Its 76,446 replay checks ground the 420 empty-section controls and save the
input receipt digest. The saved seven-cell refinement history is inherited;
literal triple regularity, matching bijections, and common-section rejection
are checked anew.

Independent root review accepted the state classification and requested the
explicit simple-support subdivision guard in the minimality proof. That
guard is now in both the general proof and all 2016 literal pin-context
tests. Normal and optimized replay perform 536810 checks each. The separate
electrical-side review confirmed the connectedness, unique-completion, and
triangle/circle subgroup arguments; its final audit is retained separately
when available.
