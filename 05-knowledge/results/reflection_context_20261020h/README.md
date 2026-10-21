# Two arithmetic reflections, marked graph operations, and lost palette support

20 October 2026. **PROVED scoped elementary statements**, with **FINITE-EXACT**
controls in `check.py` and `receipt.json`. No literature-priority claim, no
black-core connectivity theorem, and no six-colour forcing construction.

## Inheritance and live concept board

The closest proved mechanism is the [ordered-graph trijection and tree
gauge](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).
Its hostile example is the actual transitive-tournament bridge: squares in
the cut aggregate do not lift to legal black detours at a specified source.
The corrected near miss is the [marked tetrahedral
bridge](../session17_prime_shapes_20261020/PRIME_BALANCE_TETRAHEDRAL_BRIDGE.md):
it transports a particular divisor involution, not divisor complement, and
needs an orientation marking. The least-used relevant sidecar is the full
twenty-state [palette Gram kernel](../first_exit_complement_games_20261020/PALETTE-GLUE.md),
including the named input triples, not just its spectrum.

The board is: **literal source; reflection action; cut and port; palette
support; arrival clock; complete observer**. The anchor is the requested
Collatz payment frontier, pursued in the [parallel debt
packet](../debt_annulus_20261020h/README.md); the niche is the two-reflection
map; the wildcard is observer completeness in idoneal orders. These are
connections between mechanisms, not identifications of their objects.

## 1. The two reflections cannot be ordinary complement and reversal

For distinct primes p,q,r let D be the ten proper nontrivial divisors of
p^2 q r, represented by (a,b,c), with 0<=a<=2 and 0<=b,c<=1, excluding
(0,0,0) and (2,1,1). The inherited chart is

    f=a-b-c,  g=b-c.

Its horizontal rows have sizes3,4,3, and the negative/zero/positive f
regions have sizes4,2,4. Use

    Q(a,b,c)=(a,c,b),
    C(a,b,c)=(2-a,1-b,1-c),
    F=CQ.

Q reverses g, F reverses f, and C is divisor complement. They commute.
The numbers of fixed points are respectively4,2,0 for Q,F,C. The action
has two two-element orbits stabilized by Q, one two-element orbit stabilized
by F, and one free four-element orbit. This is stronger information than
the row counts.

On simple graphs with an ordered vertex set of size m>=3, let K be full
graph complementation and R the reversal of the vertex order. They generate
a four-group. K fixes no graph. KR also fixes no graph: the pair consisting
of the first and last vertices is fixed by R, so its edge bit would have
to equal its complement. Only R can have fixed points.

**Obstruction.** There is no equivariant map from D to this ordered-graph
space under any isomorphism between the two four-groups. Indeed both Q and
F fix a point of D; the image of each such point must be fixed by the
corresponding target operation. There are two distinct such operations on
the source and only one possible target operation. Injectivity is not
needed. For m=2 the target action already has only two distinct operations.

Adding decorations with an equivariant forgetful map to these same graphs
does not evade the obstruction: compose with that forgetful map. Changing
the graph operation, removing the order, or changing the objects can evade
it, but then the new preserved predicate must be specified.

Under the fixed-Hamilton-path trijection, K reverses all nonpath arcs;
R takes the converse and reverses the marked path; KR reverses only the
path arcs and then reverses the marked path. Thus the obstruction applies
to precisely these natural marked tournament operations. It does not
assert an obstruction to every possible unmarked tournament encoding.

## 2. A faithful replacement with a two-edge switch

Here is an explicit alternative, not only an obstruction. Encode (a,b,c)
as a graph G on the ordered vertices0,1,2,3 with bits

    x03=[a=0],  x12=[a=2],  x01=b,  x23=c,
    x02=x13=0.

The inverse is a=1+x12-x03, b=x01, c=x23. The mutually exclusive bits03,12
and the two excluded extreme words retain the full ten-element source.
Let rho=(03)(12), gamma=(01)(23), and M={01,23}. Then

    G(Qd)=rho(G(d)),
    G(Cd)=gamma(G(d)) symmetric-difference M,
    G(Fd)=gamma(rho(G(d))) symmetric-difference M.       (1)

These are literal edge operations. Gamma swaps03 with12, fixes01 and23,
and swaps the two pinned zero bits. It exchanges the two extreme a-states
and toggles b,c. Rho fixes03 and12 and swaps b,c. The operations commute
and are involutions. This proves (1) and reversibility for every source.

The second operation is **not full graph complementation**. Its specified
two-edge switch is exactly the change needed by this encoding. The four
ordered vertices and the mask M are part of the representation. Multiplying
or comparing divisors requires decoding and the original named primes;
graph adjacency alone is not multiplication or divisibility.

Lift G through the trijection to a five-vertex tournament with the marked
path0,1,2,3,4: an edge{i,j} of G reverses the free arc i->j+1. Its four
possible backward arcs are4->0,3->1,2->0,4->2. Rho now is exactly converse
plus marked-path reversal. The other reflection remains an explicitly
marked arc operation; it is not a path change in one fixed tournament.

**Failure boundary.** The lifted source is strong exactly when a=0 or
(a,b,c)=(1,1,1). Backward intervals cover every path cut precisely in these
four cases. The operation representing C sends every one of those four
strong sources to a nonstrong source. Thus the reversible encoding does
not preserve the strong core, let alone prove a black detour. This is a
cheap decisive hostile test, with all ten rows retained in the receipt.

The successful connection preserves the complete divisor-reflection action
and one natural tournament reversal. It supplies an information register,
not a black-core reduction. A detour must additionally retain directed
entry/exit ports, prove each path legal and black, and reach the literal
target; the [endpoint-four packet](../endpoint_four_20261020h/README.md)
performs those separate checks for its stated layer.

### Incoming resolution: keep the ordinary operations by retaining two more states

During the final sync, incoming commit098cb67c1d supplied the independent
[minimal twelve-state resolution](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md).
It retains each divisor as a projection, splitting only pq and pr into two
orientation states each. Its graph carrier uses ordinary complement K and
order reversal R. Twelve states are minimal: over each of the two F-fixed
divisors, the fixed-point-free target KR needs a positive even fibre.
The original obstruction and this resolution are compatible: the new source
is a cover of D, not an equivariant section from D into the old graph space.

The two constructions can be joined explicitly. With edge-bit order
`01,02,03,12,13,23`, the incoming graph codes map to our switch-chart codes
as follows:

| Incoming code |0|63|4|59|1|32|62|31|2|47|61|16|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Switch-chart code |0|33|8|37|5|36|40|9|1|1|32|32|

Call this quotient H and our switched reflection S. Directly on all twelve
states,

    H K = S H,       H R = R H.

It collapses precisely the two orientation pairs over pq and pr, retains
the actual divisor, and preserves the fixed-path blue predicate. It does
not preserve full graph complement as the same operation: K becomes S.
There is no equivariant reverse section retaining K,R, by section1. This
is an exact tradeoff between context and operation, rather than a numerical
coincidence between ten and twelve. Neither construction transfers arbitrary
path incidence or black-core connectivity.

`incoming_bridge.py` reproduces the37 exact checks of this new quotient;
normal and optimized output agree. The minimality and all-exponent orbit
capacity theorem belong to the incoming packet, not to this session.

## 3. Freely composing palette gates destroys their exclusions

This statement concerns exact boundary-relation composition. It does not
silently assume that a geometric graph enforces monochromatic input triples.
For three-subsets A,B of six named colours, the inherited two-gate relation
has multiplicity576 K(A,B), where K(A,B)=|A intersection B|. Let V be the
20-by-6 incidence matrix, so K=VV^t and V^t V=6I+4J. In particular K has
row sum30. Multiplication gives

    K^2=6K+36J,
    K^r=6^(r-1) K + (30^r-30*6^(r-1))/20 J, r>=1.      (2)

The second formula follows by induction using KJ=30J; it also follows from
the inherited two nonzero eigenspaces. Every entry of K^2 is positive.
Even complementary palettes, excluded by one stage, have36 unscaled
extensions after two stages when the intermediate palette is freely summed.
All longer freely summed chains retain full endpoint support. The actual
count for r copies of the inherited relation includes the factor576^r.

For P=K/30 the normalized formula is particularly informative:

    P^r(A,B)=1/20 + 5^(-(r-1)) (K(A,B)/30-1/20).

Thus the nonuniform information contracts by1/5 at each further stage,
while the exclusions have already disappeared at stage two. This is an
exact explanation of the factor5 in this particular construction, not a
universal claim about fivefold geometry. A small nonuniform signal cannot
substitute for a zero entry in a forcing argument.

The mechanism generalizes without a spectral guess. For all s-subsets of
q named colours, 1<=s<q, the analogous incidence matrix satisfies

    V^t V = binom(q-2,s-1) I + binom(q-2,s-2) J,

with the second binomial read as zero when s=1. Consequently

    K^2 = binom(q-2,s-1) K + s^2 binom(q-2,s-2) J.

For s>=2 every entry is positive after two stages. After normalization
the nonconstant nonzero modes have ratio `(q-s)/(s(q-1))`. The boundary
s=1 instead gives K=I: a singleton palette retains its exact colour.
This identifies precisely which forgotten interface causes the failure;
there is no assertion that an actual geometric network implements a freely
summed palette interface without enforcing its named triples.

The construction target is therefore sharper: a new palette network must
retain or constrain its intermediate palette, introduce a coupled relation,
or otherwise prevent this positive two-stage completion. The actual named
ports and their equality obligations are indispensable. Neither a rank-six
summary nor arbitrary independent repetitions of this gate removes all
separating colourings. This is compatible with the inherited fact that the
entire present geometric carrier is four-colourable.

### A sharper network obstruction: retain a relative twenty-core

For these twenty three-colour palettes, K(A,B)=0 exactly when B is the
complement of A. Thus every assigned neighbouring palette excludes just
one of the twenty choices at a variable. Consider a finite simple
constraint graph whose edges carry this K relation, with boundary variables
pinned to actual palettes. Count neighbours in the graph, not parallel
copies of the same relation. If every unpinned variable can be successively
deleted when its degree in the remaining graph, including all pinned
neighbours, is at most19, then every admissible boundary assignment extends.

Proof: reinsert variables in reverse deletion order. At each insertion,
at most19 assigned neighbours exclude at most19 palettes; choose one of
the remaining twenty. Each selected edge has positive multiplicity, so
the complete assignment has positive extension weight. Boundary-to-boundary
relations must already be satisfied, and the boundary is never deleted.

Consequently an obstruction created by unpinned variables in this interface
model needs a nonempty relative twenty-core, or extra relations besides K.
The threshold is sharp for externally pinned palettes: nineteen distinct
neighbour palettes force the remaining complementary choice, while all
twenty distinct neighbour palettes exclude every choice. These are abstract
boundary tests, not twenty geometrically synchronized palettes already
constructed in the plane. The elementary greedy principle is familiar;
the application identifies the exact obstruction threshold of this inherited
palette kernel. It parallels the relative six-core test for ordinary
six-colouring, with a different state space and a different degree bound.

This yields a concrete next construction target instead of another free
chain: coupled palette constraints that survive this deletion, with all
their physical synchronization and carrier obligations retained. No relative
core alone is asserted sufficient for forcing: without pinned constraints,
assigning the same palette everywhere satisfies K on every graph, however
dense it is.

## 4. What the65 idoneal numbers contribute

The repo's adopted completeness result is
[THM-4566, Hilbert class field planes](../../../01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md).
This packet uses its current scope; it does not reopen its adopted analytic
input. The useful transfer is a precise finite-group fact, independent of
the numerical length of the list.

For a finite abelian group G, all characters G->{+1,-1} jointly separate
elements exactly when2G=0. Their common kernel is2G: each character kills
doubles, while linear functionals separate nonzero vectors in G/2G.
For the proper ideal class group of the quadratic **order** of discriminant
-4n, the idoneal property is exactly this exponent-two property. Thus its
genus/sign observers are complete; a general quadratic class group can
hide a nontrivial fibre behind the same sign data.

The [conductor-growth packet](../escape_renewal_idoneal_20261017d.md) already
provides the useful hostile family. Along n=462*4^k, the first two class
groups are C2^3 and C2^4; at k=2 the group is C4 times C2^3. The sign observer
still has16 values but the class group now has32 elements. A missing
coordinate has appeared despite unchanged prime support. Quadratic order
and maximal-order field must not be interchanged here.

This gives a concrete design rule: prove an observer injective on the
actual guarded state set before replacing that state by its character
values. It does not establish that mod15/mod35 signs determine a Collatz
source, its joint residue, its terminal, or its arrival clock. Their known
hidden joint-flow coordinates remain necessary. The analogous failure in
section3 is visible directly: summing out a palette restores forbidden
boundary pairs.

Nor can the65 objects by themselves be a finite universal debt controller.
The [finite-barrier theorem](../guarded_reset_debt_20261020/FINITE_STATE_BARRIERS_ADDENDUM.md)
already proves: a65-state positive physical-height controller with fixed
nonempty valuation words of length at most L misses some source congruent
to-1 modulo2^(65L+1). Its first65L valuations are1; a repeated state would
give an expansive cycle, contrary to the barrier. This is a specialization
of an existing obstruction, not a new disproof of Collatz. Parameterized
unbounded words, additional state, or a different proved resource remain
legitimate routes.

## Reproduction and remaining obligations

Run from the repo root:

    python -X utf8 05-knowledge/results/reflection_context_20261020h/check.py
    python -O -X utf8 05-knowledge/results/reflection_context_20261020h/check.py

The explicit finite universe is all ten divisors; every ordered graph of
sizes2 through6; every marked-path operation through graph size4; every
named palette pair at chain lengths1 through5. Tests use explicit failures,
not Python asserts. There are71,746 checks. Direct arc relabelling is
independent of the compressed graph formulas; graph strongness uses actual
reachability, independently checked against interval coverage. The source
and target fixed-point counts and the complementary-palette zero are the
hostile controls. The all-order results are the proofs above, not an
extrapolation of these finite checks.

Independent session audit: `family_algebra` reviewed the reflection map,
fixed-point obstruction, failure of strongness, palette powers, and the
scopes of the inherited class-group and barrier results. Its optimized
checker replay matched all71,746 checks and receipt LF SHA256
`1738a33eead0ea52df17e821d71e01237908f3ccf02bc83f6bec24cb6d707392`.
The same reviewer separately checked the general q,s incidence formula and
the relative twenty-core extension with its pinned-boundary hypotheses.

Still OPEN: all-layer black-core connectivity; an unconditional geometric
six-colour synchronizer beyond the current carrier; a source-uniform
finite-arrival annulus estimate; and universal grounded Collatz coverage.
The shared theme is to retain a complete boundary state. No reduction
between these open problems is asserted.
