# Session51: which information survives a mathematical quotient?

24 October 2026. **PROVED scoped structural statements; FINITE-EXACT
certified computations; CITED arithmetic inputs.** No claim of literature
priority, complete all-order metagraph wiring, or a scalar recurrence from
order alone. This synthesis distinguishes new extractions from the inherited
mechanisms identified in each theorem packet.

## 1. Counting: one full new order, and the next correction

The parent-only pipeline computes

    A000568(10) = (97,336,068 - 5,508)/10 = 9,733,056.

Its input is the authenticated order-eight parent atlas and full actions.
It constructs the order-nine atlas itself:191480 classes are selected by an
intrinsically distinguished singleton, and56 exceptional classes by marked-root
transport. Only232 singleton-free pointed states require exceptional handling,
with2088 order-eight canonicalizations. The191536 output classes carry their
complete automorphism actions and labelled orbit mass2^36.

The next pass evaluates97213472 pointed order-ten extensions. Production
canonicalization never exceeds order eight; the order-ten evaluator only
extends and refines. A completely separate audit uses its own order-nine
canonicalization, verifies every representative and full group, and agrees
on every per-parent contribution using incoming-edge refinement. The target
count is not an executable input. This is not a claim that all research,
including the independent audit, avoided order-nine canonicalization.

See [construction and prediction](parent-nine/THEOREM.md),
[independent audit](audit-parent-growth/AUDIT.md),
[structural correction](order-ten/THEOREM.md), and
[its audit](audit-order-ten/AUDIT.md).

### What the correction measures

For an ordinary tournament T, let C(O) be the actual stable directed
refinement cell containing the vertex orbit O. Define

    e(T) = sum over Aut(T)-vertex orbits O of (|C(O)| - |O|).

This is lost distinction under refinement, NOT automorphism group order,
Hamiltonian path count, or a count of blue edges. A rigid regular tournament
has maximal refinement error despite having no nontrivial automorphism.

The exact positive order-ten error polynomial is

    12t^14 +46t^18 +4t^36 +84t^40 +14t^72.

It has160 classes, first moment5508, and second moment229416. Exhaustive
odd-cell partition analysis gives strong error3992; the inherited ordered
strong-component law E(z)=A(z)^2 E_s(z) transfers it to5508.

### Order eleven is not left at a conjectured error

The independently audited nonregular classification gives

    E_s,11 = E_regular11 +9166,
    E_11   = E_regular11 +19480.

The known total nonregular positive polynomial is

    48t^14 +348t^18 +8t^36 +256t^40 +28t^72.

Its688 classes are exact, not a lower-bound census. Odd automorphism orbits
and degree constraints isolate the regular sector; in particular regular
error44 is impossible. The small nonrigid classification finds exactly18
nonrigid regular-eleven classes:4 vertex-transitive,7 of error22,5 of
error66,2 of error88. Their78 pointed types give root-moving deficit120.

If M10 is the number of almost-regular ten-vertex classes, then

    R11 = (M10+120)/11,
    E_regular11 = 10 M10 -120.

A targeted parent construction gives M10=13333. Independently, its full reroot
trace is13453, hence R11=1223. Only almost-regular parents are needed, not a
complete order-ten atlas or an order-eleven canonicalizer. The resulting
regular positive polynomial is

    7t^22 +5t^66 +2t^88 +1205t^110.

Thus the complete totals are

    E_s,11=142376,   E_11=152690,

with1907 positive-error ordinary classes at order eleven. This is NOT yet
A000568(11): the unrestricted upper trace U11 remains uncomputed here.

Sources: [nonregular theorem](order-eleven/THEOREM.md),
[audit](audit-order-eleven/AUDIT.md),
[regular computation](regular-eleven/THEOREM.md),
[audit](audit-regular-eleven/AUDIT.md),
[nonrigid structural classification](regular-eleven-nonrigid/THEOREM.md),
[independent structural audit](audit-regular-eleven-nonrigid/AUDIT.md).

### An atlas-free check of the regular coefficient

The standard score generating polynomial gives the labelled regular count

    L11 = [x1^5 ... x11^5] product_(i<j)(xi+xj)
        = 48,251,508,480.

Two independently written residual-degree recursions compute this coefficient,
without using the parent atlas. The classified nonrigid full-group orders
are3 (six classes),5 (seven),9 (one),11 (three),55 (one), so their labelled
mass is

    11! (6/3 +7/5 +1/9 +3/11 +1/55) =151,764,480.

The remainder divided by11! is1205 rigid classes. Their error110 each plus
the nonrigid positive error660 gives133210 again. This is a genuine
score-coefficient/orbit-mass cross-check, not a bare-metagraph reconstruction.
The coefficient formula and orbit-stabilizer theorem are inherited standard
mechanisms; the small complete nonrigid classification makes them useful here.
The [standalone atlas-free packet](regular-eleven-score-mass/THEOREM.md) and
[its independent audit](audit-regular-eleven-score-mass/AUDIT.md) supply a
complete executable route using no target or parent atlas.
See the [first degree recursion](regular-eleven/degree_mass.py),
[independent recursion](audit-regular-eleven/degree_dp.py), and
[exact arithmetic receipt](validation/atlas_free_regular_check.json).

### Two different defects, not one vague notion of symmetry

For a tournament with q vertex orbits, the root-moving contribution is n-q.
It vanishes on rigid tournaments. By contrast regular e=n(q-1) is largest
on rigid tournaments. For the regular-eleven family, root-moving deficit120
is small while refinement error133210 is large. Keeping these two meanings
separate explains why classifying a few nonrigid types can determine a large
refinement correction. It does not supply all higher-order classifications.

## 2. A complete symmetry decoder in a declared all-size domain

Assume each actual stable cell is a singleton or a consecutive cyclic prime
tournament R_p. Distinct prime sizes have uniform crossings. For equal prime
sizes, a nonuniform crossing's full translation stabilizer is either zero
or a nonzero-slope line in F_p squared.

Line stabilizers become gain edges; zero stabilizers pin both cells. A gain
component supports one free rotation exactly when its cycle products are
consistent and no pin reaches it. This gives the complete ordinary
permutation action, not merely its cardinality, and

    e(T)=sum_p p(p-1) times the number of frozen p-cells.

At p=3 every nonuniform crossing gives a line, so only cycle inconsistency
freezes cells. At p=5 a single rigid crossing can freeze two cells on a tree.
The triangle-only cycle picture therefore needs genuine edge-rigidity data.
Arbitrary regular prime cells, including the bad regular-seven example, are
not silently substituted for R_p.

The [decoder theorem](cyclic-cell-decoder/THEOREM.md) and
[audit](audit-cyclic-decoder/AUDIT.md) state the exact domain and recognition
checks. The [loss witness](decoder-loss-witness/README.md), independently
[audited](audit-loss-witness/AUDIT.md), has two ten-vertex tournaments with
identical complete decoder output but H=15465 and15549. Symmetry data does
not determine path-fibre mass or a full metagraph carrying that mass. This
witness alone does not establish different bare unweighted adjacency.

## 3. Converse, actual ports, and actual blue edges

Converse-refining the160 ordinary defective classes gives

    6t^14 +23t^18 +2t^36 +54t^40 +7t^72

on92 merged nodes, weighted error3234. Exactly24 are self-converse, all of
error40. Direct Hamiltonian-word and independent subset-DP controls prove
that all24 are mixed:4200 literal blue presentation ports and337304 black
ports. None is an all-blue point of symmetry. The other68 defect nodes are
black-only. Self-converse, rigid under ordinary automorphisms, and all-blue
are three distinct properties.

The [converse theorem](converse-errors/THEOREM.md) and
[audit](audit-converse-errors/REPORT.md) retain that distinction. Its all-order
weighted strong-component identity is

    C(z,t) = (1+C_s(z,t))/(1-S(z^2,t^2)),
    M(z,t) = (A(z,t)+C(z,t))/2.

Here S counts strong ordinary classes, C_s strong self-converse classes,
C all self-converse classes, and t marks refinement error. Paired strong
components double BOTH size and error. These formulas do not encode the
actual partners of ports.

Actual partners were therefore reconstructed separately. Encoding each of
the4200 blue path words and flipping all nonpath bits gives90 induced blue
edges among the24 nodes:4 loops and86 nonloop edges. There are180 internal
ports and4020 boundary edges. Four degree-one crossing nodes are isolated
in the induced slice; the20 degree-two crossing nodes form one component.
There are no bridges or articulation vertices in this induced graph, and
its cycle rank is71=90-24+5. This integer has no demonstrated connection to
the Ogg prime71; it is an Euler count of a specified multigraph.

See [literal wiring](blue-defect-wiring/THEOREM.md) and
[independent audit](audit-blue-defect-wiring/REPORT.md). The [complete blue neighborhood](blue-defect-neighborhood/THEOREM.md), with
[independent audit](audit-blue-defect-neighborhood/REPORT.md), identifies all
4020 boundary targets:1031 distinct zero-error self-converse classes. The
selected incident-edge graph has1055 vertices and4110 edges, in components
containing4+238 and20+793 positive/zero nodes. It does not include edges
between two boundary vertices, and these are not claimed to be full-metagraph
components.

The [all-size score-layer rule](score-layer-separation/THEOREM.md), with
[independent audit](audit-score-layers/REPORT.md), explains a separation
without enumerating a large graph. At n=2p let a layer have p vertices of
degree a and p of degree b=2p-1-a, and put k=b-a. For a Hamiltonian path
starting at s and ending at t, the inherited nonpath-flip identity implies

    c3(C(T,path)) - c3(T) = d_T(s)-d_T(t)-1.

Thus the allowed neighboring triangle counts are c-k-1,c-1,c+k-1;
blue paths exclude the equal-degree endpoint case. These are necessary values,
not promises of realizable paths or neighbors. Adjacent layers' enclosing
intervals are separated by(p-2)(k-1). Consequently different layers share no
neighbor and have full-metagraph distance at least three, not infinite distance.
An edge inside this union of score layers is possible only in the balanced
layer k=1, with high-degree start and low-degree end. No primality or regular
internal-cell assumption is needed. The inherited necessary loop leg law
also makes all unbalanced layers loop-free; it is not a sufficient loop test.
At order ten the degree-one source layer has c3=30 and blue targets26/32,
while the degree-two layer has c3=40 and blue targets38/40. Their boundary
classes therefore cannot coincide. Finite connectivity inside each selected
incident graph still requires the actual witness computation.

## 4. Relational repair: what a small packet forgets

The inherited exact hole rule remains valid only with option-conditioned
full supports. Equal charge labels need not have the same partners.
The new native example adds an authenticated owner-value relation between
existing1457 and1394 ports, keeping every exterior record, complete kernels,
all declared dyadic guards, converse, and literal changed clocks. A shared
record is inspected by two factors, not allocated twice.

The [native packet](relational-trades/THEOREM.md) and its
[independent audit](audit-relational-trades/AUDIT.md) cover65536 states,
21696 fibres and105112 endpoint pairs. There are762 normalized primitives,
34 newly indivisible. The100 maximal trades give100 z^8(1+z)^3 and800 strict
critical pairs. Ordinary thresholds at grades9/10 rise from14/30 to16/34; the global
strict maximum remains16/34. The explicit grade8/charge165 hostile has new ordinary
optimum10/22 but strict optimum16/34. The shorter ordinary route reacquires
resources and is not a no-reuse route.

The broader [binary-interval theorem](implication-intervals/THEOREM.md),
with [independent audit](audit-filters-intervals/REPORT.md), identifies a
precise structural obstruction. Inside a fixed endpoint interval, every
Boolean unary/binary relation becomes an implication constraint. Its
relation-valid switch sets are ideals of an SCC poset. Complete additive
guards select the zero-defect ideals.

Strict routes are exactly neutral packet partitions whose contracted
dependency graph is acyclic. Small zero-sum packets can create a cycle when
contracted, even when the original graph is a path. Merging the packet-SCCs
repairs execution but must charge the whole merged event. This is why bounded
relation arity, treewidth one, bounded initial SCC size, fixed C2-squared
charges, nonzero increments and exact populations still do not bound event
size. The resource-nonreuse interpretation assumes the declared binary, globally
disjoint reservations and no unencoded acceptance predicates. The native
signed-equality case is an order-free special case; it is not an authentication
of arbitrary implication families at every order.

Joint costs also matter: the audit exhibits the exact frontier
{(6,12),(4,14)}, where coordinatewise minima(4,12) are unattainable. The
complete records, integer margins, continuing retired-resource history and
both cost coordinates cannot be erased by an arithmetic charge quotient.

## 5. Primorials: an exact comparison, not a moonshine identification

The [prime-offset note](prime-offsets/THEOREM.md) recovers the supplied
primorial-indexed Fortunate sequence, lesser-offset sequence (including its
omitted initial term), increasing support and complement. The elementary
bound F_n>p_n makes membership through71 a complete finite question, not a
prefix experiment. Two hundred recursive Lucas certificates and independent
quadratic-form and elliptic-curve computations support the comparisons.

At30030, offset17 is the first success on both sides, giving primes30013 and
30047. Any surviving nontrivial offset below17 squared must be prime. The
window is sharp:30030 plus/minus289 are both prime with composite offset289.
An allowed wheel residue still need not give a prime endpoint.

Among Ogg's15 primes, the Fortunate set misses{2,11,29,31,41}. Combining the
positive and negative offset sets still misses{2,31} and includes primes
outside Ogg's set. All these finite value decisions are global over sequence
indices by the cutoff lemma. Fortune's conjecture and its lesser analogue
remain unproved; the Ogg geometric classification is cited, not established
by the bounded arithmetic check.

The actual bridge to repair is complement logic with explicit quantifiers:
residue screens are overapproximations, whereas exact branch supports already
carry witnesses. Sieve conjunctions and alternative repairs use different
AND/OR operations. No map from supersingular points or prime certificates to
native tournament routes is supplied. See the
[independent arithmetic audit](audit-filters-intervals/REPORT.md).

## 6. What the common mechanism does and does not unify

The session does not identify all these objects. It identifies recurring
ways a projection can lose the very predicate being asked about:

- a color cell forgets vertex orbits;
- a gain kernel forgets the crossing matrix and Hamiltonian mass;
- an ordinary count forgets converse pairing and port partners;
- an arithmetic zero-sum packet forgets precedence and resource history;
- a wheel survivor forgets primality and first-success minimality.

Each positive result restores a specified sidecar or proves that it can be
omitted in a declared domain. The lesson is not that the same small numbers
must encode one hidden object. It is that exact quotients, preserved
predicates, and lost coordinates tell us when information can be transferred.

## 7. Remaining obligations

- Compute U11, or compress the parent information needed for it, before
  claiming A000568(11) by this route. E11 alone is not that count.
- Extend the cyclic-cell decoder beyond its recognized internal cell types.
- Recover black partners and global connectivity beyond the completed blue
  defect neighborhood; local bridge statements are not global ones.
- Authenticate native relational families with controlled neutral packet
  contraction. The abstract chain obstructions prohibit an unqualified
  uniform bound under arbitrary pins.
- Keep the arithmetic offset questions separate from supersingular geometry
  unless an actual predicate-preserving map is found.

All universal claims above rest on proofs with explicit hypotheses. Finite
censuses audit those proofs and identify their exact family outputs; they are
not silently promoted to all-size theorems. Publication status is announced
separately after remote verification.
