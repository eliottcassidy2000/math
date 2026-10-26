# Unnamed continuation: a two-vertex decoder and the order-eight response census

2026-10-26. **INHERITED PROVED:** U1 below; **PROVED:** U2 reset lemma.
**FINITE-EXACT:** exhaustive fixed-path presentations through order eight.
**OPEN:** all-order classification by strong-atom response histograms.
No priority claim. This packet extends the protocol of
[HYP-9416 / fresh re-marking two-port quotient](../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md).

## U1. One singleton on each side recovers every strong joint-port state

Let C be a strong tournament of order n>=3 with distinct ordered ports s,t.
Form S=e star [C;s,t] star e by the inherited rooted hinge product. Name
the new left and right vertices a,b. Thus a->s, C\{s}->a, t->b,
b->C\{t}, and b->a. No labels on S are supplied to the decoder.

**Theorem (PROVED).** The ordinary isomorphism class of S determines the
joint rooted isomorphism class [C;s,t]. Therefore this single named
two-sided attachment observation is injective on all strong states.

**Proof.** In S, a has outdegree one and b has indegree one. They are the
unique ordered pair (a,b) with these degrees whose deletion leaves a strong
core. Indeed an original vertex v has outdegree

    d_S(v)=d_C(v)+1-1_(v=s)+1_(v=t).

Strongness gives 1<=d_C(v)<=n-2. Hence the only possible other outdegree-one
vertex is s, and the only possible other indegree-one vertex is t. If the
candidate deleted low-degree vertex is s instead of a, a remains with its
only out-neighbour deleted, so the remaining core is not strong. If the
candidate high-degree vertex is t instead of b, b remains with its only
in-neighbour deleted, with the same consequence. The case of both wrong
choices also fails. Deleting a,b leaves C, which is strong. Finally s is
the unique out-neighbour of a and t the unique in-neighbour of b. These
are intrinsic tests and commute with every isomorphism. QED.

The incoming independent HYP-9416 audit at `9f2a9722ef` proved this same
result while this packet was in progress. Its degree-only decoder and our
strong-core-deletion decoder are alternative formulations; U1 is credited
to that incoming result, not claimed as a new classification. The original
nine order-at-most-two contexts can be replaced by this single
order-one/order-one context. It is the whole free
class of S, not the scalar H(S); scalar sandwiches can collide. We do not
claim one-sided contexts suffice, or that a response histogram is an
individual rooted state.

The joint class is also sufficient for fresh fixed-port re-marking and
complement: isomorphisms biject the eligible Hamilton paths and their
rooted complemented targets, with multiplicity. Thus any exactly lumpable
quotient preserving this attachment observation is discrete on strong
joint classes. Successors may be nonstrong; this is not a claim that the
strong subset is closed under the dynamics.

## Current finite result

The standard-library C++ census exhausts all 2,097,152 fixed-Hamilton-path
presentations at order eight. It finds 6,880 ordinary classes, including
6,008 strong classes. All 6,008 have distinct joint response histograms,
and all still have distinct histograms after correctly summing multiplicities
over forgotten ports. There are 2,032,504 strong presentations and 269,860
joint target classes. This is finite evidence, not an all-order proof.

Even the unweighted set of ordinary partner classes distinguishes all 6,008
strong atoms. The corresponding joint support also distinguishes them.
Thus multiplicity and joint ports are necessary for composition in general,
but this census does not show that either is necessary merely to identify
a strong source. The two classification problems must stay distinct.

The legacy `atom_responses.py` producer still used a dictionary overwrite
when forgetting ports. It now sums coefficients. The smallest mass-loss
witness is the sole strong four-vertex class: true mass five became four.
This repairs the producer, not the independently audited HYP-9416 counts,
which already used the correct pushforward. Larger erroneous loss amounts
depend on the arbitrary order of overwritten keys; the C++ and Python
diagnostics need not agree on that discarded quantity.

The separate Python checker verifies the sandwich decoder on every root
pair of the complete through-seven strong atlas, checks randomized relabelled
larger controls, and computes unnamed multiset and normalized probability
closures through seven. Those normal and optimized checks agree.

## U2. The literal unnamed reset/complement quotient has a closed formula

This is a different state space from HYP-9416's arbitrary two-port classes.
Let S=X_n be the literal presentations containing the identity Hamilton path,
and tau the deterministic complement reversing all non-path pairs. We keep
tau as a separate allowed operation, not only a composite refresh-then-tau.

**Reset lemma (PROVED).** Suppose q:S->Q is the observation and R is a
transition kernel whose whole row at x depends only on q(x). Then the
coarsest exactly lumpable partition preserving q and both R,tau is

    x equivalent y iff (q(x),q(tau x))=(q(y),q(tau y)).        (U2)

The same statement holds for counted successor multisets. **Proof:** the
two observations at time zero and after tau force this pair in every such
quotient. Tau swaps the pair. Equal pairs give equal q and hence equal
entire R rows before any quotient; their pushforwards to pair cells are
equal. Thus the pair partition is sufficient as well as necessary. The
argument works for any finite involution, without a tournament hypothesis.

For unrestricted fresh re-marking, select a uniformly random Hamilton path
of the current tournament and put it in position order. The output is uniform
over the presentations of its ordinary class A: Aut(A) acts freely on
ordered Hamilton paths, with each presentation having exactly |Aut(A)|
preimages. Thus R is a reset depending only on A. For q(x)=[x], U2 gives

    the exact minimal state ([x],[tau x]).

For q(x)=[x;0,n-1], unrestricted reset still depends only on q. Fixed-port
refresh (choose only paths from the current first to last vertex) is also
a reset depending only on q, uniform on its rooted-presentation fibre.
These paths are nonempty on X_n. Both reset choices therefore have the
same exact minimum for the joint-port observer:

    ([x;0,n-1],[tau x;0,n-1]).

This is a pair of separate rooted isotypes; it need not retain a common
vertex alignment. Fresh reset erases that alignment before it can be tested.
Named guards do not, explaining why the previous named-carrier theorem
has a different answer. If the observation also includes the literal
blue/black flag b, append b; b(tau x)=b(x).

**Composition boundary.** The joint pair is closed under uncoloured ordered
join: compute its source by rooted ordered join and its target by rooted
hinge product. It is not sufficient for blue-coloured joins. The checker
uses the inherited seven-vertex twins x,y with the same two rooted types
and colour. The joins J(x,sigma x) and J(y,sigma x) have respectively blue
and black presentations. Thus even adding each input's colour does not
restore the missing reflection matching. General prime substitutions also
retain the inherited need for all ordered run ports. U2 is an exact answer
for its declared reset/complement alphabet, not for every assembly grammar.

The independent literal refinement agrees with U2 on every presentation
through seven. The C++ census counts pair cells through eight:

| n | literal states | free pair | free pair + colour | joint pair | joint pair + colour |
|---:|---:|---:|---:|---:|---:|
|6|1024|732|740|1024|1024|
|7|32768|24654|24754|32268|32272|
|8|2097152|1821889|1822204|2067580|2067596|

Counts concern the entire X_n, not only strong sources. At seven, equality
of the joint pair identifies 500 presentations beyond singleton counting;
these are genuine unnamed-reset identifications and named-guard separations.

For HYP-9416's different V1 state space, our independent exact refinement
reproduces 11033 strong cells for counted multisets but 11031 for normalized
probability kernels at seven (out of 14162 strong joint states). Empty
successor rows remain distinguished as failure. The normalization distinction
is real, although both protocols have the same pair answer in U2's reset
setting. Bisimulation/lumpability is not silently identified with every
possible weaker trace-distribution equivalence.

## Phase, arithmetic connections and reproducibility

[CONNECTIONS.md](CONNECTIONS.md) develops the typed response-calculus proposal,
the 81/27 hostile law, Frobenius versus determinant parity, the literal
F8 x F8 / F16 cyclotomic bridge, and the limits of a Fermat/Dudeney analogy.
It connects the recovered Vitali, square-sum, graceful, Petersen, Gaussian
and Penrose mechanisms without asserting an unproved isomorphism.

**Universe and controls.** `census.cpp` exhausts every fixed-path bit string
through eight, canonicalizing by all degree-compatible permutations.
`audit.py` uses a separately implemented Python canonicalizer and DFS,
checks the complete inherited strong atlas through seven, all its ordered
root pairs, the corrected aggregation helper, both versions of the V1
closure, and 160 seeded larger sandwich controls. `literal.py` independently
refines the literal reset/complement partitions through seven and checks
the blue-join hostile. `phase.py` enumerates both 12-vertex path sets and
all their complement partners, using the pinned independent audit's
canonical labeller. That last dependency is explicit, not an independent
second canonicalization of the 12-vertex census. The eight-vertex exhaustive
extension has one implementation; the through-seven overlap has two.

Run from the repository root:

```text
python -B 05-knowledge/results/unnamed_carrier_dissection_20261026/verify.py
python -B 05-knowledge/results/unnamed_carrier_dissection_20261026/verify.py --replay
python -B 05-knowledge/results/unnamed_carrier_dissection_20261026/verify.py --full
python agents/check_docs.py
```

Default verification checks seals, pinned executable inputs, receipt agreement,
and independent overlapping results. `--replay` also reruns all Python
calculations normally and under `-O`; `--full` additionally compiles the C++
census with g++ and repeats order eight. These commands use only local inputs
and standard Python/C++ libraries. Receipts record deterministic mathematical
outputs, not timings. The manifest covers the packet except itself, with LF
normalization stated explicitly. Literature and recovery provenance are in
`SOURCES.json`; historical source hashes are not assertions that those
living files can never change. The proofs have a local mathematical audit
and computational controls, not a new independent-agent review.
