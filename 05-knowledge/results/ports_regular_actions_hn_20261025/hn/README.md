# Six, seven and eight states: an exact circle-port language, and the 30030 boundary

2026-10-25. **PROVED** elementary interface and obstruction statements;
**FINITE-EXACT** controls; **INHERITED PROVED / CITED** prime-offset and
unit-circle mechanisms. No Hadwiger–Nelson lower bound, Fortunate conjecture,
or tournament-connectivity theorem follows from these comparisons.

## 1. Inheritance and the precise change of question

The [circle-port theorem](../../circle_ports_polygon_palettes_20261021/README.md)
proves that the unit graph on the actual radius-one neighbours of a planar
point is a disjoint union of paths and six-cycles. It also gives a sufficient
extension condition: every available rim list has at least two colours.
Its unresolved small-list cases are handled exactly below, rather than
being called unsafe because that sufficient test fails.

The canonical hostile is a prescribed boundary colouring that fails to
extend on a graph of small chromatic number. The [palette robustness
theorem](../../zero_heavy_robustness_20261021a/palettes.md) even gives
unbounded conditional repair radius on four-colourable unit graphs.
The least-used sidecar is the *name* of a forced colour, together with the
cycle's initial colour. The current operations are path extension, closing
an actual six-cycle, and forgetting count multiplicities.

The arithmetic source is [Session51 prime offsets, section3](../../session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md).
Its exact 30030+/-17 example is recovered in section6; none of its finite
wheel predicates is silently promoted to a prime oracle.

## 2. A minimal q+2-state interface

Let C be a set of q>=2 named colours. A word L_1,...,L_m of subsets of C
specifies the allowed lists on an ordered path. Accept precisely when that
path has a proper colouring from its lists; the empty path is accepted.

**Theorem.** The minimal complete deterministic automaton for this language
has exactly q+2 states:

    dead; free; and forced(c) for each c in C.

Start in free and accept every state except dead. For an incoming list L,

    dead -> dead,
    free -> classify(L),
    forced(c) -> classify(L minus {c}),

where classify(empty)=dead, classify({c})=forced(c), and every set of
size at least two is classified free.

**Correctness.** Let R be the exact set of colours attainable at the last
vertex of a properly coloured prefix. The next attainable set is

    R'={b in L: there exists a in R with a!=b}.

If R is empty there is no continuation. If R={c}, this is L\{c}. If
|R|>=2, every b in L differs from at least one a in R, so R'=L. Thus the
displayed quotient loses no future path-feasibility information. The
initial free state correctly treats the first list as unconstrained.

**Minimality.** Every state is reachable in one step, using respectively
an empty, two-element, or singleton list. Dead is distinguished from
every other state by the empty continuation. Free and forced(c) are
distinguished by the one-list continuation {c}. The same continuation
distinguishes forced(c) from forced(d), d!=c. These q+2 pairwise different
future languages prove minimality without an unlabelled-colour quotient.

Consequently q=4,5,6 produce exactly **6,7,8 interface states**. This is
an explanation of those numbers in this particular carrier, not a claim
that the graph has that chromatic number. The two extra states are semantic
conditions (dead and free), not additional colours or interchangeable bits.

## 3. Complete obstruction chains on a path

Assume all lists are nonempty. A path is uncolourable iff it contains a
contiguous segment whose endpoint lists are singletons, whose internal
lists all have size two, and which propagates a forced colour from one
end into conflict with the other end. More explicitly, the internal lists
successively have the form {c_(i-1),c_i}; the final singleton equals the
last propagated colour.

**Proof.** Such a segment is plainly uncolourable. Conversely run the
automaton until its first dead transition. The last list is the singleton
containing the previous forced colour. Trace back to the last singleton
before it. Between these two singletons every transition must remain
forced. A list of size at least three, or a two-element list not containing
the previous forced colour, would make the state free. From free, only a
new singleton could restore forcedness, contradicting the chosen last
singleton. Hence all internal lists have exactly the stated form. If no
earlier singleton existed, the state could never have left free. QED.

This gives short obstruction certificates and an exact repair location:
one must break a forced chain, change one of its endpoint prescriptions,
or change the underlying graph. A free state erases all earlier path
feasibility dependence, but it does not erase counting multiplicity.

## 4. Closing the six-cycle requires the anchor

For an actual six-cycle, try each a in L_0. Set its first vertex to a,
start in forced(a), and process L_1,...,L_5. The chosen a succeeds iff
the final state is free or forced(c) with c!=a. The cycle is colourable
iff at least one a succeeds. This is exact and keeps the initial colour
as a sidecar. Paths need no closing test. The independent checker also
compares this procedure with complete proper-colouring enumeration.

A sharp hostile to forgetting the anchor or retaining only list sizes is

    {0}, {0,1}, {1,2}, {0,2}, {0,1}, {0,1}.       (A)

It forces 0,1,2,0,1,0 around C6, violating the closing edge. Replace the
fifth list by {1,2}; the sizes are unchanged, but

    0,1,2,0,2,1

is a proper list-colouring. Thus even one singleton on an even cycle
can matter, and neither the last free/forced tag alone nor list-cardinality
data decide closure.

These are realizable *prescribed-port* phenomena on planar unit graphs.
Place a unit regular hexagon around a centre o, use palette C={0,...,4},
and prescribe centre colour4. For each desired list L_i in (A), attach
one off-circle leaf to rim vertex i for each colour in
F_i=C\(L_i union {4}), and prescribe those leaf colours. Each new leaf
can be chosen on the unit circle about its parent, avoiding all previous
vertices and unintended unit contacts: every other centre excludes at
most two points, and only finitely many circles and points are present.
Also avoid the unit circle about o. This gives an induced unit graph
whose exact rim lists are (A). It has20 vertices (centre, six rim
vertices, thirteen leaves) and is only three-colourable, since it is a
six-wheel with pendant leaves. The same construction realizes the repaired
lists. The argument proves existence of placements; no explicit coordinate
certificate is asserted by the finite list checker.

For Hadwiger–Nelson forcing the remaining quantifier is still essential:
an off-circle context must defeat **every** allowed exterior colouring
and **every** candidate centre colour. One chosen bad list word does not
discharge that universal requirement. The existing size-two condition and
the exact automaton now form one complete path/C6 interface for any fixed
exterior assignment.

## 5. What this small interface cannot count

If v_c is the number of prefix colourings ending in colour c, extension
by a list L obeys

    v'_c = 1_(c in L) * (sum_d v_d - v_c).          (1)

Two prefixes can have the same free state and the same attainable colour
set while their multiplicities differ. For example at q=3, one list
{0,1} gives counts (1,1,0), whereas the two lists C,{0,1} give
(2,2,0). Appending {0} distinguishes their extension counts, one versus
two. No claim of an exact counting recurrence on q+2 Boolean states follows.

More precisely, counts of all singleton continuations recover the entire
q-vector: those counts are w_c=sum v-v_c; summing gives
sum v=(sum w)/(q-1), and v_c=sum v-w_c. Thus any exact representation
supporting all named singleton queries must distinguish these full boundary
count vectors. This is the relevant warning for weighted tournament-owner
recurrences: language compression and multiplicity compression are different.

There is also a precise action obstruction. Palette permutations fix dead
and free. At q=6 their eight-state action therefore has two globally fixed
states, whereas a regular F2^3-action on eight branch conics has no fixed
point for any nonidentity translation. These two eight-element carriers
cannot be identified equivariantly while preserving their stated actions.
The count eight alone supplies no bridge; named states and stabilizers do.

## 6. What 30030+/-17 actually supplies

The arithmetic is inherited, with independent small trial-division checks
here. We have 30030=2*3*5*7*11*13, the product of the first six primes.
Every integer d with2<=d<17 has a prime divisor in that product, making
both30030-d and30030+d composite. The next prime17 gives

    30013 and30047,

both prime; thus the first allowed prime endpoint on either side occurs
at offset17. The gap is34, not a twin-prime gap. The eighth prime19
already distinguishes the two signs:30011 is prime, while
30049=151*199. At prime offset23 both endpoints fail:
30007=37*811 and30053=41*733.

For a general primorial P and next prime r, an already-prime endpoint
P+/-d greater than every prime dividing P forces gcd(d,P)=1. If
1<d<r^2, that forces d prime: a composite d would have a prime divisor
below r. At30030 the bound r^2=289 is sharp, because both29741 and30319
are prime at the composite offset17^2. The implication is from endpoint
primality to small-offset primality, not conversely.

The exact connection to boundary interfaces is logical: a finite wheel
stores necessary local exclusions; an exact path state stores the
existence of a compatible full prefix. The latter may safely merge states
with identical future languages; the former cannot be called a witness
merely because no tested prime rejected it. The first six prime factors,
seventh prime17 and eighth prime19 do not identify six colours, seven
canonical forms or eight branch conics.

There is now a more concrete arithmetic-to-counting connection. In the
[fixed-target tournament recurrence](../counting/README.md), take past
odd-cycle histories A=(30013,30047) and B=(30029,30031). Both have actual
degree 60060 and two cycles. Their response f(d)=sum gcd(c,d) is 2 for
every positive odd d<59. At d=59 it is respectively 2 and 60, since
30031=59*509. Thus their normalized future states agree while the remaining
degree R<59, but split once a 59-cycle is affordable. The prime-offset
example has become an exact pair of bounded-continuation states: the map
takes endpoint integers to permutation-cycle lengths; it preserves degree
and all currently affordable gcd interactions, and discards distinctions
detectable only by a later query. No geometric or colouring claim follows
from that arithmetic map. The checker verifies the responses independently.

The source note proves finite value-membership cutoffs for Fortunate and
lesser-Fortunate values, but their general primality conjectures remain
open. This packet adds no conjectural prime-gap bound.

## 7. Sources, controls and next construction

Primary sources checked2026-10-25:

* [De Grey, The chromatic number of the plane is at least5](https://arxiv.org/abs/1804.02385)
  provides the established finite non-four-colourable unit-distance graphs.
* [Alon–Bucic–Davies, Hadwiger–Nelson Problem for Typical Norms](https://arxiv.org/html/2609.36220v1)
  distinguishes typical-norm results from the Euclidean problem. No result
  about typical norms is transported to the Euclidean circle here.
* [Cejchan–Krizek–Somer, On Remarkable Properties of Primes Near Factorials
  and Primorials](https://cs.uwaterloo.ca/journals/JIS/VOL25/Krizek/krizek3.pdf),
  section3, gives the classical short-window sufficient conditions and
  keeps the Fortune conjectures separate.

Run `python check.py` and `python -O check.py`. The script exhausts
transition identities through q8; independently enumerates proper
colourings for bounded path words and all q3 C6 list words whose nonempty
lists have size at most two; checks the minimal-state distinguishing words,
weighted-count hostile and physical-list examples; and uses independent
trial division for the stated offsets. No graph colouring claim is inferred
from the wheel, and no full geometric placement is inferred from list counts.

The concrete HN next target is to search *realizable exterior contexts*
whose induced list words reject for every exterior colouring and centre
choice. The q+2 interface makes each local continuation exact; a remaining
global search must preserve all contexts that can produce those lists.
For tournament counting the exact transferable demand is to retain the
weighted continuation vector, not to reuse this automaton without a map.
