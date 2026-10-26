# Independent review of the assembled endpoint-four theorem

2026-10-25. **ACCEPTED / PROVED with declared finite certificates:** the
assembled [PORT_FRONTIER.md](PORT_FRONTIER.md) establishes, for every
n>=12, every directed closed blue marked Hamiltonian path with endpoint
outdegrees(n-5,4), a three-black detour to its exact labelled K target.
The proof is an exhaustive inheritance assembly, not an extrapolation
from the new small-order controls. No additional clone, mixedness,
strong-target, port-overlap or isomorphism filter is present.

**Subsequent integration:** [HEAD_ASSEMBLY.md](HEAD_ASSEMBLY.md) closes
the finite head discussed below and gives n>=8. Section5's order11
paragraph records the boundary at the time of this n>=12 review; it
is superseded by that separately audited assembly. The layer transport
and higher-layer/global connectivity obligations remain unchanged.

This review reads the already audited sources and their exact domains.
It does not repeat a large enumeration. The next-boundary discussion
below records the necessary degree change when passing to the closed
end of an open presentation.

## 1. Quantifiers and actual operations

Take an arbitrary labelled tournament T and an arbitrary directed
Hamiltonian word P satisfying closure, its actual positional anti
condition, and the stated endpoint degrees. Labelling vertices by
their positions in P is a coordinate choice. Every selected word is
then transported back along that same labelling. K is equivariant
under this coordinate change, so the target remains the original
labelled K_P(T).

The conclusion quantifies every such marked word, not merely one
chosen word in each class or one source per port set. It asserts the
existence of a three-word ledger for each input; it does not assert
that every Hamiltonian word or one universal ledger works.

Closure itself makes the source tournament strong. There is no
additional strongness test and, especially, no restriction requiring
the target to be strong. Open or closed **middle** means the independent
orientation of b versus c, and is not confused with closure of P.
The intermediate words need not remain in endpoint layer4 or remain
closed. They are tested at their own successive source tournaments.

All ledgers retain all original vertices and arcs. Clone classes may
change across K moves, but no quotient is applied in the proof. The
identity on unordered edge sets and the separate direction tests imply
literal target equality. Consequently no assumption that a clone
quotient commutes with a K operation is needed.

## 2. Exhaustion at n>=15

The predecessor
[PORT_FRONTIER.md](../port_residual_integral_20261025/PORT_FRONTIER.md)
and its source audit were read at their exact scopes. Its base theorem
covers these mutually exhaustive alternatives outside the endpair:

1. open middle, with arbitrary incoming port set I;
2. closed middle with c outside I;
3. closed middle with c in I but y outside I.

Since |I|=3, the only complement is a closed middle with
I={h,y,c},2<=h<=n-4. The first port cannot be b because a->b is a
path edge. No condition such as disjoint I and r_P(I) is introduced.

For h=2 the second-port theorem applies from n>=11. For h=3 and h=4,
the two directions of the actual chord1 versus3 are separately covered
from n>=12. All h>=5 through n-4 are covered from n>=12. Hence every
threshold is met at n>=15 and all possible h are exhausted.

The base's arbitrary-k cyclic-cut condition is n>4k-2, which becomes
n>=15 at k=4; it has not been lowered by assertion. The independent
finite systems at n15..26, followed by the structural endpoint-row
chains for n>=27, are explicit dependencies of that base. This
computational head is part of the proof, not sampled evidence for its
unbounded tail. The later failed-prefix/chord residual descriptions
are not extra source filters: the completed h families cover their
whole port domains without those conditions.

## 3. The three low orders and the six old failures

The complete [low-order audit](../port_residual_integral_20261025/LOW_ORDERS.md)
and [independent review](../port_residual_integral_20261025/LOW_ORDERS_REVIEW.md)
partition each marked source universe by its exact three-element I
and the two bits T(b,c),T(x,c). Every cell is nonempty, and every
completion remains represented. The counts are:

| Order | All branch cells | Base cells | Closed-endpair residual cells |
|---|---:|---:|---:|
|12|336|322|14|
|13|480|464|16|
|14|660|642|18|

The counts follow from4*binom(n-3,3) and2(n-5), not from deleting
unfavourable source types. The residual cells are precisely the same
I={h,y,c} families described above, whose thresholds are all satisfied
at n12,13,14. At n13 every base cell already passes.

At n12 exactly two complete base cells remain. The new
[three-chord proof](order12/README.md) and
[independent review](order12/INDEPENDENT_REVIEW.md) cover their entire
18-dimensional source families:

    I={2,7,9} or{6,7,9},        1->10, 2->10.

If the old first word H is blue, g=r_H r_P forces the three chords
4->9,10->5,7->5. The new ledger is legal and black whenever these
three chords hold; if one fails then H is black and the old ledger
works. This is implication plus contrapositive, not an unjustified
equivalence between the three guards and H-blue. The selector starts
at the original T and chooses a complete ledger there.

At n14 the four failed cells are exactly those listed in the
[order14 proof](order14/README.md), with the original closed-middle
and x,c branch bits retained. The
[one-bit extension](order14/GUARD_EXTENSION.md) and its
[independent review](order14/INDEPENDENT_REVIEW.md) show that the
new ledger works throughout T(3,5)=0, while T(3,5)=1 implies the old
final word F is black. The inherited first two words are universally
black on each full source cell. Each bit half is a full nonempty
27-dimensional family. The proof does not restrict the source to
the four individual F-blue tournaments.

Thus all six old failures are covered on their entire original
affine cells. At every branch, legality, actual-source blackness and
the exact target are separate proved obligations. Testing the last
word at the target B is valid because its own K operation preserves
its anti status. No partially executed old ledger is spliced into
a new one.

Combining12,13,14 with the inherited n>=15 theorem proves exactly
the assembled n>=12 statement.

## 4. Open/closed transport changes the degree layer

This distinction is essential when stating what remains. Let a directed
blue word have source endpoint scores

    (n-1-k,k).

Under its K move the first endpoint gains its retained outgoing path
arc relative to full converse, and the last loses its retained incoming
path arc relative to full converse. Therefore the new endpoint scores
are

    (k+1,n-2-k).                                      (1)

The closing pair is a nonpath pair and reverses, so this exchanges
open and closed presentations and sends the layer parameter to

    k' = n-2-k.                                     (2)

This agrees with
[THM-4729, L1, one closed end per line](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
It can also be checked directly from the definition of K, without any
class quotient.

In particular an **open** endpoint-four presentation with scores
(n-5,4) has a closed partner with scores(5,n-6), in layer n-6.
It remains layer4 only at n=10. The present theorem at n>=12 does
not therefore settle every open endpoint-four presentation. Conversely,
it does cover the open partners of its proved closed layer, which
have scores(5,n-6).

A black detour reverses because every word is directed and black at
both ends of its own K edge. Thus proving redundancy for **every**
closed layer would cover every blue line by choosing its closed end.
There is no need for a separate universal open-line theorem once all
closed layers are solved. At this fixed layer, however, the degree
change in(2) cannot be ignored.

## 5. The next genuine boundaries

There are three distinct remaining obligations.

**The finite head at order11.** This assembly does not prove the full
closed endpoint-four layer at11. The closed-middle endpair subcase
y,c in I already has the196,608-presentation finite certificate in
[THM-4732, section3](../../../01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md).
That is a specified subcase, not the whole order11 layer. The genuine
remaining finite-head question is the part outside already certified
subfamilies; it should not be described as if the entire endpair
subcase were still untouched. No optimality claim for threshold12
follows, and a failed fixed ledger is not a disconnection witness.

**The next unbounded port universe.** The first next layer is endpoint
five: closed blue sources with scores(n-6,5) and an arbitrary
four-element interior incoming set I. Existing consecutive-tail and
protected/unprotected-terminal-port results cover some such sets, but
not that complete universe. Merely replacing |I|=3 by4 in the base
proof fails at a specific point: the three parity-pair demands plus
one extra forced port require at least four ports, which contradicts
three but not four. A new forcing argument, repair, or reduction
between layers is needed. Beyond endpoint five, a proof covering
only finitely many fixed layers still leaves arbitrary k, including
the large closed layers reached from open small-k presentations by(2).

**The global connectivity step.**
[THM-4732, Lemma2](../../../01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md)
gives the precise criterion: the line metagraph with pure-blue vertices
deleted must be connected, and every blue line between non-pure-blue
vertices must have black-connected ends. The present theorem settles
a specific part of the second condition. It supplies neither the first
condition nor all the remaining layers of the second.

The pure-blue deletion is load-bearing: the explicit order-seven
C3 -> 1 -> C3 target retained in the predecessor frontier has no
incident black word. That hostile remains valid even though the
new theorems repair the separate order12/14 ledger failures.
Any eventual all-layer statement must preserve the appropriate
non-pure-blue endpoint domain instead of extending an unfiltered
small-order claim across this obstruction.

**Disposition:** approve the n>=12 closed endpoint-four theorem and
its all-marked-word quantifier. Keep the degree transport(1)--(2),
the order11 finite-head distinction, arbitrary higher port sets, and
the pure-blue-deleted connectivity criterion explicit when describing
the remaining all-order problem.
