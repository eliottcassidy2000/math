# Every closed endpoint-four line is repaired from order twelve

2026-10-25. **PROVED by audited inheritance and the six new repairs:**
for every tournament T of order n>=12, every directed closed blue
Hamiltonian word P with endpoint outdegrees(n-5,4) has a three-black
detour to its exact labelled K_P(T). No strongness, mixedness, module,
port-overlap, or clone restriction is imposed. The complete proof uses
the predecessor's declared finite affine certificates as well as
structural infinite-order arguments. Full black-core connectivity is
still **OPEN**.

## Statement and literal semantics

For a directed word Q, K_Q retains its path arcs and reverses all other
arcs. Blue means that Q's positional reversal is an anti-automorphism
of its actual source; black means that it is not. Closed means that
the last vertex points to the first. We prove that three directed
words Q1,Q2,Q3, each black at its successive source, satisfy

    K_Q3(K_Q2(K_Q1(T))) = K_P(T).

The equality is on every original vertex and arc. In particular this
is valid when intermediate clone classes split or merge: the words
are applied to the full tournaments, with no quotient operation.

Label P by0,...,n-1 and set a=0,b=1,x=2,y=n-3,c=n-2,z=n-1.
The exact endpoint row is N^-(a)={z} union I, with |I|=3 and
I contained in{2,...,n-2}. Reflection determines the other endpoint
row. The middle b versus c is a separate orientation from closure
of the full path.

## Exhaustive proof

For n>=15 the conclusion is the audited
[previous exhaustion theorem](../port_residual_integral_20261025/PORT_FRONTIER.md).
Its base bank handles open middle, closed middle with c outside I,
and closed middle with c in I but y outside I. The only remaining
port sets are I={h,y,c},2<=h<=n-4. Its h2 repair applies from11,
and its two h3 branches, two h4 branches, and all h>=5 repairs apply
from12. These thresholds and full source universes were explicitly
audited in that packet. Neither failed-prefix guards nor disjoint
reflected ports are added to this exhaustion.

At n=12,13,14 the
[complete low-order affine audit](../port_residual_integral_20261025/LOW_ORDERS.md)
and its [independent review](../port_residual_integral_20261025/LOW_ORDERS_REVIEW.md)
test the same base cases on all completions of every exact port row
and relevant middle/chord branch. At13 every base cell passes. The
residual I={h,y,c} cells at each of these three orders are already
covered by the just-stated structural repairs, whose thresholds are
all satisfied.

The only failed base cells are the following six. The new theorems
cover their full source universes, not only the sources making a
particular old word blue.

| Order | Exact base cells | Exhaustive source-level selector |
|---|---|---|
|12|I={2,7,9} or{6,7,9};1->10 and2->10|If4->9,10->5 and7->5 all hold, use the new rectangle; otherwise use the inherited ledger.|
|14|I={3,5,10} or{3,7,10},2->12; I={5,9,10} or{7,9,10},12->2; all have12->1|If5->3, use the common new rectangle; otherwise use the inherited branch-specific ledger.|

The [order12 proof](order12/README.md) and
[independent review](order12/INDEPENDENT_REVIEW.md) establish the
three-bit selector. If its old first word is blue, the composition
of the two anti-reversals forces all three guards. Under those guards
the new joins are directed and three fixed ordered-pair contradictions
prove blackness. The alternate ledger works on every other completion.

The [order14 proof](order14/README.md), its
[one-bit extension](order14/GUARD_EXTENSION.md), and
[independent review](order14/INDEPENDENT_REVIEW.md) give the second
selector. The failed final word's two seven-cycles force5->3, but the
new rectangle only needs that one arc, not the full symmetry. Three
fixed ordered-pair contradictions prove blackness uniformly on its
27-dimensional guarded half-cell. The other half uses the old ledger.

Each selected ledger starts at T and has three moves. In each case
the exact unordered-edge identity is

    E(Q1) XOR E(Q2) XOR E(Q3) = E(P).

Three complements give the required target; directedness and blackness
are proved separately at their actual sources. No portions of different
ledgers are spliced. This proves the claim at12 and14, hence at every
n>=12. The [assembly review](PORT_FRONTIER_REVIEW.md) audits the
inheritance partition and its thresholds.

## What this closes, and what remains

Any walk can replace a closed blue endpoint-four edge at n>=12 by
this black segment while preserving its endpoints. This removes the
six old word failures as possible obstructions in that layer. It does
not establish that arbitrary blue words are closed, belong to this
degree layer, or can be moved into it through black steps.

The earlier pure-blue order7 target remains a real obstruction to an
unfiltered all-order replacement theorem; its exact witness is retained
in the [previous frontier](../port_residual_integral_20261025/PORT_FRONTIER.md#4-genuine-small-order-hostile-and-scope-boundary).
The threshold12 is proved, not claimed optimal. Order11 still has a head
outside its previously certified subcases. The next degree layer has
four incoming interior ports, so the older contradiction forcing at
least four ports no longer rules it out.

There is also an exact open/closed transport: if P has endpoint scores
(n-1-k,k), its K target has scores(k+1,n-2-k). The closing pair flips,
so every blue line has a closed marked end, but its layer changes by
k->n-2-k. Open-four therefore becomes closed layer n-6; the present
fixed-layer theorem does not automatically cover it. Proving all closed
layers would remove this open/closed issue. Global connectivity still
requires the pure-blue-deleted metagraph condition identified in the
[assembly review](PORT_FRONTIER_REVIEW.md).

The next structural question is whether a change of marked path can
move a general admissible blue edge into a proved layer while preserving
a usable literal route to the original target. An isomorphism, correct
edge XOR, or degree count alone does not provide that route.
