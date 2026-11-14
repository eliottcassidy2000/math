# The 67/69 arc as a two-sided programmable exact selector

**Proved:** a local arc-rooted criterion which transports a unit odd
step through arbitrarily many independently variable transitive runs;
and a sublinear upper bound for every fixed exterior satisfying the
criterion. No universal placement claim.

## The literal hidden single-arc pair

The inherited [two-visit pair](../chanpak_j11_hybrid_20261114/TWO_VISIT_WINDOWS.md)
has rows A=(2,12,17,37,75,23,47) and B=(2,12,17,37,43,71,31).
Their displayed difference is reversal of triangle 4,5,6. Relabel
B by swapping 2 with 3 and 4 with 5. In the new labels it is exactly
A with the arc 3->2 reversed. Marks 0 and 1, and vertex 6, are fixed.
This is a statement about the actual graphs, not equality of a
score or path-count summary.

Write x=2 and y=3 for the testers and C={4,5,6} for the exterior.
Originally the two vertices 0->1 are the sink-side switch tape:
both lose to every vertex of C. Vertex 0 loses to both testers,
and vertex 1 beats both. A source-side tape can be inserted before C.
Inside either tape any tournament is permitted; each vertex
independently has one contact bit, meaning it beats both testers
or loses to both. All arcs between the parts obey W -> C -> G
(and W -> G), with W the source tape and G the sink tape.
The testers have outgoing subsets U,V of C, respectively.
The two arms differ only in y->x versus x->y.

## A local skew response controls both entire tapes

Define D_C(I,O) as the number of ordered pairs of nonempty directed
paths covering C disjointly, with first path ending in I and second
starting in O. In the tape template, if an automorphism of C sends
U to V, then

    H(x->y) - H(y->x)
        = H(W) H(G) [ D_C(V,U) - D_C(U,V) ].              (1)

Empty tapes have H=1. For proof, all paths not using the toggled
edge cancel under the identity map. In paths using it, contract
the oriented edge into an ordered domino D, with separate allowed
predecessor and successor sets. The two versions have respective
exterior entry/exit masks (C minus V,U) and (C minus U,V).
For any fixed before/after allocation of W,G, the parts on either
side of D necessarily follow the order W,C,G. Only when C directly
touches D on that side does its orientation matter to that entry
or exit constraint. If at most one such constraint is active, the
automorphism (or its inverse) matches the marginal counts for each
fixed number of vertices of C before D. The tape constraints,
including which contact bit may touch D, do not see the change.

If both constraints are active, all of W is before D and all of G
after D, so they contribute H(W)H(G). The exterior contribution is
D_C(C minus U,V)-D_C(C minus V,U). The unrestricted incoming
terms D_C(C,V) and D_C(C,U) are equal by assumption. This reduces
the difference to the bracket in (1). These cases partition all
ordered domino uses; no arithmetic walk count is substituted for
a literal path count.

For the original C use a->c->b->a, U={a}, V={b}. The bracket is 2:
in a two-part cover with the second path starting at b, its old
predecessor c in the triangle can be last in the first part for
both possible split sizes; the analogous predecessor of a is b,
which does not meet the other allowed entry condition. This
explains the cancellation across all higher run / visit terms.

## The gap survives changed cores, and arbitrarily many lengths

Three portable exterior choices are below. In each case the
automorphism sends U to V and the bracket is exactly 2. U,V are
bitmasks; p is the vertex-image permutation.

| Exterior outgoing rows | U,V | p | tape-free H pair | H pair for G=01 |
|---|---|---|---|---|
| (4,1,2) | 1,2 | (1,2,0) | 11,13 | 67,69 |
| (6,12,24,17,3) | 1,16 | (4,0,1,2,3) | 89,91 | 481,483 |
| (14,28,104,112,37,67,19) | 9,68 | (2,3,5,6,4,0,1) | 2325,2327 | 13397,13399 |

The last choice is a seven-vertex exterior, not an inflation of
the triangle. The bracket and the physical counts are independently
computed by pair and path subset DP in pair_structural.py.

If both tapes are transitive, (1) is exactly +2 for **every
pair of binary words**, not just words grown by repeating a single
letter. For example W=0^a1^b and G=1^c0^d1^e, with all five lengths
positive, are five independently variable lengths in a physical
exact pair. At all lengths one the counts are 1489 and 1491.
The successive runs of opposite contact type cannot merge in the
marked construction. More alternating runs give arbitrarily many
independently chosen lengths. This does not assert distinct scalar
outputs for every distinct parameter packet. Conversely for a
fixed bracket 2, any nontransitive tape introduces a directed
triangle and has H>=3, so synchronously replacing an entire tape
by an arbitrary tournament does not retain gap two.

A transitive prefix/suffix of the exterior whose vertices have
common tester contacts may also be painted independently, but
then these pieces literally concatenate with W/G. Treating those
pieces as two further tapes is a parameterization, not an enlargement
of the family.

## Whole programmable families still cannot fill all heights

Fix any one exterior and a genuinely nonzero skew as above, so
U,V are distinct masks of equal cardinality. Choose a in U minus V
and b in V minus U. Strip leading 1s from transitive W and trailing
0s from transitive G: they are universal sources and sinks, so do
not affect H. Denote X a count bound, not a tournament.

For G alone, suppose the final two 1 contacts are at positions j,n,
necessarily with the last gate at position n after normalization.
Take disjoint odd cycles x->a->(a subset of G ending at j)->x
and y->b->(a subset of G ending at n)->y. Before j each gate
has three independent states (first cycle, second cycle, unused);
between j and n it has two (second cycle, unused). The endpoint is
one gate per cycle, so the number of extra vertices on each cycle
must be even. With k=j-1 three-way slots and l=n-j-1 two-way slots,
at least 3^k 2^l/4 assignments meet the parity. For l=0 this is
the elementary parity average (3^k+2+(-1)^k)/4; for l>=1 first
enforce even parity in the first slot class, then use half of the
two-way assignments in the middle to enforce the second. In the
[odd-cycle-packing count](https://arxiv.org/html/2412.10572v3#S3.SS4),
each such disjoint pair has weight 4. Thus H>=3^k 2^l.
There are at most 2^k choices at those k positions and the
intervening positions are fixed zeros. Summing over k and l with
that bound gives O(X^log_3(2)) possible normalized G words,
including the unique word per length with just one live 1:
using one tester cycle gives a logarithmic length bound. The same
argument in the other direction holds for W and its first two 0s,
using b as a triangle/exterior return to x and a to y. It works for
any such fixed C because those two exterior vertices are distinct.

When both tapes have at least two live contacts, pin the first two
0s of W and last two 1s of G. Let s be the number of W positions
after its second 0 plus G positions before its penultimate 1,
and p the number of positions in their two internal gaps.
Join both banks into each of two disjoint cycles, of the form
tester -> W subset -> one distinct exterior vertex -> G subset
-> tester. Each starts at one pinned 0 and ends at one pinned 1.
The s positions choose either cycle or neither. Each of the p
internal positions can choose one specified cycle (the one with
the earlier W start or later G end) or neither. Both cycles have
base even length 4. For s>=2 reserve the last two shared positions:
for any choices at the other positions, at least one assignment
of these two makes the added counts odd in both cycles. Thus there
are at least 3^s 2^p/9 valid pairs, and H>=4*3^s*2^p/9.
For a fixed s,p there are at most (s+1)(p+1)2^s normalized tape
pairs, splitting the two lengths between the two tapes.
Summing with 3^s 2^p<=9X/4 bounds them by

    O(X^log_3(2) (1+log X)).

For s<2 all variable bits outside the internal fixed runs have
only O(1) choices; one tester cycle bounds each length by O(log X).
If one tape has only one live contact its word is unique per
length and the other tape obeys the one-bank bound, again
O(X^log_3(2)(1+log X)). Empty normalized tapes fit the same bound.
The upper bound transfers to the number of *counts* or occupied
windows below X. Indeed all counted odd-cycle collections exist as
induced configurations in both orientations of the tester edge,
and do not use that edge.

A finite set of fixed exterior cores, even with arbitrarily many
free tape runs, is therefore not a cofinite selector. To get a
uniform proof from this construction one must select changes of
the exterior or a different absolute placement mechanism. Equation
(1) states the exact invariant to preserve in that search.

Replay: python3 pair_structural.py.
