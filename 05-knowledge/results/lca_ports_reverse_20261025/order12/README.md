# A three-chord selector closes both order-twelve residual cells

2026-10-25. **PROVED:** both complete order12 source families left by the
inherited base bank have literal three-black repairs. The selector reads
only three actual chords. A blue failed word forces precisely those
chords through its four three-cycles. **FINITE-EXACT:** the standalone
checker certifies the full affine families, every symmetric hostile
completion, and explicit omitted-guard controls. The
[independent review](INDEPENDENT_REVIEW.md) accepts the proof and supplies
a separate all-completion control.

## 1. Inheritance and exact source universe

The closest mechanism is the complete
[low-order affine audit](../../port_residual_integral_20261025/LOW_ORDERS.md),
with its [independent review](../../port_residual_integral_20261025/LOW_ORDERS_REVIEW.md).
It leaves exactly two order12 base cells: the old first word can be
blue, while its second and third words are always black. The canonical
hostiles are the two literal source rows preserved there. The corrected
near miss is to treat a legal word or a correct edge XOR as proof of
blackness. The least-used sidecar is the *source plus failed-word*
symmetry: two actual anti-reversals compose to an automorphism of the
same source.

Anchor: repair both complete cells. Niche: remark a failed word and
reuse the existing base bank. Wildcard: exploit the four three-cycles
to obtain a local three-chord rule. The concept board is literal ports /
two anti-reversals / three forced chords / commuting path-edge changes /
exact clone labels / source-family quantifiers.

Let P=(0,1,...,11) be a directed closed blue Hamiltonian word in T.
Here blue means that r_P(v)=11-v is an anti-automorphism. Let

    N^-_T(0)={11} union I,
    I={2,7,9} or I={6,7,9},
    1->10,       2->10.                                  (S)

The endpoint outdegrees are(7,4). The last two conditions say open
middle and x->c, with x=2,c=10. All other arcs are unrestricted subject
to P-blue. There is no source isomorphism, strongness, mixedness,
module, disjoint-port or clone condition.

For a directed word Q, K_Q retains every Q path arc and reverses every
other arc. Blackness is tested at that word's actual source. Put
B=K_P(T). A word has the same blue/black status at both ends of its own
K edge, because its positional reversal preserves its path-edge set.

The two source families each have affine dimension18. This dimension
is a finite exact check of the source universe, not a premise of the
proof below.

## 2. The old ledger and what actually failed

The inherited open-penultimate ledger is

    H  =(0,1,10,11,2,3,4,5,6,7,8,9),
    H2 =(1,10,0,11,2,3,4,5,6,7,8,9),
    H3 =(1,2,3,4,5,6,7,8,9,10,0,11).                    (O)

All three words are directed, their path-edge XOR is E(P), and H2,H3
are black throughout both complete source families. These facts are
part of the audited low-order bank; this packet independently checks
them on the full18-dimensional affine cells. The H2 obligation is
also the n12 open-penultimate relaxed affine system, which is
inconsistent without even imposing a three-port row bound.

Only H may be blue. Thus (O) already repairs every source where H
is black. It remains to repair its blue locus while preserving the
original target B.

## 3. The failed symmetry forces the replacement joins

If H is blue, both r_H and r_P are anti in T. Their product is the
literal automorphism

    g=r_H r_P=(0 6 2)(1 7 3)(4 10 8)(5 11 9).           (1)

This is obtained by directly reading the two word reversals.
In the two source families, the entire outgoing row of11 is respectively

    {0,2,4,9}       or       {0,2,4,5}.                 (2)

Consequently 8->11,3->11 and6->11 are true in both families. Transport
the first two arcs through (1):

    T(4,9)=T(10,5)=T(8,11)=1,
    T(7,5)=T(3,11)=1.                                  (3)

Hence H-blue implies the three-chord condition

    G:  4->9, 10->5, 7->5.                             (4)

Nothing here uses a regularity count or an isomorphism representative.
The precise positions of all four cycles and both endpoint rows matter.

## 4. A complete guarded three-black repair

In either full source family (S), assume only G; H itself need not be
blue. Use

    Q1=(0,1,2,3,4,9,10,5,6,7,8,11),
    Q2=(0,1,2,3,4,9,10,5,7,8,11,6),
    Q3=(0,1,2,3,4,5,7,8,9,10,11,6).                    (N)

### Legality at the actual intermediate tournaments

Q1's new joins are4->9,10->5 and8->11. The first two are G, and the
third follows from (2). Every other join is a P edge.

At S1=K_Q1(T), Q2 removes6 from between5 and7 and places it after11.
Its new joins are5->7 and11->6. Both pairs were omitted by Q1, so
they are the complements of the original7->5 and6->11. The first
is G and the second follows from (2). All other joins are retained.

Q3 restores the original P joins4->5,8->9,10->11. Each was omitted
by both Q1 and Q2, so two complements restore it. The new joins5->7
and11->6 are retained from Q2; its other joins are retained as well.

This also shows that each of the three guards in G is necessary for
this particular new ledger: omitting the first or second can break a
Q1 join, while omitting the third can break the Q2 join5->7.

### Exact target by two independent path changes

Define unordered-pair sets

    X={45,89,10-11,49,5-10,8-11},
    Y={56,67,57,6-11}.

The words satisfy

    E(Q1)=E(P) XOR X,
    E(Q2)=E(P) XOR X XOR Y,
    E(Q3)=E(P) XOR Y.

Thus their XOR is E(P). Three full complements together with that
identity give exactly K_P(T), on the original labelled vertices.
The edge calculation does not substitute for the legality proof.

### Three short blackness witnesses

For Q1 at T, the actual reversal sends the ordered pair(0,9) to the
reversed pair(10,11). The source port and path bits are

    T(0,9)=0,        T(10,11)=1.

For Q2 at S1, the actual reversal sends(0,4) to(5,6). The pair04 was
omitted by Q1 while56 was retained. Since4 is not a port,

    S1(0,4)=0,       S1(5,6)=1.

For Q3 test at B, the other end of its own K edge. Its reversal sends
(1,10) to(2,11). Both pairs are nonpath pairs of P. Open middle gives
T(1,10)=1; the port9 in I gives11->2. Therefore

    B(1,10)=0,       B(2,11)=1.

Every word has a literal arc contradiction to its anti-reversal.
All three are black. This proves (N) throughout the guarded source
families; the affine check finds dimension15 for each such family.

## 5. Three-bit decision and complete coverage

Read the three chords in G.

- If all three hold, use (N).
- Otherwise H cannot be blue, by (3). Use the old ledger (O).

These alternatives cover every completion of both original source
families, with three black moves in either case. No search or full
anti-reversal test is required to choose the repair.

Together with the independently audited order12 base cells and the
already proved I={h,9,10} residual repairs, this closes **every order12
closed blue endpoint-(7,4) marked path** by three black moves. The
inheritance partition is the one used in
[PORT_FRONTIER.md](../../port_residual_integral_20261025/PORT_FRONTIER.md):
the base handles every open middle and every closed middle except the
closed endpair; the h2/h3/h4/h>=5 families handle that entire endpair
from n>=12 (h2 already from n>=11). No new source filter is introduced.
This is an order12 consequence, not an all-layer connectivity theorem.

## 6. Independent five-black fallback retained

There is another exact route that helped isolate the new symmetry.
When H is blue, remark the source using H instead of P. This is a
coordinate change, not a K move. H is closed because9->0. Its first
vertex still has outdegree7, so anti-reversal makes its last vertex
have outdegree4. Its new incoming interior port set is

    {3,4,9}   for original I={2,7,9},
    {3,8,9}   for original I={6,7,9}.

Both are outside the order12 failed-cell list. The existing base bank
therefore replaces the H edge by three black moves in the original T.
Appending H2,H3 gives five black moves to the original B. The necessary
remarked middle and chord are exactly T(1,8) and T(10,8).

The checker retains the full literal replacement words in FALLBACK:
one branch for1->8, and two for8->1 split by10 versus8. It checks all
six resulting affine families, including their exact intermediate
targets. Their dimensions are2,1,1 for each original I, covering its
whole3-dimensional H-blue family. This fallback is independent of the
new three-chord repair, but is weaker in path length.

## 7. Controls, reproduction, and scope

[check.py](check.py) is standalone standard-library Python. It imports
no earlier producer or solver. It uses one independent forward-arc bit
per unordered pair, imposes exactly (S), and checks consistency before
any implication. The complete affine checks cover:

- each original18-dimensional source cell, including legality and
  target of (O), and universal blackness of H2,H3;
- each15-dimensional G cell, including all three new blackness
  equations and the displayed individual witnesses;
- each3-dimensional H-blue cell, with all16 literal completions
  independently replayed using integer row bits;
- all six five-black fallback cells;
- six omitted-guard controls, one per guard in each source cell.

The omitted-guard controls keep the other two guards true. They exhibit
the actual illegal Q1 or Q2 join in (N), and the old ledger repairs those
same sources with three black words. The two inherited literal blue-H
hostiles remain blue under (O) and are repaired by (N). Thus neither
the hostile nor the guard is hidden by changing the source universe.

The dimensions give eight H-blue labelled sources per original cell.
This explains why the failed symmetry is small; it was not assumed
as a restriction on the full source families. The total two-cell
universe has2*2^18 marked labelled sources. A short block-menu search
suggested (N): it considered1,992 words formed by at most four retained
forward intervals, and the first eligible middle edge-path produced
the repair. The proof uses only (1)--(4), the displayed words and the
arc witnesses, not that discovery search.

From the repository root:

    python 05-knowledge/results/lca_ports_reverse_20261025/order12/check.py
    python -O 05-knowledge/results/lca_ports_reverse_20261025/order12/check.py

The [normal](normal.json) and [optimized](optimized.json) receipts are
byte-identical; all gates are explicit exceptions. All12 vertices and
all intermediate arcs remain literal. No clone merger is approximated
or discarded. Neither a general all-order three-chord theorem nor full
black-core connectivity is claimed.
