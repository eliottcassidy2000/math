# A recurrence on involution pairs for reflection-fixed joint covers

**2026-10-29. PROVED all-size bijection and arithmetic-operation bound;
FINITE-EXACT controls. Same-agent audit. This counts fixed covers, not all
native covers, and does not supply a free pairing proof of periodic parity.**

## J1. Source, ports and inherited mechanism

Let T be a tournament on 2m vertices with a supplied fixed-point-free
anti-involution rho. Thus rho^2=1 and u->v iff rho(v)->rho(u).
Choose three sources s0,s1,s2 in distinct rho-orbits and prescribe ends

    (rho(s2), rho(s1), rho(s0)).                          (J1)

These six ports are distinct. On ordered spanning three-path covers let
J reverse every path, apply rho, and reverse the path list. J preserves
these exact ports. We count its fixed covers.

The inherited mechanism is the reversal-fixed cover coordinate d in
[Session24, PATHCOVERS section 2](../exception_pathcover_transfers_20261021/session24-pathcovers/PATHCOVERS.md).
The refinement here keeps prescribed ports and gives an explicit subset
recurrence on rho-orbits for even sources. The motivating hostile is
[TERMINAL_SELECTOR S6](TERMINAL_SELECTOR.md#s6-the-precise-parity-boundary-and-the-reflection-near-miss):
reflection has fixed covers, so those covers cannot be discarded when
trying to prove evenness.

An older comparison is [the October6 chessboard note, Proposition7.1](../chessboard_weave_20261006.md).
Its Klein-group parity proof for knight paths uses stabilizers of size at
most two. Our one reflection has genuine fixed points. The shared
operation is a group action on paths; the missing hypothesis is a free
or even-sized orbit, not a similar-looking picture. No knight-graph
symmetry is transported to the tournament.

## J2. Exact bijection

Remove the three rho-orbits of the sources. There remain h=m-3 free
two-element orbits, indexed by a set F. For A subset F define:

- O(A): the number of directed paths from s0 to rho(s2), using exactly
  one representative of each orbit in A and no other internal vertices.
- M(A): the number of directed paths starting at s1, using exactly one
  representative of each orbit in A and no others, whose last vertex v
  satisfies v->rho(v). The path of one vertex s1 is allowed.

Then the fixed joint response is

    Fix_J(T;s0,s1,s2) = sum_(A subset F) O(A) M(F\A).     (J2)

**Proof.** In a fixed cover, the third path is the rho-reverse of the first.
Their disjointness means the first contains exactly one vertex from each
rho-orbit used by this pair. Its source and end force the two prescribed
outer source orbits. The middle path is its own rho-reverse. It has even
length because rho has no fixed vertex. Its first half Q contains one
representative of every remaining orbit, starts at s1, and ends at v with
v->rho(v). The whole middle path is Q followed by rho(reverse(Q)).

Thus splitting a fixed cover gives one term of J2. Conversely join an
O(A) path, an M(F\A) half-path and their forced reflected pieces. The
anti-involution supplies all reflected arcs; the middle-edge condition
supplies the only extra join. The orbit partition gives disjointness and
full coverage. These constructions are inverse, preserving all six ports.
No division by two, unlabelled orbit identification or connectivity guess
occurs. Each reconstructed object consists of exactly three paths. QED.

## J3. Explicit recurrence and complexity

For either start s in {s0,s1}, let D_s(A,v) count paths from s to v using
one representative of each orbit in A. Initialize D_s(empty,s)=1. For
every unused orbit i and its two vertices w, update

    D_s(A union {i},w) += D_s(A,v) * [v->w].             (J3)

Then

    O(A)=sum_v D_s0(A,v)[v->rho(s2)],
    M(A)=sum_v D_s1(A,v)[v->rho(v)].                      (J4)

There are O(2^h m) states and O(2^h m^2) updates. Combining complementary
subsets costs O(2^h), so a supplied guarded source has an exact query
algorithm using O(2^(m-3)m^2) arithmetic operations and O(2^(m-3)m)
stored integers. Bit costs depend on count size. This is an all-size
bound for fixed covers. It is not an improvement on the linear bulk
transfer for the special periodic family; its usefulness is that it
applies to any supplied anti-involution without a bounded-width order.

The saved binary choice is literal: selecting one representative fixes
the other half by rho. It is available only after the anti-involution,
the paired source ports, and the central-edge condition are retained.
Forgetting that last condition would invent paths.

## J4. Controls and stopping boundary

Forty-eight seeded tournaments with supplied anti-involutions at orders
6,8,10,12, with independently varied source-orbit choices, agree with
full ordered-cover enumeration. That enumeration also checks

    number of all covers = number of fixed covers mod 2.

For the periodic actual targets and original endpoint assignment, J2
gives the following fixed counts:

| n | Fixed covers |
|---|---:|
| 8 | 2 |
| 10 | 6 |
| 12 | 26 |
| 14 | 106 |
| 16 | 428 |
| 18 | 1634 |
| 20 | 6534 |
| 22 | 25028 |
| 24 | 99892 |
| 26 | 383384 |
| 28 | 1529820 |
| 30 | 5872408 |
| 32 | 23433410 |

The first three agree with the previous independent full enumeration.
The recurrence J2--J4 is proved at all sizes; the displayed values are
FINITE-EXACT. No short scalar recurrence is inferred from this list.
The inherited full-cover parity theorem and J imply that these periodic
fixed counts are even at every size, but that implication is not a new
combinatorial pairing proof.

The six-vertex transitive tournament with rho(v)=5-v has exactly one
fixed cover for sources (0,1,2): ((0,3),(1,4),(2,5)). Thus even order plus
an anti-involution alone does not force evenness. Six is the smallest
possible order with these six distinct ports.

Neither twelve-vertex 81/27 hostile admits an anti-isomorphism to itself:
its score multiset differs from its converse score multiset. The checker
therefore rejects application of this symmetry reduction to either
source. Their genuine joint responses remain those in DEFECT_TRANSFER;
the guard is a domain boundary, not a new separation theorem for them.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/reflection_query.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/reflection_query.py
```

The precise remaining target is a local pairing or a finite closed transfer
for the fixed periodic covers that proves their evenness directly. General
actual-complement closure remains open and cannot be obtained merely by
folding sources that fail the anti-involution guard.
