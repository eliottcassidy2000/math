# An actual141 source reset retaining seven vertices

**PROVED literal source replacement and uniform-family formulas;
FINITE-EXACT local reset boundaries.** This repairs the141 source obstruction
at that particular count. It is not an arbitrary-target or recursive reset.

## 1. The retained source, vertex map, and complete continuation

Use T_a(2) from NOTE and BIRTHS, with vertices0,1 the11 prefix, vertices2,3,4
the tail01,11,10, ports a=5,b=6, and root r=7. It has H141 and, by BIRTHS,
no induced supertournament with H425 at any order.

There is nevertheless the following explicit replacement, keeping seven of
its eight named vertices and every arc between those seven:

1. Delete port b=6. The remaining seven-vertex graph W has H33.
2. Add x beating exactly old vertices3,4,5. Its new index is7 after deletion;
   the outmask is56. The count becomes69.
3. Add y with the same old out-set3,4,5, and with x->y. The outmask is again56.
   The count returns to141.

The new graph is the history lane's alternate one-apex source U, with rows

    U=(254,508,504,240,480,448,384,256,9),

under this exact persistent-label map:

|Old T_a vertex|0|1|2|3|4|5|7|
|---|---:|---:|---:|---:|---:|---:|---:|
|U vertex|1|2|3|6|7|8|0|

The two new vertices are U4 and U5. Thus the original root r and port a
survive, while port b is deliberately replaced. The old last tail vertex4
maps to U7, the converse1-ready root used by the history construction.
In the actual append ordering, add z beating only vertex4 (outmask16), then
w beating only z (outmask512). These are the two converse ready births.
The complete exact trace is

    141 -> 33 -> 69 -> 141 -> 215 -> 425.

The history lane owns the positive alternate-source/425 witness; this note
adds the explicit map from the obstructed source to that source. Every arc
between the seven persistent vertices is unchanged throughout. Counts at
the temporary graphs are recorded rather than suppressed.

For T_b(2), first apply the count-neutral edge flip on {3,5} from NOTE.
This gives T_a(2) on the identical vertex set. Then use the same replacement.
It does not contradict the all-additions obstruction: one original port
has been deleted, and the count temporarily falls.

## 2. Exact family extension and the unique equality scale

The same map exists for every k>=1. Delete port b=k+4 from T_a(k), retaining
the original root and port a. Since r precedes the entire chain and a is
the one remaining apex, the literal one-apex word is

    1 0^k 1 0^2.

The two insertions beat exactly the middle tail vertex, last tail vertex,
and port a; each outmask is `7*2^(k+1)`, and the first new vertex beats the
second. They give the words `1 0^k 1 0^3` and `1 0^k 1 0^4`, up to the
explicit relabelling of the transitive spine.

For any r>=1 the inherited g recurrence gives exactly

    H(1 0^k 1 0^r)=(2^(r+1)-1)*2^k+2^r+1.

To see it, the binary integer is `2^(k+r+1)+2^r`; peel its r terminal
zeros. The weighted quotient sum contributes
`(2^(r+1)-2)*2^k+2^r-1`, the final odd quotient contributes one parity unit,
and `g(2^(k+1)+1)=g(2^k)=2^k+1` supplies the remainder.

Consequently this same physical replacement has count trace

    34*2^k+5 -> 7*2^k+5 -> 15*2^k+9 -> 31*2^k+17.

It preserves the original count iff `12-3*2^k=0`, namely **only at k2**.
This is the exact first failed implication in attempting to turn the141
reset into a uniform return machine. The operation and named vertices
generalize; the equality of the old and new Hamilton counts does not.

## 3. Two cheap broader reset tests

The all-order source obstruction tells us an old vertex or arc must change.
The replacement above changes a port, so the following nearby possibilities
were tested without imposing readiness as a search filter.

* At k1 delete port b, leaving H19 on six vertices. Exhaust all8192 ordered
  pairs of new-vertex attachments which could restore H73. There are34
  exact count hits. Only4 have the required complete uniform-prefix response
  `34*2^k+5`; every one is merely the old T_a(k) with a neutral singleton.
  In particular no strong new source implements this two-birth return.
* At k2 the deleted-port graph has H33 on seven vertices. Exhaust all32768
  ordered pairs restoring H141. There are22 hits, including the positive
  reset above. None of the strong hits has the required response when one
  prefix vertex stays fixed and the other is inflated to an arbitrary
  transitive block. This also tests the two possible single-exception prefix
  templates. It does not exhaust arbitrary long moving-cut words.

The complete module-placeholder response from BIRTHS provides the identity
test, not equality at a few sampled prefix lengths. The finite universes
and every uniform hit are retained in `reset.json`.

## 4. The fixed-order neutral edge component

At order8, consider literal single-edge reversals that leave H141 unchanged,
starting from T_a(2). Exhausting the connected component gives exactly eight
labelled graphs and ten edges. Every graph is isomorphic to T_a(2) or T_b(2).
The operation graph is K2 times P4: the swap of prefix labels0,1 commutes
with the three-step inner path

    T_a -- T_b -- (T_b with tail labels2,3 exchanged)
        -- (T_a with tail labels2,3 exchanged).

The actual flipped edges in that path are {3,5}, {2,3}, {2,5}; the independent
K2 edge is {0,1}. This is a real operation graph, not a tournament imposed
on eight objects. It preserves H and the vertex set but retains no1-ready
source and cannot reach a425-compatible source, because all eight graphs
belong to the two all-order obstructed isomorphism types from BIRTHS.

This is a **fixed-order, single-edge-neutral** boundary. It does not rule out
simultaneous multiple-edge changes, temporary count excursions, vertex
replacement, or other source families. The positive replacement in section1
uses precisely the kind of temporary excursion that this restricted graph
does not permit.

## 5. A broader uniform replacement bank: positive operators, failed branch

The next test allows deletion of **any** one of the six nonprefix vertices
from either k1 source, followed by two arbitrary births. That is twelve
universes of8192 labelled attachment pairs. Require both strongness and the
complete uniform-prefix response `34*2^k+5`; readiness is not a search filter.

There are exactly eight labelled strong operators meeting these conditions:
four from T_a and four from T_b. They all yield the same **pointed**
isomorphism type, with the free prefix vertex fixed. One representative at
k1 has rows

    (126,236,120,16,226,72,8,109).

An actual all-k operator from T_a deletes the original root r, then uses
the two base outmasks8 and109; every comparison with the old prefix vertex
is extended uniformly to the entire prefix. The response coefficients are
`A1=73,A2=34,A_r=0` for r>=3, proving the resulting H=`34*2^k+5` for every
k>=1. The result is strong at the base, and replacing one vertex by a
transitive block preserves strongness. Thus this is a genuine uniform
count-preserving source replacement, though it replaces the marked root.
Every other hit has an explicit prefix-fixing isomorphism to this one in
`reset.json`.

The positive operator still fails the required branch at the first hostile
scale. At k2 its new nine-vertex H141 source has complete strong frontier
under425 of sizes `27,11,0`. None has count425. The low-ratio ear theorem
applies since `(9+2)*141=1551>425`; the proper-divisor condition is unchanged.
Therefore **none** of the eight operators yields a source admitting a425
extension at any order. This is an all-order source failure for the new
pointed type, not merely a failed two-birth test. It is a different actual
H141 source from the successful one in section1.

This completes the broader test with both its positive structure and its
precise stopping reason: bounded vertex replacement can preserve the entire
count family, yet the retained marked response still does not support the
next desired ternary value. A more permissive replacement or additional
retained source coordinate is needed.

Reproduction: `python -B reset_check.py` and `python -B -O reset_check.py`.
The checker builds every intermediate graph, verifies the seven-vertex arc
ledger and target isomorphism, checks the formula on k1..5, recomputes the
two complete replacement universes, and closes the neutral edge component.
The all-k identity and uniqueness of k2 follow from section2's calculation.
It additionally exhausts the twelve broader replacement universes and
closes the new pointed type's425 frontier; every retained frontier count
is independently checked by direct subset DP.

Normal and optimized receipts are byte-identical, SHA256
`e40cf5da7ca78c0b2661321319bbbb0a66f959b7adcf64efc822323ea2e82f35`.
The history lane independently verified the positive reset, all persistent
arcs, the continuation425, and the all-k formula boundary in
[`../history/reset_audit.json`](../history/reset_audit.json). That independent
audit covers sections1–2; the broader finite bank in section5 has the two
counting paths inside the local verifier, not a separate-agent audit.
