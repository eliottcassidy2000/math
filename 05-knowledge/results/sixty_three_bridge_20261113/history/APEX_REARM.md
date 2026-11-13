# Same-count apex reflection, a three-arc rearm, and the63 fibre

**Status: PROVED** for the all-word operation and infinite family below.
**FINITE-EXACT**, using the inherited complete two-class classification,
for the minimum-width statement at order8. This is a new retained-port
use of an old converse symmetry; neither the existence of H=63 nor the
one-apex counting formula is new.

## 1. Inheritance and the changed question

The closest proved source is the one-apex binary signature in
[single-core-cycle-count.md](../../../variables/single-core-cycle-count.md),
already listing1001100 and1100110 in May2026.
[THM-344](../../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md)
classifies the two eight-vertex H=63 tournaments; both have complete
odd-cycle conflict graph K31. The same signatures are credited again in
[THM-4748](../../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md).

The hostile example is arithmetic division63/3=21: a factor in the
integer count need not be a removable triangle SCC. The corrected near
miss is to preserve H but forget which root is ready. The least-used
sidecar here is the pair of endpoint deletion counts of a binary word.

The concrete new question is: can the SAME named root, its sole incoming
arc, every other root incidence, and the entire transitive carrier all
remain fixed while a same-H modification changes readiness from0 to1?
The answer is yes, uniformly; at63 three arc changes are necessary and
sufficient within the eight-vertex H=63 fibre.

## 2. Typed state and the all-word operation

Let w=w0...w_(m−1) be a binary word, m≥2. Its apex tournament A(w) has
named vertices a,v0,...,v_(m−1), with vi→vj for i<j and

    a→vi iff wi=1.

Write h(w)=H(A(w)). For a named root r, use the state

    (h,s,q,F)=(H(T), #HPs starting at r, H(T−r), E_r(2)),

where E_r(2) counts ordered two-slot path covers whose first path is
nonempty and ends at r; the second slot may be empty. A root is
**unit-ready** when its indegree is1 and q−s=1.

Define w* by w*_i=1−w_(m−1−i). On the SAME vertex names replace A(w)
by A(w*). This changes only arcs incident to a; it preserves every arc
of the named transitive carrier. The map a↦a, vi↦v_(m−1−i) identifies
A(w*) with the converse of A(w), so

    h(w*)=h(w).                                             (1)

Assume w0=1 and w_(m−1)=0, and retain r=v0. Both tournaments have
exactly the SAME root star: a→r and r→vi for i>0. In particular the
same named predecessor a is retained. They are both strong: the path
a→v0→...→v_(m−1)→a is a spanning directed cycle.

Their complete rooted states are

    A(w):  (h, h(w[1:])−w1, h(w[1:]), h+w1),
    A(w*): (h, h(w[:-1])−(1−w_(m−2)),
               h(w[:-1]), h+1−w_(m−2)).                    (2)

Here h of a one-letter word is1, so m=2 is included. To prove(2), delete
r. The remaining carrier is transitive, hence an HP beginning at a is
possible in exactly one way and exists iff w1=1. Since r beats every
other remaining vertex except a, deleting r from an HP that starts at r
gives every HP of T−r except these w1 paths. Thus q−s=w1.

For F, concatenate the first slot ending at r with the second slot.
The successful concatenations are in bijection with all HPs, by cutting
after r. A failed concatenation must have its second slot begin at a.
Its first slot must then be the singleton r: a is the only possible
predecessor of r and is already in the other slot. Therefore the
number of failures is precisely the number w1 just counted, giving
F=h+w1. Apply the same argument to w*, and use converse symmetry after
deleting the last letter of w to identify its deletion count.

Consequently every word beginning10 and ending00 admits the exact
same-H rearm

    q−s: 0 → 1,                                            (3)

while retaining the root, its complete star, its predecessor, and all
carrier arcs. This is an involutive move with an explicitly checked
entry condition. It does not claim every H has such a source word.

The mutation count is the number of indices i for which
w_i=w_(m−1−i). It can grow with m. Starting10 and ending10 is a genuine
boundary: the mirror retains defect0. For example1010 has H=9 and is
fixed by the mirror. An empty word, a different endpoint pattern, or a
different retained root is not covered by(2).

## 3. A constant-three-arc infinite family

For b≥1 put B=2^b and choose

    w_b = 1 0^b 1^b 00,
    w_b*= 11 0^b 1^b 0.

These words differ at exactly the three zero-based positions

    1, b+1, 2b+1.                                         (4)

Thus the rearm always changes exactly three apex arcs, for arbitrarily
large carrier size, while fixing the root star. The common count and
states are

    h_b = 3B²+4B−1,
    before: (h_b, 3B−1, 3B−1, h_b),
    after:  (h_b, B²+2B, (B+1)², h_b+1).                  (5)

To obtain the count, the inherited odd-cycle formula gives, for positive
a,b,c,d and A=2^a,B=2^b,C=2^c,D=2^d,

    h(1^a0^b1^c0^d)
      =3+(A−1)(B−1)+(C−1)(D−1)+(A−1)(D−1)BC.             (6)

Indeed all directed odd cycles pass through the apex. Those meeting
only the first adjacent pair of runs contribute
((A−1)(B−1)+1)/2; the second adjacent pair contributes
((C−1)(D−1)+1)/2. Those from the first1-run to the final0-run contribute
(A−1)(D−1)BC/2, since the nonempty middle carrier supplies equal even
and odd subset counts. The inherited identity H=1+2·#odd cycles proves
(6). Setting a=1,c=b,d=2 proves h_b. Deleting the first letter leaves
leading zeros followed by1^b00, with count3B−1. Deleting the last
letter in(6) gives(B+1)². Equation(2) proves the states.

The first values are19,63,223,831,... at orders6,8,10,12,... . Both
states already have gcd(h,q)=1: gcd(h_b,3B−1) divides2, while
gcd(h_b,(B+1)²)=1 follows from h_b≡−2(B+2) modulo(B+1)².
Thus this repairs readiness, not an absent gcd property. The inherited
unit-ready NR operation can then produce the exact target3h_b+2; no
claim that its output is again unit-ready is made.

## 4. Literal63 repair

Rows list outgoing-neighbor bit masks, with labels0,...,7. Take

    T=(16,5,17,7,74,31,175,63), root r=6, apex a=4.

Deleting4 leaves the transitive order

    [6,7,5,3,1,2,0],

and its apex word is1001100. Reverse exactly the three pairs

    {2,4}, {3,4}, {4,7},

to obtain

    U=(16,5,1,23,198,31,175,47).

Its apex word on the SAME carrier order is1100110. The entire star at6
is unchanged, with sole predecessor4, and

    state(T,6)=(63,11,11,63),
    state(U,6)=(63,24,25,64).                              (7)

These are the two inherited converse isomorphism classes, now connected
by a concrete same-name operation preserving the required port.

## 5. The exact width≤2 obstruction at63

The finite universe is all subsets of one, two, or three unordered pairs
of each of the two literal eight-vertex cores. A move flips its selected
pairs simultaneously and is retained only when its final count is63.
This permits arbitrary pair positions; no root-star filter is imposed
while establishing the obstruction.

For T in§4 the ONLY neutral single-pair flips are

    p={0,2}, q={1,3}, z={5,7}.

The ONLY neutral two-pair flips are {p,q}, {p,z}, {q,z}. Each single
move just exchanges the names in a transitive two-vertex module. They
commute, so their closure is an eight-state cube. Each fixes root6 and
its rooted state. By equivariance, no sequence of neutral width≤2
moves from T can leave this cube or perform the rearm(7).

For the other representative

    V=(2,16,3,23,165,79,31,111),

the corresponding neutral single pairs are {0,2},{3,6},{5,7}, and again
the only neutral doubles are their pairwise combinations. Therefore,
using the inherited two-class completeness from THM-344, the entire
labelled eight-vertex H=63 fibre under neutral width≤2 moves splits into
80640/8=10080 eight-state cubes. In particular the two converse classes
cannot be connected by any such sequence. The three-pair move(7) is
therefore sharp, even allowing arbitrarily many narrower neutral steps.

The label count uses rigidity of both cores. Each score class has size
at most two; an automorphism fixes every singleton class, and cannot
exchange a two-vertex class because its directed internal edge would
reverse. Thus each class contributes 8! labelled tournaments, giving
2·8!=80640 in total.

For each representative the complete width-three probe finds seven
neutral triples: four in the same isomorphism class and three in the
opposite class. This last number is supplementary; the sharp rearm
requires only the two-width exclusion and explicit witness.

This is a fixed-order, fixed-count statement. It does not prohibit
passing through a different count or adding vertices. The small twin
swaps do not refute unrestricted reconfiguration, just as a disconnected
fixed-color graceful fibre does not refute global graceful connectivity.

## 6. The exact feedback-rank cost of the NR exit

Let τ(T) denote the minimum number of vertices whose deletion makes T
acyclic. For ANY nonempty rooted tournament(T,r), perform the ordinary
root-child birth, retain its new child y (N), then perform a birth at y
and retain y (R). Call the second new vertex z. Then

    τ(NR(T,r)) = τ(T−r)+1.                                (8)

The induced graph on r,y,z is the directed triangle r→y→z→r, disjoint
from the induced copy T−r. Every feedback set must therefore contain at
least one vertex of this triangle and a feedback set of T−r. Conversely,
delete r and a minimum feedback set of T−r. The remaining graph is the
ordinal join of y→z with the remaining acyclic old graph. This proves
both inequalities. Equivalently, the maximum transitive-subtournament
order obeys β(NR(T,r))=β(T−r)+2. These are properties of the actual
two births, independent of their Hamilton count.

For a unit-ready source in§2, w begins11 and ends0. Its deletion T−r
contains the triangle a→v1→v_(m−1)→a, and deleting a makes it
transitive. Thus τ(T−r)=1 and the NR output has EXACT feedback rank2.
In particular no choice of an apex in that output restores the one-apex
chart. This gives the precise reuse boundary: the same-H repair is
uniform across all eligible sources, but its accepted NR exit requires
a larger construction class for the next repair.

The rank-one wrapper principle in the supplied Session23 synthesis is
an older special-case inspiration. Equation(8) is the rooted-deletion
version for these literal births, with its entry object T−r explicit.

## 7. Reproduction and scope

Run `python3 probe.py` for the two literal one/two/three-pair slices,
with direct-permutation checks of the source H values and subset DP
for all flips. Run `python3 reflection_check.py` for all-word rooted
controls, the family, the mirror failure1010, and the actual63 NR target.
It also checks(8) on every labelled tournament of orders1 through4 at
every possible root.
The saved receipts state their finite universes. Neither experiment
re-enumerates all tournaments or substitutes for the all-word proof.

The operation preserves H and a physical port but may change other
endpoint counts and future contexts. It does not show that every odd
target other than7,21 has an entry source; that remains the Chan–Pak
coverage problem. The new gain is an exact repair that a width-two
neutral search provably cannot find at63.
