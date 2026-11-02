# A minimal carrier for square guarded actual complements

Status: **PROVED** for the specified operation language; **FINITE-EXACT**
square-sum instances. General strong-atom composition and unrestricted
re-marking closure remain **OPEN**.

## Sources and operations

Fix labeled vertices, an undirected guard graph G with m edges, and a bank
P of distinct directed spanning paths using G. A source is any full
tournament T on these vertices. Observe which named bank paths are legal;
their number is the guarded Hamilton response H_G when the bank is the
entire path set of G. For each legal path P, the command K_P is the
ACTUAL complement: retain its directed path edges and reverse every
other tournament edge. The next command may choose any currently legal
bank path. Thus legal re-marking within the bank is allowed.

Names and the common vertex frame are retained. We count unmarked-source
observer states with a named selected-path command, not rooted metagraph
vertices carrying a current path, and not states modulo vertex relabeling.

Restriction to G is a closed projection: legality uses only its edges,
and K_P on those edges depends only on those same edges. Each fibre has
2^(binom(n,2)-m) full tournaments. This is a literal relation map, unlike
identifying sources merely because their current counts happen to agree.

Encode guard orientations by x in F2^m. For a path P let E(P) be its
undirected-edge mask and let A_P be its legality cube. Then

    K_P(x)=x xor (1_G xor E(P)),   domain A_P.              (1)

Every A_P has 2^b points, b=m-n+1. K_P is an involution preserving A_P.
Two cubes intersect exactly when their directed arc requirements agree
on every common edge. The inherited
[path kernel](../joint_incidence_clocked_frames_20261102/PATH_KERNEL.md)
already counts these joint intersections. We now close them under (1).

## A family theorem when path compatibility is a matching

Assume the compatibility graph of the bank has only isolated vertices
and disjoint edges. Write i for its isolated paths. For each compatible
pair {P,Q}, put d=|E(Q) minus E(P)|=|E(P) minus E(Q)|>=1 and
c=2^(b-d). The smallest deterministic carrier preserving every named
legality guard and every future legal K command has

    i + sum_pairs [3+2*indicator(d>=2)] + indicator(dead sources exist) (2)

states. This theorem applies to every finite bank satisfying the stated
hypothesis; it is not confined to square sums.

Proof. Isolated P has its entire cube as one state with a K_P loop.
For a compatible pair let I=A_P intersect A_Q. Its c points have both
guards. K_P(I) and K_Q(I) are disjoint sets of c points with only the P
and Q guard, respectively. A command applied at the corresponding leaf
returns to I. The three types form a star:

    P returning leaf <-- K_P --> both <-- K_Q --> Q returning leaf.

The remaining P-only sources number 2^b-2c and are closed under K_P;
the analogous Q-only sources form another closed type. These two types
exist exactly when d>=2. No third bank path can appear, by the matching
hypothesis. All sources with no legal bank path form one dead state.

The resulting quotient is closed. It is minimal because different path
guard sets distinguish all types except a returning versus closed leaf
with the same name; one application of that named K command distinguishes
those two. Thus one refinement round beyond joint legality is sufficient,
and sometimes necessary. This proves both bounds in (2).

All four sets outside I are essential in the argument. A present response
H_G=1 alone cannot tell whether complementing its legal path exposes a
second path. This is the exact extra bit on each paired leaf, rather than
a claim that a global free bit solves general attachment memory.

## The square sum graph on 25 integers

For vertices 1..25 and edges i+j a square, independent path enumeration
gives m=32 and 20 directed Hamilton paths. Their compatibility graph is
4K2 disjoint union 12K1, with d=2 on each pair. Therefore b=8,c=64.
The minimum in (2) is

    12 + 4*5 + 1 = 33 states.                           (3)

There are 4,864 active guard orientations. Among the 32 active quotient
states, 12 have weight 64, eight have weight 128, and 12 have weight 256.
The dead state has weight 2^32-4864=4,294,962,432. All weights multiply
by 2^268 for full labeled tournaments. The resulting H_G distribution is

    4,294,962,432 + 4,608 z + 256 z^2.

The 24 initial joint-guard classes split into exactly 32 active classes
after one complement refinement, and no further split occurs. The exact
program [audit.py](audit.py) enumerates every active source, constructs
every legal literal target, and verifies every quotient fibre. The
theorem above separately explains and proves its minimality.

Other full square-path banks give 3 states at n=15,16,17 and 7 states at
n=23, including the dead state. These are declared finite orders, not a
claim about all square-sum graphs from 25 onward.

## Why the scope cannot silently grow

At n=15 the single extra guard edge is {1,3}, at distance four along the
selected square-sum path. Complete that path transitively, and then reverse
only this extra chord. These two full tournaments belong to the SAME
restricted carrier state, but their ordinary Hamilton values are 1 and 9.
Both remain legal sources for the named square-sum path. This is a direct
hostile against promoting (3) to a full tournament or full-metagraph
carrier. The source and reversed-chord calculation are reproduced exactly.

The all-size boundary is already proved in
[THM-4807, untraced-arcs and legal-path reconstruction](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md):
the FULL named Hamilton-path legality language separates every pair of
distinct labeled tournaments. Consequently a carrier preserving all those
guards already needs 2^binom(n,2) states before complement transitions are
considered. The present compression is possible precisely because the
arithmetic path bank is smaller. Unnamed observations and selected-frame
transport are separate quotient questions.

The twelve-vertex 81/27 hostile remains inherited from
[DEFECT_TRANSFER.md](../overnight_rank_return_20261029/DEFECT_TRANSFER.md):
equal two-run mass 9,011,898 with ranks 24/22 modulo three, and equal
three-run mass 80,614,764 with ranks 430/435. Those tensors observe a
different, richer language. Neither they nor actual target relations are
identified by the present square-guard quotient.

## Q8 transport and the test for a legitimate quotient

The inherited
[four-marker Q8 remainder](../joint_incidence_clocked_frames_20261102/TRANSPORT_AND_CLOCKS.md)
belongs to the M12 hole-move model: four named tracked tokens leave eight
possible source frames; five tracked tokens determine the frame. That
result supplies no action of Q8 on the square-guard sources.

For any proposed group quotient, the operational test is precise: each
guard and output must be constant on quotient fibres, and every enabled
command must map an entire fibre into one target fibre. A Q8 orbit can
be discarded only after those conditions are established for the actual
operations. Arbitrary re-marking can fail them even when a fixed set of
markers succeeds.

There is also a useful algebraic boundary here. If all formulas (1) are
extended to TOTAL maps on F2^m, they are commuting translations of order
two. Every group generated by them is elementary abelian, so no faithful
Q8 transport can be obtained from these fixed-frame translations alone.
This statement concerns the total extensions, not commutativity of partial
commands: their legality domains must still be checked. Source-dependent
frame permutations could create a larger affine action, but would require
an explicit new map and guard audit.

Thus the 33-state reduction comes from joint incidence plus exact
transition refinement. It neither assumes nor needs quaternion symmetry.
