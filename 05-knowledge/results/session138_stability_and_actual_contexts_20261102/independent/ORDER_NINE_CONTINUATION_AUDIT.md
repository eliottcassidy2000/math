# Order-nine continuation: independent theorem audit

**Verdict: ACCEPT — PROVED, with a finite exact graph certificate.** This is
an addition for checkpoint four; it does not alter earlier frozen packets.
For every weighted tournament of order nine,

    Delta >= (175/256)V + (105/256)tau + C4/16 + (5/8)C5 + C6/4 + C7 + C8.

Consequently the sharp variance coefficient is **175/256**. For V>0,
equality in the unreserved inequality Delta >= 175V/256 holds exactly for
an ordinal sum of four fair two-vertex blocks and one singleton, in any
block order. This does not assert equality in the stronger reserved bound.
Together with earlier results, the still-unsettled exact-coefficient orders
are 10 through 1499. No order-ten computation was performed.

## Analytic audit and dependency boundary

The reviewed final manuscript is `stability/SHARP_NINE_GRAPH_CERTIFICATE.md`;
its precise bytes and those of its supporting files are pinned in
`ORDER_NINE_REVIEWED_SOURCES.json`. The earlier candidate filename is not
part of this packet. Status/name promotion made no mathematical change.

The auditor independently reconstructed all 30 cycle-collection types,
including the empty type, from partitions with parts at least two and total
at most nine. Determinant signs are (-1)^(sum(j-1)); the signed full
cover v=(-1)^9 det(A) has coefficient (-1)^(number of cycles). The inherited
odd-cycle expansion gives coefficient 2^(number of cycles) when every
cycle is odd. These rules reproduce both full signed-gap ledgers, including
3+3+3 and 2+2+2+3. No choice of a presumed determinant sign is made.

For each selected four-cycle, the five-vertex complementary margins are
u_R+v_R and u_R-v_R. Adding C4_R/2 gives exactly the relevant signed
four-cycle groups: a pair of four-cycles is counted twice, while each other
collection has a unique selected four-cycle. Thus the actual coefficient
of M44 is -1, not -2. The inherited order-five bounds u_R+v_R >= 1/16 and
u_R-v_R >= 5/16 give a common C4/16 reserve. All these remainders have order
five, so no empty-remainder or low-order exception is being used.

The other negative couplings satisfy C6-M26 >= C6/4 and
C5-2M225 >= 5C5/8. For a six-vertex complement, weighted degrees give

    W_R = (sum_v d_v^2)/2 - sum_e q_e^2 >= V_R^2/3 - V_R/4,

using sum_v d_v=2V_R and 0<=q_e<=1/4. Together with the previously proved
local triangle inequality tau_R>=W_R, this implies

    M22_R-(3/2)tau_R <= 75/256-(V_R-15/16)^2/3 <= 75/256.

Summation over the selected triangle counts M33 twice and M223 once.
It therefore gives tau+3M33-M223 >= 181tau/256. The remaining unused
products in each signed ledger are nonnegative. The independent symbolic
checker records their exact lists rather than assuming that signs cancel.

Write P=V-M22+M222-M2222. The resulting reduction is

    Delta >= P+(181/256)tau+M23+C4/16+(5/8)C5+C6/4+C7+C8.

The triangle inequality also gives M23>=P21(q), where P21 counts a wedge
and a vertex-disjoint edge, without an induced-subgraph condition. The
separately affine polynomial

    P+(19/64)W(q)+P21(q)-(175/256)V

has endpoint residual

    1024F = 81m-64m2+16m3-4m4+19w+16p21.

The graph audit below proves this residual nonnegative. Since tau>=W,
the surviving triangle reserve is 181/256-19/64=105/256. The coefficient
19/64 is sharp for this graph reduction, witnessed by P3 disjoint-union
three edges; no optimality of the actual triangle reserve is asserted.

## Independent finite certificate and implementation controls

`n9_audit/AUDIT.md` contains the complete coverage proof and resource record.
The independent implementation does not import producer representatives or
pattern-counting code. It uses adjacent-transposition BFS for seven-vertex
orbits, equitable color refinement followed by exhaustive within-cell
permutations for exact eight-vertex canonicalization, and vertex-addition
matching/wedge recurrences rather than the producer's fixed pattern masks.
Refinement alone is not treated as a complete isomorphism test.

The exact universe is:

- all 2,097,152 labelled seven-vertex graphs, partitioned into 1,044 orbits;
- 133,632 eight-vertex augmentations, canonically reduced to 12,346 classes;
- all 256 ninth-vertex neighborhoods for each class: 3,160,576 cases.

Minimum residual is zero, maximum 60,372, with 14 zero augmentations.
A separate full vertex-subset DP and wedge-first count checks 344 snapshots,
including every zero. Every augmentation checks m2+w=binom(m,2).
These are full-universe coverage checks, not a random tournament sample.
Signed count intermediates fit safely in 32 bits; the unsigned digest wrap
is intentional. Checks are explicit failures, not removable assertions.

Normal, optimized, and completed UBSan runs have byte-identical receipts.
The initial two UBSan compilation attempts failed under a 24 MiB data cap;
those attempts remain recorded as failures. Increasing only the compilation
budget to 64 MiB allowed a terminal successful build. Measured compiler
peak was 50,468 KiB; successful certificate runtime peaks were 6,132,
6,144, and 6,164 KiB. There were no sanitizer diagnostics; stderr contains
only the expected representative-count line. The 40 MiB runtime target
was met. No claim is made that the failed attempts passed.

`n9_audit/cycle_ledger.py` independently reconstructs the partition ledgers,
selected-cycle multiplicities, nonnegative discarded products, rational
square completion, endpoint scaling, and sharp fair-pair value. Normal and
`python -O` outputs agree. This small symbolic control supplements, rather
than replaces, the analytic inequalities and full graph enumeration.

## Equality and firewall

The zero endpoint graphs have three isomorphism types: empty,
4K2 disjoint-union K1, and P3 disjoint-union 3K2. The last two can differ
by one edge. Hence the continuous multiaffine residual has genuine zero
segments; it is false to infer endpoint-only zeros from this certificate.

Instead, equality in the unreserved tournament bound forces tau=0 by the
strictly positive 105/256 reserve. Every supported binary completion is then
triangle-free, hence transitive; multilinearity gives H=1. The previously
audited H=1 classification and scalar minimization force four fair ordinal
pairs and one singleton. Conversely those blocks have H=1, u=81/256,
v=0, V=1, and Delta=175/256.

The earlier Pfaffian termwise-square counterexample is a separate packet
and is neither used nor duplicated here. This proof uses only valid total
margin bounds, never the refuted unit-coefficient termwise assertion.
