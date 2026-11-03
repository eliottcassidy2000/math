# Concurrent seam recurrence: exact transport and the universal extension

The independent [golden seam packet](../golden_seam_recurrence_20261103/README.md),
published at59fc9abf05 while this checkpoint was being prepared, derives the
same golden recurrence in different coordinates. It also supplies a larger
finite alphabet. The overlap is explicitly credited; neither packet's
chronology is used as a priority claim.

## I1. Exact coordinate and word-order transport

Our processing word is the reverse of its source-order word. Its singleton
a is our s, its triangle b is our c, and the Q4 masks agree. Its coordinates
(h,x,y) are related to ours by

    (h,x,y)=(A+B,(C-B)/2,(B+C)/2),
    (A,B,C)=(h+x-y,y-x,x+y).

The halves are integral on actual words, since successive singleton H
responses are odd. With row vectors and

    L=[[1,0,0],[1,-1/2,1/2],[0,1/2,1/2]],

the transfer matrices obey M_l L=L P_l. `integrate.py` independently
reconstructs M using actual path/endpoint/cover DP and (G4), then checks
this equality for **all960 supplied charts**:1 singleton,1 triangle,
5 strong4,50 strong5 and903 strong6 charts. It does not copy the incoming
six-state forest compiler to derive M. Their full source rows are checked
to be strong and to contain the declared Hamilton spine.

The golden reachability totals and all ten strong5 misses agree exactly.
Our bound10001 includes one more attained value than its bound10000.
Both packets independently obtain the same order12 aggregate rational
series. Our source-aware55->111 exhaustion is an additional control.

## I2. The finite six-vertex extension is already complete

The incoming packet supplies eight six-vertex charts repairing all eight
extra strong5 misses. Its fixed65-symbol bank then passes a holdout through
30000, with only7 and21 missing. The full30000 saturation is **INHERITED
FINITE-EXACT**, not independently rerun here. Its actual eight repair words
are independently checked in this packet by literal graph DP and strongness:

| H | Source-order word in incoming notation |
|---:|---|
|47|a s14|
|83|a s268|
|103|a s330|
|491|a a s248 a|
|709|a s261 a a a|
|833|a s85 b|
|3227|a r35 s262|
|4619|a s453 r29|

Here rj,sj are the actual strong5/strong6 chart masks, not integer scalar
placeholders. `integration.json` preserves reversed processing words,
literal target rows and their three-coordinate states. It also pins the
incoming evidence by LF-normalized SHA256.

## I3. What the all-size proof now resolves

Incoming THEORY T4 explicitly leaves arbitrary block-size three-coordinate
closure unproved. Our general (G1) run bijection, (G3) reconstruction and
the universal forest intertwining in RECURRENCE G4 resolve that scoped
question for every finite supplied block, with nonnegative integer matrices
in (A,B,C) coordinates. The reason is structural: a single exceptional
cross-edge permits at most the runs T1 B1 T2 B2. It does not depend on a
finite atlas, bounded block size, or polynomial-time compilation.

This does not prove that any finite chosen alphabet reaches every allowed
odd value. A uniform block formula and a uniform target inverse are
different obligations;37's absent golden state is a permanent counterexample
to inferring the latter from the former.

## I4. Retain which end may grow

Incoming THEORY T6 proves minimum rational linear dimension9 for both-end
word attachments, retaining the full3-by-3 word matrix. Its three-versus-nine
hostile independently transports as

    v(s c)=v(q5)=(4,1,9),
    H(s s c)=9,   H(s q5)=13.

The equal states tolerate every future extension on the processing right;
prepending on the processing left separates them. This is a necessary
qualification of our coarsest-observer theorem, not a contradiction. The
incoming nine-dimensional minimum is inherited; this packet replays its
small actual-source hostile rather than its full9-by-9 rank certificate.

The next proof target is a guarded inverse for the fixed65-symbol bank, or
a counterexample beyond its holdout. A source-dependent larger block may
be compiled by the universal formula, but its actual existence and choice
must still be proved. The incoming fixed-wrapper obstruction also rules
out treating fixed word concatenation as uniform3H+2 or3H+4 gadgets on the
entire golden family.
