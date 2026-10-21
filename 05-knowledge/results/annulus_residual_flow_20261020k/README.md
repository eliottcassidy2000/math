# Source-owned temporal flow: 2,764 residual intervals become an optimal 94

2026-10-20. **PROVED:** the source-owned horizon-capacity criterion,
handoff certificates, weak-dual obstruction, and nested residual inequality.
**FINITE-EXACT:** optimal whole-interval completion at N4096 and N8192.
**OPEN:** a uniform153N bound and universal Collatz entry. Optimality here
is only inside the explicitly frozen sibling-interaction bank.

Files: [checker](check.py), [complete receipt](receipt.json).

## 1. Inheritance and the changed operation

The immediate source is [the preceding block packet](../annulus_block_completion_20261020j/README.md).
It retained21,152 native three-pair blocks at N4096 and completed6,088 of
8,852 positive interaction intervals by putting each inside a single
negative interval, without overlapping uses of a container. It left
2,764 positive intervals. Its sharp timing-blind norm27 theorem and the
112-source, budget336/full-bill156 control are unchanged.

The present operation permits a positive interval to move between
negative owners at an actual horizon. For example, negative intervals
[0,1),[1,2) jointly pay positive[0,2), although neither contains it.
The time coordinate is essential: two negative copies each of[0,1)
and[2,3) have total area4, but cannot pay[0,3) at horizon1.

Source objects and conventions are unchanged: primitive positive sources
s in[4N,256N), enlarged charges A_(64N)(s), fixed627-point frontier V,
actual first-arrival terminal and depth, and the full canonical-P braid
clock. These are not literal fully odd U-word clocks. The signed charge
is v_s=(-1)^clock A_(64N)(s), with Q(a,b)=a²-ab+b² and
beta(v,w)=Q(v+w)-Q(v)-Q(w). Compatible mod15/mod35 registers jointly
encode mod105. Sources stopped in V are excluded before forming sibling
runs, rather than assigned an invented later arrival.

The exact original ledger remains

    Delta E=A_d+I_d,   A_d+M_d=153N,
    Delta E<=153N iff I_d<=M_d.                        (1)

The mechanism is inherited from [native sibling guards](../sibling_annulus_guards_20261020i/README.md)
and [source-owned first-arrival accounting](../session18_native_annulus/README.md).
The relevant correction is the six-source overlapping-exposure hostile:
separate safe blocks cannot discard their mutual interactions. Here
every original pair address remains owned, and every such interaction
outside the selected sibling bank stays in R.

Board: original pair address / half-open clock interval / reusable capacity
at different times / fixed-source refinement / dual horizon obstruction /
unmatched physical interaction. The anchor is a better residual certificate;
the niche is a finite flow with an exact dual; the boundary case is a
temporal gap that total capacity cannot fill.

## 2. Exact capacity and source-preserving handoffs

Each nonzero native three-pair block has selected debit either
-3*1_J(d) or +3*1_I(d), where I,J are half-open integer depth intervals.
The zero blocks retain their identities but supply no capacity. Let P
be the positive block bank and G the negative bank. Define

    p(d)=#{I in P:d in I},    b(d)=#{J in G:d in J}.

**Horizon-capacity theorem.** A selected set S of whole positive blocks
is jointly paid by the full negative bank at every horizon if and only if

    #{I in S:d in I} <= b(d) for every d.              (2)

This is the exact signed sum -3b(d)+3 active_S(d)<=0. It is not an estimate
of the other pair interactions. It neither spends future cancellation
nor withdraws any original source norm from153N.

This is a global scalar payment, so the two matched owners may belong to
different actual terminals. No charge vector is moved between terminals,
and no new physical collision is asserted. Every selected term keeps its
original source pair, terminal and clock. All other pair terms, including
old/new and cross-block interactions, stay in the unchanged remainder R.

There is an equivalent source-indexed handoff certificate: partition each
selected I into half-open pieces and assign each piece to a negative J
containing it, with no two pieces assigned to one J overlapping in time.
Such pieces imply(2). Conversely at each integer cell[d,d+1), match the
active positive owners to distinct active negative owners using(2),
then concatenate consecutive cells with the same two owners. This gives
a finite certificate. A persistent matching makes it shorter by keeping
an assignment until either owner ends.

This refinement does not create new physical arrival events at a handoff.
It partitions the indicator of an already authenticated interaction.
The artificial endpoint terms cancel on adjacent pieces. One negative
pair address can serve different positive owners at disjoint horizons;
it cannot supply two units at the same horizon. The entire new allocation
replaces the old one-container allocation. They must not be stacked as
independent uses of the same negative bank.

A useful network representation makes the ownership explicit. Use one
node for each horizon endpoint. A selected positive interval[a,b) is a
forward arc a->b of capacity1; the backward edge d+1->d has capacity b(d).
Sending one unit around the corresponding cycle consumes one unit at
each covered horizon. Hence integer circulations and feasible whole-job
subsets correspond exactly, with cost-1 on every job arc. The implementation
does not trust a solver or a fractional relaxation: it supplies and
independently checks actual integer jobs, source-owned segments and a
matching obstruction certificate.

## 3. A small exact obstruction certifies the optimum

For any finite set D of horizons, let Hit(D) be the positive intervals
meeting at least one d in D. Every feasible selected set satisfies

    number of omitted intervals
       >= |Hit(D)| - sum_{d in D} b(d).               (3)

Proof: assign each retained interval in Hit(D) to its first horizon in D.
At a given horizon at most b(d) retained intervals can be assigned. An
interval meeting several horizons is counted once in Hit(D), not several
times. Empty D gives the bound0.

Thus equality between(3) and the number actually omitted proves optimality
without assuming that the discovery algorithm is complete. The producer
sweeps horizons, processes ends before starts, and discards a latest-ending
mutable job at any overload. The verifier only needs feasibility plus(3).
Previously completed jobs can be frozen by subtracting their load from b.

The best dual in(3) is also computed exactly by a separate dynamic program.
If h<d is the preceding selected horizon, the newly hit intervals at d
are exactly those covering d with start>h. Indeed an interval covering d
and any earlier selected horizon must cover the latest one h. Therefore

    F(d)=max_{h<d}[F(h)+#{I:h<start(I)<=d<end(I)}-b(d)], (4)

with the alternative of starting at d and the empty solution0. Endpoints
are literal integer horizons, so no boundary interval is silently included.
This proves the dual computation; matching primal and dual values prove
the reported finite optima. No general greedy-optimality assumption is used.

## 4. The N4096 refinement is nested and globally optimal

The source compiler reconstructs the same21,152 blocks and126,912 physical
source rows and matches the previous bank hash. There are8,852 positive,
8,906 negative, and3,394 zero blocks. Freeze all6,088 formerly completed
positive block IDs. The residual problem has2,764 jobs.

The new allocation completes2,670 more, for8,758 completed intervals and
94 remaining. It is certified by the same two horizons in both the
nested and unrestricted problems:

| Source problem | Dual horizons | Distinct hit jobs | Available capacity | Forced omissions |
|---|---|---:|---:|---:|
| All8,852 jobs | 9,115 | 3,782 | 3,688+0 | 94 |
| Residual2,764; earlier6,088 fixed | 9,115 | 1,856 | 1,762+0 | 94 |

At depth9,3,779 original positive intervals meet capacity3,688. At115,
four meet capacity0; one of these already met depth9. Thus the union has
3,782, not3,783, members. The obstruction adds separated temporal demands
while retaining their shared owner. A maximum single-cut shortage91 alone
would not prove the necessary94 whole-job omissions.

The receipt includes all94 omitted block IDs and their564 physical source
rows, with actual terminals, clocks, mod15/mod35 and charges. It also
contains all16,629 ownership segments. One positive block, ID20119 on
[21,107), uses13 different negative owners. Its complete source rows and
the thirteen seam intervals are retained. This is a genuine improvement
over one-container completion, not a change of the positive source set.

For compactness the complete126,912-row bank is specified by its original
N, canonical run construction and source-row hash. Reproduction recomputes
every actual route and charge, independently checks the charge coordinates,
then authenticates every block with the inherited native auditor. The
flow addresses are therefore tied to that source bank, not transferable
to an abstract equal-size interval list.

## 5. Stronger residual guards; whole jobs and marked pieces differ

Define R_d to be the sum of every original old/new and unordered new/new
interaction not among the selected sibling pair addresses. It includes
all interactions between different blocks. It is unchanged from the prior
packet. The original sibling debit is exactly3[p(d)-b(d)]. For the94
unmatched whole jobs, let U_new(d) count those active at d. Then

    I_d <= R_d+3U_new(d) <= R_d+3U_old(d),              (5)

where U_old was the previous2,764-job residual. The second inequality
holds pointwise because the old6,088 selected IDs remain selected. The
new upper bound improves by as much as7,401, exactly at depth18. Its
remaining positive debit peaks at282, at depth17. Optimal omitted-job
cardinality does not claim optimal peak debit or a unique optimizer.

If partial positive intervals are allowed, match min(p(d),b(d)) actual
positive and negative owners independently at each horizon. This leaves
marked, source-owned pieces, not completed whole intervals. The sharper
pointwise statement is

    I_d <= R_d+3 max(p(d)-b(d),0) <= R_d+3U_new(d).     (6)

The left residual peaks at273 at depth9 in this bank. Every cut's actual
canonical owner assignment and leftover IDs are recomputed and hashed.
No such partial match removes its unmatched pieces later or certifies the
entire source interval. Consequently(6) does not contradict the94 minimum
for whole-job completion.

All257 inherited physical horizons are checked at N4096, including those
after the last sibling interval. R is read from the unchanged full-source
ledger; its original I is never changed. The maximal excess of(5) over
M is-201674, reached after all arrivals, agreeing with the inherited
complete control. The original N7 unit-rate hostile remains a failure;
this packet changes no sources or interaction values there.

For every N and every supplied native bank, either

    R_d+3U_new(d)<=M_d,
    or the stronger R_d+3max(p(d)-b(d),0)<=M_d,         (7)

at all horizons suffices for the153N increment bound. This is a quantified
certificate implication. Neither a uniform estimate on R nor sufficiently
small native deficiencies for every N is proved.

## 6. One higher-scale control and the remaining obstruction

A single additional scale N8192 is checked, without a broader census.
The bank has42,321 triples and253,926 physical source rows, with17,142
positive and16,998 negative intervals. The exact optimum retains16,983
and omits159. Its dual uses horizons16 and126:

    |Hit({16,126})|=9385,    b(16)+b(126)=9225+1,
    9385-9226=159.

The receipt contains32,251 ownership segments, and a positive interval
requiring11 successive owners in this deterministic allocation. Its
marked-segment debit peaks at474 at depth16. Thus94 is not a universal
scale-independent residual count. N8192 is a native sibling-bank control;
no new full annulus R/I census or full153N conclusion is claimed there.

The precise boundary is now smaller: at N4096 the remaining94 whole
intervals cannot all be paid by these negative sibling terms, even with
arbitrary lawful handoffs. That does not imply that their original
sources, or the complete annulus, are unpaid. Additional owned interaction
families, moving norm or a bound on R can still pay them. Increasing the
flow solver's search budget cannot remove this certified bank-relative
obstruction; the available source-owned capacity or the guarded target
must change.

## 7. Reproduction and validation

The independent finite universe is every subset of the six nonempty
integer intervals on[0,3], for every capacity vector in{0,1,2}³:1,728
problems. Every primal subset and every dual horizon set is independently
enumerated. Controls include duplicate intervals with different owners,
handoff across a seam, an uncovered horizon despite excess total area,
wrong types, duplicated ownership segments and out-of-owner claims.
The native verifier keeps full source/clock guards and authenticates both
primal allocation and dual obstruction. No external solver is needed.

```text
python -B 05-knowledge/results/annulus_residual_flow_20261020k/check.py --verify
python -B -O 05-knowledge/results/annulus_residual_flow_20261020k/check.py --verify
```

Generated receipt:882,379 checks. Verification adds one saved-receipt
comparison. LF SHA256:
`35974096cbf893fb209ce5c47d915a9837b4a048ac7579b1d0806b46e106ed4d`.

Independent peer proof/API audit: **PASS**, including the exact dual,
whole-job/marked-segment distinction and unchanged-R scope. The peer
independently verified the repaired horizon guards under optimized Python
with nine bool/float hostiles and the valid handoff control. Normal and
optimized complete runs match the saved receipt; the parent also replayed
the complete pre-hardening receipt (the final delta adds four type checks).
