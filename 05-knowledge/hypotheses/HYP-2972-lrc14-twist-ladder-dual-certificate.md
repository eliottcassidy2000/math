---
id: HYP-2972
title: LRC14 twist-ladder dual certificate
status: PROOF-INTERFACE / finite-bank global twist witness certificate; not a proof of LRC14
source: codex-2026-06-24-S155
script: 04-computation/lrc14_twist_ladder_dual_certificate_codex_s155.py
result: 05-knowledge/results/lrc14_twist_ladder_dual_certificate_codex_s155.out
related:
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2956
  - HYP-2953
  - HYP-2934
  - HYP-2932
  - HYP-2928
  - HYP-2901
  - HYP-2908
  - THM-566
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2972: LRC14 Twist-Ladder Dual Certificate

## Claim

The HYP-2963 labelled-packet bank has a surprisingly small **global rational
twist certificate**:

```text
for every audited row S,
there exists reduced a/q with q <= 42 such that
  ||v a/q|| >= 1/14 for every v in S.
```

This is a different proof route from HYP-2965/HYP-2968.  It does not start
from safe components, boundary endpoints, or lift intervals.  It starts from a
finite set of rational phases and treats a failed ladder as a dual blocker
hypergraph:

```text
vertices = rational twists a/q
hyperedge for speed v = twists blocked by v
```

The theorem target is dynamic, not fixed finite.  HYP-2901 already shows that
no fixed finite denominator ladder can prove all of LRC14, because committed
lcm speeds kill every denominator below the committed wall.  HYP-2972 says the
current reduced packet bank is small enough for `q<=42`, and the only
`q<=27` obstruction is visibly a resource packet.

## Computation

Script:

```text
04-computation/lrc14_twist_ladder_dual_certificate_codex_s155.py
```

Stored output:

```text
05-knowledge/results/lrc14_twist_ladder_dual_certificate_codex_s155.out
```

Default bank:

```text
single_limit=180
two_swap_limit=36
alias_depth=4
lcm_tail_max=5
qmax=64
q27_m=27
audited rows=21913
```

## Ladder Census

The main finite-bank result:

```text
ladder                  twists  certified  missed  max first q
q<=27                      229      21908       5            27
q<=42                      541      21913       0            41
q<=64                     1259      21913       0            41
qdiv+Q27(m<=27)           3219      21913       0            77
q<=64+Q27(m<=27)          3893      21913       0            41
```

The full ladder's first-witness denominators begin:

```text
14:6050
8:2306
7:2195
9:2172
10:2075
11:2059
13:1863
12:1682
19:349
17:325
16:132
15:129
```

For `qdiv>14` rows, the first-witness denominators are:

```text
19:349
17:325
16:132
15:129
21:92
22:57
23:49
18:23
20:18
41:5
25:4
27:3
24:1
```

## The q<=27 Dual Obstruction

The smaller `q<=27` ladder misses exactly the divisor-loaded lcm-tail packet
family in this bank:

```text
{1,2,3,4,5,6,7,8,9,10,11,13,84m},  m=1..5.
```

Every one of these rows has the same first full-ladder rescue:

```text
t = 17/41.
```

The first slack numerator is `1`, with tight speed `5`.  Thus the rescue is
barely strict at the integer inequality layer:

```text
14 * dist_num(v * 17 mod 41) - 41 >= 0,
minimum = 1.
```

This matters because `41 = 3*14 - 1` is already the Farey/K33 denominator from
the `3/41` near-miss lane.  The lcm-tail covering packet and the K33-unit
excess packet are seeing the same first non-apex denominator.

The dual blocker profile for `q<=27` is not a single magic speed.  The
committed tail speed has the largest load, but many failed twists have private
blockers among the AP-core speeds.  This is useful: the blocker hypergraph is
a resource ledger, not a one-coordinate scalar.

## Proof Program

The proposed route is:

```text
Moon-core / packet reductions
  -> attempt a bounded twist ladder
  -> if a twist survives, it is an exact LRC14 certificate
  -> if none survives, read the speed-vs-twist blocker hypergraph
  -> private blocker or high-load resource gives descent, state-lift, or next rung
```

This turns the fixed finite-denominator obstruction from HYP-2901 into a
recursion rule:

```text
committed wall at X
  -> no fixed q<=X proof
  -> blocker resource must explain the next rung
```

In the present bank, the explanation is concrete:

```text
q<=27 obstruction = lcm-tail packet
next rung         = 41
uniform witness   = 17/41
```

## Relationship To Earlier Routes

This is not a replacement for HYP-2965 or HYP-2968.  It is a different
coordinate system.

HYP-2965 asks:

```text
Which exact boundary bridge stays open?
```

HYP-2968 asks:

```text
Which of fourteen lift packets has positive interval mass?
```

HYP-2972 asks:

```text
Which rational twist already witnesses the row, and if none in a ladder does,
which speeds form the finite blocker cover?
```

The quotient preserves the threshold witness predicate for selected twists.
It destroys open intervals between twists and fine endpoint ownership.  That is
the cost, but also why it is meaningfully different from the previous angle.

## Tournament Analysis

Vertices are proof coordinates, not runners:

```text
q_threshold_gate
twist_ladder_primal
q41_rescue_rung
q27_fiber_ladder
blocker_hypergraph_dual
committed_denominator_wall
source_spectrum_packet
exact_haar_interval
boundary_gap_packet
raw_apex_residue_tournament
raw_runner_set
```

Pair observable:

```text
which coordinate preserves an exact LRC witness or a finite dual blocker
certificate while surviving denominator-wall adaptation.
```

Switch/gauge:

```text
lexicographic retention of primal twist witness, blocker dual,
denominator-wall adaptability, packet compatibility, and anti-scalarization.
```

The default tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
c3=0
scc=(1,1,1,1,1,1,1,1,1,1,1)
hp=1
```

Hamiltonian path:

```text
source_spectrum_packet
> twist_ladder_primal
> blocker_hypergraph_dual
> q41_rescue_rung
> q27_fiber_ladder
> exact_haar_interval
> boundary_gap_packet
> committed_denominator_wall
> q_threshold_gate
> raw_apex_residue_tournament
> raw_runner_set
```

## Status

Confirmed as a finite-bank certificate for the HYP-2963 bank.  Not a full proof
of LRC14.  The next theorem is the dynamic ladder theorem:

```text
after standard LRC14 reductions, every primitive packet either has a bounded
twist witness, or its finite blocker hypergraph yields a descent/state-lift/
next-rung certificate compatible with HYP-2901.
```
