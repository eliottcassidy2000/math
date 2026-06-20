---
id: HYP-2677
title: LRC(14) packet-sign tournament atlas
status: OPEN; exact atlas, classifier not proof
source: codex-2026-06-20-S49
depends_on:
  - HYP-2676
  - HYP-2674
  - HYP-2675
  - THM-546
  - HYP-2639
  - HYP-2459
related:
  - HYP-2638
  - HYP-2648
  - HYP-2450
  - OPEN-Q-108
---

# HYP-2677 - LRC(14) Packet-Sign Tournament Atlas

## Claim Being Tested

The six one-missed-sector packet signs from HYP-2674/HYP-2676 are useful, but
only as chamber data with an address attached.  The exact atlas tests three
quotients:

```text
six signs as K4 edge orientations,
six missed sectors as tournament vertices with pair-pressure gauges,
three opposite K4 edge-pair balances retaining exact packet mass.
```

The guardrail comes from the coefficient-tiling bridge: unmarked sign
tournaments are too lossy.  The atlas therefore keeps exact packet magnitudes,
additive profile, squarefree profile, run count `V`, and the incoming THM-546
S2 rational Abel pressure

```text
|Delta_w| / ((6/49) V(E')/w).
```

## Computation

Script:

- `04-computation/lrc14_packet_tournament_atlas_codex_s49.py`
- output: `05-knowledge/results/lrc14_packet_tournament_atlas_codex_s49.out`

The script imports the exact packet telescope from HYP-2676, then adds:

- K4 edge-sign tournaments under four sector-to-edge gauges;
- all sector-to-K4-edge permutation fingerprint distributions;
- opposite K4 edge-pair sums such as `(01,23)`, `(02,13)`, `(03,12)`;
- six-sector vertex tournaments using declared pairwise gauges:
  `scalar_signed_rank`, `cyclic_pair_sum`, `opposite_compensated`, and
  `absolute_opposite`;
- score histograms, directed triangle counts, SCC sizes, Hamiltonian-path
  counts, and tie counts.

All displayed LRC quantities are exact Fractions.

## Exact Findings

### B14 near-speed leaders

The top twelve B14 near-speed `Delta_w` rows all have packet sign word
`++++++`.  As K4 edge signs, every one collapses to the same transitive
tournament:

```text
K4 lex scores=(3,2,1,0), directed_3cycles=0, scc=(1,1,1,1), HP=1.
```

The six-sector cyclic-pair tournament is also identical across the twelve:

```text
scores=(3,3,3,2,2,2), directed_3cycles=8, scc=(6,), HP=41.
```

Thus the K4 sign quotient correctly locates the same-sign pocket, but it does
not separate the dangerous finite leaders.  The missing address is magnitude
and additive model data.  For example the B14 top row is

```text
E'=(0,2,4,6,7,8,9,10), w=12
Delta_w=5347/30870
sign=++++++
V=58
AbelPressure=5347/18270
sumset excess=3
K2=9/4
```

### KPS third-pocket contrast

The KPS third-pocket row remains the clean contrast:

```text
E'=(0,3,5,16,28,30,33), w=35
Delta_w=1171/452760
sign=++-+--
V=84
AbelPressure=1171/133056
sumset excess=13
K2=26/7
```

As a K4 sign word, `++-+--` is not rigid.  Over all sector-to-edge
permutations, its exact K4 type distribution is:

```text
10x scores=(2,2,1,1), c3=2, scc=(4,), HP=5
 6x scores=(3,2,1,0), c3=0, scc=(1,1,1,1), HP=1
 2x scores=(2,2,2,0), c3=1, scc=(3,1), HP=3
 2x scores=(3,1,1,1), c3=1, scc=(3,1), HP=3
```

Under the lex K4 gauge it happens to look transitive, so the K4 sign quotient
alone can miss the cancellation.  The opposite-pair balances keep the address:

```text
lex opposite-pair balances:
  1+6 = 27/784
  2+5 = 1543/25872
  3+4 = -23/6468
```

The negative `3+4` pair is the useful "weighted positive / weighted negative"
signal: cancellation appears before absolute values.

The six-sector `cyclic_pair_sum` quotient sees the topology directly:

```text
KPS cyclic_pair_sum:
  scores=(4,3,3,3,2,0), c3=4, scc=(5,1), HP=15.

B14 leader cyclic_pair_sum baseline:
  scores=(3,3,3,2,2,2), c3=8, scc=(6,), HP=41.
```

So KPS is not merely "same K4 sign class with a few negatives."  In the cyclic
sector-pair gauge, negative pair mass flips enough arcs to break the fully
strong six-sector shape into a `5+1` SCC topology.

### HYP-2675 boundary versus true-wide

The HYP-2675 rows split before magnitudes:

```text
boundary row:
  E'=(0,2,4,6,8,10,12,14), w=15
  sign=+++-+-
  Delta_w=391/41160
  AbelPressure=391/18816

true-wide base:
  E'=(0,4,6,8,10,12,14,15), w=16
  sign=++++++
  Delta_w=943/5880
  AbelPressure=943/2970
```

The boundary row has lex K4 type `scores=(3,1,1,1), c3=1, scc=(3,1), HP=3`,
while the true-wide row is transitive in the K4 sign quotient.  This supports
HYP-2675's split: boundary collar and true-wide base are different proof
obligations, not one `span>14` branch.

## Hypothesis Triage

The atlas records:

- **PASS H1:** B14 near-speed leaders collapse to one transitive K4 sign type.
  This locates the `++++++` pocket but does not separate it.
- **PASS H2:** KPS `++-+--` has multiple possible K4 types under
  sector-edge relabelling.
- **PASS H3:** HYP-2675 boundary and true-wide rows split by sign before
  magnitude.
- **WARN H4:** the `opposite_compensated` six-sector quotient is not a KPS
  separator; KPS shares a type with some B14 leaders.  The better separator is
  `cyclic_pair_sum`.
- **PASS H5:** exact magnitudes inside the same `++++++` chamber create
  multiple six-sector opposite-compensated tournament types.

## Integration With Incoming THM-546 S2

During this session, incoming mac-mini S2 sharpened THM-546 to the rational
signed Abel bound:

```text
|Delta_w(E',w)| <= (6/49) V(E') / w.
```

Here `V(E')` is exactly the one-missed-sector run count used by the packet
atlas.  Therefore HYP-2677 is not trying to replace the analytic estimate.  Its
role is narrower:

```text
THM-546 S2 handles signed Abel size.
HYP-2677 classifies the tight ungapped/same-sign chambers where size alone
does not explain the finite leaders.
```

A later rebase also brought in THM-547, which closes the boundary-collar branch
`second-largest <= 14` modulo a feasible finite check.  In this atlas the
HYP-2675 boundary row is therefore a calibration point for sign/topology, not
the live proof frontier.  The remaining structural target is the true-wide
branch (`second-largest > 14`) together with the finite `++++++` packet models.

## Proof Route

1. Treat K4 edge signs as a finite chamber coordinate.  `++++++` is the same
   transitive K4 type for every B14 near-speed leader, so this quotient is a
   locator, not a proof.
2. Attach exact opposite-pair balances.  Negative or small opposite-pair sums
   preserve the weighted positive/negative cancellation that sign words alone
   erase.
3. Use the six-sector `cyclic_pair_sum` gauge as a topology-of-information
   diagnostic: negative pair mass flips arcs in a fixed mod-7 cyclic
   tournament and can expose cancellation as SCC/triangle changes.
4. Classify `++++++` rows by finite Ruzsa/Freiman model plus opposite-pair
   packet balances; prove the complement either has a cyclic-pair arc flip,
   small exact pair mass, or enough signed Abel cancellation under THM-546 S2.
5. Use THM-547 to demote the boundary collar to finite-check/calibration status;
   focus new structural work on true-wide rows and same-sign finite models.
6. Reattach HYP-2648 state-word and HYP-2639 relation-shell addresses before
   any final LRC step.

## Tournament Analysis

Alternate vertex sets considered: runners, packet signs, K4 edges, missed
sectors, opposite K4 edge pairs, additive profiles, squarefree profiles,
relation-shell labels, and proof obligations.

Primary quotient A:

- vertices: K4 vertices;
- edge observable: packet sign assigned to a K4 edge;
- switch/gauge: sector-to-edge assignment, with `+` using the gauge direction
  and `-` flipping it;
- tie path: none for nonzero packet signs.

Primary quotient B:

- vertices: missed sectors `1..6`;
- pairwise observable: exact packet-pair pressure, especially
  `cyclic_pair_sum(s,t)=cyclic_base(s,t)*(P_s+P_t)`;
- switch/gauge: the fixed mod-7 cyclic orientation on sectors;
- tie Hamiltonian path: lower sector wins exact zero ties.

The quotients preserve packet sign, exact packet-pair mass, and finite
tournament topology.  They destroy full wall endpoints, measured state words,
and relation-shell visibility, so they are not final proof objects.

## Honest Status

No LRC(14) proof is claimed.  HYP-2677 turns the user's tournament prompt into
an exact classifier and gives one useful next theorem: finite `++++++`
Ruzsa/Freiman packets plus opposite-pair balances should cover the B14 leaders,
while the high-growth complement should be discharged by cyclic-pair arc flips
or the signed Abel estimate.
