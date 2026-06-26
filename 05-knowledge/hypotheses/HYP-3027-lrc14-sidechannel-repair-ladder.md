---
id: HYP-3027
title: LRC14 side-channel repair ladder for automatic quotient failures
status: EVIDENCE / full-bank quotient-repair audit and zipper-theorem target; not a proof
source: codex-2026-06-26-S190
script: 04-computation/lrc14_sidechannel_repair_ladder_codex_s190.py
result: 05-knowledge/results/lrc14_sidechannel_repair_ladder_codex_s190.out
related:
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3022
  - HYP-3021
  - HYP-3020
  - HYP-3018
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-2963
  - HYP-2997
  - HYP-2995
  - HYP-2992
  - THM-572
  - OPEN-Q-108
---

# HYP-3027: LRC14 Side-Channel Repair Ladder

## Claim

The HYP-3016/HYP-3017 automatic quotient failure, HYP-3018 normal-fan support
repair, HYP-3020 discrepancy-height trident, HYP-3021 decoy-generator
classification, HYP-3022 barcode/normal-fan refinement, HYP-3023 automatic
fiber zipper, HYP-3024 convergence audit, and HYP-3025 arc-Cech topology can be
sharpened into a repair-ladder theorem target:

```text
inside each automatic-word fiber,
the first nonzero repair cochain among exact M/q, boundary topology,
C27/K33 packet labels, tail magnitude, and zipper/sidecar coordinates
must open a strict component, descend to a known family,
be dual-annihilated, or emit F7/THM-572 residual debt.
```

The full HYP-2963 bank supports this as an exact finite guardrail.  Automatic
words are not safe quotients; exact `M` repairs open/boundary status but still
does not repair theorem-route purity; packet labels or the maximal guarded
non-route packet signature repair route purity in the audited bank.

## Computation

The S190 script imports the regenerated HYP-2963 classifier from HYP-3017 and
places the HYP-3023 magnitude-cocycle splitter beside the repair joins on the
same bank:

```text
packets=21913
Q-WITNESS                 14676
BOUNDARY-AP-GW                2
BOUNDARY-PETAL-SPORADIC       4
K33-STATE-LIFT                3
COVERING-MOMENT            7228
```

It then compares quotient joins:

```text
automatic_word:                       143 mixed-route fibers, 1 mixed-status fiber
word_plus_q_threshold:                  5 mixed-route fibers, 1 mixed-status fiber
word_plus_q_factor_power:             366 mixed-route fibers, 0 mixed-status fibers
word_plus_tail_exact:                2175 mixed-route fibers, 2 mixed-status fibers
word_plus_M:                          366 mixed-route fibers, 0 mixed-status fibers
word_plus_M_q:                          1 mixed-route fiber, 0 mixed-status fibers
word_plus_boundary_topology:            1 mixed-route fiber, 0 mixed-status fibers
word_plus_packet_labels:                0 mixed-route fibers, 0 mixed-status fibers
word_plus_guarded_packet_no_route:      0 mixed-route fibers, 0 mixed-status fibers
```

The result separates two notions that had been blurred:

- `M` is enough to separate strict-open from boundary in this bank.
- `M` is not enough to classify theorem route; it leaves `366` mixed-route
  fibers.
- `M+q_threshold` almost repairs route, leaving one two-row mixed fiber:
  `two drop(10,13)->add(20,26)` versus
  `two drop(8,12)->add(16,24)`.
- Safe-component topology similarly almost repairs route, leaving one
  K33/covering mixed pair:
  `two drop(12,13)->add(26,36)` versus `single swap 12->72`.
- C27/K33/transfer packet labels, or the full guarded non-route signature,
  are route-pure on the full audited bank.

## AP/Goddyn-Wong Word Repair

The AP/GW word is still the canonical stress fiber:

```text
MFCMMCCFFFCCC
raw rows=639
routes: Q-WITNESS=603, BOUNDARY-AP-GW=2,
        BOUNDARY-PETAL-SPORADIC=1, COVERING-MOMENT=33
```

Within that word:

- `q_threshold=14` reduces the fiber to `32` rows but still mixes AP/GW,
  petal, and covering routes.
- `q_factorization + power_lift_guard` isolates AP/GW as a two-row fiber.
- `M` isolates AP/GW as a two-row boundary fiber.
- packet labels split AP and GW into singleton fibers.

## Tournament Analysis

Tournament Analysis uses quotient-repair carriers as vertices, not runners.
The pairwise observable is route purity, open/boundary purity, family purity,
retained magnitude, retained topology, retained packet labels, noncircularity,
proof cost, and discriminating fiber count.

The S190 repair-carrier tournament is transitive:

```text
score_hist={0:1,1:1,...,13:1}
c3=0
scc=(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
hp=1
```

High-retention path:

```text
word_plus_guarded_packet_no_route
> word_plus_packet_labels
> word_plus_route
> word_plus_filter_exit
> word_plus_M_q
> word_plus_boundary_topology
> word_plus_q_threshold
> word_plus_q_class
> automatic_counts
> automatic_word
> word_plus_tail_bucket
> word_plus_q_factor_power
> word_plus_M
> word_plus_tail_exact
```

The ordering is a proof-carrier ranking, not a conceptual beauty contest:
global route-purity dominates status-purity, so exact `M` is demoted despite
being the open/boundary repair coordinate.

## Proof Target

The next theorem should not say "automatic quotients work."  It should say:

```text
For every primitive LRC14 packet in a fixed automatic sidecar fiber,
the first nonzero repair cochain in the ordered ladder
  exact scale/q, boundary topology, packet transfer labels, harmonic/dual
  certificate, residual F7/THM-572
is forced to fire.
```

This would turn automatic language into a scheduler for proof obligations,
not a scalar invariant.  The two remaining non-route mixed pairs from
`word_plus_M_q` and `word_plus_boundary_topology` are the smallest concrete
test cases for the local zipper step.
