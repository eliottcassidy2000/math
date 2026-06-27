---
id: HYP-3023
title: LRC14 automatic fiber zipper and magnitude-cocycle splitter
status: EVIDENCE / full-bank quotient splitter and proof-interface target; not a proof
source: codex-2026-06-26-S187
script: 04-computation/lrc14_automatic_fiber_zipper_codex_s187.py
result: 05-knowledge/results/lrc14_automatic_fiber_zipper_codex_s187.out
related:
  - HYP-3022
  - HYP-3021
  - HYP-3020
  - HYP-3019
  - HYP-3018
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3009
  - HYP-3008
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3023: LRC14 Automatic Fiber Zipper

## Claim

The full HYP-2963 labelled-packet bank supports the S175/S181 guardrail in a
sharper form:

```text
automatic word              = warning side channel
residue-terminal automaton  = strong but still mixed side channel
magnitude cocycle           = first tested non-route splitter with zero
                              mixed theorem-route fibers on this bank
barcode and packet zippers  = certificate-anchor refinements
```

This is not a proof of LRC14.  It is a proof-interface result: any finite
automaton / Moser / fibbinary / Fermat-Catalan quotient must prove a
familywise magnitude-cocycle lemma before it can be used as a theorem
coordinate.

## Computation

Script:

```text
04-computation/lrc14_automatic_fiber_zipper_codex_s187.py
```

Stored output:

```text
05-knowledge/results/lrc14_automatic_fiber_zipper_codex_s187.out
```

The S187 script imports the HYP-2963 classifier and audits the default
`21913`-packet bank with a zipper ladder:

```text
automatic_word
residue_terminal_fiber
magnitude_cocycle
barcode_shadow
packet_zipper
route_labelled_packet
```

The full-bank split table:

```text
split                     fibers    mixR    mixF    maxF  maxMix  purity
automatic_word               225     143     178    1179    1179   36.4%
residue_terminal_fiber     16555     265     265      30      30   98.3%
magnitude_cocycle          21909       0       0       2       0  100.0%
barcode_shadow             21913       0       0       1       0  100.0%
packet_zipper              21913       0       0       1       0  100.0%
route_labelled_packet      21913       0       0       1       0  100.0%
```

The largest mixed fibers before the magnitude split are:

```text
automatic_word:
  max_mixed=1179
  routes={'Q-WITNESS': 1124, 'COVERING-MOMENT': 55}

residue_terminal_fiber:
  max_mixed=30
  routes={'Q-WITNESS': 27, 'COVERING-MOMENT': 2, 'BOUNDARY-AP-GW': 1}
```

Thus automatic words and residue-terminal fields are useful search telemetry
but unsafe quotients.  The exact magnitude cocycle

```text
(M, q_threshold, farey_excess, lacunary_tail_ratio)
```

has no mixed theorem-route fibers in this full-bank run.

## Target Word

The AP/Goddyn-Wong automatic word from HYP-3017 remains the first target:

```text
MFCMMCCFFFCCC
```

Inside that word fiber:

```text
rows=639
routes={'Q-WITNESS': 603, 'BOUNDARY-AP-GW': 2,
        'BOUNDARY-PETAL-SPORADIC': 1, 'COVERING-MOMENT': 33}
families={'infinite q-witness family': 603,
          'AP/Goddyn-Wong boundary family': 2,
          'unit-petal sporadic family': 1,
          'open-Haar loose family': 29,
          'covering boundary-moment family': 4}
distinct M values=33
distinct residue-MFC fibers=181
```

Residue-terminal fibers are still mixed inside this target word:

```text
residue_terminal_fiber: fibers=181, mixed_route=27
largest_mixed_rows=30
routes={'Q-WITNESS': 27, 'BOUNDARY-AP-GW': 1, 'COVERING-MOMENT': 2}
```

The sample mixed fiber contains AP, covering rows, and q-witness rows:

```text
AP                                  M=1/14  route=BOUNDARY-AP-GW
magnitude liar 12->96               M=8/101 route=COVERING-MOMENT
single swap 12->180                 M=3/37  route=COVERING-MOMENT
single swap 11->53                  M=1/11  route=Q-WITNESS
two drop(11,12)->add(25,26)         M=1/11  route=Q-WITNESS
two drop(11,13)->add(25,27)         M=1/11  route=Q-WITNESS
```

After the magnitude split:

```text
magnitude_cocycle: fibers=638, mixed_route=0
barcode_shadow:    fibers=639, mixed_route=0
packet_zipper:     fibers=639, mixed_route=0
```

This makes `MFCMMCCFFFCCC` a concrete family-template problem: prove the
`33` exact magnitude values in this automatic word route to AP/GW equality,
q-witness, petal, covering, or existing dual-certificate anchors.

## Tournament Analysis

Vertices were quotient / sidecar bundles over HYP-2963 packets, not runners.
The pairwise observable was route purity, max mixed-fiber size, magnitude
retention, safe-topology retention, packet-label retention, finite-state
checkability, and proof cost.  The switch was majority comparison of those
observable vectors, with the zipper order used as the tie Hamiltonian path.

Fingerprint:

```text
score_hist={0: 1, 1: 1, 3: 3, 5: 1}
directed_3cycles=1
scc_sizes=[3, 1, 1, 1]
hamiltonian_path_count=3
score_order=packet_zipper > barcode_shadow > magnitude_cocycle >
            route_labelled_packet > residue_terminal_fiber >
            automatic_word
```

The nontrivial SCC is a useful warning: the exact packet zipper, barcode
shadow, and magnitude cocycle should be treated as a linked proof bundle, not
as a scalar ranking.

## Assumption Challenge

The candidate vertex sets considered were runners, gaps, fixed circle
sections, endpoint owners, residues, residue-automaton fibers, persistence
bars, magnitude cocycles, C27/K33 packet labels, and proof obligations.

The chosen vertices are quotient/sidecar bundles because they preserve the
predicate being tested: theorem-route purity and strict-open versus boundary
equality at threshold `1/14`.  They intentionally destroy raw runner identity
and raw sequence-name information.  The challenged assumption is that an
automatic word can become a proof coordinate without a fiber-purity theorem.

## Proof Target

The theorem-shaped next step is:

```text
Inside every automatic/residue packet fiber, the magnitude cocycle either
  (a) is an AP/Goddyn-Wong boundary equality,
  (b) opens a strict safe component,
  (c) descends to a q-witness or known source-family formula,
  (d) is annihilated by a barcode / Fejer / Ramanujan / Haar certificate, or
  (e) emits K33/F7/THM-572 residual debt.
```

The full-bank evidence says this is the right coordinate to attack.  The
first local proof target is the `MFCMMCCFFFCCC` fiber and its `33` exact
magnitude values.
