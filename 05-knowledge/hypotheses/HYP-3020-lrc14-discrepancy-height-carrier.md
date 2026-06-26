---
id: HYP-3020
title: LRC14 discrepancy-height trident carrier
status: EVIDENCE / bounded carrier scout and proof-interface target; not a proof
source: codex-2026-06-26-S184
tangent: T1102
script: 04-computation/lrc14_discrepancy_height_carrier_codex_s184.py
result: 05-knowledge/results/lrc14_discrepancy_height_carrier_codex_s184.out
related:
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3011
  - HYP-3009
  - HYP-3008
  - HYP-2997
  - HYP-2995
  - HYP-2991
  - HYP-2989
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3020: LRC14 Discrepancy-Height Trident Carrier

## Claim

The HYP-3016/HYP-3017 automaton route-purity failure suggests a better
carrier than a raw finite-state word:

```text
residue discrepancy  +  height/Mahler scale  +  Hensel singular-lift status
```

This "trident" is a proof-interface coordinate, not a certificate.  Its job is
to split the mixed automaton fibers where AP/Goddyn-Wong boundary atoms share a
finite-state shadow with strictly open rows, while still pointing toward
Fejer/Ramanujan/Haar or p-adic exits.

## Computation

Script:

```text
04-computation/lrc14_discrepancy_height_carrier_codex_s184.py
```

Stored output:

```text
05-knowledge/results/lrc14_discrepancy_height_carrier_codex_s184.out
```

The bounded bank contains named rows plus the AP single-swap atlas through tail
`180`:

```text
row_bank=2173 distinct primitive rows
status_counts={'boundary': 2, 'open': 2171}
min_positive_safe_measure=1/1260 at K33_12_to_36
```

The negative automaton result persists:

```text
automatic_word                       mixed_status=1
residue_mfc_pairs_mod14              mixed_status=2
residue_l1_14_27_41                  mixed_status=2
hensel_2_3_7                         mixed_status=1
mahler_height_bucket                 mixed_status=1
residue_discrepancy_height           mixed_status=1
```

The combined trident removes boundary/open mixing in this bounded bank:

```text
discrepancy_height_hensel_trident    fibers=2167
                                     mixed_status=0
                                     mixed_route=6
                                     largest_fiber=2
```

This is useful evidence but also a warning: the trident is very close to exact
packet identity in the bounded scout.  The next theorem must compress it
without reintroducing the automaton fiber leak.

## Collision Readout

The AP residue+MFC fiber contains AP and many open tails:

```text
AP13               boundary  height=(4,13,M<0.90)
loose_12_to_26     open      height=(5,26,M<0.90)
loose_12_to_96     open      height=(7,96,M<0.90)
```

The GW one-dipole fiber also mixes:

```text
GW_12_to_24        boundary  height=(5,24,M<0.90)
GW_fiber_12_to_38  open      height=(6,38,M<0.90)
GW_fiber_12_to_52  open      height=(6,52,M<0.90)
GW_fiber_12_to_150 open      height=(8,150,M<0.90)
```

The height coordinate already sees the lost magnitude cocycle, while residue
discrepancy and Hensel counts tell which denominator/no-lift channels should
receive the next certificate.

## Beck-Fiala Incidence Readout

The script also treats packet fields as a feature hypergraph:

```text
feature_token_count=353
row_feature_arity_range=(12,12)
max_token_support=2005
```

Read as a proof heuristic: each packet belongs to boundedly many feature
families, so a Beck-Fiala-style signed-incidence argument is plausible only
after high-support tokens are split by another coordinate.  The high-support
tokens are exactly the leaky coarse fields such as common Hensel status and
common residue discrepancy.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
exact_labelled_packet
safe_topology_barcode
discrepancy_height_hensel_trident
residue_discrepancy_height
hensel_singular_lift_guard
erdos_turan_residue_discrepancy
mahler_height_proxy
automatic_word_sidecar
raw_scalar_family_name
```

Pairwise observable:

```text
predicate purity
route purity
magnitude retention
discrepancy retention
local lift stability
certificate handoff
finite cost
anti-scalar guard
```

Fingerprint:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1),(8,1)]
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_order=discrepancy_height_hensel_trident > exact_labelled_packet >
  residue_discrepancy_height > safe_topology_barcode >
  erdos_turan_residue_discrepancy > mahler_height_proxy >
  hensel_singular_lift_guard > automatic_word_sidecar >
  raw_scalar_family_name
```

The tournament rank is not a beauty contest.  It says the trident is the best
next exploratory carrier because it is cheaper than exact packet identity but
retains the coordinates that the automaton quotient destroyed.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, p-adic
roots, feature tokens, persistence bars, and proof obligations.

The chosen vertices are proof carriers.  The trident preserves the
boundary-vs-open predicate in the bounded bank, route-purity telemetry,
magnitude/Farey scale, residue-discrepancy side data, and p-adic singular-lift
flags.  It destroys exact endpoint-owner geometry, exact Fejer atom banks, and
full safe-component topology unless HYP-3015 barcode and HYP-2981 interval
certificate fields are reattached.

## Theorem Target

Discrepancy-height compression theorem:

```text
Inside each mixed residue/automatic fiber, either
  (a) the compressed discrepancy-height-Hensel signature is AP/GW boundary,
  (b) nonzero height cocycle opens a strict safe component,
  (c) residue discrepancy routes to a Ramanujan/Haar/Fejer certificate,
  (d) Hensel singular status routes to finite p-adic lift debt, or
  (e) the packet is named K33/F7/THM-572 residual debt.
```

The immediate engineering pull is to add the trident as a sidecar to the full
HYP-2963 classifier, then try to coarsen the sharp S184 signature until it is
as small as possible while keeping `mixed_status=0` on the full packet bank.
