---
id: HYP-3016
title: LRC14 automaton fiber-mixing and magnitude-cocycle guardrail
status: SYNTHESIS / exact bounded fiber-mixing audit and proof-interface guardrail; not a proof
source: codex-2026-06-26-S175
script: 04-computation/lrc14_automaton_fiber_mixing_codex_s175.py
result: 05-knowledge/results/lrc14_automaton_fiber_mixing_codex_s175.out
related:
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3011
  - HYP-3010
  - HYP-3009
  - HYP-3008
  - HYP-3002
  - HYP-2997
  - HYP-2963
  - HYP-2928
  - HYP-2927
  - THM-572
  - OPEN-Q-108
---

# HYP-3016: LRC14 Automaton Fiber-Mixing

## Claim

The Moser/fibbinary finite-automaton fields are useful packet side channels,
but residue plus automaton terminal state is not theorem-safe for the LRC14
tight-locus problem.  Exact bounded fiber mixing shows why:

```text
fixed residue-language fiber can contain both
  AP/Goddyn-Wong zero-open boundary atoms
and
  strictly open perturbations with positive safe measure.
```

Therefore any HYP-2963 quotient using automatic language data must also retain
or discharge a magnitude/Farey coordinate inside each residue-automaton fiber.
The proposed name for that retained lost coordinate is the
`magnitude_cocycle`.

## Computation

Script:

```text
04-computation/lrc14_automaton_fiber_mixing_codex_s175.py
```

Output:

```text
05-knowledge/results/lrc14_automaton_fiber_mixing_codex_s175.out
```

The S175 row bank contains named LRC14 rows and the AP single-swap atlas

```text
{1,...,13}\{drop} union {add}, 14 <= add <= 180
```

with exact rational interval-union safe-measure computation at threshold
`1/14`.

Summary:

```text
row_bank=2172 distinct primitive rows
status_counts={'boundary': 2, 'open': 2170}
min_positive_safe_measure=1/1260 at K33_12_to_36
```

The only zero-open rows in this bank are AP and Goddyn-Wong:

```text
AP13          mu=0
GW_12_to_24   mu=0
K33_12_to_36  mu=1/1260
loose_12_to_26 mu=426/35035
loose_12_to_96 mu=5219/840840
```

## Fiber-Mixing Evidence

The key table:

```text
mfc_word                 mixed_boundary_open=1
terminal_word            mixed_boundary_open=1
residue_multiset_mod14   mixed_boundary_open=2
residue_mfc_pairs        mixed_boundary_open=2
residue_terminal_pairs   mixed_boundary_open=2
residue_mfc_bitlen_pairs mixed_boundary_open=0 in this bounded bank
exact_speed_tuple        mixed_boundary_open=0
```

The important collision is concrete.  Under `residue_mfc_pairs`, AP shares its
fiber with strictly open rows:

```text
AP13             mu=0
loose_12_to_96   mu=5219/840840
loose_12_to_26   mu=426/35035
```

The Goddyn-Wong one-dipole fiber also mixes:

```text
GW_12_to_24      mu=0
swap_12_to_38    mu=426/35035
swap_12_to_52    mu=426/35035
swap_12_to_150   mu=426/35035
```

So coarse automatic state does not merely fail to distinguish AP from GW.  It
fails to distinguish boundary equality from strict openness inside the same
residue-language fiber.

Adding bit length removes mixed boundary/open fibers in this bounded atlas, but
that is a proxy rather than a theorem invariant.  The theorem-facing coordinate
is exact magnitude/Farey scale, together with endpoint owners and safe
components.

## Tournament Analysis

Vertices are quotient/proof carriers, not runners:

```text
exact_labelled_packet
farey_magnitude_height
residue_mfc_bitlen_pairs
residue_terminal_pairs
residue_mfc_pairs
terminal_word
mfc_word
perfect_power_word
gap_ratio_bucket
raw_counts
```

Pairwise observable:

```text
fiber_purity
finite_state_checkability
magnitude_retention
residue_endpoint_retention
route_compatibility
anti_scalar_guard
proof_cost
```

Fingerprint:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1),(8,1),(9,1)]
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_order=exact_labelled_packet > farey_magnitude_height >
  residue_mfc_bitlen_pairs > residue_terminal_pairs > residue_mfc_pairs >
  terminal_word > mfc_word > perfect_power_word > gap_ratio_bucket > raw_counts
```

The transitive order is not a ranking of mathematical beauty.  It is a proof
admissibility order: exact packet data and Farey magnitude dominate because
the finite-state residue shadows are mixed.

## Proposed Packet Fields

Add these fields to HYP-2963 sidecars where automaton data are used:

```text
residue_automaton_fiber_id
automaton_terminal_word
magnitude_cocycle
farey_magnitude_height
fiber_anchor_row
safe_component_measure
safe_component_count
largest_safe_component
bitlength_proxy_status
fiber_mixing_exit
```

## Theorem Target

Magnitude-cocycle packet theorem:

```text
For a primitive LRC14 packet in a fixed residue-automaton fiber, either
  (a) the magnitude cocycle is zero and the packet is AP/GW boundary equality,
  (b) the magnitude cocycle opens a strict safe component,
  (c) the packet descends to a proved single-swap / family formula,
  (d) a Fejer/Ramanujan/Haar/endpoint certificate annihilates the cocycle, or
  (e) the cocycle emits named K33/F7/THM-572 residual debt.
```

This is the automaton version of the no-free-slider rule.  The finite automaton
can label the packet, but it cannot replace the magnitude coordinate.
