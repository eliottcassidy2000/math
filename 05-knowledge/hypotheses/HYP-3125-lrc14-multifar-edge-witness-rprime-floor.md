---
id: HYP-3125
title: LRC14 multi-far edge-witness Rprime floor
status: EVIDENCE / executable edge-floor synthesis scout; not a proof
source: codex-2026-06-27-S271
tangent: T1199
technique: LTI-260
tournament_technique: LTT-158
script: 04-computation/lrc14_multifar_edge_witness_floor_codex_s271.py
result: 05-knowledge/results/lrc14_multifar_edge_witness_floor_codex_s271.out
related:
  - HYP-3127
  - HYP-3124
  - HYP-3123
  - HYP-3122
  - HYP-3121
  - HYP-3120
  - HYP-3118
  - HYP-3116
  - HYP-3112
  - HYP-3108
  - HYP-3106
  - HYP-3101
  - HYP-2968
  - HYP-2963
  - THM-573
  - THM-572
  - THM-082
  - OPEN-Q-108
---

# HYP-3125: LRC14 Multi-Far Edge-Witness Rprime Floor

## Claim

The open `r=2..6` covering branch from HYP-3121 should be treated as an
edge-witness floor problem.  After the `u=14t` lift, write
`S = R union 14Q`.  The target is

```text
Rprime(R,Q) =
  mu(R-safe intersect Q-safe) / (mu(R-safe) * mu(Q-safe)) >= c_floor > 0
```

on the multi-far packets, with finite named exceptions routed to a terminal
packet.  HYP-3124 says a directed edge is not a scalar relation; it is a
two-ended witness carrying a four-sector deck and both endpoint-deletion
children.  Here the directed edge is:

```text
R-safe packet -> Q-safe packet
```

and `Rprime` is the diagonal sector mass of that edge after normalization.
The proposed proof object is:

```text
edge_floor_packet =
  edge_witness_id
  + tail_R_safe_packet
  + tip_Q_safe_packet
  + outside_four_sector_deck
  + edge_tail_tip_sector_word
  + tail_deletion_child_Rprime
  + tip_deletion_child_Rprime
  + Rprime_lower_bound_candidate
  + multifar_r_bucket
  + EH_level_distribution_proxy
  + major_arc_residue_exception_set
  + gaussian_heat_kernel_width
  + finite_ruler_desmoothing_threshold
  + asano_contraction_status
  + lee_yang_zero_free_after_contraction
  + phi4_kappa4_sign
  + normal_fan_chamber_id
  + chiral_guard_word
  + terminal_exit_or_named_debt
```

This is not a proof yet.  It is a sharper carrier for the proof obligation:
either a Gaussian/EH-style distribution estimate gives a positive uniform
edge diagonal after finite-ruler desmoothing, or the bad packet is forced to
carry a named Asano/Lee-Yang/phi4/Cech/H7 exit.

## Evidence From S271

The script
`04-computation/lrc14_multifar_edge_witness_floor_codex_s271.py` performs a
bounded midpoint-grid audit in the lifted variable `u in [0,14)`, measures
`R_safe`, `Q_safe`, their joint mass, `Rprime`, tail-deletion children,
tip-deletion children, individual edge-sector ratios, and a residue-bin
level-of-distribution proxy.  Output is stored in
`05-knowledge/results/lrc14_multifar_edge_witness_floor_codex_s271.out`.

Main numerical readout:

| row | r | Bonferroni | Rprime | tail-child Rprime range | tip-child Rprime range |
|---|---:|---:|---:|---:|---:|
| `{1..11,13,84}` | 1 | `-0.130833` | `0.513784` | `[0.730195,0.976879]` | `[1.000000,1.000000]` |
| `{1..10,13,84,154}` | 2 | `-0.182512` | `0.925326` | `[0.881755,1.018485]` | `[0.953575,0.972941]` |
| probe r=3 | 3 | `-0.264655` | `0.935073` | `[0.875850,1.029144]` | `[0.928880,1.007078]` |
| probe r=4 | 4 | `-0.287333` | `0.966976` | `[0.930200,1.050907]` | `[0.943605,1.013613]` |

The first two rows reproduce the HYP-3121 lesson: Bonferroni is negative, but
quasi-independence is positive.  The new S271 point is recursive.  Tail
deletion usually improves the floor, while tip deletion exposes the
multiple-of-14 side as the sharper child.  That is exactly the shape predicted
by HYP-3124: one-sided safe mass is not a legal witness; the floor must carry
both endpoint children.

The individual R-speed / Q-speed edge ratios are close to `1`, with worst
edge `(7,6)` in the sampled rows.  That suggests the obstruction is not a
single bad pair edge but a packet-level interaction among many sector cuts.
The proof should therefore look for an average level-of-distribution theorem
over edge-sector residues plus a finite major-arc exception ledger, not a
local pairwise Bonferroni repair.

## EH, Gaussian, And Asano Roles

The Elliott-Halberstam connection is only a structural analogy.  No prime
distribution theorem is being imported into LRC.  The transferable idea is a
level-of-distribution sidecar:

```text
EH_level_distribution_proxy =
  L2/BV average cancellation of edge-sector residue errors
  outside a finite major-arc exception set.
```

The Gaussian role is a positivity engine.  Replace hard safe indicators by a
Gaussian/heat-kernel or Selberg-type minorant, prove a smoothed diagonal lower
bound, then pay a finite-ruler desmoothing debt using HYP-3123 normal-fan/Cech
fields.

The Asano role is a legality check for contraction.  The tail and tip
variables of an edge-sector partition packet should only be contracted if the
contraction preserves a Lee-Yang-style zero-free sidecar.  HYP-3122 supplies
the phi4/quartic cumulant stabilizer; THM-082 supplies the deletion-contraction
warning that the tail-to-tip convention matters.

## Tournament Analysis

Tournament vertices are repair operators and proof carriers, not runners,
primes, Gaussian samples, or raw ratios.  Pairwise observable: majority
retention over LRC predicate retention, attack on the multi-far floor,
tail/tip recursion, decorrelation power, zero-free contraction, finite packet
fields, formalization readiness, and failure guardrails.

S271 fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,5:2,6:2,8:1,9:1}
directed_3cycles=2
scc_sizes=[4,1,1,1,1,1,1]
hamiltonian_path_count=5
selected_path =
  edge_witness_two_ended_RQ_packet
  -> uniform_multifar_Rprime_floor
  -> gaussian_heat_minorant_smoothing
  -> normal_fan_cech_finite_ruler
  -> asano_lee_yang_edge_contraction
  -> chiral_cross_sector_guard
  -> phi4_quartic_cumulant_stabilizer
  -> H7_state_lift_zero_kernel
  -> EH_level_distribution_proxy
  -> raw_Bonferroni_or_scalar_p0
```

The non-transitive component is useful rather than a flaw: it says the next
proof step should not choose between edge recursion, the raw floor, Gaussian
smoothing, and residue distribution.  It should join them in one packet.

## Post-Rebase Integration With HYP-3127

After this S271 scout was committed locally, incoming HYP-3127/S68 arrived on
`origin/main` with a stronger formulation of the same route: the multi-far
floor is proposed as an Asano contraction of single-far Lee-Yang factors.
That does not obsolete this packet.  It changes its role.

HYP-3127 promotes Asano from a legality sidecar to a candidate main engine:
coverage is multi-affine by inclusion-exclusion, each far element is a tip
factor, and the HYP-3124 tail/tip recursion is the contraction order.  It also
sharpens the Gaussian and EH roles: Gaussian is the free-field/decoupled
limit, while EH is a level-of-distribution/SPEC-bound ideal rather than a
load-bearing theorem.  The S271 `edge_floor_packet` is therefore the row schema
needed to audit HYP-3127's three obligations:

1. single-far zero-free polydisk status,
2. `SPEC` or signed-cancellation bound for the floor constant,
3. recursion termination/monotonicity through tail and tip deletion children.

Read this HYP-3125 artifact as the local diagnostic and packet ledger for the
HYP-3127 reduction, not as a competing route.

## Next Theorem Target

For every `r=2..6` multi-far covering packet:

1. Construct the HYP-3124 two-ended edge packet
   `R-safe packet -> Q-safe packet`.
2. Prove an LRC-specific level-of-distribution bound for the edge-sector
   residue errors after removing a finite major-arc exception set.
3. Use Gaussian/heat smoothing to obtain a positive diagonal lower bound.
4. Desmooth through the HYP-3123 finite-ruler/Cech component payload.
5. If the desmoothing or distribution step fails, require an
   Asano/Lee-Yang/phi4 sidecar or route the packet to H7 state-lift debt.

Challenged assumption: `Rprime` is not a scalar property of the set `S`.
It is information in a directed proof edge.  The tail knows how much safe mass
survives after deleting an `R` constraint; the tip knows how the multiple-of-14
lonely packet changes after deleting a `Q` constraint.  A proof must preserve
both recursions or name the coordinate it destroys.
