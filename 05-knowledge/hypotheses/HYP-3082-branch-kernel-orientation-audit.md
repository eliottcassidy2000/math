---
id: HYP-3082
title: Branch-kernel orientation audit for the HYP-2963 LRC14 packet bank
status: EVIDENCE / bounded proof-interface audit; not a proof
source: codex-2026-06-26-S250
tangent: T1166
technique: LTI-231
tournament_technique: LTT-129
script: 04-computation/lrc14_branch_kernel_orientation_codex_s250.py
result: 05-knowledge/results/lrc14_branch_kernel_orientation_codex_s250.out
forum_posts:
  - poke-forum/posts/20260626-lrc14-branch-kernel-orientation-audit/post.md
related:
  - HYP-3081
  - HYP-3079
  - HYP-3078
  - HYP-3077
  - HYP-3076
  - HYP-3075
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3068
  - HYP-3067
  - HYP-3066
  - HYP-3060
  - HYP-3058
  - HYP-3057
  - HYP-3056
  - HYP-3054
  - HYP-3009
  - HYP-2996
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3082: Branch-Kernel Orientation Audit For The HYP-2963 LRC14 Packet Bank

## Claim

The S250 audit makes HYP-3081 executable after the S249 branch-tournament
guardrail, the HYP-3079 Lean q-cusp finite-tail ledger, the HYP-3078 q-cusp
principal-part gate, the HYP-3077 median-hull scheduler, the HYP-3076/HYP-3075
arithmetic sidecar discipline, the HYP-3074 route-state closure layer, the
HYP-3073 polymer/Dirichlet stub, the HYP-3072 cross-carrier stub, and the
HYP-3071/HYP-3070 observability/control layers on the existing HYP-2963/S168 packet
bank.  In that bounded bank, every hard non-AP/GW packet has a declared
section/grid exit or named K33 state-lift debt before the Robbins bridge test
is applied.  After known sidecars are attached, the protected branch graph has
no naked bridge.

This is not a proof of LRC14.  It is a proof-interface checkpoint: the local
packet-bank branch graph is Robbins-safe, so the next theorem obligation is no
longer "find a missing scalar" but "prove every primitive residual reaches
this branch graph, then discharge K33/THM-572 and covering-family exits without
creating a new naked bridge."

## Evidence

The script audits the default packet bank:

```text
audited packets = 21913
hard non-AP/GW packets = 7235
candidate F7 packets = 0
```

Route census:

```text
Q-WITNESS                14676
COVERING-MOMENT          7228
BOUNDARY-PETAL-SPORADIC     4
K33-STATE-LIFT              3
BOUNDARY-AP-GW              2
```

The section/grid ledger exactly matches the residual-section packet-grid
vocabulary:

```text
direct_q_witness_section             14676  orthogonal_zero
open_haar_witness_section             6040  vertical_owner_strip
covering_boundary_moment_section      1188  nested_refinement
unit_petal_c27_strip_section             4  horizontal_owner_strip
k33_state_lift_section                   3  cross_handoff
ap_gw_boundary_section                   2  same_tile_indicator
```

The raw scalar collapse is unsafe:

```text
raw scalar star:
  nodes = 6
  edges = 5
  bridges = 5
  naked_bridges = 5
```

After attaching branch/section/grid/finalizer/no-lift sidecars, the protected
branch graph has no naked bridge:

```text
protected branch graph:
  nodes = 80
  edges = 83
  bridges = 69
  naked_bridges = 0

contracted protected core:
  nodes = 1
  edges = 0
  bridges = 0
  strong_orientation_exists = True
```

The bridge certificates are not a single reason.  They split across direct
q-witnesses, AP/GW stops, C27 owner-strip discharge, K33 state-lift debt,
section/grid exits, positive open or nested-refinement exits, exact schema,
power-lift sidecars, route kernels, Desargues/Beal finalizer gate, and named
residual debt.

## Branch Ledger

The theorem-facing ledger rows are:

```text
Q-WITNESS                direct_q_witness_section           orthogonal_zero        rows=14676 exit=direct-q-witness
COVERING-MOMENT          open_haar_witness_section          vertical_owner_strip   rows= 6040 exit=strict-safe-component
COVERING-MOMENT          covering_boundary_moment_section   nested_refinement      rows= 1188 exit=strict-safe-component
BOUNDARY-PETAL-SPORADIC  unit_petal_c27_strip_section       horizontal_owner_strip rows=    4 exit=petal-discharge
K33-STATE-LIFT           k33_state_lift_section             cross_handoff          rows=    3 exit=state-lift-debt
BOUNDARY-AP-GW           ap_gw_boundary_section             same_tile_indicator    rows=    2 exit=boundary-equality
```

The low-frontier representatives are AP, GW `12->24`, K33 `12->36`, P10+GW,
petal `10->20`, petal `13->26`, P10+K33, and
`two drop(12, 13)->add(26, 36)`.  The four C27/petal representatives carry
`q=27=3^3` power stress and are therefore routed through the HYP-3058/HYP-3081
no-lift discipline rather than treated as scalar coincidences.

## Coarse Quotient Warning

Several coarse quotients still create mixed route fibers:

```text
automatic_word                                      mixed_fibers=143
branch                                              mixed_fibers=33
q_factorization                                     mixed_fibers=36
power_lift_guard                                    mixed_fibers=14
automatic_word+q_factorization                      mixed_fibers=375
automatic_word+q_factorization+power_lift_guard     mixed_fibers=366
automatic_filter_exit                               mixed_fibers=0
```

The important readout is not that `automatic_filter_exit` is magically final.
It is that coarse scalar or sequence shadows remain theorem-unsafe until the
branch sidecars decide whether their route bridges are protected exits or real
proof obligations.

## Tournament Analysis

Vertices are branch carriers, residual sections, Haar grid exits, proof gates,
no-lift guards, and named residual debts, not runners.  The pairwise observable
is retained exact scale, boundary/open status, route handoff, owner current,
observer-cut payload, bridge protection, no-lift discipline, and named
residual exit.  The switch is the lexicographic retained-payload vector.

The S250 tournament is transitive:

```text
score_hist = {0:1, 1:1, ..., 12:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
tie_path = observer_cut_payload_orbit
         > Robbins_no_bridge_assembly
         > labelled_packet_bank
         > residual_section_exit
         > Haar_grid_exit
         > endpoint_owner_strip
         > C27_petal_branch
         > K33_state_lift_branch
         > covering_moment_branch
         > Fermat_Catalan_no_lift_guard
         > Desargues_Beal_finalizer_gate
         > named_residual_debt
         > raw_scalar_shadow
```

## Assumption Challenge

The audit explicitly did not assume tournament vertices must be runners or
arcs.  It considered runners, gaps, fixed circle sections, section boundaries,
wall crossings, residues, cover arcs, Fourier modes, matroid circuits, branch
ears, endpoint interfaces, observer-cut payloads, and proof obligations.  The
chosen quotient preserves exact `M>=1/14`, q-threshold, open/boundary status,
endpoint owner, C27/K33 route handoff, section exit, power-lift guard, and
residual name.  It destroys raw runner identity, scalar `M` bucket, automatic
word alone, and raw power or graph numerology after their sidecars have been
named.

## Limits

HYP-3082 does not prove:

```text
global reduction to the HYP-2963 packet bank
THM-572
the family theorem for all covering tails
K33 state-lift impossibility
```

It gives a cleaner closure target: prove that every primitive residual reaches
the S250 protected branch graph, then discharge the remaining K33/THM-572 and
covering-family exits without producing a new naked proof bridge.
