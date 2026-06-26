---
id: HYP-3034
title: LRC14 arc-boundary path lift and owner-deletion persistence
status: EVIDENCE / targeted topology carrier scout; not a proof
source: codex-2026-06-26-S198
tangent: T1115
script: 04-computation/lrc14_arc_boundary_path_lift_codex_s198.py
result: 05-knowledge/results/lrc14_arc_boundary_path_lift_codex_s198.out
related:
  - HYP-3031
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3018
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3034: LRC14 Arc-Boundary Path Lift

## Claim

The HYP-3030 status-topology gate should be strengthened from a Betti scalar
to a boundary-operator carrier:

```text
closed danger-arc Cech nerve
  -> GF(2) boundary ranks d1,d2
  -> explicit H1 representative
  -> owner-deletion persistence
```

This mines the older tournament path-homology and deletion-persistence work,
but changes the vertex set.  The useful vertices are not runners and not a
runner tournament.  They are the individual closed danger arcs and the proof
boundary operators on their Cech nerve.

## Computation

Script:

```text
04-computation/lrc14_arc_boundary_path_lift_codex_s198.py
```

Stored output:

```text
05-knowledge/results/lrc14_arc_boundary_path_lift_codex_s198.out
```

The script uses the HYP-2963 default row bank:

```text
packets=21913
single_limit=180
two_swap_limit=36
alias_depth=4
lcm_tail_max=5
```

It deliberately avoids recomputing exact `M` for every row.  It first runs a
lightweight exact interval status scan over all rows, then computes the
expensive GF(2) path-boundary representative only on the residue-terminal
boundary/open collision surface and named controls:

```text
path_lift_rows=41
```

## Readout

The targeted path-lift surface has exactly two closed-H1 rows:

```text
closed_h1_rows=2
status_counts={'boundary': 2}
source_family_counts={'named-frontier': 2}
apgw_path_cycle_rows=2
```

They are AP and Goddyn-Wong:

```text
AP:
  C0/C1=1/1
  rank_d1=90
  rank_d2=84
  representative_edges=58
  cycle_owner_speeds=(1,2,3,4,5,6,7,8,9,10,11,12,13)

GW 12->24:
  C0/C1=1/1
  rank_d1=102
  rank_d2=90
  representative_edges=58
  cycle_owner_speeds=(1,2,3,4,5,6,7,8,9,10,11,13,24)
```

Deleting any one owner speed from either displayed representative kills the
closed-H1 signal:

```text
AP deletion beta1: each owner -> 0
GW deletion beta1: each owner -> 0
```

The two residue-terminal boundary/open collisions are path-lift separated:

```text
collision 0: size=30 status={'boundary':1,'open':29} h1={1:1,0:29}
collision 1: size=11 status={'boundary':1,'open':10} h1={1:1,0:10}
```

The lightweight full-bank status scan reconfirms the coarse gate fact:

```text
coarse_et_unit_gate fibers=21702
mixed_status_fibers=0
```

## Interpretation

HYP-3030 already showed that AP/GW are the only boundary rows inside the
residue-terminal status collisions and that they carry closed arc beta `(1,1)`.
HYP-3034 adds two pieces of structure:

1. The equality atom has an explicit GF(2) boundary representative, not just a
   scalar `beta1=1`.
2. That representative is owner-essential on the collision surface: deleting
   any owner speed from the displayed representative kills the closed cycle.

This gives a cleaner proof target for the topology-first gate.  The theorem
should not merely say "closed beta1 appears."  It should say that any zero-open
packet either carries an AP/GW-type owner-essential closed boundary cycle, or
emits named F7/THM-572 residual debt.

## Tournament Analysis

Vertices are proof carriers / boundary operators, not runners:

```text
arc_boundary_path_lift
owner_deletion_persistence
arc_cech_beta
coarse_et_unit_status
residue_terminal_word
```

Pairwise observable:

```text
status_purity
H1_exactness
deletion_persistence
topology_retention
endpoint_ownership
quotient_defect_visibility
proof_cost
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1]
hamiltonian_path_count=1
score_order=arc_boundary_path_lift >
            arc_cech_beta >
            owner_deletion_persistence >
            coarse_et_unit_status >
            residue_terminal_word
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, individual danger arcs, endpoint
cells, boundary cocircuits, residues, Fourier modes, matroid circuits, path
chains, and proof obligations.

The chosen vertices are closed danger arcs and boundary operators because the
preserved predicate is boundary/open status at threshold `1/14`, plus the exact
closed-cover obstruction when no strict safe interval exists.

This quotient destroys raw runner order and scalar route labels.  Owner
deletion records only the speeds supporting the chosen H1 representative.
Therefore this is not a route-purity classifier by itself; it is the topology
front gate before HYP-3030/HYP-3028 route scheduling.

## Proof Target

Arc-boundary path-lift theorem:

```text
Every primitive zero-open LRC14 packet either
  (a) has an owner-essential closed danger-arc H1 representative
      of AP/GW type,
  (b) is discharged by a named K33/F7/THM-572 or harmonic residual,
  (c) or is impossible.
```

Once (a)-(c) is established, coarse ET+unit status convergence may forget
route labels while preserving the LRC yes/no predicate, and open route-mixed
fibers can be scheduled by magnitude, barcode, Fejer/Ramanujan/Haar,
q-witness, covering, petal, or state-lift certificates.

## Next Work

1. Run a full closed-H1 scan only after adding a cached Cech sidecar.  The S198
   script intentionally limits expensive representatives to collision/control
   rows.
2. Prove owner-essentiality for the AP/GW closed cycle without enumerating
   every arc.
3. Compare owner-deletion persistence with HYP-3029 largest-stalk descent and
   HYP-3018 active-bottleneck normal-fan supports.
4. Add sidecar fields:

```text
closed_arc_boundary_rank_d1
closed_arc_boundary_rank_d2
closed_arc_h1_rep_edge_count
closed_arc_h1_owner_support
owner_deletion_beta1_word
arc_boundary_path_lift_exit
```
