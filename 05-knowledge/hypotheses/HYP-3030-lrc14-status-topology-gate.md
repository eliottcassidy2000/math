---
id: HYP-3030
title: LRC14 status-topology gate for coarse zipper fibers
status: EVIDENCE / full-bank bridge between topology and arithmetic status gates; not a proof
source: codex-2026-06-26-S194
tangent: T1111
script: 04-computation/lrc14_status_topology_gate_codex_s194.py
result: 05-knowledge/results/lrc14_status_topology_gate_codex_s194.out
related:
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-3018
  - HYP-3016
  - HYP-3015
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3030: LRC14 Status-Topology Gate

## Claim

The proof target after HYP-3024 and HYP-3025 should split into two gates:

```text
arc-Cech boundary cycle gate
  -> coarse ET + Henselian-unit status gate
  -> route/certificate scheduler
```

The first gate handles equality atoms.  The second gate may forget exact route
labels as long as it preserves boundary/open status.

This is not a proof of LRC14.  It is a sharper proof-interface claim: the
tested full HYP-2963 bank says boundary/open status is controlled before full
theorem-route purity is needed.

## Computation

Script:

```text
04-computation/lrc14_status_topology_gate_codex_s194.py
```

Stored output:

```text
05-knowledge/results/lrc14_status_topology_gate_codex_s194.out
```

The script audits the HYP-2963 default `21913` packets, focusing on:

1. The residue-terminal fibers that still mix boundary and open packets.
2. The route-mixed fibers that survive HYP-3024's coarse ET+unit gate.

## Readout

Residue-terminal fibers:

```text
fibers=16555
mixed_status=2
mixed_route=265
max_route_mixed=30
boundary_apgw_cycle=2
boundary_noncycle=0
open_with_closed_beta1=0
open_min_safe_topes=4
```

The two boundary/open collisions are exactly AP and GW against open impostors.
The boundary rows have:

```text
closed_arc_beta=(1,1)
open_arc_beta=(6,0)
safe_topes=0
boundary_owner_sums_mod_14=(0,0,0,0,0,0)
```

Every open cohabitant in those two mixed status fibers has closed arc
`beta1=0` and at least four strict safe topes.

Coarse ET+unit gate:

```text
fibers=21702
mixed_status=0
mixed_route=15
max_route_mixed=4
aggregate_status={'open': 38}
aggregate_routes={'COVERING-MOMENT': 15, 'Q-WITNESS': 23}
open_with_closed_beta1=0
open_min_safe_topes=4
```

Thus the remaining route-mixed coarse fibers are harmless for the direct LRC
predicate: all `38` packets in those fibers are strictly open and have closed
arc `beta1=0`.

## Interpretation

HYP-3024 made the arithmetic observation: coarse Erdos-Turan clocks plus
Henselian unit-root counts have `0` mixed boundary/open fibers on the full
bank.

HYP-3025 made the topological observation: AP/Goddyn-Wong are closed arc cover
cycles, while named positive rows have closed arc `beta1=0`.

HYP-3030 says these observations should be ordered, not merged into one bulky
signature.  First prove that zero-open packets must be AP/GW arc-boundary
cycles or a named residual.  Then use the coarse ET+unit gate only as a
status-preserving quotient and let magnitude/barcode/Fejer/K33 machinery route
the remaining open fibers.

This also clarifies the incoming switchboard stack: HYP-3029 tests safe-component stalk descent, HYP-3026 supplies the
large labelled packet fusion carrier, HYP-3027 tests repair ladders for
automatic quotient failures, and HYP-3028 names the residual status-gate target.
HYP-3030 is the topology witness inside that stack: the equality atoms are
exactly the closed arc-cycle rows before the residual open-route scheduler even
starts.

## Tournament Analysis

Vertices are proof gates, not runners:

```text
arc_boundary_cycle_gate
coarse_et_unit_status_gate
magnitude_route_splitter
barcode_packet_scheduler
raw_residue_terminal_word
```

Fingerprint:

```text
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
directed_3cycles=0
hamiltonian_path_count=1
score_order=arc_boundary_cycle_gate >
            coarse_et_unit_status_gate >
            magnitude_route_splitter >
            barcode_packet_scheduler >
            raw_residue_terminal_word
```

## Assumption Challenge

Alternate vertex sets considered: runners, individual danger arcs, endpoint
cells, boundary cocircuits, residues, Fourier/ET clocks, p-adic unit roots,
magnitude fibers, barcode packets, and proof obligations.

The chosen vertices are proof gates because the preserved predicate is not the
full route label.  It is the yes/no LRC predicate: whether the row has a strict
safe interval at threshold `1/14`.

The quotient destroys raw runner identity and exact route labels.  This is
allowed only after the arc topology has separated equality atoms and the
coarse ET+unit gate has certified status constancy.

## Proof Target

Status-topology gate theorem:

```text
Every primitive zero-open LRC14 packet either
  (a) carries the AP/GW arc-boundary cycle
      closed beta=(1,1), open beta=(6,0), six zero owner sums,
  (b) emits named K33/F7/THM-572 residual debt, or
  (c) is impossible.

Once (a)-(c) is established, every remaining coarse ET+unit route collision
is strict-open and can be discharged by q-witness, covering, magnitude,
barcode, Fejer/Ramanujan/Haar, or state-lift certificates.
```
