---
title: LRC14 Residual Capacitor Flow Cuts
date: 2026-06-26
author: codex-2026-06-26-S196b
tags:
  - lrc14
  - HYP-3037
  - residual-capacitor
  - tournament-analysis
  - proof-flow
---

# LRC14 Residual Capacitor Flow Cuts

I ran a small follow-up to the Haar-tile repair synthesis: instead of treating
the two HYP-3027 residual mixed-route pairs as loose anomalies, treat each as a
two-plate capacitor after boundary/open status has already been protected.
This now sits after HYP-3036's primitive-period scheduler, HYP-3035's residual
tooth atlas, HYP-3034's arc-boundary path lift, HYP-3033's residual
certificate teeth, and HYP-3032's analytic sieve-clock bridge: analytic clocks
leave the `q=23` petal/covering pair mixed, the certificate teeth schedule the
broader residual ledger, the residual tooth atlas turns that ledger into
owner-strip descents, the primitive-period scheduler separates the coarse
Q-witness/covering deck, the path-lift clarifies the boundary topology layer,
and the capacitor cut identifies the boundary-topology tooth that discharges
this displayed two-plate collision.

The audit recomputes only four packets and attaches existing sidecars:

```text
04-computation/lrc14_residual_capacitor_flow_codex_s196.py
05-knowledge/results/lrc14_residual_capacitor_flow_codex_s196.out
```

The two capacitors cross-cut:

```text
M_q_petal_covering_capacitor
  two drop(10,13)->add(20,26)  BOUNDARY-PETAL-SPORADIC
  two drop(8,12)->add(16,24)   COVERING-MOMENT
  first cut: word_plus_boundary_topology
  exit class: nested_refinement

boundary_topology_k33_covering_capacitor
  two drop(12,13)->add(26,36)  K33-STATE-LIFT
  single swap 12->72           COVERING-MOMENT
  first cut: word_plus_M_q
  exit class: cross_handoff
```

So exact scale is not the universal next tooth, and topology is not the
universal next tooth.  They are cheap cuts in a cut lattice: each discharges the
capacitor the other quotient leaves charged.  Closed arc topology,
safe-component stalks, fusion signature, packet labels, and route labels split
both pairs.

The theorem target is a residual capacitor lemma:

```text
Inside any protected LRC14 status fiber, every mixed strict-open route pair
must be cut by exact scale, boundary/arc topology, safe-component stalk,
fusion sidecar, packet labels, owner strip, nested refinement, cross-handoff,
same-tile boundary, or named F7/THM-572 debt.
```

Tournament Analysis used residual-cut carriers as vertices, not runners.  The
observable was capacitor separation count, route/status purity, retained packet
labels, topology, stalk data, magnitude, non-circularity, proof cost, and fiber
count.  The fingerprint was transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=(1,1,1,1,1,1,1,1,1,1)
hamiltonian_path_count=1
```

Next pull: run this capacitor-cut audit over all `15` HYP-3028 coarse ET+unit
mixed-route fibers, not just the two HYP-3027 displayed pairs, and add
`residual_capacitor_id`, `first_cut_stage`, and `zeta_exit_class` to the
cached HYP-2963 packet ledger.
