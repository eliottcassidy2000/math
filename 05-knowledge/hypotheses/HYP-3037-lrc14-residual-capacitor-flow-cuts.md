---
id: HYP-3037
title: LRC14 residual capacitor flow cuts
status: EVIDENCE / finite residual-cut scout and theorem target; not a proof
source: codex-2026-06-26-S196b
tangent: T1118
related:
  - HYP-3036
  - HYP-3035
  - HYP-3034
  - HYP-3033
  - HYP-3032
  - HYP-3031
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3024
  - HYP-3023
  - HYP-2992
  - HYP-2991
  - HYP-2990
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3037: Residual Capacitor Flow Cuts

## Claim

After the status gates have removed boundary/open ambiguity, a mixed open-route
fiber should be treated as a residual capacitor: two plates carry the same
coarse charge under a quotient, and the proof obligation is to find the first
named cut that separates them.

The two residual mixed-route pairs from the HYP-3027 repair ladder cross-cut
each other.  This is downstream of HYP-3036, HYP-3035, HYP-3034, HYP-3033, and
HYP-3032: the primitive-period scheduler separates the coarse residual
Q-witness/covering deck, the residual tooth atlas turns the full coarse ledger
into first owner-strip descents, the arc-boundary path-lift pass makes the
AP/GW closed-H1 cycle owner-essential, the residual certificate-teeth pass
supplies the broad topology-bucket-plus-unit-scale scheduler for all 15 coarse
residual fibers, and the analytic sieve-clock bridge keeps the `q=23`
petal/covering pair mixed as a squarefree-clean residual.  The capacitor cut
here explains which non-analytic sidecar cuts the displayed two-plate
collisions.

```text
M_q_petal_covering_capacitor
  two drop(10, 13)->add(20, 26)  BOUNDARY-PETAL-SPORADIC
  two drop(8, 12)->add(16, 24)   COVERING-MOMENT
  first cut: word_plus_boundary_topology
  exit class: nested_refinement

boundary_topology_k33_covering_capacitor
  two drop(12, 13)->add(26, 36)  K33-STATE-LIFT
  single swap 12->72             COVERING-MOMENT
  first cut: word_plus_M_q
  exit class: cross_handoff
```

Thus exact scale is not the universal next tooth, and topology is not the
universal next tooth.  They form a min-cut pair: each cuts the capacitor that
the other quotient leaves charged.

## Computation

Script:

```text
04-computation/lrc14_residual_capacitor_flow_codex_s196.py
```

Output:

```text
05-knowledge/results/lrc14_residual_capacitor_flow_codex_s196.out
```

The script recomputes only the four target packets from the labelled
counterexample classifier, attaches the existing arc-Cech topology,
safe-component stalk, and fusion sidecars, and evaluates these stages:

```text
raw_residual_pair
automatic_word
word_plus_M_q
word_plus_boundary_topology
closed_arc_topology
safe_component_owner_stalk
safe_component_exact_stalk
fusion_signature
packet_label_sink
route_label_sink
```

On the four target packets:

```text
raw_residual_pair              fibers=2 pair_splits=0 mixed_route=2
automatic_word                 fibers=2 pair_splits=0 mixed_route=2
word_plus_M_q                  fibers=3 pair_splits=1 mixed_route=1
word_plus_boundary_topology    fibers=3 pair_splits=1 mixed_route=1
closed_arc_topology            fibers=4 pair_splits=2 mixed_route=0
safe_component_owner_stalk     fibers=4 pair_splits=2 mixed_route=0
safe_component_exact_stalk     fibers=4 pair_splits=2 mixed_route=0
fusion_signature               fibers=4 pair_splits=2 mixed_route=0
packet_label_sink              fibers=4 pair_splits=2 mixed_route=0
route_label_sink               fibers=4 pair_splits=2 mixed_route=0
```

The petal/covering pair shares automatic word and exact `M+q`; it is split
when boundary topology includes the strict safe-mass bucket.  The K33/covering
pair shares automatic word and coarse boundary topology; it is split already
by exact `M+q`.

## Theorem Target

Residual capacitor lemma:

```text
Fix a protected LRC14 status fiber after automatic/residue and coarse status
gates.  If a fiber contains two strict-open packets with different theorem
routes, then the first nonzero cut among exact scale, q threshold, boundary
topology, closed arc topology, safe-component stalk, fusion sidecar, and packet
label must classify the pair as owner_strip, nested_refinement, cross_handoff,
same_tile_boundary, or named residual debt.
```

In max-flow/min-cut language, the proof ladder is not a scalar ranking.  Status
gates remove boundary/open capacity; route capacitors are the remaining finite
cut obligations.  A putative counterexample family would have to survive every
cut without entering AP/Goddyn-Wong boundary equality, C27/petal descent,
covering moment descent, K33/THM-572 cross-handoff, Fejer/Haar/Ramanujan dual
annihilation, or named F7 debt.

## Tournament Analysis

Vertices are residual-cut carriers, not runners or arcs.

Pairwise observable:

```text
capacitor separation count,
route/status purity,
retained packet labels,
topology,
stalk data,
magnitude,
non-circularity,
proof cost,
fiber count
```

Switch/gauge:

orient `A -> B` when `A` has lexicographically larger observable vector.  The
tie Hamiltonian path is the stage order listed above.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=(1,1,1,1,1,1,1,1,1,1)
hamiltonian_path_count=1
```

High-retention path:

```text
packet_label_sink
> fusion_signature
> route_label_sink
> closed_arc_topology
> safe_component_exact_stalk
> safe_component_owner_stalk
> word_plus_boundary_topology
> word_plus_M_q
> raw_residual_pair
> automatic_word
```

The ranking is theorem-facing: packet labels and fusion signatures retain most
information, but exact `M+q` and boundary topology are the cheaper min-cuts for
the two known capacitors.

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, endpoint walls,
wall-crossing events, residue words, cover arcs, safe-component bars, dyadic
Haar tiles, packet rows, route labels, residual pairs, and proof obligations.

The chosen vertices are residual-cut carriers because the preserved predicate
is route purity after boundary/open status is already protected.  This choice
destroys raw runner identity and exact route labels until a sidecar separates
the capacitor or declares named residual debt.

The challenged assumption is that the next proof tooth should be globally
ordered.  The audit says the right abstraction is a cut lattice: exact scale
and boundary topology are incomparable cheap cuts, and each discharges a
different residual capacitor.

## Next Pull

1. Run the capacitor-cut audit over all `15` HYP-3028 coarse ET+unit
   mixed-route fibers, not just the two HYP-3027 displayed pairs.
2. Add `residual_capacitor_id`, `first_cut_stage`, and `zeta_exit_class` to a
   cached HYP-2963 packet sidecar.
3. Prove family lemmas for the two observed exit classes:
   `nested_refinement` for the petal/covering exact-scale collision and
   `cross_handoff` for the K33/covering topology collision.
