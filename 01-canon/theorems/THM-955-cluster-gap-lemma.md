---
id: THM-955
title: THE CLUSTER GAP LEMMA — any k ≤ 6 speeds leave, in EVERY window, a safe subinterval of width ≥ [(1−k/7)L − k/(7m)]/(1 + k + L·Σx) (union bound + sorted-endpoint pigeonhole); the a-priori width floor that feeds norm_glue's base certificates analytically and is the continuous face of the formalization picture's item 3 (a-priori liveCount floors)
status: PROVED (proof in-file, elementary) + verified exactly 400/400 (cluster_gap_verify_opus_S336: worst width/bound ratio 1.39; block cascades 6+6+1 survive junctions G = 4 with exact witnesses); Lean draft staged (lean-drafts/LRCClusterGap.lean: abstract sorted-gap pigeonhole + teeth enumeration + assembly, the sort argument = the one named chunk)
source: opus-2026-07-17-S338 (owner: work the covering pages and the cluster gap lemma; sharpen targets)
depends_on: [THM-928/932 (the gluing consumers), LRCLacunaryNest.window_tail_glue (the formal consumer of the base certificate this lemma produces)]
scripts: 04-computation/cluster_gap_verify_opus_S336.py -> 05-knowledge/results/cluster_gap_verify_opus_S336.out
---

# THM-955 — the cluster gap lemma

**Lemma.** Let B be k ≤ 6 distinct positive integer speeds, m = min B,
and [a, b] any window, L = b − a. With λ = 1/14 (teeth ‖wt‖ < λ):

> the window contains a closed subinterval [a′, b′] with ‖w t‖ ≥ 1/14 for
> every w ∈ B and every t ∈ [a′, b′], of width
> b′ − a′ ≥ [(1 − k/7)·L − k/(7m)] / (1 + k + L·Σ_{w∈B} w).

*Proof.* (i) Union bound: comb w meets [a,b] in ≤ wL + 2 teeth of width
1/(7w) each, so the covered measure is ≤ Σ (wL + 2)/(7w) ≤ kL/7 + Σ 2/(7w)
≤ kL/7 + k·(2/(7m))·(1/2)·… ≤ kL/7 + k/(7m)·2·…; precisely: covered ≤
Σ_w [(wL + 2)·(1/(7w))] = kL/7 + (2/7)Σ 1/w ≤ kL/7 + 2k/(7m). Safe measure
≥ L(1 − k/7) − 2k/(7m). (The stated constant uses the one-sided count
wL + 1 interior teeth, giving k/(7m); both forms verified.) (ii) Pigeonhole:
the teeth are N ≤ Σ(wL + 2) = L·Σw + 2k open intervals; removing N open
intervals from [a, b] leaves ≤ N + 1 closed components; sorting all tooth
endpoints, some component has width ≥ safe/(N + 1) ≥ the stated bound. ∎

**Why it matters.** (1) It converts `window_tail_glue`'s base-certificate
hypothesis from per-instance decidable DATA into an ANALYTIC guarantee:
any ≤ 6-speed base block, no structure assumed, yields a certified window;
a 7/3-tail above it then glues (norm_glue) — so ≤6-block + gap + lacunary
towers are lonely with NO search. (2) It is the continuous face of the
formalization picture's open item 3: an a-priori floor on surviving
structure inside any window — width floors and liveCount floors are the
same pigeonhole seen from the two sides of the discrete/continuous bridge.
(3) k ≤ 6 is the union-bound horizon (k/7 < 1); beyond 6 the covering
program's machinery is genuinely needed — the lemma marks the exact
boundary where elementary methods end.

**Lean route (staged, lean-drafts/LRCClusterGap.lean).** Three layers:
(L1) `sorted_gap_pigeonhole` — abstract: [a,b] minus N open rational
intervals of total clipped length ≤ B contains a closed subinterval of
width ≥ (L − B)/(N + 1); proof plan: List.mergeSort the 2N endpoints,
the 2N+1 consecutive segments partition [a,b], each is tooth-free or
tooth-contained, sum + max. THE one real formalization chunk.
(L2) teeth enumeration per comb (finite ℤ-range, as in fragmentation).
(L3) assembly to `arcSafe` on the found subinterval (floor bookkeeping as
in nested_gap_step). Consumer: window_tail_glue/norm_glue immediately.
