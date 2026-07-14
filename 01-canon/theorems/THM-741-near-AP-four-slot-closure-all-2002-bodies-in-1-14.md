---
id: THM-741
title: NEAR-AP FOUR-SLOT CLOSURE — every 13-speed family with AT LEAST 9 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 9-element body E ⊆ {1,…,14} (all C(14,9)=2002) and all v₁<v₂<v₃<v₄ not in E, {E,v₁..v₄} is lonely. Proof = the THM-735 Bonferroni tree at j=4: legs J4 (one inequality, all four ≥ V₁(E)) / J3 (per-v₁ exact bodies) / J2 (per-(v₁,v₂)) / J1 (per-(v₁,v₂,v₃) tail) / bottom (exact-ℚ sweeps of covering quadruples via lcm-multiples) — with PROVED P1/P2 LEMMA-SKIPS at every level (subtrees where the next Bonferroni threshold already fires from the parent's exact data close without computing the child body; sound because P1/P2 are one-level bounds off exact data)
status: CLAIMED (kind-pasteur-2026-07-13-S128 cont.5) — the OVERNIGHT run (owner-directed): multiprocessing over bodies, resume-safe per-body JSONL checkpoint (scratchpad; harvested to 05-knowledge/results on completion), probe-calibrated. Regression built in: the subtree (E={1..9}, v₁=10) must reproduce THM-738's body-{1..10} numbers (V=154, 143, 7537, 27) exactly. Upgrades to PROVED when all 2002 bodies close clean
source: kind-pasteur-2026-07-13-S128 (cont.5)
depends_on:
  - THM-735   # the simultaneous multi-peel lemma (j=4,3,2,1 legs) + P1/P2 peel lemmas (THM-733)
  - THM-731 / THM-732 / THM-366
related:
  - THM-734 (j=2, 364 bodies), THM-738 (j=3, 1001 bodies) — this is rung j=4 of the same ladder
  - THM-737/739/740 (opus: pack-clock / cluster-clock / two-cluster — tiling the j≥7 seam from the coherent side)
  - MISTAKE-122 (j≤6), MISTAKE-141 (exact thresholds), HYP-6540 (calibration)
---

# THM-741 — the near-AP four-slot closure (2002 bodies, overnight run)

Statement and method as in the title. New at this rung: (a) the general-Q bottom now has THREE
slots available to cover Q(E) before v₄; (b) lemma-skips at every level — at the v-level with the
parent's exact (r,m), the child's data is bounded by P1 (r' ≤ v·m + (15/7)r) and P2
(m' ≥ (6/7)m − 8r/(49v)); if the next-level Bonferroni threshold already fires with these bounds,
the whole subtree closes with zero further computation (one-level applications only — the crude
multi-level composition is known to fail, HYP-6540).

## Evidence log

- [ ] probe calibration (sample bodies, serial timing) + extrapolation, worker count chosen
- [ ] overnight run launched (resume-safe); regression subtree check passes
- [ ] all 2002 bodies clean; tight census; verdict
