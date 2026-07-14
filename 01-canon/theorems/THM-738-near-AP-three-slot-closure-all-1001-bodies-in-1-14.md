---
id: THM-738
title: NEAR-AP THREE-SLOT CLOSURE — every 13-speed family with AT LEAST 10 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 10-element body E ⊆ {1,…,14} (all C(14,10)=1001) and all c<a<b not in E, {E,c,a,b} is lonely. Proof = THM-735's Bonferroni tree body-by-body: leg J3 (one inequality, all three slots ≥ V₁(E), clustering-immune), leg J2 (per-c exact bodies, j=2), leg J1 (per-(c,a) exact bodies, THM-732(iii) tail), bottom = exact-ℚ sweeps of covering triples (enumerated as lcm-multiples of the per-(c,a) missing-divisor set) + THM-366 for non-covering. Strictly contains THM-734 (two-slot, 364 bodies) and THM-735(iii) (body {1..10})
status: CLAIMED (kind-pasteur-2026-07-13-S128 cont.4) — the 1001-body run is THIS SESSION's computation (script lrc14_thm738_1001_body_tree_kps_S128c4.py); upgrades to PROVED when all bodies close clean. NOTE the general-body covering condition: a 10-element body E ⊆ {1..14} need not cover q ≤ 10 (e.g. {1..6,11..14} misses q=8,9,10 multiples...), so the per-body missing set Q(E) ⊆ {2,…,14} is computed exactly and the bottom covering test is Q-general (unlike THM-735(iii) where body {1..10} covers q ≤ 10 automatically)
source: kind-pasteur-2026-07-13-S128 (cont.4)
depends_on:
  - THM-735   # the simultaneous multi-peel lemma (j=3,2,1 legs) — this is its 1001-fold body sweep
  - THM-731 / THM-732   # per-peel covariance + crude disc
  - THM-366   # non-covering dispatch
related:
  - THM-734 (the j=2 analogue, 364 bodies), THM-736 (mac-mini far-peel Farey), THM-737 (opus pack-clock — closes the coherent-pack sector complementarily)
  - HYP-6540 (calibration), klein HYP-6505 (isolation wall — bypassed by simultaneity), MISTAKE-122 (j≤6), MISTAKE-141 (all thresholds exact)
---

# THM-738 — the near-AP three-slot closure (1001 bodies)

Statement and method as in the title; per body E: m_E > 0 exact (automatic from LRC(≤13) but computed),
V₁(E) = minimal integer with 3/V₁ < 4·m_E/((99/70)·r_E); legs J2/J1 with exact per-c and per-(c,a)
bodies; bottom triples enumerated via the exact missing-divisor lcm (or fully when the pair already
covers Q(E)), swept in exact ℚ.

## Evidence log

- [ ] 1001-body run: per-body V₁ distribution, level counts, bottom sweeps, timing
- [ ] tight census (expect only non-covering tights, never swept; any covering L=0 flagged loudly)
- [ ] verdict
