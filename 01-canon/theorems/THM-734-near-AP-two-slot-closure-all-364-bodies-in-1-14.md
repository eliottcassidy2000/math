---
id: THM-734
title: NEAR-AP TWO-SLOT CLOSURE — every 13-speed family with AT LEAST 11 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 11-element body E ⊆ {1,…,14} (all C(14,11)=364 of them) and every pair a<b of positive integers not in E, {E,a,b} is 1/14-lonely. Proof = THM-733's method run body-by-body: per-body exact (r_E,m_E) (m_E>0 automatic from settled LRC(≤13): 11 speeds have max-min ≥ 1/12 > 1/14), monotone-certified A₀(E), per-a exact thresholds, finite exact-ℚ box sweeps; every L=0 pair must be non-covering (THM-366) — the box census doubles as a TIGHT-FAMILY CENSUS of the near-AP two-slot region
status: PROVED (kind-pasteur-2026-07-13-S128 cont.2; Lean transcription open). Sweep complete: 364 bodies, 245,994 box pairs swept in exact ℚ, 58.5 s total. Every body: m_E > 0 (computed exactly; min m = 313/9702 ≈ 0.0323 at {1,2,3,4,5,7,8,9,11,12,13}), A₀ certified exactly (max A₀ = 455 at {1,2,3,4,8,9,10,11,12,13,14}, which also has the worst box, 6055 pairs). TIGHT CENSUS: exactly TWO L=0 families in the entire region — the AP {1..13} and the Goddyn–Wong doubling {1..11,13,24}, both missing q=14 (THM-366, t=1/14). ZERO covering L=0 pairs. Scope: contains THM-733 (body {1..11}) and ALL of kps cont.58's enumerated multi-killer extremals — k=10's {1..10,13,22,84} (L = 2227/105105 exactly) in body {1..10,13}; k=9's in {1..9,13,14}. k≤8 multi-killer families have ≥3 outliers >14, outside the two-slot format (M ≥ 2/23, loose, per cont.58); unbounded-outlier bodies need the multi-scale induction (P1/P2 compose) — klein-S289's isolation wall applies there, not here
source: kind-pasteur-2026-07-13-S128 (cont.2)
depends_on:
  - THM-731   # certificate chain (klein-S287)
  - THM-732   # exact Bernoulli/Dedekind disc_v + far-element tail
  - THM-733   # the method (peel lemmas P1/P2, monotone f, box sweep) — this theorem is its 364-fold body-by-body extension
  - THM-366   # non-covering dispatch (t=1/q witness) for the tight pairs
related:
  - HYP-6495/6505/6510  # exposure/collapse structure; klein-S289 isolation classifier (my boxes handle the non-isolated regime EXACTLY because bodies are bounded; klein's infinite non-isolated wall {1,90..101}-type has ≤1 speed in {1..14}, outside this class — consistent)
  - kps cont.58/59/68 (the multi-killer skeleton this executes body-by-body)
---

# THM-734 — the near-AP two-slot closure (364 bodies)

**Statement.** Every 13-speed family (distinct positive integers) with at least 11 speeds in
{1,…,14} satisfies LRC(14): it has a time t with all ‖wt‖ ≥ 1/14. Strictly positive measure of
such times (L > 0) except at exactly-tight families, each of which is non-covering and dispatched
by THM-366's witness t = 1/q; the sweep enumerates every tight family in the region as a census.

**Method (per body E, all 364).** Exactly THM-733's three legs with body-dependent data:
1. m_E > 0 for every 11-element body — from LRC(≤13) (settled): max-min ≥ 1/12 > 1/14 gives an
   open neighborhood of good times. So A₀(E) is finite for every body.
2. Uniform leg a ≥ A₀(E): peel lemmas P1/P2 + THM-732(iii), one exact monotone inequality.
3. Per-a leg (a ∉ E, a < A₀): exact (r_a, m_a); box {(a,b): a < b ≤ v₀(a)} swept in exact ℚ.

## Evidence log (all done — script lrc14_thm734_body_sweep_kps_S128c2.py, output in 05-knowledge/results/)

- [x] 364-body sweep complete: 245,994 exact box pairs, 58.5 s. Per-body A₀ ∈ [34, 455]; boxes ∈ [3, 6055].
- [x] tight census: EXACTLY the two predicted families — {1..13} (AP) and {1..11,13,24} (GW doubling),
      both missing q=14 → THM-366. No new tight family anywhere in the region: the census is an
      independent computational confirmation that the Goddyn–Wong list is complete on this domain.
- [x] zero covering L=0 pairs.
- [x] spot theorem-content check: cont.58's k=10 covering-min extremal {1..10,13,22,84} is inside body
      {1..10,13}'s sweep; its exact L = 2227/105105 ≈ 0.0212 > 0.

## Remark — why this region was the right target

The whole program's danger zone (near-AP, where M approaches 1/14) lives here: every known extremal,
near-extremal, and tight family has ≥11 speeds in {1..14}. The theorem says: with two speeds free to
be ANYTHING, the near-AP region cannot produce a counterexample to LRC(14). What escapes the format is
only the k≤8 multi-killer configurations (≥3 large outliers — loose by cont.58's enumeration, M ≥ 2/23)
and genuinely multi-scale families (klein-S289's non-isolated wall) — the compact core is now certified
around, and the surviving analytic content is exactly the multi-scale/non-isolated stratum.
