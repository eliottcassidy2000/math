---
id: THM-738
title: NEAR-AP THREE-SLOT CLOSURE — every 13-speed family with AT LEAST 10 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 10-element body E ⊆ {1,…,14} (all C(14,10)=1001) and all c<a<b not in E, {E,c,a,b} is lonely. Proof = THM-735's Bonferroni tree body-by-body: leg J3 (one inequality, all three slots ≥ V₁(E), clustering-immune), leg J2 (per-c exact bodies, j=2), leg J1 (per-(c,a) exact bodies, THM-732(iii) tail), bottom = exact-ℚ sweeps of covering triples (enumerated as lcm-multiples of the per-(c,a) missing-divisor set) + THM-366 for non-covering. Strictly contains THM-734 (two-slot, 364 bodies) and THM-735(iii) (body {1..10})
status: PROVED (kind-pasteur-2026-07-13-S128 cont.4; Lean transcription open). Run complete: 1001 bodies in 1659 s (27.7 min). Totals: 6,056,745 exact per-(c,a) bodies (leg J1), 4,677,712 bottom exact-ℚ sweeps (dominated by the 35 "flood" bodies with Q(E)=∅ — bodies containing {8,…,14}-covering sets where EVERY triple is covering; worst body {1,2,3,8,9,10,11,12,13,14}: V₁=368, 62.7 s). Per-body missing sets Q(E) ⊆ {2,…,14} computed exactly (general bodies DO miss small q — e.g. Q={4,8,12,14} occurred — handled by the lcm-multiples bottom enumeration). REGRESSION: body {1..10} reproduces THM-735(iii) exactly (V₁=154, 143 J2 bodies, 7537 J1 pairs, 27 sweeps). RESULT: ZERO tights among the 4.68M swept covering triples, ZERO covering L=0 — every swept family has L > 0 strictly; non-covering triples dispatched by THM-366. The theorem contains THM-734 (≥11 in {1..14}) and THM-735(iii) (body {1..10}) as strict special cases
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

## Evidence log (all done — script lrc14_thm738_1001_body_tree_kps_S128c4.py, output in 05-knowledge/results/)

- [x] 1001-body run: 1659 s; V₁ ∈ [~30, 368]; 6,056,745 J1 bodies; 4,677,712 bottom sweeps
      (35 flood bodies with Q=∅ dominate: up to 183k sweeps each; all clean)
- [x] tight census: NONE among swept covering triples (the known tights {1..13}, {1..11,13,24} are
      non-covering, correctly never enumerated — dispatched by THM-366); ZERO covering L=0
- [x] verdict: ESTABLISHED. Regression vs THM-735(iii) exact.

## Scoping the next rung (j=4)

Sampled 9-element bodies give j=4 thresholds V₁ ≈ 180–230 (e.g. {1..9}: r=20, m=0.181, V₁=209).
C(14,9) = 2002 bodies with one more tree level ⟹ estimated 6–30 h — an overnight/cron target.
j=5,6 will likely need the per-peel exact-disc tightening (THM-735's CS form) to keep trees tame.
