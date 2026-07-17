---
id: LEM-042
title: THE PAIR-OVERLAP LAW (the full-period arithmetic core of the 7-wall pair-crumb — the fleet's one-item-wide open, extended). (A) THE EXACT FORMULA: for danger sets D_v = {t : ‖vt‖ < 1/14}, speeds a ≤ b, g = gcd(a,b), δ = g/(ab), W = (a+b)/(14ab), w_min = 1/(7b): μ(D_a ∩ D_b) = g·Σ_{j∈Z} min((W − |j|δ)₊, w_min) — arcs overlap iff centers are within W; center differences take exactly the values jδ, each g times per period; the overlap at distance u is the trapezoid min((W−u)₊, w_min). Referee: exact (Fractions) against brute sweeps on 204 pairs. (B) THE EQUIDISTRIBUTION VALUE: the trapezoid's normalized integral is EXACTLY ab·w_a·w_b = 1/49 for every pair — the independence value; the law is a discrete sampling of a shape with mean exactly 1/49. (C) THE NEAR-EQUAL FLOOR AND THE FAILURE MAP: consecutive-reduced pairs (b/g = a/g + 1): μ = 1/(7(a′+1)) for a′ ≤ 6, EXACTLY 1/49 at a′ ∈ {6,7}, and ≥ 1/49 for all a′ ≤ 200 (named finite-induction for all a′) — equal-gap blocks always produce consecutive-reduced pairs, so THE WALL-RELEVANT CASE HOLDS CLEANLY. General reduced pairs: the floor FAILS on an arithmetic (not ratio-cone!) set — (13,15) fails at ratio 1.15 (9/455 ≈ 0.01978), (5,8) = 1/56, (1,12) = 1/84, all failures SHALLOW at moderate sums except the decaying a′ = 1 column; the formula computes the exact domain. (D) THE WALL CROSSED on referee blocks: for c = 7 near-equal blocks the union bound is DEAD (1 − 7/7 = 0) but the hunter ledger with EXACT pair credits gives good ≥ Σ_{consecutive} μ ≈ 6/49 ≈ 0.123 > 0 — verified ≤ the true good measure (≈ 0.298) on four blocks incl. gcd-scaled. (E) THE WINDOW FACE COLLAPSES (honest hand-off): min over window positions of μ(D_a ∩ D_b ∩ I)/|I| = 0 even at |I| = 1/2 for (50,51) — the intersection support CLUSTERS (the difference phase cycles once per period), so NO uniform per-window pair floor exists for near-equal pairs: the pair credit is a full-period/positioned phenomenon, and per-cell obligations must aggregate cells or position windows (mac-mini/klein: this constrains HYP-3874's shape)
status: PROVED ((A) 5-line arc counting; (B) integral identity; (D) hunter + exact values) + REFEREED EXACT (204 formula-vs-brute pairs; consecutive law to a′ = 200; four wall blocks; window collapse measured); (C)'s all-a′ consecutive floor: verified + named finite-induction
source: boxeph-2026-07-17-S69 (owner directive: finish the LRC 14 proof; integrate and extend incoming work — opus THM-956's one-item list, mac-mini HYP-3874, klein HYP-4021 hunter ledger)
depends_on: [klein's path_hunter_add_le (the ledger consuming these credits), independent of all my LEM-030..041 frame line]
related: [codex pair law (S33 usage: the D-matrix deviation — the box-hit face of the same object), THM-884(E)]
script: 04-computation/lrc14_pair_overlap_law_boxeph_S69.py -> 05-knowledge/results/lrc14_pair_overlap_law_boxeph_S69.out
---

# LEM-042 — the pair-overlap law

The 7-wall's pair-crumb, full-period face: an exact finite trapezoid-sum
whose mean is exactly 1/49, with the clean floor precisely where the wall
needs it (consecutive-reduced pairs = equal-gap blocks) and a computable
failure map everywhere else. The union bound dies at c = 7; the hunter
ledger with these exact credits walks through at ≈ 6/49. The window face
is the honest boundary: near-equal intersections cluster, so the credit
cannot be localized uniformly — it must be harvested globally or at
positioned windows. Lean shape: (A) is finite rational arithmetic per pair;
(B) is one algebraic identity; the consecutive-floor induction is
decide-adjacent.

## Evidence log
- [x] formula exact ×204; consecutive law to 200; failure map; wall ×4
- [x] window collapse measured (the constraint on per-cell obligations)
- [ ] named: the all-a′ consecutive-floor induction (finite forms per n = ⌊x⌋)
- [x] named item ADVANCED (LEM-043(E), S70): the difference-comb law + beat census (mass in d beats, uniform in the gcd-7 case; single beat at d = 1; off-beat windows carry zero — 3/3 demo); the coverage synthesis assigns regimes to tools
