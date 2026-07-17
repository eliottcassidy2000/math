---
id: LEM-043
title: THE CONE FLOOR AND THE c = 7 WALL DICHOTOMY (the pair-crumb assembly layer, full period). (A) THE DEFECT BOUNDS: 1/49 − g/(7b) ≤ μ(D_a ∩ D_b) ≤ 1/49 + g/(7b) — the overlap is the δ-sampling of an even unimodal trapezoid with integral exactly w_a·w_b, and sampling a nonincreasing tail loses/gains at most one cell (one line each side). (B) THE CONE FLOOR, exact: c₀(7/3) := min{μ over reduced pairs with b′/a′ ≤ 7/3} = 1/63, attained at (4,9) — PROVED by exhaustive finite check (b′ ≤ 700, integer arithmetic 7gΣ vs 2ab) + the (A)-tail (1/49 − 1/4900 > 1/63 for b′ > 700). Also c₀(13) = 1/91 at (1,13) for the all-comparable ≤ 13× convention. (C) THE c = 7 WALL DICHOTOMY (full period): every 7-block either has all consecutive ratios ≥ 7/3 — the LACUNARY branch, lonely universally by opus-S336 — or has some consecutive ratio < 7/3, and then klein's path_hunter_add_le with the ONE cone-floor credit beats the union bound's death (7 × 1/7 = 1): good = 1 − μ(∪ D_i) ≥ c₀ = 1/63 > 0. Referee: structured blocks both branches, credits ≥ 1/63 and ≤ true good throughout. (D) THE c = 8 BOUNDARY, sharper than expected: pair credits need Σ > (c−7)/7 = 1/7; the uniform floor 7·(1/63) < 1/7 FAILS, but consecutive-type blocks have every credit STRICTLY above 1/49 (a′ ≥ 8), and the demonstrated block [100..107] crosses with Σ = 0.14293 > 1/7 — c ≥ 8 is per-family decidable (exact credits), uniformly NAMED not claimed. (E) THE DIFFERENCE-COMB LAW: D_a ∩ D_b ⊆ {t : ‖(b−a)t‖ < 1/7} (triangle inequality); the mass distributes over d = b − a beats (exactly uniform when the reduced pair has period structure, e.g. gcd-7 pairs; single beat at d = 1); off-beat sweep windows of width 2/a carry ZERO in the clustered regime (2464, 2465 demo) — positioning is essential exactly where cite_cluster7_lonely already operates by the C0-citation route: THE COVERAGE SYNTHESIS: clustered (ratio ≤ 1 + 1/1232) → kps's closer; loose-comparable (some ratio < 7/3) → this dichotomy; lacunary → opus's branch; the middle → exact credits (decide)
status: PROVED ((A) one-line sampling bounds; (B) finite check + tail; (C) hunter + (B); (E) subset law one line) + REFEREED EXACT (500 random defect pairs; exhaustive cone scans; five 7-blocks; the c = 8 instance; beat censuses; the clustered beat-miss 3/3)
source: boxeph-2026-07-17-S70 (owner directive: work the remaining proof surface; integrate incoming — opus lacunary branch + THM-956, klein hunter ledger, kps clustered closer, LEM-042)
depends_on: [LEM-042 (the exact pair formula), klein path_hunter_add_le (kernel-pure), opus S336 lacunary branch, kps cite_cluster7_lonely (the clustered face)]
script: 04-computation/lrc14_cone_floor_wall_boxeph_S70.py -> 05-knowledge/results/lrc14_cone_floor_wall_boxeph_S70.out
---

# LEM-043 — the cone floor; the wall dichotomy

At c = 7 the union budget is exactly zero, so ONE positive pair credit
crosses the wall. The credit exists uniformly on the lacunary-complement
cone: c₀(7/3) = 1/63, exact. Every 7-block therefore falls to one of three
PROVED tools — lacunary (opus), one-credit hunter (here), clustered closer
(kps) — and the coverage synthesis names which tool owns which regime. The
c = 8 boundary is decidable per family and crossable in the consecutive
regime; its uniform version is the named residue.

## Evidence log
- [x] defect bounds ±g/(7b) (500 exact pairs); c₀(7/3) = 1/63, c₀(13) = 1/91
- [x] dichotomy on five blocks; c = 8 instance crosses (0.14293 > 1/7)
- [x] difference-comb subset law; beat censuses; clustered beat-miss 3/3
- [ ] named: uniform c ≥ 8 (needs strict-excess accounting or higher-order credits)
- [ ] named: Lean rendering — (A)/(B) are integer decides; (C) composes with
      path_hunter_add_le (already kernel-pure)
