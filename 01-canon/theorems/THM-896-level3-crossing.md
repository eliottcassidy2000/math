---
id: THM-896  # renumbered from 895 (helmholtz-perron first-pushed)
title: THE LEVEL-3 CROSSING — the m′ = 8 wall falls to order-3 Bonferroni: the truncation uncovered ≥ μE − S₁ + S₂ − S₃ is a valid LOWER bound (classical Bonferroni, odd truncation — NO tree needed, unlike level 2 where the full S₂ has the wrong sign and Hunter's tree was required); with equidistributed densities (2/13, 4/169, 8/2197) the coefficient wall 2197 − 338m′ + 52C(m′,2) − 8C(m′,3) is COERCIVE FOR m′ ≤ 11 (margins /2197: 643, 501, 355, 197, 19 — another knife-edge — at m′ = 7..11; dead at 12): THE LEVEL LADDER'S PRICE SCHEDULE FOR LRC(14) IS 6.5 (union), 7.5 (Hunter tree), 11.5 (Bonferroni-3); REFEREED: on three real radius-8 packets (8 lifted residues over prefix {1,2,3,4}) the EXACT bound is positive (+0.146, +0.149, +0.126) and within 8–16% of the true uncovered mass — the level-3 certificate is not just coercive but TIGHT; (T) the triple sawtooth: generic triples sit at (2/13)³ = 0.003641 (exact hit at (77,143,169)) with RESONANCE ENHANCEMENTS (near-equal (150,151,152): 3.3×; near-1:2:3: 2×) — the error budget's danger is enhancement-side on S₃ (upper bounds needed: the triple beat lemma, named); (F) THE FERTILITY FIELD (the toothpick probe, n = 5 → 6): the metagraph growth automaton is SYMMETRY-SUPPRESSED, not boundary-localized — fertility ≈ the Burnside shadow 32/|Aut| (regular/Paley-5: 7 ≈ 32/5, the minimum; transitive: 21; max 24 at mid-level x = 16): the anti-toothpick verdict, with the extension bipartite graph 12 × 56 carrying 215 incidences (mean parent-count 3.84 of 6)
status: A PROVED (classical Bonferroni direction + exact coefficient arithmetic) + REFEREED EXACT (3/3 radius-8 packets coercive, tightness 8–16%); T empirical (4 exact triples; the law's shape confirmed, the triple beat lemma named open); F computed exact (12-class extension census)
source: opus-2026-07-16-S327 (owner: compute the level-3 tensor and cross the m′ = 8 wall; the toothpick/automaton probe on A139250/arXiv:1004.3036)
depends_on:
  - THM-856/863/864 (the level-2 machinery this ladder extends)
  - the Kikuchi-ladder reflection (S326: this is its first new rung built)
related: [THM-855 F6 (closure at order 2 — the moment ladder's other face), A139250 (the toothpick automaton), A000568/A000570 (the paired-diagram reading)]
verification: 05-knowledge/results/level3_crossing_toothpick_opus_S327.out
---

# THM-896 — the level-3 crossing

## (A) The wall falls

Bonferroni: for events A₁..A_m, the inclusion–exclusion truncation ending at
−S₃ LOWER-bounds the uncovered measure — no spanning-tree selection needed
(the parity of the truncation supplies the sign that level 2 lacked; Hunter's
tree was the level-2 workaround, and the ladder alternates: full sums work at
odd truncations). With densities at equidistribution:

> uncovered/μE ≥ 1 − m′(2/13) + C(m′,2)(4/169) − C(m′,3)(8/2197),
> i.e. (2197 − 338m′ + 52C(m′,2) − 8C(m′,3))/2197:
> **643, 501, 355, 197, 19 at m′ = 7, 8, 9, 10, 11; −187 at 12.**

The level ladder's price schedule: union 6.5 → tree 7.5 → **Bonferroni-3
11.5** (the m′ = 11 margin 19/2197 is the new knife-edge). Referee, exact
(prefix {1,2,3,4}, μE = 13/21, three random radius-8 packets, S₁/S₂/S₃ in
exact ℚ): bounds +0.14613, +0.14850, +0.12563 vs actual uncovered 0.15987,
0.16161, 0.14919 — **coercive AND tight (within 8–16%)**. The m′ = 8 wall is
crossed with exact certificates; radii 8–11 are asymptotically closed by the
same schema once the triple error budget is proved.

## (T) The triple sawtooth (the next error budget)

ρ₃ at generic triples hits (2/13)³ = 0.003641 exactly-ish ((77,143,169):
0.003641); resonances ENHANCE: (150,151,152): 0.0118 (3.3×), (99,199,300):
0.0072 (2×). Since S₃ enters NEGATIVELY, enhancements hurt: the needed lemma
is the TRIPLE BEAT UPPER BOUND (the THM-864 analogue at triples: overlap
localizes on the joint beat; enhancement bounded by the pairwise data) —
named as the level-3 quantification, with the same sub-orbit machinery
expected to work.

## (F) The fertility field (the toothpick probe)

One-vertex extensions n = 5 → 6, per class: fertility (distinct 6-classes
reached) = 7 (regular/Paley-5, x = 0 — the MINIMUM), 13–23 (x = 8..32
levels), 24 (max, at x = 16), 21 (transitive). The pattern is the **Burnside
shadow ≈ 32/|Aut|**: symmetry suppresses growth; the maximum sits mid-level.
**The metagraph's growth automaton is symmetry-suppressed, not
boundary-localized — the anti-toothpick verdict** (A139250's growth lives at
exposed ends; the metagraph's lives at asymmetric middles). The owner's
paired-diagram reading stands upstream: A000568 is even for n ≥ 3 because
reversal acts freely off SC — the metagraph IS a diagram of reversal-pair
"toothpicks" with SC classes as the unpaired centers; the growth rule that
generates 2, 4, 12, 56 is the extension bipartite graph (215 incidences at
5 → 6, mean parent-count 3.84), now measured.

## Named next

- The triple beat lemma (the level-3 THM-864): sub-orbit APs at the joint
  relation; enhancement caps.
- Level-5 (ending −S₅): the next odd rung's wall (the quartic/quintic
  coefficient arithmetic — does it reach m′ = 13, the whole problem?).
- The fertility field at 6 → 7 (does the Burnside shadow law persist? the
  exact fertility ≈ f(|Aut|, level) regression).
