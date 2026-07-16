---
id: THM-897
title: THE LEVEL-5 WALL AND THE TRIPLE BEAT LEMMA — (W) the order-5 Bonferroni wall W₅(m′) = Σ_{k≤5}(−1)^k C(m′,k) 2^k 13^{5−k} is COERCIVE FOR m′ ≤ 15 (dead at 16), and THE FULL 13-COMB SCALE-ONE PROBLEM (m′ = 13, no prefix) sits INSIDE with margin exactly 35035/371293 = 9.44% — and 35035 = 5·7²·11·13: the ladder price schedule completes as 6.5 (union) → 7.5 (tree) → 11.5 (B3) → 15.5 (B5) ⊇ 13; (T1) THE TRIPLE BEAT LEMMA, scale-separated case, PROVED in three lines from THM-815's per-interval discrepancy: |μ(W ∩ D₃) − (2/13)μ(W)| ≤ 22·κ_W/(169·x₃) (κ_W = #components), refereed 4/4 with 5–50× slack; (T2) the same-scale cascade: μ₃/((2/13)μ₂) ∈ [0.93, 3.25] on the battery (consecutive triples ≈ 3.25 ≈ the joint beat; (77,143,169) EXACTLY 1.000) — enhancement bounded small, the Fourier/sub-orbit proof with the triple resonance integer named; (R) THE ADMISSIBILITY DISCOVERY (three exact 13-packet referees, all honestly NEGATIVE with identified causes): packet 1 (residues 1..13, moderate lifts): S₅ 7× enhanced — the 13-MULTIPLE (residue-0 = the deep fibre) aligns with the 1/13 grid; packet 2 (no 13-multiple, larger lifts): S₂ exact to 4 digits but S₄/S₅ still 1.4×/3.7× — ADDITIVE QUADRUPLES (x_i + x_j = x_k + x_l: S321's energy axis) which lift size cannot remove; packet 3 (greedy Sidon): additive resonances killed but S₅ 12× — the construction accidentally created a 7-DILATE cluster (the THM-863 ratio axis): THE LEVEL-5 REGIME IS DELIMITED JOINTLY BY THE ADDITIVE AND MULTIPLICATIVE NON-RESONANCE AXES — the repo's two coordinate systems are exactly the two admissibility conditions; whether any 13-packet at accessible scales satisfies both is the named open question (if empty at small scales, the full-problem certificate is genuinely asymptotic and the prefix recursions remain essential)
status: W PROVED (exact coefficient arithmetic); T1 PROVED + REFEREED (4/4); T2 empirical with the proof route named; R = three exact negative referees with structural diagnoses (the honest boundary of the level-5 schema)
source: opus-2026-07-16-S328 (owner: prove the triple beat lemma and compute the level-5 wall)
depends_on:
  - THM-896 (level 3), THM-815 (the per-interval discrepancy lemma T1 uses)
  - THM-863/864 (the ratio axis), THM-882/S321 (the additive/coherence axis)
related: [THM-769 (the deep fibre — packet 1's diagnosis), the Kikuchi ladder reflection (S326), A139250-adjacent S327 fertility]
verification: 05-knowledge/results/level5_wall_triple_beat_opus_S328.out
---

# THM-897 — the level-5 wall and the triple beat lemma

## (W) The wall

W₅(m′)/13⁵ at m′ = 11..17: 57067, 46013, **35035**, 23309, 9867, −6435,
−26949. Coercive through m′ = 15; the FULL scale-one problem (13 combs, no
prefix) sits inside at margin 35035/371293 ≈ 9.44%, with the numerator
factoring as 5·7²·11·13. The Bonferroni ladder's price schedule for LRC(14):
union 6.5 → Hunter tree 7.5 → B3 11.5 → **B5 15.5**, and 13 < 15.5: the
whole problem is, at the level of equidistributed densities, a level-5
statement.

## (T1) The triple beat lemma (scale-separated) — proved

For W a finite union of κ_W intervals and any comb D₃:
|μ(W ∩ D₃) − (2/13)μ(W)| ≤ 22κ_W/(169x₃). *Proof:* THM-815's discrepancy
lemma (6) applied per component, both sides. ∎ Applied to W = D₁ ∩ D₂ ∩ E
(κ_W ≤ 2(x₁+x₂)/13 + κ_E): the triple overlap factors through the pair
overlap with conditional density 2/13, sharply when x₃ ≫ x₁ + x₂. Referee:
4/4 scale-ordered triples, errors 5–50× under the cap.

## (T2) The same-scale cascade (empirical, proof route named)

μ₃/((2/13)μ₂) on the battery: 0.93–3.25; the consecutive triple
(150,151,152) peaks at 3.254 (the joint beat), and (77,143,169) sits at
1.000 exactly. The proof is THM-864's sub-orbit engine run on the
pair-window orbit against the third comb, with resonance parameter the
TRIPLE RESONANCE INTEGER |h·x₃·c − p·A| (c the pair's Bezout inverse) —
one session at THM-864 rigor.

## (R) The admissibility boundary (three honest negatives)

Exact S₁..S₅ on three 13-comb packets:
1. residues 1..13, lifts ≤ 25: BONF5 = −0.319 (S₅ 7× equid) — the
   13-multiple (deep fibre, THM-769) aligns with the grid;
2. no 13-multiple, lifts 20–60: BONF5 = −0.095 (S₂ EXACT to 4 digits;
   S₄/S₅ at 1.4×/3.7×) — ADDITIVE QUADRUPLES persist at any lift size;
3. greedy Sidon (additive-clean): BONF5 = −0.642 (S₅ 12×) — the
   construction created a 7-dilate cluster: MULTIPLICATIVE resonance.

> **The level-5 equidistribution regime is delimited jointly by the additive
> axis (S321 energy: no repeated pair-sums) and the multiplicative axis
> (THM-863: no small-q ratios/dilates).** The two coordinate systems this
> program built are exactly the two admissibility conditions of its top
> rung. Named open: does a 13-packet at accessible scale satisfy both? A
> both-clean construction (Sidon + pairwise ratio-non-resonant) is the next
> computation; if the joint regime is empty at small scales, the
> full-problem certificate is genuinely asymptotic and the prefix
> decompositions (THM-863's trichotomy, the radius recursions) remain the
> effective route — the ladder tells us WHY: each Bonferroni level trades
> comb count for purer non-resonance demands.

## Named next

- The both-clean 13-packet search (Sidon × ratio-non-resonant × no
  13-multiple): one targeted construction session.
- T2 referee-grade (the triple resonance integer sub-orbit proof).
- The quadruple LOWER bound (S₄ enters positively at level 5: the
  four-fold beat FLOOR — the harder direction, as Hunter's tree was at
  level 2: perhaps a hypertree selection at level 4).
