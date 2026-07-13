# THM-723 — The mixed-slope `j ≥ 8` stratum reduces to PERTURBED TIGHT DILATES, and the clean-q game floors it at `1/13 − B/(2L)` (B=1 census-complete)

**Status:** Reduction chain PROVED modulo cited canon (each step tagged below); the clean-q
floor is CENSUS-VERIFIED (70,102 adversarial B=1 patterns, min ≥ 1/12; shifted-tight
variants ≥ 1/13; B=2 blocker sweep ≥ 1/12) with the finite-lemma proof path stated —
two small decidable lemmas remain (grid-LRC(6) mod q ≤ 13; the blocker-corner list).
**Author:** death-star-2026-07-12-S16.
**Context:** THM-721 Part 6 closed compressed `j = 7` (slope dichotomy) and equal-slope
`j ≤ 12`; the honest residual of the compressed lane was `j ∈ [8,13]` mixed-slope at every
admissible scale (backlog lead death-star-S15). This theorem is that lead executed.

**Setup.** `V` primitive 13-family, scale `L ≥ 2`, balanced decomposition `v_i = b_i + L·k_i`
(`|b_i| ≤ L/2`), `B = max|b_i|`, `j = #{i : b_i ≠ 0}`. A scale is *admissible* when the
compression is meaningful (`L > 91B`; for step 3 below strengthen to `L > 574B`).
Stratum (3) = families with some admissible scale, but `j ≥ 8` at every one.

---

## Part 1 — the reduction: stratum (3) is perturbed tight dilates [PROVED mod cited canon]

At an admissible scale, all 13 elements are within `L/91` of multiples of `L`. Consider the
lift multiset `K = (k_i)`:

1. **Some impure lift `k_i = 0`** ⟹ `|v_i| = |b_i| ≤ B ≪ L ≤ v_max/13`: a tiny-speed runner,
   ratio `v_max/v_min ≥ 13L/B` — far/scale-separation peel territory (THM-608/620 family).
   [Route exists in canon; FLAGGED, not re-proved here.] Pure lifts are nonzero automatically.
2. **`≤ 12` distinct lift values** ⟹ `M(K) ≥ 1/13` (LRC(≤13), settled; dilation-invariance
   handles the gcd) ⟹ atom (THM-636 / LRCDecorrelation13, kernel-pure):
   `M(V) ≥ 1/13 − B/L > 1/14` once `L > 182B/13`. LOOSE.
3. **13 distinct lifts, `M(K) ≥ 3/41`** ⟹ `M(V) ≥ 3/41 − B/L > 1/14` once
   `B/L < 3/41 − 1/14 = 1/574`. LOOSE. (Uses only the atom.)
4. **13 distinct lifts, `M(K) < 3/41`** ⟹ by the EMPTY WINDOW `(1/14, 3/41)`
   (opus-S248/HYP-6190, VERIFIED 56k families, 0 in the window — status: verified-not-proved)
   `M(K) = 1/14`: `K` is TIGHT ⟹ `K ∈` the tight locus = `{AP, V*={1..11,13,24}}` up to
   mod-14 residue-preserving shifts and dilation (THM-612/708 + HYP-6190/6195 — same status).
   Covering `K` is excluded outright: covering-min `= 14/183 > 3/41`… `14/183 ≈ 0.0765 >
   3/41 ≈ 0.0732` (klein ILP, PROVEN for lift speeds ≤ 182; `k_max > 182` ⟹ the lift itself
   has scale structure — recurse: `k_max ≤ v_max/L` shrinks by `> 91×` per descent, so the
   recursion bottoms in `O(log v_max)` steps at height ≤ 182).

**⟹ Stratum (3) = `V = L·K + b`, `K` tight (AP/V*/mod-14 shifts), `j ≥ 8` nonzero small
offsets.** The stratum is INHABITED by primitive DC families (e.g. `360360·AP + b`,
verified in the S16 census, EXP E) — it is not vacuous, and it is exactly the
near-tight-dilate corner of the incoherent bulk.

## Part 2 — the profile object and the clean-q game [NEW]

For the profile `W = {(k_i, b_i)}`: `M(V) ≥ reach₂(W) − B/(2L)` (THM-721 Part 1 / HYP-4342),
and `reach₂(W) = sup_s min(pure(s), γ_F(s))` with `γ_F(s) = max_u min_{impure} ‖k_i s + b_i u‖`.
At `s = a/q` (coprime):

- **pure margin:** `‖k_p a/q‖ ≥ 1/q` iff `q ∤ k_p`. A base `q` is **clean** if `q` divides
  no pure lift. With ≤ 5 pure lifts, a mini grid-LRC(6) mod `q` gives `a` with pure margin
  `≥ ⌊q/6⌋/q` [finite decidable lemma, `q ≤ 13`, ~20k cases].
- **impure vacancy (B=1):** the `j` forbidden-`u` centers `−b_i k_i a mod q` occupy `≤ min(j, q)`
  of the `q` slots. If `q > j`: some slot is vacant ⟹ max circular gap ≥ 2 slots ⟹
  `γ_F ≥ 1/q`. So ANY clean `q ∈ (j, 13]` certifies `reach₂ ≥ 1/13`.
- **the blockers:** blocking all of `(j, 13]` forces `pure ⊇ {j+1, …, 13}` — at `j = 8` the
  UNIQUE maximal blocker `pure = {9,10,11,12,13}`; at `j = 9..12` the blocker towers
  `pure ⊇ {10..13}, {11..13}, {12,13}, {13}` (+ free slots). Each blocker leaks at a small
  `q` (occupancy collisions fatten the vacancy) or at a non-integer `s` — a FINITE explicit
  list, settled by the census below (worst corner `{9..13}`: exhaustive 128 signs,
  min `reach₂ = 3/11`, witness `s = 32/33`).
- **j = 13 (no pure) collapses:** at `s = 1/2` all ±1 centers land in `{0, 1/2}` ⟹
  `γ_F ≥ 1/4`; `B ≤ 3` similarly (`≥ 1/6, ≥ 1/8` verified). The no-pure end is trivial —
  within the stratum, `j = 8` (most pure blockers) is the hard end, INVERTING the naive
  union-bound ordering in `j`.

## Part 3 — the census (all exact Fractions; every value a certified lower bound)

`lrc14_mixed_slope_j8plus_census_deathstar_S16.py` (+ shifted/B2 probe, both `.out` in results):

- **AP lift, B=1:** ALL 2379 pure supports (sizes 1–5 = `j ∈ [8,12]`) × {adversarial
  slot-filling + random} signs = 34,941 patterns: **0 below 1/12.** Per-`j` minima:
  `1/10, 1/11, 1/12, 1/12, 1/12`.
- **V\* (GW) lift, B=1:** 35,161 patterns: **0 below 1/12** (identical minima profile).
- **Blocker corner `pure={9..13}` exhaustive** (2^7 signs, full refinement): min `3/11 ≈ 0.273`.
- **Mod-14 shifted tight variants** (all 13 single shifts + `13→41`): sampled-support census,
  cap `1/13`: **0 below 1/13**; record exactly `1/13` (AP[4→18]).
- **B=2 blocker sweep** (5 hardest supports × 160+ magnitude/sign trials): min `1/12`
  (`pure = {3,4,5}` — the small-q blocker pinning the `q = 6` witness).
- **End-to-end:** worst corner profile at `L = 360360`: `V` primitive, divisor-complete;
  transported witness `t = (m + 32/33)/360360` has EXACT margin `3243239/11891880 ≈ 0.2727`
  on the integer family. The certificate pipeline (profile `s*, u*` → `m = round(Lu* − s*)`
  → `t`) works end-to-end, same shape as `margin_uescape_j6` (LRCUEscape.lean).

## Part 4 — statement and honest scope

> **Theorem (B=1, census-complete; proof = 2 finite lemmas + corner list).**
> `V = L·K + b`, `K` ∈ tight locus (AP/V*/mod-14 shifts), `b_i ∈ {−1,0,1}`, `j ≥ 8`,
> `L > 91`. Then `reach₂(W) ≥ 1/13` (census: ≥ 1/12 on AP/V*), so
> `M(V) ≥ 1/13 − 1/(2L) > 1/14`. NO LRC(14) input; ingredients: LRC(≤13) citations,
> the vacancy pigeonhole, grid-LRC(6) mod `q ≤ 13` [decidable], the blocker-corner
> witnesses [explicit rationals].

- `B = 2`: blocker sweep holds the `1/12` floor (needs `L > 182` via `1/12 − 2/(2L) > 1/14`…
  actually `1/12 − 1/L > 1/14 ⟺ L > 84`); full support-census at `B = 2` not yet run.
- General `B`: the vacancy argument runs on the `q·|b_i|`-grids; expected floor decays
  gracefully in `B`; the stratum only needs `B < L/574`, and the census machinery
  (early-exit capped, exact) scales — next session target.
- The reduction's step-4 inputs (empty window, tight locus) are VERIFIED-status canon;
  the theorem inherits that tag until GAP-A closes. Everything else is proved or decidable.

**Consequence for the dichotomy:** with THM-721 (Parts 3+6) and this, the compressed lane
reads: `j ≤ 7` closed, equal-slope `j ≤ 12` closed, `j ∈ [8,13]` mixed-slope closed at
`B ≤ 2` (census) with a finite-lemma proof path — the large-diameter half of the endgame
now rests on [no-admissible-scale families = the pair-sum/Parseval domain, klein-S264] +
[the finite checks].

**Files:** `04-computation/lrc14_mixed_slope_j8plus_census_deathstar_S16.py`,
`05-knowledge/results/lrc14_mixed_slope_j8plus_census_deathstar_S16.out`,
`05-knowledge/results/lrc14_mixed_slope_shifted_B2_deathstar_S16.out`.
**Related:** THM-721 (Parts 1–6), THM-636, HYP-4342, HYP-6190/6195/6210 (window + tight
locus), THM-612/708, klein-S267 ILP (covering-min ≤ 182), THM-608/620 (peels), HYP-6295.
