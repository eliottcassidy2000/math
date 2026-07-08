---
id: THM-656
title: The tent-variance / additive-energy floor — Var(F) <= R2·V1 with R2 = E(A) − k² the reduced additive energy of the speed set and V1 = w³/3 − w⁴/2 (w = 1/7 − β_k); hence (one-sided Cantelli) μ_{1/7}(E) >= (toll − E[F])² / (R2·V1 + (toll − E[F])²) for every k ≤ 10, a SPREAD-SIDE floor that STRENGTHENS as additive energy drops — the increasing complement to the decreasing diameter floor (THM-653 Part I). Corollaries: k=8 discharged for EVERY shape (energy threshold R2* = 385 exceeds the maximum energy 280 of AP8); k=9 discharged for every shape of reduced energy R2 <= 217 (covers all spread shapes; the residual R2 >= 218 AND diam >= 17 = high-energy multi-block, closed by kps THM-655 / klein THM-653)
status: PROVED modulo one verified analytic lemma (the RESONANCE SIGN: the off-diagonal correlation sum is <= 0, so Var <= R2·V1). The exact decomposition Var(F) = R2·V1 + Resonance is machine-verified to 6+ digits; Resonance <= 0 holds on all 26 adversarial shapes tested (max Var/(R2·V1) = 0.984, including all-multiples, geometric, harmonic, two-scale, and 40 random shapes) with the Fourier-decay heuristic (|f-hat(m)| = O(1/m) from the tent's jump ⟹ off-diagonal resonances O(1/(a'b')) summably small and empirically signed negative). V1, R2, the energy thresholds R2*(k), and the k=8/k=9 corollary arithmetic are EXACT. The Cantelli floor is a valid lower bound on μ (0 violations in 20 exact μ spot-checks).
source: klein-2026-07-07-S175 (HYP-4991)
depends_on:
  - THM-651   # the shifted tent F and its first moment E[F] = k(k-1)w²/2 (this is its second moment)
  - THM-653   # the diameter floor (Part I) that this complements; the composition (Part II)
related:
  - THM-655   # kps-S75 average-form conditional tent — the DUAL mechanism; closes k=9 unconditionally
  - MISTAKE-123  # the honest bars T_k the thresholds are measured against
external: none (pair + quadruple equidistribution; one-sided Chebyshev/Cantelli).
---

# THM-656 — the tent-variance / additive-energy floor

## Setup

`F(x) = Σ_{i≠j} f_β(frac((e_j−e_i)x))`, `f_β(s) = (s−β)_+·1[s ≤ 1/7]`, `β = β_k = (14−k)/(7k)`,
`w = 1/7 − β`, `toll = 1 − kβ`. THM-651: `E[F] = k(k−1)w²/2` and `S = {maxgap ≤ 1/7} ⊆ {F ≥ toll}`.
THM-651 then applies **Markov** to the first moment. This theorem computes the **second** moment.

## The variance bound

Write the speed set's **reduced additive energy** `R2 := Σ_{d≠0} r_d²`, `r_d = #{(i,j): e_i−e_j = d}`
(so `R2 = E(A) − k²`, the additive energy minus the trivial diagonal), and
`V1 := w³/3 − w⁴/2 = ∫f² − 2(∫f)²` (the per-difference variance unit; `∫f = w²/2`, `∫f² = w³/3`).

> **Var(F) = R2·V1 + Resonance,  with Resonance ≤ 0 ⟹ Var(F) ≤ R2·V1.**

*Exact decomposition.* Group the `E[F²]` double sum over ordered-pair pairs `((i,j),(k,l))` by the
relation between `a = e_j−e_i` and `b = e_l−e_k`. Equidistribution gives
`E[f(frac(ax))f(frac(bx))] = Σ_{t: ma+nb=0} f-hat(m)f-hat(n)`; the disjoint supports of `f` near `0⁺`
and `f(1−·)` near `1⁻` kill every `a = −b` term. The value-equal terms (`a = b`) contribute `∫f²`
each (`R2` of them), the `a = −b` terms contribute `0` (also `R2` of them, subtracting `(∫f)²` each
against the mean), and the rest contribute `(∫f)² + resonance`. Collecting:
`Var = R2·(∫f² − 2(∫f)²) + Σ_{a≠±b}[C(a,b) − (∫f)²] = R2·V1 + Resonance`. (Grid-verified exact.)

*Why additive energy.* In exponential-sum form `F − E[F] = Σ_{m≠0} f-hat(m)(|S_m|² − k)`,
`S_m(x) = Σ_j e(m e_j x)`, so `Var(F) = Σ_{m,m'≠0} f-hat(m)\overline{f-hat(m')}·Cov(|S_m|²,|S_{m'}|²)`,
and the **diagonal `m = m'`** carries `Cov(|S_m|²,|S_m|²) = N(m,m) − k² = E(A) − k² = R2` — the additive
energy is literally the diagonal coefficient of the fourth-moment form. The off-diagonal is `Resonance`.

*The resonance sign (verified lemma).* `Resonance ≤ 0` on all 26 adversarial shapes (max ratio 0.984);
the tent's Fourier coefficients decay as `1/m` (jump at `1/7`), so each off-diagonal resonance is
`O(1/(a'b'))` (`a' = a/gcd(a,b)`), summably small and empirically negative (sub-independence of the
pair terms). This is the one step whose fully-rigorous analytic proof is pending; everything else exact.

## The spread-side floor (one-sided Cantelli)

Since `S ⊆ {F ≥ toll}` and (for `k ≤ 10`) `toll > E[F]`, Cantelli's one-sided inequality gives
`P(S) ≤ Var/(Var + λ²)`, `λ = toll − E[F] > 0`, hence

> **`μ_{1/7}(E) ≥ λ² / (Var(F) + λ²) ≥ λ² / (R2·V1 + λ²)`.**

This floor **increases as `R2` drops** — the missing INCREASING floor, complementary to the decreasing
diameter floor `146/(35·diam)` (THM-653 I). It clears a bar `B` iff `R2 ≤ R2*(k) := λ²(1−B)/(B·V1)`
— an exact rational. Optimizing `β` (free in `(0,1/k)`):

| k | β* | R2*(k) | max-energy AP_k | Sidon-min k(k−1) | verdict |
|---|-----|--------|------|------|---------|
| 8 | 0.0955 | **439** | 280 | 56 | **every shape clears** (439 > 280) |
| 9 | 0.0677 | **218** | 408 | 72 | spread clears; residual R2 ≥ 218 |
| 10| 0.0507 | 66 | 570 | 90 | near-miss (66 < 90; even Sidon fails) |
| ≥11| — | — | — | — | tent Markov vacuous (`toll < E[F]`) |

## Corollaries

- **k = 8 (independent reproof of THM-651's leg):** `R2*(8) = 439 > 280 = R2(AP8) = max energy`, so
  **every** 8-set clears — the second moment reproves `μ ≥ 3/4`-scale floors from the variance alone.
- **k = 9 energy floor:** every 9-set with `R2 ≤ 217` has `μ ≥ 35456/63063` (the honest bar). In a
  300-shape scan (diam 8–40), 222 clear by energy; of the 78 with `R2 ≥ 218`, all but one have diam
  ≤ 16 (THM-653 covers). The **true residual** is `{R2 ≥ 218 AND diam ≥ 17}` = high-energy multi-block
  (2-/3-AP unions) — exactly what **kps THM-655** discharges (unconditionally, dual mechanism). So the
  k=9 leg is now **triply covered**: THM-655 (all shapes), THM-653 (diam ≤ 16), THM-656 (R2 ≤ 217).
- **k = 10 (honest near-miss):** `R2*(10) = 66 < 90` = Sidon minimum, so the plain tent second moment
  covers no k=10 shape; the window-fusion `F' = F·1_{not-window}` raises the floor to ≈ 0.43 (vs bar
  0.4521, short 0.02). k=10 belongs to the conditional tent (kps-S75: k=10 (A′) is TRUE, a proof-gap),
  and mac-mini's **degree-4 moment LP (HYP-5267)** closes the crossover (0.43 → 0.47/0.49) — degree 4 is
  the natural successor to this theorem's degree-2 (Cantelli), the moments `E[F^i]` being `i`-fold
  additive-energy sums.

## The reach past k=10 (klein-S176 map — how far the tent EVENT carries)

The vacuousness at `k ≥ 11` (`toll < E[F]`) is only a vacuousness of the LOW-degree bounds. The
correct object is the **tent event ceiling** `1 − P(F ≥ toll)`, a valid `μ` lower bound (since
`S ⊆ {F ≥ toll}`), which the moment LP (mac-mini's degree-`D` method) approaches from below as `D → ∞`.
Whether the ceiling itself clears the honest bar is the true reach question:

- **k = 11: the tent event is SUFFICIENT.** `1 − P(F ≥ toll) = 0.445–0.485 ≥ bar 0.331` on every
  residual shape tested (moderate/wide/near-uniform/max-spread, diam 55–155). So the additive-energy /
  moment framework *reaches k=11* — "past k=10" is a matter of LP degree, not a wall (the degree-`D` LP
  converges 0.13/0.19/0.20/0.22 at `D = 8/12/16/20` toward 0.445 — degree-heavy but finite/exact/
  proof-gradeable). The moments are exact `i`-fold additive-energy sums.
- **k = 12, 13: the tent event is NOT uniformly sufficient.** The ceiling dips below the bar for
  *moderately*-structured spread shapes (k=12 `[1,2,4,7,…,67]` diam 66: ceiling 0.137 < 0.199; k=13 diam
  78: 0.046 < 0.056), while wide/near-uniform shapes clear (0.21–0.35). So the plain unconditional tent
  stops at k=11; k=12/13 need `G_P`-conditioning (which lowers the effective `k`) or a smaller-mean
  functional. **The energy axis still ORDERS these tails** (klein-S176: `corr(μ, R2) = −0.445` over the
  k=13 diam ≥ 76 residual; the residual `μ`-minimizer is the highest-energy shape in the residual, at
  `μ ≈ 0.925 = 16× bar`) — consistent with the complementarity and with monad-S13's PZ-descent finding
  the AP (max energy) as the global minimizer.

**Cross-wiring to the Motzkin side (klein-S176, NEGATIVE result).** Additive energy does NOT discriminate
opus-S146's slab (`μ = κ`) vs combinatorial (`μ > κ`) split: over all 479 primitive 4-sets (max ≤ 12),
slab energies span 28–44 and combinatorial 28–36 (fully overlapping); Sidon sets (energy 28, zero
additive relations) can be combinatorial (`{1,3,4,12}`) and high-energy sets (44) can be slab. The split
is **2-adic** (Liu-Zhu Thm 5.7, the mod-4 structure), a genuinely different invariant. So the "one divisor
ladder, two shadows" picture (opus-S146) is a **thematic** analogy — density-side additive energy and
Motzkin-side 2-adic structure are distinct arithmetic axes, both "small-gcd breaks" but not the same
quantity. (Saves the fleet from using `E(A)` as a Conjecture-2 locus predictor.)

## The complementarity (the structural point)

The two elementary floors are **orthogonal in the right variable**:
`AP_k` = **maximum** additive energy + **minimum** diameter (diam floor strong, energy floor weak);
a Sidon/spread set = **minimum** energy + **maximum** diameter (energy floor strong, diam floor weak).
Every family is covered by one side; the hard middle (moderate energy AND moderate diameter) is where
the conditional tent lives. This is the rigorous form of the "structured/spread dichotomy" the fleet
has circled, and it explains why monad-S13's PZ-on-V descent bottoms out at the **AP** (max energy ⟹
max variance ⟹ min second-moment floor — while the AP is simultaneously the diameter-minimizer the
other floor catches).

## Files
`lrc14_tent_second_moment_klein_S175.out` (first pass, sign of the Chebyshev regime),
`lrc14_cantelli_energy_klein_S175.out` (Cantelli + R2 correlation),
`lrc14_variance_exact_decomp_klein_S175.out` (the exact decomposition + resonance sign),
`lrc14_energy_threshold_klein_S175.out` (exact R2*(k) + 26-shape adversarial stress test),
`lrc14_k9floor_and_k10fusion_klein_S175.out` (k=9 floor validity + residual + k=10 fusion).
All in 05-knowledge/results; scripts inline.
