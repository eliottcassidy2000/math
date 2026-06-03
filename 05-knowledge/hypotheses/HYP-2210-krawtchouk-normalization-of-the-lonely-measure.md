# HYP-2210 — Krawtchouk normalization of the lonely measure: the depth distribution is a weight enumerator

**Session:** claudebox-2026-06-03-S620. **Frame:** the flow-shell / twisted-involution n=14 line (HYP-2205) +
the four-lenses overlap hierarchy (HYP-2200), now with **Krawtchouk normalization**. **Threads:** HYP-2160 (χ =
resonance-energy signature), HYP-2155 (resonance energy), Krawtchouk/Basic.lean.

## The depth distribution is a weight enumerator
At clock `t` the **codeword** is `X(t) = (1_{A₁}(t),…,1_{Aₙ}(t)) ∈ {0,1}ⁿ` (which runners are forbidden), and
`depth(t)` is its Hamming weight. So the covering-depth distribution `p_k = meas{depth = k}` is the **weight
enumerator** of the time-indexed code `{X(t)}`, and `p₀` (lonely measure) is its weight-0 coefficient. The natural
normalization of a weight enumerator is the **Krawtchouk / MacWilliams transform** — and it makes the LRC resonance
structure explicit.

## The Krawtchouk normalization (verified exact + formalized)
With the `±1` characters `sᵢ = 1 − 2·1_{Aᵢ}`, Fourier inversion on the cube gives
`2ⁿ · ∏ᵢ(1 − xᵢ) = ∑_{S⊆[n]} ∏_{i∈S}(1 − 2xᵢ)` (formalized `two_pow_mul_prod_one_sub`). Integrating:

> **`p₀ = (1/2ⁿ) ∑_{k=0}^n ρ_k`,  where  `ρ_k = ∑_w K_k(n,w) p_w`** is the Krawtchouk transform of the depth
> distribution (`K` = `Math.Krawtchouk.K`), and `ρ_k = ∑_{|S|=k} R_S`, `R_S = E_t ∏_{i∈S}(1−2·1_{Aᵢ})` the
> level-`k` resonance correlation.

Verified EXACT for n=14 (13 runners, δ=1/14): `p₀ = (1/2¹³)∑ρ_k` matches the direct measure for the wall config
`{1..11,13,14}` (0.01216), the AP (0.00000), and random (0.13423).

## Where the resonance lives (verified)
- **Levels 0,1 are the independent baseline, exactly.** `ρ₀ = 1`, `ρ₁ = n(1−4δ) = 13·5/7` — zero excess over the
  binomial baseline `C(n,k)(1−4δ)^k` for every config. (Formalized baseline: `K_at_zero` — `K_k(n,0) = C(n,k)`.)
- **All resonance lives in `ρ_{k≥2}`.** The wall and AP carry large excess from `k=2` up (`+2.1`/`+2.4` at k=2),
  growing through the high levels; random is near baseline (`+1.0`). The Krawtchouk normalization isolates exactly
  the correlation the flow-shell crosses (HYP-2205) create.
- **Independent model recovered:** `R_S = (1−4δ)^{|S|}` ⟹ `∑ρ_k = (2−4δ)ⁿ` ⟹ `p₀ = (1−2δ)ⁿ = (6/7)¹³ ≈ 0.1348`.
  The normalized lonely ratio `p₀/(6/7)¹³` reads the resonance penalty: wall 0.090, AP 0.000, random 0.996.

## The improvement: loneliness = Krawtchouk positivity (Delsarte LP)
`p₀ > 0 ⟺ ∑_k ρ_k > 0` is exactly a **Delsarte linear-programming / Krawtchouk-positivity** condition on the depth
weight enumerator. The flow-shell cross/apex structure (HYP-2205) determines the `ρ_k`; a Krawtchouk-LP lower bound
on `∑ρ_k` would prove LRC(14). Because levels 0,1 are pinned at baseline and `ρ_k` is `σ`-invariant (the
twisted-involution symmetry, HYP-2205, makes the weight enumerator even), the LP is over the `k≥2` resonance only —
a small, structured program. This connects the n=14 shell-width inequality to the coding-theory LP machine, and the
`ρ_k` are the χ-signature of HYP-2160 in the Krawtchouk basis.

## Formalized (math-lean, sorry-free)
- `Math/LonelyRunner/KrawtchoukNormalization.lean`: `two_pow_mul_prod_one_sub` (the `2ⁿ` normalization /
  weight-enumerator Krawtchouk transform), `char_sum_level_zero`.
- `Math/Krawtchouk/Basic.lean`: `K_at_zero` (`K_k(n,0) = C(n,k)`, the resonance baseline) — first content beyond the
  bare definition in the repo's Krawtchouk module.

## Open
- The Krawtchouk-LP lower bound: find a feasible dual polynomial certifying `∑ρ_k ≥ c > 0` for every
  multiple-of-14 config (⟹ LRC(14)).
- Formalize `ρ_k = ∑_w K_k(n,w) p_w` (the depth-to-Krawtchouk transform) and the level-1 baseline `ρ₁ = n(1−4δ)`.
- The MacWilliams dual: is the dual weight enumerator's positivity the apex/2q seam in the Krawtchouk basis?
