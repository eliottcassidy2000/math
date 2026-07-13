# The two-scale constant is small and adversarial-robust — the rigorous gap is cross-sector cancellation, and dilation reframes the tail as primitive-spread

*klein-2026-07-12-S273. Owner: work on the explicit two-scale constant. This is the single
quantitative item left on both density-route base-row tails (k=8 `Φ ≤ cap₉`, k=9 `J ≥ 432/91`; see
[[lrc14-density-route-base-rows-tails-are-explicit-thm710-transfers]], HYP-6260/6272). The single-peel
bound is already PROVED (THM-700, kind-pasteur): `|p₀(E'∪{w}) − Plat(E')| ≤ V(E')/(6w)`, `V(E') ≤
14·Σe'`. Its own "Scope" flags two open pieces — the crude constant overshoots empirical by ~10², and
it only bites when `w ≫ Σe'`. This session pins the real situation for the actual functionals `Φ`, `J`.*

---

## Three findings that reframe "the constant"

### 1. Dilation-invariance ⟹ the tail is the *primitive-spread* regime (margin 0.146, the loosest)

`Φ` and `J` are **exactly dilation-invariant**: `Φ(cE) = Φ(E)` for integer `c ≥ 1` (substitute
`y = frac(cx)`, uniform as `x` is uniform, so the sector-occupancy law of `{frac(c e_i x)}` equals that
of `{frac(e_i y)}`). Verified to machine precision: `Φ(7·{0..7}) = Φ({0..7})` exactly.

This kills a red herring. A first sweep found "band max `Φ = 0.438`" at large diameter — but the
maximizers were `4·{0..7}` (`d=28`) and `5·{0..7}` (`d=35`): **dilations of consec-8, imprimitive, =
the compact extremal rescaled**, not new tail cases. THM-719 already restricts to *primitive* 8-cores,
so they are correctly excluded. Re-running the band over **primitive** cores only (`gcd = 1`), the
genuinely-spread tail has

> `max Φ over primitive 8-cores, d = 26..300  =  0.348`  (at `{0..6, 41}`), margin **`+0.146`** to `cap₉`.

So the two-scale constant governs only the primitive-spread tail, where `Φ ≈ 0.34` — the *loosest*
regime, with margin `0.146` (larger than the compact regime's `+0.056`). The mirror holds for `J`: its
tail moves the *safe* direction (up, toward `≈5.7`), margin `0.945` above the threshold.

### 2. The true constant is small and adversarial-robust — `C_Φ ≈ 0.9`, `C_J ≈ 7`

Measuring `err·w = |Φ(C∪{w}) − Φ_∞(C)|·w` with the grid scaled to `w` (`N_g = 400w`, so `frac(wx)` is
always resolved — the first run's `err·w ≈ 13–15` at `w ≈ N_g` were pure **aliasing artifacts**):

| | clean (prime) `w` | adversarial `w` (multiples of `lcm(C)`) | worst |
|---|---|---|---|
| `Φ`, `C = {0..6}` | 0.02–0.17 | 0.09–0.44 | **0.44** |
| `Φ`, `C = {0..5,7}` | 0.07–0.42 | 0.05–0.72 | **0.72** |
| `Φ`, `C = 2·{0..6}` | 0.01–0.05 | 0.18–0.89 | **0.89** |
| `J` (k=9 8-clusters) | — | — | **≈ 7** |

The decisive point: **adversarial `w` (sharing all of `C`'s factors, where the per-offset geometric-sum
cancellation fails) give `err·w` only ~2–3× the clean value — still `O(1)`, not `O(Σe')`.** So the
sharp constant is *uniform* (`C_Φ ≈ 0.9`, `C_J ≈ 7`, numerical estimates ±20%), and the cross-sector
cancellation that produces it is real and robust, not a clean-`w` accident. With these,

> `D₀ = C / margin`:  `Φ`: `0.9 / 0.146 ≈ 6`;  `J`: `7 / 0.945 ≈ 8`  — **both inside their exhaustive
> boxes** (k=8 `d ≤ 25`, k=9 `d ≤ 18`). If the empirical constants were rigorous, there is **no
> intermediate band** and both rows close outright.

### 3. The precise rigorous gap: a ~14× (Φ) / ~26× (J) constant improvement, all in cross-sector cancellation

The crude THM-700 constant is `V/(6w) ≈ 49/w` for `Φ` on `C={0..6}` (`Σe'=21`), giving `D₀_crude ≈
504`; for `J` the larger coefficients push `C_J^crude ≈ 440`, `D₀_crude ≈ 500`. Closure needs `D₀ ≤
box`, i.e.

> **`C_Φ ≤ 0.146·25 ≈ 3.65`  and  `C_J ≤ 0.945·18 ≈ 17`.**

Empirical (`0.9`, `7`) clears these with `4×` / `2.4×` slack; the crude (`49`, `440`) misses by `14×` /
`26×`. The clean, provable improvements do **not** span the gap:

- **`sin(πℓ/7)` factor.** Keeping `|ĝ_s(ℓ)| = |sin(πℓ/7)|/(π|ℓ|)` instead of `≤ 1/(π|ℓ|)` replaces
  `Σ_{ℓ≠0} 1/ℓ² = 2ζ(2) = 3.29` by `Σ_{ℓ≠0}|sin(πℓ/7)|/ℓ² ≈ 1.80` — a factor **`1.82×`** (the zeros on
  `7∣ℓ` contribute little; it is the small-`ℓ` magnitudes that matter). Best clean rigorous constant:
  `≈ 27/w` for `Φ`, `D₀ ≈ 185` — still an order over.
- **Cauchy–Schwarz + Parseval** on the `w`-multiples gives an *absolute* bound `|error| ≤ √6/2 ≈ 1.2`
  with **no `1/w` decay** (the energy `Σ_ℓ|f̂_s(ℓw)|²` has no automatic `1/w` gain). Useless for the
  crossover.

Everything past `1.82×` is the **cross-sector cancellation** THM-700 named as open: the seven signed
correlations `∫ f_s(x)·g_s(wx)\,dx` (with `f_s = 1\{E'` misses exactly sector `s\}`, disjoint) cancel
below their individual total-variation bounds. The adversarial-`w` robustness (finding 2) is direct
evidence this cancellation is genuine and uniform, but the crude `Σ_s V(f_s)` bound is blind to it.

## The remaining lemma, stated crisply

> **Lemma (cross-sector two-scale constant, OPEN).** There is an absolute `C₀` (independent of `E'` and
> `w`) with `|Σ_{s=0}^{6} ∫₀¹ f_s(x) g_s(wx)\,dx| ≤ C₀/w`, where `f_s = 1\{E'` misses exactly `s\}`,
> `g_s(y) = 1\{y∈[s/7,(s+1)/7)\} − 1/7`. Equivalently: the joint law of `(`missed sector of `E'`,
> sector of `frac(wx))` decouples at rate `1/w` with a constant not scaling in `Σe'`.

The data says `C₀ ≈ 1` (and its `Φ`/`J`-weighted forms `≈ 0.9`, `≈ 7`). This is a **2-D
joint-equidistribution / discrepancy** statement — the natural tool is a van der Corput / Weyl
second-moment estimate on the coupled system `(v_i x, wx)`, exploiting that the "missed-sector" process
`Σ_s s·f_s` has bounded complexity. That is the sharp, self-contained target for a focused session.

## Honest scope and the fallback

- **Not closed rigorously this session.** The uniform constant `C_Φ ≤ 3.65`, `C_J ≤ 17` remains
  conjectural (measured, not proved); the clean bounds reach only `≈ 27` / `≈ 240`.
- **But the tail is TRUE with wide margin, and pinned.** Dilation reduces the rows to primitive cores;
  the primitive-spread band is safe to `d = 300` (max `Φ = 0.348`, margin `0.146`); the constant is
  uniform and small even against adversarial `w`. The gap is a single clean lemma, not a structural
  unknown.
- **Fallback if the lemma resists:** the `sin`-improved rigorous constant (`≈ 27/w`) gives `D₀ ≈ 185`,
  and the primitive band `26 ≤ d ≤ 185` is a *finite* check — feasible only if THM-719's "max-Φ-per-
  diameter at near-consec" structure extends (a structured sweep, not brute `C(184,7)`). Monotonicity
  of `max_d Φ` (owner's original "tail-monotonicity") would also close it.

*Files: `04-computation/lrc14_two_scale_constant_klein_S273.py` (+out),
`lrc14_two_scale_constant_controlled_klein_S273.py` (+out), `lrc14_J_two_scale_constant_klein_S273.py`.
HYP-6285. Sharpens THM-699/700's crude constant toward the true crossover; consumes THM-687/688/710.
Twin of [[the-far-element-J-recursion-splits-the-threshold-from-the-exact-min-klein-S271]] and
[[the-k8-deg3-row-tail-is-an-explicit-Phi-transfer-not-a-cited-limit-klein-S272]] (the row closures
this constant would make fully rigorous). Hands the cross-sector lemma to kps (THM-700 author).*
