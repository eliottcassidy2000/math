# Audit 4 (independent, adversarial) of section 2d (S22) of `collatz_five_mirrors_20260929.md`

**Object.** Section 2d ("the renewal structure of the resonant coefficient"), the S22 rows of the section-7 status
table, and the S22 clause of the status paragraph. Supporting material read: sections 2, 2b, 2c, 4;
`collatz_five_mirrors_gap/renewal/renewal_bound/renewal_cert_20260929.py` and their `.out` files;
`collatz_five_mirrors_coherent_20260929.py` (for `closed_family`).
**Method.** Every claim was re-derived before the note's proof was read, and every number was recomputed with my own
code, `04-computation/experiments/collatz_five_mirrors_20260929_audit4.py` (output
`05-knowledge/results/collatz_five_mirrors_20260929_audit4.out`, 13 s), which imports none of the repo's S22 scripts:
the law `mu_h` itself on `Z/3^h` for `h <= 10` (forward DP), a closed recursion written from the definition (exact
shrinking window, truncations `a <= 40, 50, 60`, extended to `n = 160`), a walk DP over the words from the top level
with state `(T_j, counter)` and no window in the total cost (the repo's DP covers `94%` of the mass at `h = 80` and
`63%` at `h = 120`; mine `1 - 10^(-10)`), the exact cost laws, the bound `B_h(s) = Σ_J mass_h(J) Ñ_(J-1)` at every
`h <= 150` and every `0 <= s <= floor(h log_2 3)`, and a 30-digit `mpmath` rerun of the closed recursion to isolate
float64 rounding. Notation as in the note: `m = log_2 3`, `θ* = ln(2(1 - log_3 2))`, `I = θ* m - ln(m - 1)`,
`e^(-I) = 3^(h*-1)`.

## 1. Claims and verdicts

1. **Constants** (`0.73814, -0.30362, 0.054979, 0.946505, 1.70951, 0.98699`; the identity `e^(-I) = 3^(h*-1)`;
   `Z(θ*) = m - 1`; tilted mean `m`; bracket `1 + 131.4 C`; `(3/2)/|ln 0.987| ≈ 115`; `≈ 6000` levels; window width
   `0.61 √h`). **HOLDS.** All reproduced (output A): `e^(-I)` and `3^(h*-1)` agree to `3·10^(-16)`; the term-by-term
   identity `-I = -(1/p) ln(2q) + ln q - ln p = -ln p - (q/p) ln q - ln 3` holds (both `-0.0549794728`); the bracket is
   `1 + 131.4 C` (`= 206` at the note's `C = 1.56`); the tilted standard deviation over `m` is `0.608`.

2. **Lemma G, the 2-adic reading** `(t 2^(-r) mod 3^n)/3^n = t/(2^r 3^n) + m_r/2^r`, `m_r = (-t 3^(-n)) mod 2^r`.
   **HOLDS**, and is an integer identity: `2^r x = t + m 3^n` with `0 <= x < 3^n` forces `0 <= m < 2^r` (because
   `0 < t < 3^n`) and `m ≡ -t 3^(-n) mod 2^r`; checked as an integer identity over all units `n <= 6`, `r <= 20`
   (`14560` cases, `0` failures) and as floats to `5·10^(-16)` (the note says `4·10^(-16)`; the `n = 5` maximum is
   `5.0·10^(-16)`). The unit hypothesis is not needed (any `0 < t < 3^n`).

3. **Lemma G, the two-term bound and the classes.** **HOLDS.** With `α = (θ + ξ_1)/2`, `β = (θ + ξ_1 + 2ξ_2)/4`,
   `|(1/2)e(α) + (1/4)e(β)|^2 = 5/16 + cos(2π(α - β))/4` and `α - β = -(2ξ_2 - ξ_1 - θ)/4 = -φ`; tail `<= 1/4`.
   `ξ_1 ≠ ξ_2` gives `cos(2πφ) ∈ {-cos(πθ/2), -sin(πθ/2)} < 0`, hence `|G| < (1 + √5)/4`. Class identification:
   `3^(-n) ≡ (-1)^n mod 4`, so `ξ ≡ (-1)^(n+1) t mod 4`; for odd `t`, `ξ ≡ 1 mod 4` iff `t ≡ (-1)^(n+1) mod 4`
   (checked for every unit, `n <= 9`, `0` failures). The bound `|G| <= 1/4 + (5/16 + cos(2πφ)/4)^(1/2)` was checked
   against every unit `n <= 9`: never violated (worst slack `2.7·10^(-3)` at `n = 9`). The corollary "`2^(-m) mod 3^n`
   has a gap iff the digits `m+1, m+2` of `-3^(-n)` differ" **HOLDS** (`ξ(2^(-m)) = (-3^(-n)) >> m` digit-wise;
   `270` cases, `0` failures; the gap holds at `32, 29, 30` of the `60` steps `m <= 60` at `n = 20, 40, 80`).

4. **Lemma G, the numerical prose.** Three statements are wrong or loose:
   (a) "the class `4 | t` with `θ <= 1/2` never exceeds `0.846`" — **FAILS**, the inequality is reversed: the class
   `4 | t, θ <= 1/2` *contains the resonance* (`0.9436, 0.9714, 0.9859, 0.9937, 0.9967` at `n = 5..9`); it is the
   class `4 | t, θ > 1/2` (equivalently `t` odd, `ξ ≡ 3 mod 4`, `θ <= 1/2`) that never exceeds `0.8462` (output B).
   (b) "the resonant set named: `t = 2^v u` small or `t = -(odd) small`" — **FAILS** for the second half: an odd `t`
   with `θ -> 1` is `t = 3^n - t'` with `t'` *even*, and the resonant class `t ≡ (-1)^n mod 4` forces `4 | t'`; the
   involution `t -> 3^n - t` maps `{v_2(t) >= 2}` onto `{t odd, ξ ≡ 3 mod 4}`. The units with `|G_9(t)| > 0.99` are
   exactly `±2^8, ±2^9, ±2^10, ±2^11, ±2^8·5`: the resonant set is `t = ±2^v u` with `v >= 2` and `2^v u` small (the
   section-7 table row already says `±2^v u`; only the body text is wrong).
   (c) "The gap classes are exactly half the units (`6562` of `13122` at `n = 9`)" — **HOLDS WITH CORRECTION**:
   `6562` gap against `6560` non-gap; half up to one, not exactly half (the two gap classes have `3281` elements each,
   the two others `3280`). "the two gap classes never exceed `0.672`" **HOLDS** (`0.6718`).

5. **The exponent-walk identity** `e(2^s Y_h/3^h) = Π_j e((2^(s - T_j) mod 3^j)/3^j)` (modular inverses for
   `s < T_j`) and **the real-phase identity** `e(2^s Y_h/3^h) = e(2^(s-A) C_w(3,2)/3^h)` for `A <= s`. **HOLDS.**
   Proof: `2^s Y_h = Σ_j 3^(h-j) 2^(s-T_j)` in `Z_3`, and `e(3^(h-j) x/3^h) = e(x/3^j)` for `x ∈ Z_3` (well defined
   on `Z_3/3^j`); for `A <= s` the same sum is the integer `2^(s-A) Σ_j 3^(h-j) 2^(P_(j-1))`. Brute force at `h = 5`
   (all `3125` words with letters `<= 5`, `s ∈ [-8, 20]`, including negative `s`): `2.6·10^(-15)`; the rational
   formula `Y_h = Σ_j 3^(h-j) 2^(-T_j) mod 3^h` agrees with the forward recursion exactly; the real-phase identity to
   `0` (output C). **Ramanujan average** `μ(3^n)/L_n` **HOLDS** (`2` is a primitive root mod `3^n`; `-1/2` at `n = 1`,
   `0` for `n = 2..7`).

6. **Lemma R' (statement and proof).** **HOLDS WITH CORRECTION: `s >= 0` is required.** All proof steps are correct:
   `T_j` is strictly decreasing in `j`, so `{j : T_j > s} = {1..J}` and `J(w) = J` iff `T_(J+1) <= s < T_J`; with
   `w = (v, a_J, u)` the exponents `κ_j` for `j > J` depend on `u` only, `κ_J = κ - a`, and for `j < J`,
   `κ_j = (κ - a) - T_j^(v)` is the exponent walk of `v` at level `J - 1` started at `κ - a < 0`, so
   `E_v Π_(j<J) ω_j(κ_j) = mu_hat_(J-1)(2^(κ-a))` by item 5; `v, a_J, u` are independent; the top-part phases
   (modulus `1`) are dropped by the triangle inequality; `Σ_(m>=1) 2^(-κ-m)|mu_hat_(J-1)(2^(-m))| = 2^(-κ) Ñ_(J-1)`
   and `Σ_κ P(T_(J+1) = s - κ) 2^(-κ) = mass_h(J)`; `J = 1` uses `Ñ_0 = 1` (`mu_hat_0 ≡ 1`), `J = h` uses the empty
   `u` with `κ = s` and `mass_h(h) = 2^(-s)`. But for `s < 0` the stated definition gives `mass_h(h) = P(0 <= s <
   a_h) = 0` while `c_h = mu_hat_h(2^s) ≠ 0`, so "for every integer `s`" is false; the lemma holds for `s >= 0`
   (which is all Theorem C needs).

7. **Lemma R' term by term.** **HOLDS**, and more strongly than stated; one number needs correction. With my
   complete walk DP (output F): `h = 40, s* = 57`: the ratios `0.016, 0.029, 0.053, 0.136, 0.273, 0.439, 0.171, 0.512,
   0.732, 0.484, 0.468, 0.206, 0.695` and the worst `0.840` (at `J = 20`) are reproduced exactly; `h = 80, s* = 120`:
   the worst ratio over *all* `J` is `0.904` (`J = 64`), the note's `0.810` being the worst over `J <= 30` (the repo's
   DP only counted `J <= 30`, and its `c_J` for `J >= 21` are incomplete at `h = 80`); `h = 120`: `0.913` (`J = 64`);
   `h = 20`: `0.884`. `Σ_J c_J` agrees with the closed recursion to `10^(-12)` relative at `h = 20..120`, and with the
   law itself to `10^(-16)` at `h = 6, 8, 10` (output D). Beyond the note: the inequality `|m_h(s)| <= B_h(s)` was
   checked at **every** `h <= 81` and **every** `0 <= s <= floor(h m)`: worst ratio `0.918` (`h = 64, s = 0`); for
   `82 <= h <= 150` worst `0.109` (output H, H2). The exact masses agree with the DP masses to `5·10^(-12)`, and the
   note's "`J <= 20` complete to `< 1%`" for the repo's DP holds.

8. **The Chernoff mass law** `mass_h(J) <= P(T_(J+1) <= s) <= e^(-θ* s) Z(θ*)^(h-J) = e^(-hI) e^(θ* δ) (m-1)^(-J)`.
   **HOLDS.** Standard Chernoff at the fixed tilt `θ* < 0` (`P(S_N <= s) <= e^(-θ s) Z(θ)^N`), `Z(θ*) = q/(1-q) =
   m - 1`, tilted mean `m`; the algebra `e^(-θ* s) Z^(h-J)` with `s = hm - δ` gives `e^(-hI) e^(θ*δ) (m-1)^(-J)`
   (checked). Numerically never violated (output G: `max_J (mass - Chernoff) < 0` at `h = 40..1000`). The normalised
   growth factors `1.35, 2.05, 1.90, 1.76, ...` (`h = 40`) and `1.72, 1.71, 1.71, 1.71, 1.70, ...` (`h = 1000`,
   after the first ratio `0.73`) are reproduced; the Gaussian-window reading (centre `≈ δ/m ≈ 4`, width `0.61 √h`)
   is a heuristic consistent with the data (the ratios cross `1.7095` near `J = 4` at `h = 40`, near `J = 2..3` at
   `h = 1000`). The mass itself peaks at `J ≈ 0.21 h + 3` (`11, 19, 28, 65, 210` at `h = 40, 80, 120, 300, 1000`).

9. **Theorem C (conditional).** **HOLDS WITH CORRECTION** (`0 <= s`). The derivation is correct:
   `Σ_J mass_h(J) Ñ_(J-1) <= e^(-hI) e^(θ*δ) [1 + C Σ_(J>=1) (m-1)^(-J) ρ^(J-1)] = e^(-hI) e^(θ*δ) [1 + C (m-1)^(-1)
   /(1 - ρ/(m-1))]`, geometric series convergent iff `ρ < m - 1`; `1 + 131.4 C` at `ρ = 3^(-1/2)` is right. Two
   remarks: H with `n = 0` forces `C >= 1` (harmless); the statement must read `0 <= s = hm - δ` (item 6), and the
   "Result in one paragraph" sentence "for every `h` and every `s <= hm`" likewise. The constant is large: with the
   observed `C = 1.56` (or `1.3`), `C' ≈ 206` (`172`), a thousand times the computed finite-`h` bounds — this is the
   Gaussian window, not an error, but "with a *constant* prefactor" should not be read as "with the constant `0.18`".

10. **`Ñ_n` to level 80 and its rate.** **HOLDS WITH CORRECTION** of the stated ranges. My recursion (window `[-60,
    0]`, `a <= 50`) reproduces the note's values (`1.000, 0.927, 0.920, 0.952, 0.896, 0.822, 0.562, 0.555, 0.445,
    0.275, 0.476, 0.667` for `n = 1..12`; `0.5736` least-squares rate over `20..80`, `0.5574` over `40..80`;
    `<= 1.292` for `20 <= n <= 80`, max `1.562` at `n = 16`); the truncations `a <= 40 / 50 / 60` agree to
    `1.3·10^(-11)` / `1.3·10^(-14)` (the note's `1.3·10^(-11)` for `40` vs `60` is right); float64 rounding is
    `1.1·10^(-14)` relative at `n = 40` and `9.7·10^(-15)` at `n = 60` against 30-digit `mpmath` with the same
    truncation (output J). But the note read only the printed levels: the true range of `Ñ_n 3^(n/2)` over `n <= 80`
    is `[0.055, 1.56]` (minimum at `n = 78`: `0.056, 0.055, 0.113` at `n = 77, 78, 79`), not `[0.12, 1.56]`; and the
    sup-norm range `N_n 3^(n/2) ∈ [0.77, 3.5]` is `[0.67, 3.77]` (`n = 77`, `n = 29`). Note also that the truncation
    is not self-certifying a priori (`2^(-50) = 9·10^(-16)` per level against `Ñ_80 = 1.6·10^(-20)`); the
    certification rests on the agreement of the three truncations plus the `mpmath` check, which is adequate.

11. **The negative family beyond level 80 (obligation (3) of the note, done here to `n = 160`).** **NEW FINDING,
    material to the status of H.** `Ñ_n 3^(n/2)` stays in `[0.029, 0.31]` for `81 <= n <= 122`, then rises to
    `0.90, 3.37, 7.97, 8.67, 2.38, 1.80` at `n = 124, 126, 128, 130, 132, 134` and returns to `0.03–0.35` at the
    even levels `136..160` (output E2; the maximum over all `100 <= n <= 160` is `8.67` at `n = 130`). The cause is a single coefficient wave: the sup `sup_m |m_n(-m)| 3^(n/2)` moves
    from `m = 60` at `n = 84` to `m = 1` at `n = 130` (about `-1.3` in `m` per level) while growing from `1.0` to
    `13.4` (about `×1.06` per level, i.e. these coefficients decay at rate `0.61`, not `0.577`), then dies; a second,
    weaker wave (`0.2–0.8`) runs from `(136, 52)` to `(160, 21)`. The wave is truncation-stable (`a <= 50` against
    `a <= 60`: `3.6·10^(-13)` relative at `n = 84..120`, `m <= 20`) and cannot be float64 rounding (absolute rounding
    error is bounded by `n · 10^(-16) · 3^(-n/2)`). The one-step Gauss sums along the wave are generic (`|G_n(2^(-m*))|
    = 0.17–0.74`) except at `n = 128` (`0.967`, where `3^128 ≡ 1 mod 2^9` makes `-3^(-128) ≡ -1 mod 2^9`, a run of
    nine equal digits), so this is *not* a Lemma-G one-step resonance but a multi-level coherence whose mechanism is
    OPEN. Consequences: (i) H is **not refuted** — the least-squares rate of `Ñ_n` over `20..160` is `0.5716`, over
    `100..160` `0.5685`, both below `0.585`; (ii) the note's description "the negative-power coefficients sit at the
    Parseval scale" and the extrapolation constant `1.3` of `renewal_cert` (iv) **fail** between `n = 124` and `134`
    (weighted norm up to `8.7` times the Parseval scale); the "`1.3%` margin" is therefore not the right measure of
    the risk — the risk is a wave that grows without saturating; (iii) the bound values at `h <= 300` are unaffected
    (the mass at `J - 1 ∈ [124, 134]` is `< 10^(-7)` at `h = 300`); at `h ≈ 600`, where the mass peak reaches
    `J ≈ 126`, those terms would enter with a factor up to `7` above the extrapolation, still a constant.

12. **"For `h <= 81` every `Ñ_(J-1)` in the bound is a computed number, so the inequality `M_res(h) <= 0.18 e^(-hI)`
    there rests on Lemma R' and float64 arithmetic only."** **FAILS as stated.** `M_res(h)` is not defined anywhere.
    If it means `|mu_hat_h(2^(s*))|` at the observed argmax: the computed bound `B_h(s*)/e^(-hI)` over `40 <= h <= 81`
    is at most `0.213` (`h = 44`; `0.202` at `h = 41`), over `20 <= h <= 39` at most `0.283`, over `h <= 19` up to
    `1.06` (output H) — so "`0.18`" is right only at the three levels `40, 60, 80` the note computed, and the sentence
    silently extrapolates across `h`. If it means the maximum over `0 <= s <= hm` (the quantity Theorem C is about),
    the bound is much larger away from the argmax: `max_s B_h(s)/e^(-hI) = 0.68–0.95` for `20 <= h <= 81` (attained
    at `s = floor(hm)`, `δ < 1`, where `e^(θ*δ) ≈ 1`), maximum `0.872` over `40 <= h <= 81` and `1.057` over all
    `h <= 81`. What is true and computed (from `Ñ_n`, `n <= 80`, the exact masses, and Lemma R'): for
    `40 <= h <= 81`, `|mu_hat_h(2^(s*))| <= B_h(s*) <= 0.213 e^(-hI)` and `max_(0<=s<=hm) |mu_hat_h(2^s)| <= 0.872
    e^(-hI)`; with `Ñ_n` computed to `n = 120`, `B_h(s*) <= 0.164 e^(-hI)` for `82 <= h <= 150` (output H2). The
    actual coefficients are far below the bound (`|m_h(s*)|/B_h(s*) = 0.07 -> 0.019` over `40..150`).

13. **The bound values `0.18, 0.14, 0.11 e^(-hI)` at `h = 40, 100, 300` and the row `0.180, 0.178, 0.136, 0.141,
    0.148, 0.120, 0.105, 0.124, 0.110`.** **HOLDS numerically WITH CORRECTION of status.** Reproduced (`0.1795,
    0.1783, 0.1364, 0.1414, 0.1477, 0.1204, 0.1052, 0.1235, 0.1102`). But for `h >= 200` the value depends on `Ñ_n`
    beyond the computed range: without `Ñ_n` for `n > 80` the bound carries the unconditional extra term
    `Σ_(J>=82) mass_h(J)` (`Ñ <= 1`), which is `6.7·10^(-7) e^(-hI)` at `h = 150` (negligible) but `0.93 e^(-hI)` at
    `h = 200`, `3.4·10^3 e^(-hI)` at `h = 250` and `10^6 e^(-hI)` at `h = 300`; with my `Ñ_n` to `n = 120` the tail
    `Σ_(J>=122) mass_h(J)` is `2·10^(-14) e^(-hI)` at `h = 200`, `4.8·10^(-6) e^(-hI)` at `250`, `1.48 e^(-hI)` at
    `300`. So the `h = 200, 250, 300` entries (and the `h = 2000` one, which the note does label) are "under the
    hypothesis with the observed constant" — and after item 11 that constant is not `1.3`. The parenthesis "(`0.08
    e^(-hI)` at `h = 2000` under the hypothesis with the observed constant)" should cover `h >= 200` (with `n <= 80`)
    or `h >= 300` (with `n <= 120`). The `renewal_bound` script's docstring says the extrapolation constant is `3.6`
    while the code uses `1.292`; harmless but should be fixed.

14. **The profile (C1).** Ceiling side **HOLDS**: rms `2.23·10^(-10)` (`n = 40`, against `3^(-20) = 2.87·10^(-10)`)
    and `2.22·10^(-20)` (`n = 80`, against `8.2·10^(-20)`) over `k ∈ [-60, -1]`; argmax `k = nm - 6.40, -6.80`.
    Floor side "`≈ 2^(-d) M(n)`" **FAILS quantitatively**: `|m_n(k_c + d)|/M(n)` equals `0.06–0.14 · 2^(-d)` for
    every `-2 <= d <= 12` at both `n = 40` and `n = 80` (output I) — the *decay* `2^(-d)` per step is right, the
    prefactor is `≈ 0.1`, not `1`; the note's listed values (`7.7·10^(-6), 2.5·10^(-6), ..., 4.6·10^(-9)` "for `d =
    1..9`") are the values at `k = floor(nm), ..., floor(nm) + 8`, i.e. the note's `d = 1` is `k - nm = -0.8`, below
    the floor; and "`4·10^(-12)` at `d ∈ [21, 60]`" is the maximum (rms `7.4·10^(-13)`). The low (corridor) side of
    the peak decays at `0.71–0.94` per step (`n = 40`) and `0.76–0.94` (`n = 80`), steepening away from the peak,
    consistent with the note's "`≈ 0.85^δ` near the peak". The three-regime Gauss-sum means (C4) were not recomputed.

15. **`J_eff = 10, 14, 19` at `h = 40, 80, 120` (`≈ 1.6 √h`).** **HOLDS**, with two caveats. Reproduced as the first
    `J` with remainder `< 10%`; the remainder is not monotone (at `h = 80` it is `0.111` at `J = 7`, `0.074` at `14`,
    back to `0.104` at `16`), so the stable-10% criterion gives `10, 17, 19` and the stable-5% `10, 18, 19`; at
    `h = 20`, `J_eff = 8` (all criteria). Adding `h = 20`, `1.6 √h = 7.2, 10.1, 14.3, 17.5` against `8, 10, 14, 19`
    does fit better than `h/8` and than a two-parameter log law (`8.2 ln h - 20 = 4.6, 10.2, 15.9, 19.3`), so the
    `√h` reading is supported by four points; it remains OBSERVED (the status paragraph lists it under VERIFIED).

16. **The per-mass weights `|c_J|/mass_h(J) · h`** (`1.94, 1.89, 1.48` at `J = 4`, etc.) **HOLD** exactly (output F;
    at `h = 120` from a complete DP, the repo's `63%`-mass DP being complete for `J <= 12`). The `h = 20` row
    (`1.59, 1.75, 1.71, 1.44, 0.28, 0.40, 0.21, 0.08, 0.03` for `J = 2..10`) continues the pattern.

17. **The `J`-series table** (`|c_J|/|full| = 0.02, 0.05, 0.10, 0.26, 0.52, 0.81, 0.26, 0.59, 0.44, 0.20, 0.10`, the
    arguments, `0.65, 0.64` and `0.59, 0.64, 0.62`). **HOLDS** (reproduced to the printed digits).

18. **The floor counter `F` (words that leave the corridor, `K = 3`).** "carry `2.4%, 7.8%` of the mass" —
    `2.4%` **HOLDS**, `7.8%` **FAILS**: the `F >= 1` mass at `h = 80` is `1.97%`; the note's `7.8% = 1 - 0.9218`
    counts the `5.9%` of the mass outside the repo DP's cost window (`A > s + 60`), which are deep-ceiling words with
    `F = 0`. "`< 2.3%, < 1.3%` of the coefficient" is the *largest single term* (`|c_(F=1)|/|full| = 2.27%, 1.32%`),
    not the total: the vector remainder `|full - c_(F=0)|/|full|` is `0.42%, 0.38%` (so the result paragraph's
    "`< 1%`" **HOLDS**, and "`100.4%, 99.65%`" **HOLDS**), and the sum of moduli `Σ_(F>=1)|c_F|/|full|` is `6.5%,
    3.4%`. The two sentences of the note are mutually inconsistent as written.

19. **The three regimes / drift sentence.** "starts `δ` bits below the critical line and drifts *up* relative to the
    line at `2 - m = 0.415` bits per level" — **wording/sign**: `κ_(j-1) - (j-1)m = (κ_j - jm) + (m - a_(j-1))`, so
    `κ_j - jm` *decreases* by `2 - m` per level on the way down: the walk moves away from the floor toward the ceiling
    (which is what the rest of the paragraph uses: "stays in the corridor until the exponent turns negative"; `J ≈
    0.21 h` **HOLDS**, the mass peaks at `J ≈ 0.21h + 3`).

20. **Statuses and framing.**
    (a) Status paragraph: "the resonant window decays at the rate `3^(h*-1)` with a constant prefactor" — under H
    Theorem C gives an *upper* bound at that rate; "decays at least at the rate" (the body says so; the status
    paragraph does not). (b) "VERIFIED — ... `J_eff ≈ 1.6 √h`" — OBSERVED/DIRECTION (item 15). (c) "(T2) is reduced
    to one constant" — the reduction `M ≍ P_h iff Σ_J γ_J ≠ 0` presupposes that the limits `γ_J` exist, which the
    body lists as obligation (2) and labels DIRECTION; in the status paragraph it reads as a result. (d) "the S21
    'levels `180..300` lean against `M ≍ P_h`' of section 2c is withdrawn as evidence: ... Theorem C forbids any rate
    above `0.9465` under H" — the withdrawal is justified **only under H** (CONJECTURAL) and **only for the
    exponential rate**: `M ≍ P_h` also asserts the prefactor (`h^(-3/2)`), which Theorem C does not address (it
    allows `M/P_h` to grow polynomially), and the observed rise of `M/P_h` by `1.26` over `180..300` is exactly a
    prefactor statement; moreover section 2c and the section-7 rows still say "lean slightly against it", so the
    note is internally inconsistent. (e) "doubly exponentially better than the `C* h^(-6409)` of the Fourier–renewal
    method" — compares a *conditional* bound on *one family* (`t = 2^s`, `0 <= s <= hm`) with the cited
    *unconditional* bound on all of `M(h)`; "doubly exponentially" is rhetoric (the honest comparison: geometric with
    an explicit constant versus polynomial with a tower constant, on a sub-family and under H). (f) "The
    Fourier–renewal method's structured/unstructured dichotomy gives `(1-c)^n` at these frequencies with a `c` far
    too small" — **not supported by the cited material**: the note cites Tao only from memory via Mazur's (2.3) and
    Mazur's paragraph on the polynomial bound `C* h^(-6409)`; nothing cited says what the method yields on the
    negative powers of two, and the note should mark this as an unverified recollection or drop it. (g) The section-7
    row "the negative family at the Parseval scale ... | VERIFIED" and the sentence "The hypothesis is VERIFIED to
    level `80`" must be re-read after item 11: VERIFIED describes `n <= 80` only, and the family leaves the Parseval
    scale by a factor `8.7` at `n = 130`. (h) The row "the bound `= 0.11–0.18 e^(-hI)` at `h = 40..300` | VERIFIED"
    conflates computed (`h <= 150`) and hypothesis-dependent (`h >= 200`) values (item 13). (i) "the profile
    asymmetry: ... floor side `≈ 2^(-d) M(n)` | OBSERVED" needs the prefactor `0.1` (item 14). (j) The Lemma G row
    "`<= 0.809` on half the units" is fine; the body's three prose errors (item 4) are not in the row. (k) The
    reading "H is a square-root-cancellation statement for one explicit character sum ... of the `×2 ×3` kind" and
    "Lemma G gives a gap at about half the steps but does not by itself control the product" are fair and correctly
    labelled CONJECTURAL/DIRECTION; item 11 adds that the product *does* develop coherent waves that one-step gaps do
    not see. No Collatz proof step is claimed, correctly.

## 2. Textual corrections to apply (exact phrase -> replacement)

1. Lemma G paragraph: "Checked against the modular sum to `4·10^(-16)` (`n <= 5`)" -> "Checked against the modular
   sum to `5·10^(-16)` (`n <= 5`)".
2. Lemma G paragraph: "the class `4 | t` with `θ <= 1/2` never exceeds `0.846`" -> "the class `4 | t` with
   `θ > 1/2` (equivalently `t` odd, `ξ ≡ 3 mod 4`, `θ <= 1/2`) never exceeds `0.846`".
3. Lemma G paragraph: "with the resonant set named: `t = 2^v u` small or `t = -(odd) small`" -> "with the resonant
   set named: `t = ±2^v u` with `v >= 2` and `2^v u` small (an odd resonant `t` is `t = 3^n - 4u'` with `4u'`
   small); at `n = 9` the units with `|G_9| > 0.99` are `±2^8, ±2^9, ±2^10, ±2^11, ±5·2^8`".
4. Lemma G paragraph: "The gap classes are exactly half the units (`6562` of `13122` at `n = 9`)" -> "The gap
   classes are half the units up to one (`6562` against `6560` at `n = 9`; the involution `t -> 3^n - t` pairs
   `{v_2(t) = 1}` with `{t odd, ξ ≡ 1 mod 4}` and `{4 | t}` with `{t odd, ξ ≡ 3 mod 4}`)".
5. Lemma R' statement: "For every `h >= 1` and every integer `s`" -> "For every `h >= 1` and every integer `s >= 0`"
   (and add after the definition of `mass_h(J)`: "for `s < 0` the identity `c_h = mu_hat_h(2^s)` holds but
   `mass_h(h)` as defined is `0`, so the lemma is stated for `s >= 0`").
6. Lemma R' checks: "at `h = 80`, `s = 120`, worst `0.810`" -> "at `h = 80`, `s = 120`, worst `0.810` over `J <= 30`
   (the DP's range; over all `J` the worst is `0.904` at `J = 64`, and `0.913` at `J = 64` for `h = 120`); the
   inequality was also checked at every `h <= 81` and every `0 <= s <= floor(h log_2 3)`, worst ratio `0.918`
   (audit 4)".
7. Theorem C statement: "for all `h >= 1` and all integers `s = hm - δ`, `δ >= 0`" -> "for all `h >= 1` and all
   integers `0 <= s = hm - δ`, `0 <= δ <= hm`"; result paragraph: "for every `h` and every `s <= hm`" -> "for every
   `h` and every `0 <= s <= hm`".
8. Theorem C paragraph: "`Ñ_n 3^(n/2) = 1.00, ..., 0.67` for `n = 1..12`, then in `[0.12, 1.56]` to `n = 80` (max
   `1.56` at `n = 16`, `1.29` for `n >= 20`)" -> "... then in `[0.055, 1.56]` to `n = 80` (max `1.56` at `n = 16`,
   min `0.055` at `n = 78`; `<= 1.29` for `20 <= n <= 80`)"; "the sup over `m <= 60` of the single coefficients has
   `N_n 3^(n/2) ∈ [0.77, 3.5]` for `n <= 80`" -> "... `∈ [0.67, 3.77]` for `n <= 80`". Result paragraph: "`Ñ_n
   3^(n/2) ∈ [0.12, 1.56]` for `n <= 80`" -> "`∈ [0.055, 1.56]`".
9. Theorem C paragraph, replace "For `h <= 81` every `Ñ_(J-1)` in the bound is a computed number, so the inequality
   `M_res(h) <= 0.18 e^(-hI)` there rests on Lemma R' and float64 arithmetic only." with "For `h <= 150` every
   `Ñ_(J-1)` that matters in the bound is a computed number (the terms with `J - 1 > 80` carry at most `Σ_(J>=82)
   mass_h(J) <= 6.7·10^(-7) e^(-hI)` at `h <= 150` even with `Ñ <= 1`), so the computed inequalities
   `|mu_hat_h(2^(s*))| <= B_h(s*) <= 0.213 e^(-hI)` (`40 <= h <= 81`; `<= 0.164 e^(-hI)` for `82 <= h <= 150` with
   `Ñ_n` to `n = 120`) and `max_(0<=s<=hm) |mu_hat_h(2^s)| <= max_s B_h(s) <= 0.872 e^(-hI)` (`40 <= h <= 81`; the
   maximum of `B_h(s)` sits at `s = floor(hm)`, where `e^(θ*δ) ≈ 1`) rest on Lemma R', the exact masses and
   float64 arithmetic only."
10. Result paragraph: "and the bound `Σ_J mass_h(J) Ñ_(J-1)` equals `0.18, 0.14, 0.11 e^(-hI)` at `h = 40, 100, 300`
    (`0.08 e^(-hI)` at `h = 2000` under the hypothesis with the observed constant)" -> "and the bound `Σ_J mass_h(J)
    Ñ_(J-1)` equals `0.18, 0.14 e^(-hI)` at `h = 40, 100` (computed) and `0.11, 0.08 e^(-hI)` at `h = 300, 2000`
    under the hypothesis with an extrapolated constant (for `h >= 200` the terms with `J - 1 > 80` are not
    negligible without H)". Section-7 row: "the bound `Σ_J mass_h(J) Ñ_(J-1) = 0.11–0.18 e^(-hI)` at `h = 40..300`
    ... | VERIFIED" -> "... `= 0.14–0.18 e^(-hI)` at `h = 40..150` | VERIFIED; `0.11–0.12 e^(-hI)` at `h = 200..300`
    | under H".
11. Corridor paragraph: "on the *floor side* `k = nm + d` they are much larger, `≈ 2^(-d) M(n)`" -> "on the *floor
    side* `k = floor(nm) + d` they are much larger, `≈ 0.1 · 2^(-d) M(n)` (the ratio to `2^(-d) M(n)` is `0.06–0.14`
    for `-2 <= d <= 12` at `n = 40, 80`)"; "(`7.7·10^(-6), 2.5·10^(-6), ..., 4.6·10^(-9)` for `d = 1..9` at `n =
    80`; `4·10^(-12)` at `d ∈ [21, 60]`" -> "(`7.7·10^(-6), 2.5·10^(-6), ..., 4.6·10^(-9)` for `d = 0..8` at `n =
    80`; max `4·10^(-12)`, rms `7·10^(-13)` at `d ∈ [21, 60]`". Section-7 row: "floor side `≈ 2^(-d) M(n)`" ->
    "floor side `≈ 0.1·2^(-d) M(n)`".
12. Corridor paragraph: "drifts *up* relative to the line at `2 - m = 0.415` bits per level" -> "moves away from the
    critical line toward the ceiling `κ = 0` at `2 - m = 0.415` bits per level".
13. `J`-series paragraph: "carry `2.4%, 7.8%` of the mass and `< 2.3%, < 1.3%` of the coefficient at `h = 40, 80`
    (the `F = 0` words give `100.4%, 99.65%` of it" -> "carry `2.4%, 2.0%` of the mass; their largest single term is
    `2.3%, 1.3%` of the coefficient, their vector total `0.4%` (the `F = 0` words give `100.4%, 99.65%` of it; the
    sum of their moduli is `6.5%, 3.4%`".
14. Status paragraph: "the resonant window decays at the rate `3^(h*-1)` with a constant prefactor" -> "the resonant
    window decays at least at the rate `3^(h*-1)` (an upper bound) with a constant prefactor"; "VERIFIED — `Ñ_n` at
    the Parseval scale to level 80 (rate `0.574`), the bound term by term, `J_eff ≈ 1.6 √h`; the hypothesis H is
    CONJECTURAL and (T2) is reduced to one constant" -> "VERIFIED — `Ñ_n` at the Parseval scale to level 80 (rate
    `0.574`; extended to level 160 by audit 4: rate `0.572`, but `Ñ_n 3^(n/2)` reaches `8.7` at `n = 130`), the
    bound term by term; OBSERVED — `J_eff ≈ 1.6 √h`; the hypothesis H is CONJECTURAL, and (T2) reduces to one
    constant `Σ_J γ_J` provided the limits `γ_J` exist (DIRECTION)".
15. Theorem C paragraph, after "The hypothesis is what the data say (VERIFIED, not proved)": add "to level `80`;
    audit 4 extended `Ñ_n` to `n = 160`: `Ñ_n 3^(n/2) ∈ [0.029, 0.31]` for `81 <= n <= 122`, then `0.90, 3.37, 7.97,
    8.67, 2.38, 1.80` at `n = 124, 126, ..., 134` (a single coefficient wave travelling from `(n, m) = (84, 60)` to
    `(130, 1)` and growing from `1.0` to `13.4` times the Parseval scale, with generic one-step Gauss sums along it),
    then `0.03–0.35` at the even levels to `n = 160`; the least-squares rate over `20..160` is `0.5716`. H is not refuted, but 'at the
    Parseval scale' is a finite-range description and the constant `1.3` of the `h = 2000` extrapolation is not
    uniform." Section-7 row "the negative family at the Parseval scale ... | VERIFIED" -> append "to `n = 80`; not
    beyond `n = 122` (audit 4)".
16. "What changed" paragraph: "and the 'levels `180..300` lean against `M ≍ P_h`' of section 2c is withdrawn as
    evidence: the drift is the widening `J`-window, and Theorem C forbids any rate above `0.9465` under H" -> "and
    the 'levels `180..300` lean against `M ≍ P_h`' of section 2c is, under H, no longer evidence against the
    *exponential rate* `3^(h*-1)` (Theorem C forbids any rate above `0.9465` under H, so the fitted `0.947–0.950`
    must be a prefactor effect); it remains evidence about the *prefactor* (the rise of `M/P_h` by `1.26` over
    `180..300`), which Theorem C does not address, and without H the S21 reading stands". Make section 2c and the
    two section-7 rows consistent with this.
17. Result paragraph: "doubly exponentially better than the `C* h^(-6409)` of the Fourier–renewal method" ->
    "geometric with an explicit constant where the cited unconditional bound on all of `M(h)` is polynomial with a
    tower constant (the comparison is between a conditional bound on one family and an unconditional bound on the
    full level)".
18. "What the hypothesis is" paragraph: "The Fourier–renewal method's structured/unstructured dichotomy gives
    `(1-c)^n` at these frequencies with a `c` far too small; H needs the sharp constant." -> either delete, or
    "(From memory, not verified against the source:) the Fourier–renewal method would give at best a geometric
    bound with an unusably small `c` at these frequencies; H needs the sharp constant."
19. Result paragraph: "Term by term the inequality `|c_J| <= mass_h(J) Ñ_(J-1)` holds with worst ratio `0.84` (`h =
    40`) and `0.81` (`h = 80`)" -> "... worst ratio `0.84` (`h = 40`), `0.90` (`h = 80`), `0.91` (`h = 120`), and
    `0.92` over all `h <= 81` and all `0 <= s <= h log_2 3`".
20. Scripts: `collatz_five_mirrors_renewal_bound_20260929.py` docstring "(Ntilde beyond 80 extrapolated as 3.6 *
    3^(-n/2), the largest observed constant)" -> "1.292 * 3^(-n/2)" (the code's value).

## 3. Overall verdict

**SOUND WITH CORRECTIONS.** The mathematics of section 2d is correct: the 2-adic reading and the two-term bound of
Lemma G (an integer identity plus elementary trigonometry), the exponent-walk and real-phase identities, Lemma R'
(for `s >= 0`), the Chernoff mass law with the identity `e^(-I) = 3^(h*-1)`, and the derivation of Theorem C from
them all re-derive cleanly, and every number I recomputed independently (the law itself to `h = 10`, the closed
recursion to `n = 160`, the `c_J` by a complete walk DP to `h = 120`, the exact masses, the bound at every `h <= 150`
and every `s`) reproduces the note's values where the note computed them — and the bound `|m_h(s)| <= B_h(s)` holds
at every level and every exponent tested. The corrections are: three factual errors in Lemma G's prose (reversed
inequality, the resonant set, "exactly half"); the missing `s >= 0`; the false sentence "for `h <= 81` ... `M_res(h)
<= 0.18 e^(-hI)`" (undefined quantity, wrong constant across `h`, and only true at the argmax `s`); the
hypothesis-dependence of the bound values at `h >= 200`; the floor-side prefactor (`0.1`, not `1`); the `7.8%` and the
"`< 2.3%, < 1.3%`" mislabels; the ranges of `Ñ_n 3^(n/2)` and `N_n 3^(n/2)`; the overstatements in the status
paragraph (`J_eff` as VERIFIED, "(T2) reduced to one constant", "decays at the rate"); the unconditional wording of
the S21 withdrawal; the rhetorical and unsupported comparisons with the Fourier–renewal method. The one substantive
new fact is item 11: the negative family, whose Parseval-scale behaviour to `n = 80` is the entire evidence for H,
develops a coherent wave that lifts the weighted norm to `8.7` times the Parseval scale at `n = 130` before dying
out; H survives (rate `0.5716` over `20..160`), but the note's "`1.3%` margin" is not the relevant measure of its
risk, and the status of H should be recorded as "consistent with the data to `n = 160`, with one wave of amplitude
`13` at `n = 128..130`; mechanism OPEN".
