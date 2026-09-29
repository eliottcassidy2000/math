# Audit 4 (independent, adversarial) of section 2d (S22) of `collatz_five_mirrors_20260929.md`

**Object.** Section 2d ("the renewal structure of the resonant coefficient") as it stands at commit `19047a3ee`
(including the paragraphs added during this audit: the extension of `Ñ_n` to level 320, the ridges of the negative
family, the `γ_J` test, the floor mechanism, and the restated `h <= 81` sentence), the S22 rows of the section-7
status table, and the S22 clause of the status paragraph. Background read: sections 2, 2b, 2c, 4;
`collatz_five_mirrors_gap/renewal/renewal_bound/renewal_cert_20260929.py` and their `.out` files;
`collatz_five_mirrors_coherent_20260929.py`.
**Method.** Every claim was re-derived before the note's proof was read, and every number was recomputed with my own
code, `04-computation/experiments/collatz_five_mirrors_20260929_audit4.py` (output
`05-knowledge/results/collatz_five_mirrors_20260929_audit4.out`, 37 s), which imports none of the repo's S22
scripts: the law `mu_h` itself on `Z/3^h` for `h <= 10` (forward DP); a closed recursion written from the
definition with an exact shrinking window (truncations `a <= 40, 50, 60`; negative family on `[-60, 0]` to
`n = 160` and on `[-600, 0]` to `n = 320`); a walk DP over the words from the top level with state `(T_j, counter)`
and no window in the total cost (the repo's DP covers `94%` of the mass at `h = 80` and `63%` at `h = 120`, mine
`1 - 10^(-10)`), run at `h = 20, 40, 80, 120, 160, 200`, plus a joint `(J, floor-flag)` version with all phases and
with the top-part phases only; the exact cost laws; the bound `B_h(s) = Σ_J mass_h(J) Ñ_(J-1)` at every `h <= 150`
and every `0 <= s <= floor(h log_2 3)`; the no-descent probability `P_h`; exact integer checks of the ridge seeds;
and a 30-digit `mpmath` rerun of the closed recursion to isolate float64 rounding. Notation as in the note:
`m = log_2 3`, `θ* = ln(2(1 - log_3 2))`, `I = θ* m - ln(m - 1)`, `e^(-I) = 3^(h*-1)`.

## 1. Claims and verdicts

1. **Constants** (`0.73814, -0.30362, 0.054979, 0.946505, 1.70951, 0.98699`; `e^(-I) = 3^(h*-1)`; `Z(θ*) = m - 1`;
   tilted mean `m`; bracket `1 + 131.4 C`; `≈ 115 ln h`; `≈ 6000` levels; window width `0.61 √h`). **HOLDS** (output
   A): `e^(-I)` and `3^(h*-1)` agree to `3·10^(-16)`; the term-by-term identity `-I = -(1/p) ln(2q) + ln q - ln p =
   -ln p - (q/p) ln q - ln 3` holds (both `-0.0549794728`); the bracket is `1 + 131.4 C` (`= 206` at `C = 1.56`);
   the tilted standard deviation over `m` is `0.608`.

2. **Lemma G, the 2-adic reading** `(t 2^(-r) mod 3^n)/3^n = t/(2^r 3^n) + m_r/2^r`, `m_r = (-t 3^(-n)) mod 2^r`.
   **HOLDS**, as an integer identity: `2^r x = t + m 3^n` with `0 <= x < 3^n` forces `0 <= m < 2^r` (because
   `0 < t < 3^n`) and `m ≡ -t 3^(-n) mod 2^r`; checked over all units `n <= 6`, `r <= 20` (`14560` cases, `0`
   failures) and as floats to `5.0·10^(-16)` (the note says `4·10^(-16)`). The unit hypothesis is not needed.

3. **Lemma G, the two-term bound and the classes.** **HOLDS.** With `α = (θ + ξ_1)/2`, `β = (θ + ξ_1 + 2ξ_2)/4`:
   `|(1/2)e(α) + (1/4)e(β)|^2 = 5/16 + cos(2π(α - β))/4`, `α - β = -φ`, tail `<= 1/4`. `ξ_1 ≠ ξ_2` gives
   `cos(2πφ) ∈ {-cos(πθ/2), -sin(πθ/2)} < 0`, hence `|G| < (1 + √5)/4`. Class identification: `3^(-n) ≡ (-1)^n mod
   4`, so `ξ ≡ (-1)^(n+1) t mod 4`, and for odd `t`, `ξ ≡ 1 mod 4` iff `t ≡ (-1)^(n+1) mod 4` (checked for every
   unit, `n <= 9`, `0` failures). The bound was checked against every unit `n <= 9`: never violated (worst slack
   `2.7·10^(-3)`). The corollary "`2^(-m) mod 3^n` has a gap iff the digits `m+1, m+2` of `-3^(-n)` differ" **HOLDS**
   (`ξ(2^(-m)) = (-3^(-n)) >> m` digit-wise; `270` cases, `0` failures; a gap at `32, 29, 30` of the steps
   `m <= 60` at `n = 20, 40, 80`).

4. **Lemma G, the numerical prose.** (a) "the class `4 | t` with `θ <= 1/2` never exceeds `0.846`" — **FAILS**, the
   inequality is reversed: the class `4 | t, θ <= 1/2` *contains the resonance* (`0.9436, 0.9714, 0.9859, 0.9937,
   0.9967` at `n = 5..9`); the class `4 | t, θ > 1/2` (equivalently `t` odd, `ξ ≡ 3 mod 4`, `θ <= 1/2`) is the one
   that never exceeds `0.8462` (output B). (b) "the resonant set named: `t = 2^v u` small or `t = -(odd) small`" —
   **FAILS** for the second half: an odd `t` with `θ -> 1` is `t = 3^n - t'` with `t'` *even*, and the resonant
   class `t ≡ (-1)^n mod 4` forces `4 | t'`; the involution `t -> 3^n - t` maps `{v_2(t) >= 2}` onto `{t odd, ξ ≡
   3 mod 4}`. The units with `|G_9(t)| > 0.99` are exactly `±2^8, ±2^9, ±2^10, ±2^11, ±5·2^8`: the resonant set is
   `t = ±2^v u`, `v >= 2`, `2^v u` small (the section-7 row already says `±2^v u`; only the body text is wrong).
   (c) "exactly half the units (`6562` of `13122`)" — **HOLDS WITH CORRECTION**: `6562` gap against `6560` non-gap
   (the two gap classes have `3281` elements each, the two others `3280`); half up to one. "the two gap classes never
   exceed `0.672`" **HOLDS** (`0.6718`).

5. **The exponent-walk identity, the real-phase identity (`A <= s`), the Ramanujan average.** **HOLDS.**
   `2^s Y_h = Σ_j 3^(h-j) 2^(s-T_j)` in `Z_3` and `e(3^(h-j) x/3^h) = e(x/3^j)` for `x ∈ Z_3`; for `A <= s` the sum
   is the integer `2^(s-A) Σ_j 3^(h-j) 2^(P_(j-1))`. Brute force at `h = 5` (all `3125` words with letters `<= 5`,
   `s ∈ [-8, 20]`, negative `s` included): `2.6·10^(-15)`; `Y_h = Σ_j 3^(h-j) 2^(-T_j) mod 3^h` agrees with the
   forward recursion exactly; the real-phase identity to `0` (output C). Ramanujan: `2` is a primitive root mod
   `3^n`, so the average is `μ(3^n)/L_n` (`-1/2, 0, ..., 0` for `n = 1..7`).

6. **Lemma R' (statement and proof).** **HOLDS WITH CORRECTION: `s >= 0` is required.** All steps are correct: `T_j`
   is strictly decreasing, so `{j : T_j > s} = {1..J}` and `J(w) = J` iff `T_(J+1) <= s < T_J`; with `w = (v, a_J,
   u)` the exponents for `j > J` depend on `u` only, `κ_J = κ - a`, and for `j < J`, `κ_j = (κ - a) - T_j^(v)` is
   the exponent walk of `v` at level `J - 1` started at `κ - a < 0`, so `E_v Π_(j<J) ω_j = mu_hat_(J-1)(2^(κ-a))`
   by item 5; `v, a_J, u` are independent; the top-part phases (modulus `1`) are dropped by the triangle inequality;
   `Σ_(m>=1) 2^(-κ-m)|mu_hat_(J-1)(2^(-m))| = 2^(-κ) Ñ_(J-1)`, `Σ_κ P(T_(J+1) = s - κ) 2^(-κ) = mass_h(J)`; `J = 1`
   uses `Ñ_0 = 1`, `J = h` the empty `u` with `κ = s` and `mass_h(h) = 2^(-s)`. For `s < 0` the definition gives
   `mass_h(h) = P(0 <= s < a_h) = 0` while `c_h = mu_hat_h(2^s) ≠ 0`, so "for every integer `s`" is false; the
   lemma holds for `s >= 0`, which is all Theorem C uses.

7. **Lemma R' term by term.** **HOLDS**, more strongly than stated; one number corrected. With the complete walk DP
   (output F): `h = 40`: the ratios `0.016, 0.029, 0.053, 0.136, 0.273, 0.439, 0.171, 0.512, 0.732, 0.484, 0.468,
   0.206, 0.695` and the worst `0.840` (at `J = 20`) are reproduced exactly; `h = 80`: the worst over *all* `J` is
   `0.904` (`J = 64`), the note's `0.810` being the worst over `J <= 30` (the repo's DP counts `J <= 30` and its
   `c_J` for `J >= 21` are incomplete at `h = 80`); `h = 120`: `0.913` (`J = 64`); `h = 20`: `0.884`. `Σ_J c_J`
   agrees with the closed recursion to `4·10^(-12)` relative at `h = 20..200`, and with the law itself to
   `10^(-16)` at `h = 6, 8, 10` (output D). Beyond the note: `|m_h(s)| <= B_h(s)` was checked at **every**
   `h <= 81` and **every** `0 <= s <= floor(h m)`: worst ratio `0.918` (`h = 64, s = 0`); for `82 <= h <= 150`
   worst `0.109` (output H, H2). The exact masses agree with the DP masses to `5·10^(-12)`; the note's "`J <= 20`
   complete to `< 1%`" for the repo's DP holds.

8. **The Chernoff mass law** `mass_h(J) <= P(T_(J+1) <= s) <= e^(-θ* s) Z(θ*)^(h-J) = e^(-hI) e^(θ* δ) (m-1)^(-J)`.
   **HOLDS** (Chernoff at the fixed tilt `θ* < 0`; `Z(θ*) = q/(1-q) = m - 1`; tilted mean `m`; the algebra with
   `s = hm - δ` checked). Never violated numerically (output G, `h = 40..1000`). The normalised growth factors
   `1.35, 2.05, 1.90, 1.76, ...` (`h = 40`) and `1.72, 1.71, 1.71, 1.71, 1.70, ...` (`h = 1000`, after the first
   ratio `0.73`) are reproduced; the Gaussian-window reading (centre `≈ δ/m ≈ 4`, width `0.61 √h`) is a heuristic
   consistent with the data. The mass itself peaks at `J ≈ 0.21 h + 3` (`11, 19, 28, 65, 210` at `h = 40, 80, 120,
   300, 1000`).

9. **Theorem C.** **HOLDS WITH CORRECTION** (`0 <= s`). The derivation is correct (`Σ_J mass_h(J) Ñ_(J-1) <=
   e^(-hI) e^(θ*δ) [1 + C (m-1)^(-1)/(1 - ρ/(m-1))]`, convergent iff `ρ < m - 1`; `1 + 131.4 C` at `ρ = 3^(-1/2)`).
   H with `n = 0` forces `C >= 1`. The statement must read `0 <= s = hm - δ` (item 6); the "Result in one paragraph"
   sentence "for every `h` and every `s <= hm`" likewise. With the observed constants `C' ≈ 200`, a thousand times
   the computed finite-`h` bounds (the Gaussian window), so "a *constant* prefactor" must not be read as "`0.18`".

10. **`Ñ_n` to level 80.** **HOLDS WITH CORRECTION** of the stated ranges. My recursion reproduces the note's values
    (`1.000, 0.927, ..., 0.667` for `n = 1..12`; rate `0.5736` over `20..80`, `0.5574` over `40..80`; `<= 1.292` for
    `20 <= n <= 80`, max `1.562` at `n = 16`); the truncations `a <= 40 / 50 / 60` agree to `1.3·10^(-11)` /
    `1.3·10^(-14)`; float64 rounding is `1.1·10^(-14)` relative at `n = 40` and `9.7·10^(-15)` at `n = 60` against
    30-digit `mpmath` with the same truncation (output J). But the note read only the printed levels: the range of
    `Ñ_n 3^(n/2)` over `n <= 80` is `[0.055, 1.56]` (minimum at `n = 78`: `0.056, 0.055, 0.113` at `n = 77, 78,
    79`), not `[0.12, 1.56]`; and `N_n 3^(n/2) ∈ [0.77, 3.5]` is `[0.67, 3.77]` (`n = 77`, `n = 29`). The
    truncation is not self-certifying a priori (`2^(-50)` per level against `Ñ_80 = 1.6·10^(-20)`); the agreement
    of three truncations plus the `mpmath` check is adequate certification.

11. **`Ñ_n` to level 320 and the spike.** **HOLDS WITH CORRECTIONS.** On the window `[-600, 0]` (`a <= 50`,
    output E3; identical to the `[-60, 0]` run where they overlap) the least-squares rates are `0.5678` (`20..320`),
    `0.5710` (`160..320`), `0.5744` (`200..320`) — exactly the note's numbers — and `max_n (m-1)^(-n) Ñ_n = 1.580`
    at `n = 130` as stated. But the spike is `Ñ_n 3^(n/2) = 0.90, 2.06, 3.37, 4.44, 7.96, 7.47, 8.67, 3.76, 2.38,
    2.38, 1.80` at `n = 124..134`: the maximum is **`8.67` at `n = 130`**, not "`7.97` at `n = 128`" (the note's
    number is the `n = 128` value; its own "`1.58` at `n = 130`" is the same spike). Past `n = 140` the values lie
    in `[0.0014, 0.17]`, `97%` of the levels below `0.1` (the note's "`0.002–0.1` on most levels" is right). **The
    sentence "so H holds numerically with `C = 1.6`, `ρ = 0.585`, to level `320`" is a logical slip**: `ρ = 0.585 =
    m - 1` is the boundary value that Theorem C excludes (its series diverges there); the bound `Ñ_n <= 1.58
    (m-1)^n` gives nothing. What the data give is `Ñ_n <= 8.7 · 3^(-n/2)` (`ρ = 0.5774`, `C = 8.7`) or, with the
    fitted rate, `ρ = 0.5744` and `C ≈ 17`, for `n <= 320`. H is not refuted; it is consistent with the data to
    `n = 320` with those constants. The phrase "The hypothesis is VERIFIED to level `320`" should read "consistent
    with the data to level `320`" (a statement about all `n` cannot be verified), and the row "the negative family
    ... a spike `7.97` at `n = 128` | VERIFIED" needs `8.67` at `n = 130`.

12. **The ridges of the negative family.** **HOLDS as OBSERVED; the seeds VERIFIED exactly; two readings need
    qualification.** Exact integer checks (output E3): `v_3(55·2^423 + 1) = 15`, `v_3(13·2^154 - 1) = 7`,
    `2^(-423) ≡ 3^n - 55 mod 3^n` for `n = 9..15` and `≡ 3^15 - 55 mod 3^16`, `2^(-154) ≡ 13 mod 3^n` for
    `n <= 7`, `2^(-480) ≡ 2^6 mod 3^6`, `2^(-154) ≡ 2^8 mod 3^5`: all true. On my own surface `v_n(m) =
    3^(n/2)|m_n(-m)|`, `m <= 600`, `n <= 320`: the dominant ridge is born at `m* = 414, 413, 412, 410, 409` for
    `n = 11..15` with `v = 3.10, 4.40, 6.11, 8.37, 12.20`, peaks at `29.4` (`n = 29`), is `25.6, 14.0, 10.3` at
    `n = 34, 50, 70`, `6–8` at `n = 96`, `5.9, 4.0, 4.0` at `n = 115, 120, 125`, `1.8, 1.6, 1.2` at `n = 130, 140,
    160`, `0.65, 0.61, 0.27, 0.24` at `n = 165, 200, 250, 300` (the note's "`0.3–0.9` at `165..300`" is `0.2–0.7`
    on my track), and reaches `m = 43` at `n = 300` with `v = 0.24`; its mean slope is `-1.28` per level and the
    remnant decays at `0.5666` per level over `34..300` (note `0.5675`). At the birth, `|G_n(2^(-m*))| = 0.981,
    0.992, 0.996, 0.994, 0.997`, all in the class (`t` odd, `ξ ≡ 3 mod 4`, `θ = 0.84–0.94`), with runs of `11, 13,
    14, 16, 17` equal digits of `-3^(-n)` from positions `413, 411, 410, 408, 407` — every number of the note
    reproduced — and at `n = 16` the run and the gap failure are gone (`|G| = 0.61`), as the mechanism predicts.
    Two qualifications. (a) The attribution of the `n = 128..130` spike: the note says it is "the level-`5` ridge
    arriving at `m = 1..10` ... while the levels `112, 116, 120, 124, 128` have coherent small-exponent phases". The
    path is consistent with mine (global argmax `m = 21, 16, 10, 1` at `n = 115, 120, 125, 130`, on top of the whole
    600-window), but the *amplification* is not explained by those levels: along the path the normalised amplitude
    grows steadily from `1.0` (`n = 84, m = 60`) through `3.6, 5.0, 5.0, 4.1, 3.8, 7.5, 9.7, 10.9` (`n = 88..120`)
    to `13.4` (`n = 128`), i.e. `×1.06` per level over `44` levels, while the one-step Gauss sums at the ridge
    position are generic (`|G_n(2^(-m*))| = 0.17–0.74` for `n = 84..124`; `0.967` only at `n = 128`, where
    `3^128 ≡ 1 mod 2^9`). A remnant "decays at `0.5675` per level" by the note's own reading, so a remnant that
    *grows* against the Parseval scale for forty levels is an unexplained object; its slope `-1.3` is neither the
    resonance line of a fixed multiplier family (`-log_2 3 = -1.585`) nor the mean drift of the walk (`-2`). The
    mechanism of the growth should be marked OPEN. (b) The size heuristic (`Ñ_n <= C n^(1/2) 3^(-n/2)` "or so") is
    correctly labelled unproved; note that it concerns the seed amplitude, whereas the observed maxima of `Ñ_n
    3^(n/2)` (`1.56, 8.67` at `n = 16, 130`) come from ridges long after their seeds. The multiplier-family
    amplitude law (`u^(-0.438)`, scatter of four) was not recomputed.

13. **"For `h <= 81` every `Ñ_(J-1)` in the bound is a computed number ... (fourth audit)"** (the restated
    sentence). **HOLDS WITH TWO CORRECTIONS.** Reproduced: `B_h(s*)/e^(-hI) ∈ [0.135, 0.213]` for `40 <= h <= 81`
    (max at `h = 44`), `0.253` at `h = 20`, `max_s B_h(s)/e^(-hI) = 0.68–0.95` for `20 <= h <= 81` (max `0.950` at
    `h = 24`; `0.872` over `40..81`), attained at `s = floor(hm)`; so "`|mu_hat_h(2^s)| <= 0.95 e^(-hI)` for `20 <=
    h <= 81`" is right. But (a) "above `e^(-hI)` for `h <= 19`" is wrong: `B_h(s*)/e^(-hI) = 1.057, 0.998, 0.89,
    0.77, 0.65, 0.55, 0.60, 0.51, 0.43, 0.47, 0.40, 0.33` for `h = 1..12` and `0.29` at `h = 16` — above `e^(-hI)`
    only at `h = 1`; (b) "`<= 0.22 e^(-hI)` at the resonant exponents" holds for `40 <= h <= 81`; for `20 <= h <=
    39` the maximum is `0.283`. Also "all `s <= h log_2 3`" must be "all `0 <= s <= h log_2 3`" (item 6). Added
    here: with `Ñ_n` computed to `n = 120`, `B_h(s*) <= 0.164 e^(-hI)` for `82 <= h <= 150` (output H2), and the
    actual coefficients are far below (`|m_h(s*)|/B_h(s*) = 0.07 -> 0.019` over `40..150`).

14. **The bound values `0.18, 0.14, 0.11 e^(-hI)` at `h = 40, 100, 300` and the row `0.180, ..., 0.110`.** **HOLDS
    numerically WITH CORRECTION of status.** Reproduced (`0.1795, 0.1783, 0.1364, 0.1414, 0.1477, 0.1204, 0.1052,
    0.1235, 0.1102`). For `h >= 200` the value depends on `Ñ_n` beyond the computed range: without `Ñ_n` for
    `n > 80` the unconditional extra term `Σ_(J>=82) mass_h(J)` (`Ñ <= 1`) is `6.7·10^(-7) e^(-hI)` at `h = 150`
    (negligible) but `0.93 e^(-hI)` at `h = 200`, `3.4·10^3 e^(-hI)` at `250`, `10^6 e^(-hI)` at `300`; with `Ñ_n`
    to `120` the tail `Σ_(J>=122) mass_h(J)` is `2·10^(-14), 4.8·10^(-6), 1.48 e^(-hI)` at `h = 200, 250, 300`.
    So the `h = 200, 250, 300` entries are under the hypothesis (through the extrapolated `Ñ_n`), like the `h =
    2000` one the note does label; the spike of item 11 does not change them (the mass at `J - 1 ∈ [124, 134]` is
    `< 10^(-7)` at `h = 300`) but would enter at `h ≈ 600` with a factor up to `7` above the extrapolation `1.3
    3^(-n/2)`. The `renewal_bound` docstring says the extrapolation constant is `3.6`; the code uses `1.292`.

15. **The profile (C1).** Ceiling side **HOLDS** (rms `2.23·10^(-10)` at `n = 40` against `2.87·10^(-10)`;
    `2.22·10^(-20)` at `n = 80` against `8.2·10^(-20)`; argmax `k = nm - 6.40, -6.80`). Floor side "`≈ 2^(-d)
    M(n)`" **FAILS quantitatively**: `|m_n(k_c + d)|/M(n) = (0.06–0.14) · 2^(-d)` for every `-2 <= d <= 12` at
    `n = 40` and `80` (output I) — the decay `2^(-d)` is right, the prefactor is `≈ 0.1`; the listed values
    (`7.7·10^(-6), ..., 4.6·10^(-9)` "for `d = 1..9`") are those at `k = floor(nm) + 0..8`, and "`4·10^(-12)` at
    `d ∈ [21, 60]`" is the maximum (rms `7.4·10^(-13)`). The low side of the peak decays at `0.71–0.94` per step,
    steepening away from it. The three-regime Gauss-sum means (C4) were not recomputed.

16. **`J_eff = 10, 14, 19` (`≈ 1.6 √h`).** **HOLDS** with caveats: reproduced as the first `J` with remainder
    `< 10%`; the remainder is not monotone (`h = 80`: `0.111` at `J = 7`, `0.074` at `14`, `0.104` at `16`), so the
    stable-10% criterion gives `10, 17, 19` and the stable-5% `10, 18, 19`; at `h = 20`, `J_eff = 8`. With four
    points `1.6 √h = 7.2, 10.1, 14.3, 17.5` fits better than `h/8` and than a two-parameter log law (`8.2 ln h -
    20 = 4.6, 10.2, 15.9, 19.3`); the reading remains OBSERVED (the status paragraph lists it under VERIFIED).

17. **The per-mass weights `|c_J|/mass_h(J) · h`, the `J`-series table, the arguments.** **HOLD** exactly (output
    F; at `h = 120` from a complete DP). The `h = 20` row `1.59, 1.75, 1.71, 1.44, 0.28, 0.40, 0.21, 0.08, 0.03`
    (`J = 2..10`) continues the pattern.

18. **The `γ_J` test (`P_h/(h^(-3/2) e^(-hI)) = 6.8, 8.1, 9.0, 9.2, 9.2`; `|mu_hat_h(2^(s*))|/P_h = 0.463, 0.462,
    0.443, 0.462, 0.500`; `3.16, 3.76, 3.97, 4.26, 4.60`; the moduli and arguments of `γ_J` at `J = 4, 5, 8, 12,
    17`).** **HOLDS**: every number reproduced (output F2: `6.81, 8.12, 8.97, 9.22, 9.20`; `3.16, 3.75, 3.97,
    4.26, 4.60`; `|γ_5| = 2.563, 2.433, 2.337, 2.279, 2.224`; `|γ_17| = 0.015, 0.573, 1.804, 3.235, 4.528`; `arg
    γ_5 = -0.86, -0.67, -0.47, -0.50, -0.54`; `arg γ_12 = +0.82, +1.04, +1.28, +1.28, +1.27`). One wording:
    "the no-descent probability has *exactly* the `h^(-3/2)` prefactor, constant `≈ 9.2`" is a five-point
    observation still rising (`6.8 -> 9.2`); "consistent with an `h^(-3/2)` prefactor" is the supportable form. The
    conclusion "the constant `Σ_J γ_J` is not computable from `h <= 200`" is right and honest.

19. **The floor mechanism (top-part coherence, ballot fractions).** **HOLDS WITH CORRECTION of the `h = 120`
    numbers.** My joint DP (output F3) gives the top-only coherence `|W|/mass = 0.882, 0.900, 0.864` for `F = 0`
    and `0.0770, 0.0623, 0.0626` for `F >= 1` at `h = 40, 80, 120` (note: `0.88, 0.90, 0.83` and `0.077, 0.062,
    0.055` — the `h = 120` values of the repo come from its `63%`-mass DP); `|W(F>=1)|/|W(F=0)| = 0.0022, 0.0014,
    0.0022` (note `0.002, 0.0015, 0.003`); every `J = 0` word carries a floor level (mass without one `= 0`) with
    coherence `0.0162, 0.0039, 0.0015` (note `0.016, 0.004, 0.002`); the `F = 0` fractions of the classes `J = 1..4`
    times `h` are `8.10, 8.61, 7.65 / 19.57, 23.31, 21.80 / 27.71, 36.25, 35.54 / 32.49, 46.04, 47.13` (note `8.1,
    8.6, 7.7 / 19.6, 23.3, 21.8 / 27.7, 36.3, 35.5 / 32.5, 46.0, 47.1`); the `(J = 4, F = 0)` coherence is `0.500,
    0.343, 0.262` (note `0.50, 0.34, 0.26`). "the `F >= 1` coherence is `<= 0.06, <= 0.013, <= 0.005`" holds for
    `J <= 4` only (`0.16, 0.027, 0.010` at `J = 8`). The mass with a floor level is `2.42%, 1.97%, 3.00%`. The
    independence statement "conditional on the crossing exponent the top and bottom parts are independent" is the
    factorisation of Lemma R' and is correct.

20. **The floor counter in the `J`-series paragraph.** "carry `2.4%, 7.8%` of the mass and `< 2.3%, < 1.3%` of the
    coefficient at `h = 40, 80`" — `2.4%` **HOLDS**; `7.8%` **FAILS** (the `F >= 1` mass at `h = 80` is `1.97%`;
    the note's `7.8% = 1 - 0.9218` counts the `5.9%` of the mass outside the repo DP's cost window, deep-ceiling
    words with `F = 0`); "`< 2.3%, < 1.3%`" is the *largest single term* (`|c_(F=1)|/|full| = 2.27%, 1.32%`), not
    the total: the vector remainder `|full - c_(F=0)|/|full|` is `0.42%, 0.38%` (so the result paragraph's "`< 1%`"
    and "`100.4%, 99.65%`" **HOLD**) and the sum of moduli is `6.5%, 3.4%`.

21. **The drift sentence** "starts `δ` bits below the critical line and drifts *up* relative to the line at `2 - m =
    0.415` bits per level" — **wording/sign**: `κ_(j-1) - (j-1)m = (κ_j - jm) + (m - a_(j-1))`, so `κ_j - jm`
    *decreases* by `2 - m` per level on the way down (the walk moves away from the floor toward the ceiling, which
    is what the rest of the paragraph uses; `J ≈ 0.21 h` **HOLDS**).

22. **Statuses and framing.** (a) Status paragraph: "the resonant window decays at the rate `3^(h*-1)` with a
    constant prefactor" — Theorem C gives an *upper* bound at that rate under H; "at least at the rate". (b)
    "VERIFIED — ... `J_eff ≈ 1.6 √h`" — OBSERVED (item 16). (c) "(T2) is reduced to one constant" — the reduction
    presupposes that the limits `γ_J` exist (obligation (2), DIRECTION in the body; item 18 shows only `J <= 5`
    converge by `h = 200`); in the status paragraph it reads as a result. (d) "the 'levels `180..300` lean against
    `M ≍ P_h`' of section 2c is withdrawn as evidence: ... Theorem C forbids any rate above `0.9465` under H" — the
    withdrawal is justified **only under H** (CONJECTURAL) and **only for the exponential rate**; `M ≍ P_h` also
    asserts the prefactor (`h^(-3/2)`), which Theorem C does not address (it allows `M/P_h` to grow polynomially),
    and the note's own `γ` test shows `|mu_hat|/(h^(-3/2) e^(-hI))` still rising (`3.16 -> 4.60`) — exactly a
    prefactor statement; moreover section 2c and the two section-7 rows still say "lean slightly against it", so
    the note is internally inconsistent. (e) "doubly exponentially better than the `C* h^(-6409)` of the
    Fourier–renewal method" compares a *conditional* bound on *one family* with the cited *unconditional* bound on
    all of `M(h)`; "doubly exponentially" is rhetoric. (f) "The Fourier–renewal method's structured/unstructured
    dichotomy gives `(1-c)^n` at these frequencies with a `c` far too small" is **not supported by the cited
    material** (Tao is cited from memory via Mazur's (2.3); Mazur's paragraph gives the polynomial bound on the
    primitive coefficients; nothing cited says what the method yields on the negative powers of two). (g) The
    result paragraph's "`Ñ_n 3^(n/2) ∈ [0.12, 1.56]` for `n <= 80`" and "worst ratio `0.84` (`h = 40`) and `0.81`
    (`h = 80`)" carry items 10 and 7. (h) The row "the bound `= 0.11–0.18 e^(-hI)` at `h = 40..300` | VERIFIED"
    conflates computed (`h <= 150`) and hypothesis-dependent (`h >= 200`) values (item 14); the row "the profile
    asymmetry ... floor side `≈ 2^(-d) M(n)`" needs the prefactor `0.1` (item 15); the row "the floor mechanism:
    top-part coherence `0.88–0.90` without a floor level, `0.055–0.077` with one" needs `0.86–0.90` and
    `0.062–0.077` (item 19); the row "Lemma R' ... term by term VERIFIED at `h = 40, 80` (worst ratio `0.84`)"
    should record `0.90` at `h = 80` and `0.92` over all `h <= 81`, all `s`. (i) The reading "H is a
    square-root-cancellation statement ... of the `×2 ×3` kind" and "Lemma G ... does not by itself control the
    product" are fair and correctly labelled; item 12(a) shows the product does develop coherent structures that
    one-step gaps do not see. No Collatz proof step is claimed, correctly.

## 2. Textual corrections to apply (exact current phrase -> replacement)

1. "Checked against the modular sum to `4·10^(-16)` (`n <= 5`)" -> "Checked against the modular sum to
   `5·10^(-16)` (`n <= 5`)".
2. "the class `4 | t` with `θ <= 1/2` never exceeds `0.846`" -> "the class `4 | t` with `θ > 1/2` (equivalently
   `t` odd, `ξ ≡ 3 mod 4`, `θ <= 1/2`) never exceeds `0.846`".
3. "resonant set named: `t = 2^v u` small or `t = -(odd) small`" -> "resonant set named: `t = ±2^v u` with
   `v >= 2` and `2^v u` small (an odd resonant `t` is `t = 3^n - 4u'` with `4u'` small); at `n = 9` the units with
   `|G_9| > 0.99` are `±2^8, ±2^9, ±2^10, ±2^11, ±5·2^8`".
4. "The gap classes are exactly half the units (`6562` of `13122` at `n = 9`)" -> "The gap classes are half the
   units up to one (`6562` against `6560` at `n = 9`; the involution `t -> 3^n - t` pairs `{v_2(t) = 1}` with
   `{t odd, ξ ≡ 1 mod 4}` and `{4 | t}` with `{t odd, ξ ≡ 3 mod 4}`)".
5. Lemma R': "For every `h >= 1` and every integer `s`" -> "For every `h >= 1` and every integer `s >= 0`" (for
   `s < 0` the identity `c_h = mu_hat_h(2^s)` holds but `mass_h(h)` as defined is `0`).
6. "at `h = 80`, `s = 120`, worst `0.810`" -> "at `h = 80`, `s = 120`, worst `0.810` over `J <= 30` (`0.904` at
   `J = 64` over all `J`; `0.913` at `h = 120`); checked at every `h <= 81` and every `0 <= s <= floor(h log_2 3)`
   with worst ratio `0.918` (audit 4)". Result paragraph: "worst ratio `0.84` (`h = 40`) and `0.81` (`h = 80`)" ->
   "worst ratio `0.84` (`h = 40`), `0.90` (`h = 80`), `0.92` over all `h <= 81` and all `0 <= s <= h log_2 3`".
   Section-7 row: "(worst ratio `0.84`)" -> "(worst ratio `0.84` at `h = 40`, `0.90` at `h = 80`; `0.92` over all
   `h <= 81` and all exponents)".
7. Theorem C: "for all `h >= 1` and all integers `s = hm - δ`, `δ >= 0`" -> "for all `h >= 1` and all integers
   `0 <= s = hm - δ`"; result paragraph "for every `h` and every `s <= hm`" -> "for every `h` and every
   `0 <= s <= hm`"; the restated sentence "for `20 <= h <= 81` and all `s <= h log_2 3`" -> "for `20 <= h <= 81`
   and all `0 <= s <= h log_2 3`".
8. Three places: "`Ñ_n 3^(n/2) ∈ [0.12, 1.56]` for `n <= 80`" (result paragraph and section-7 row) and "then in
   `[0.12, 1.56]` to `n = 80` (max `1.56` at `n = 16`, `1.29` for `n >= 20`)" -> "`[0.055, 1.56]`" with "(max
   `1.56` at `n = 16`, min `0.055` at `n = 78`; `<= 1.29` for `20 <= n <= 80`)"; "`N_n 3^(n/2) ∈ [0.77, 3.5]` for
   `n <= 80`" -> "`∈ [0.67, 3.77]`".
9. "spikes to `7.97` at `n = 128` (`N_n 3^(n/2) = 13.4` there" -> "spikes to `8.67` at `n = 130` (`7.96` at
   `n = 128`, `N_n 3^(n/2) = 13.4` there"; section-7 row "a spike `7.97` at `n = 128`" -> "a spike `8.67` at
   `n = 130`".
10. "so H holds numerically with `C = 1.6`, `ρ = 0.585`, to level `320`" -> "so `Ñ_n <= 1.58 (m-1)^n` to level
    `320`; H needs `ρ < m - 1`, and the data give `Ñ_n <= 8.7 · 3^(-n/2)` (`ρ = 0.5774`, `C = 8.7`) or `ρ =
    0.5744` with `C ≈ 17` to level `320`".
11. "The hypothesis is VERIFIED to level `320`:" -> "The hypothesis is consistent with the data to level `320`:";
    status paragraph "VERIFIED — `Ñ_n` below the critical rate to level 320 (rate `0.568–0.574`, constant `1.58`)"
    -> "VERIFIED — `Ñ_n` below the critical rate to level 320 (least-squares rate `0.568–0.574`; `Ñ_n <= 8.7 ·
    3^(-n/2)`, the maximum `8.67` at `n = 130`)".
12. Ridge paragraph: "The spike of `Ñ_n` at `n = 128` is the level-`5` ridge arriving at `m = 1..10` (its path
    `m ≈ 155 - 1.3 (n - 5)`: `m* = 63, 48, 24, 11` at `n = 84, 96, 112, 124`) while the levels `112, 116, 120, 124,
    128` have coherent small-exponent phases, since ... Lemma G's gap fails there." -> append: "; but this does
    not account for the amplification: along the path the normalised amplitude grows from `1.0` at `n = 84` to
    `13.4` at `n = 128` (`×1.06` per level over forty levels) with generic one-step Gauss sums at the ridge
    position (`0.17–0.74` for `n = 84..124`), i.e. a remnant that grows instead of decaying at `0.5675`; the
    mechanism of this growth is OPEN (audit 4)". Also "`0.3–0.9` at `n = 165..300`" -> "`0.2–0.9` at `n =
    165..300`".
13. Restated `h <= 81` sentence: "(`0.25` at `h = 20`, above `e^(-hI)` for `h <= 19`)" -> "(`0.25` at `h = 20`,
    `0.28` for `20 <= h <= 39`, rising to `1.06 e^(-hI)` at `h = 1`, above `e^(-hI)` only there)"; "and `<= 0.22
    e^(-hI)` at the resonant exponents (fourth audit)" -> "and `<= 0.22 e^(-hI)` at the resonant exponents for
    `40 <= h <= 81` (`<= 0.29` for `20 <= h <= 39`; `<= 0.17` for `82 <= h <= 150` with `Ñ_n` to `n = 120`)
    (fourth audit)".
14. Result paragraph: "and the bound `Σ_J mass_h(J) Ñ_(J-1)` equals `0.18, 0.14, 0.11 e^(-hI)` at `h = 40, 100,
    300` (`0.08 e^(-hI)` at `h = 2000` under the hypothesis with the observed constant)" -> "and the bound `Σ_J
    mass_h(J) Ñ_(J-1)` equals `0.18, 0.14 e^(-hI)` at `h = 40, 100` (computed) and `0.11, 0.08 e^(-hI)` at `h =
    300, 2000` under the hypothesis with an extrapolated constant (for `h >= 200` the terms with `J - 1 > 80` are
    not negligible without H)". Section-7 row: "the bound `Σ_J mass_h(J) Ñ_(J-1) = 0.11–0.18 e^(-hI)` at `h =
    40..300` ... | VERIFIED" -> "... `= 0.14–0.18 e^(-hI)` at `h = 40..150` (computed), `0.11–0.12 e^(-hI)` at `h =
    200..300` (under H) ... | VERIFIED / under H".
15. "on the *floor side* `k = nm + d` they are much larger, `≈ 2^(-d) M(n)` (`7.7·10^(-6), 2.5·10^(-6), ...,
    4.6·10^(-9)` for `d = 1..9` at `n = 80`; `4·10^(-12)` at `d ∈ [21, 60]`" -> "on the *floor side* `k = floor(nm)
    + d` they are much larger, `≈ 0.1 · 2^(-d) M(n)` (the ratio to `2^(-d) M(n)` is `0.06–0.14` for `-2 <= d <= 12`
    at `n = 40, 80`; `7.7·10^(-6), 2.5·10^(-6), ..., 4.6·10^(-9)` for `d = 0..8` at `n = 80`; max `4·10^(-12)`, rms
    `7·10^(-13)` at `d ∈ [21, 60]`"; section-7 row "floor side `≈ 2^(-d) M(n)`" -> "floor side `≈ 0.1 · 2^(-d)
    M(n)`".
16. "drifts *up* relative to the line at `2 - m = 0.415` bits per level" -> "moves away from the critical line
    toward the ceiling `κ = 0` at `2 - m = 0.415` bits per level".
17. "carry `2.4%, 7.8%` of the mass and `< 2.3%, < 1.3%` of the coefficient at `h = 40, 80` (the `F = 0` words give
    `100.4%, 99.65%` of it" -> "carry `2.4%, 2.0%` of the mass; their largest single term is `2.3%, 1.3%` of the
    coefficient and their vector total `0.4%` (the `F = 0` words give `100.4%, 99.65%` of it; the sum of their
    moduli is `6.5%, 3.4%`".
18. Floor-mechanism paragraph: "`0.88, 0.90, 0.83` for `F = 0` and `0.077, 0.062, 0.055` for `F >= 1` at `h = 40,
    80, 120`" -> "`0.88, 0.90, 0.86` for `F = 0` and `0.077, 0.062, 0.063` for `F >= 1` at `h = 40, 80, 120`
    (complete DP)"; "and the `F >= 1` coherence is `<= 0.06`, `<= 0.013`, `<= 0.005`" -> "and the `F >= 1`
    coherence is `<= 0.06`, `<= 0.013`, `<= 0.005` for `J <= 4` (`0.16, 0.027, 0.010` at `J = 8`)"; section-7 row
    "`0.88–0.90` without a floor level, `0.055–0.077` with one" -> "`0.86–0.90` without a floor level, `0.062–0.077`
    with one".
19. "`(the no-descent probability has exactly the `h^(-3/2)` prefactor, constant `≈ 9.2`)`" -> "(consistent with an
    `h^(-3/2)` prefactor, the constant still rising slowly, `6.8 -> 9.2`)".
20. Status paragraph: "the resonant window decays at the rate `3^(h*-1)` with a constant prefactor" -> "the resonant
    window decays at least at the rate `3^(h*-1)` (an upper bound) with a constant prefactor"; "the bound term by
    term, `J_eff ≈ 1.6 √h`; the hypothesis H is CONJECTURAL and (T2) is reduced to one constant" -> "the bound term
    by term; OBSERVED — `J_eff ≈ 1.6 √h`; the hypothesis H is CONJECTURAL, and (T2) reduces to one constant `Σ_J
    γ_J` provided the limits `γ_J` exist (DIRECTION; only `J <= 5` have converged by `h = 200`)".
21. "and the 'levels `180..300` lean against `M ≍ P_h`' of section 2c is withdrawn as evidence: the drift is the
    widening `J`-window, and Theorem C forbids any rate above `0.9465` under H" -> "and the 'levels `180..300` lean
    against `M ≍ P_h`' of section 2c is, under H, no longer evidence against the *exponential rate* `3^(h*-1)`
    (Theorem C forbids any rate above `0.9465` under H, so the fitted `0.947–0.950` must be a prefactor effect); it
    remains evidence about the *prefactor* (`M/P_h` rising by `1.26` over `180..300`, `|mu_hat|/(h^(-3/2) e^(-hI))`
    rising `3.16 -> 4.60` to `h = 200`), which Theorem C does not address, and without H the S21 reading stands".
    Make section 2c and the two section-7 rows consistent with this.
22. "doubly exponentially better than the `C* h^(-6409)` of the Fourier–renewal method" -> "geometric with an
    explicit constant where the cited unconditional bound on all of `M(h)` is polynomial with a tower constant (a
    conditional bound on one family against an unconditional bound on the full level)".
23. "The Fourier–renewal method's structured/unstructured dichotomy gives `(1-c)^n` at these frequencies with a `c`
    far too small; H needs the sharp constant." -> delete, or "(From memory, not verified against the source:) the
    Fourier–renewal method would give at best a geometric bound with an unusably small `c` at these frequencies; H
    needs the sharp constant."
24. `collatz_five_mirrors_renewal_bound_20260929.py` docstring "(Ntilde beyond 80 extrapolated as 3.6 * 3^(-n/2),
    the largest observed constant)" -> "1.292 * 3^(-n/2)" (the code's value).

## 3. Overall verdict

**SOUND WITH CORRECTIONS.** The mathematics of section 2d is correct: the 2-adic reading and the two-term bound of
Lemma G (an integer identity plus elementary trigonometry), the exponent-walk and real-phase identities, Lemma R'
(for `s >= 0`), the Chernoff mass law with `e^(-I) = 3^(h*-1)`, and the derivation of Theorem C all re-derive
cleanly; and every number recomputed independently — the law itself to `h = 10`, the closed recursion to `n = 320`
on a 600-wide window, the `c_J` by a complete walk DP to `h = 200`, the exact masses, the bound at every `h <= 150`
and every exponent, the ridge seeds, the `γ_J` and floor-mechanism tables — reproduces the note's values where the
note computed them (the ridge birth data and the `γ_J` table to the last digit), and the inequality `|m_h(s)| <=
B_h(s)` holds at every level and exponent tested. The corrections are: three factual errors in Lemma G's prose
(reversed inequality, the resonant set, "exactly half"); the missing `s >= 0`; the boundary value `ρ = 0.585` in
"H holds numerically" (Theorem C needs `ρ < m - 1`); the spike (`8.67` at `n = 130`, not `7.97` at `128`); the two
misreadings in the restated `h <= 81` sentence ("above `e^(-hI)` for `h <= 19`", "`<= 0.22`" beyond `40 <= h <=
81`); the hypothesis-dependence of the bound values at `h >= 200`; the floor-side prefactor (`0.1`); the `7.8%`
and "`< 2.3%, < 1.3%`" mislabels; the `h = 120` floor-mechanism numbers from an incomplete DP; the ranges of
`Ñ_n 3^(n/2)` and `N_n 3^(n/2)`; the overstatements in the status paragraph (`J_eff` as VERIFIED, "(T2) reduced to
one constant", "decays at the rate", "VERIFIED to level 320"); the unconditional wording of the S21 withdrawal; the
rhetorical and unsupported comparisons with the Fourier–renewal method. The substantive open point is item 12(a):
the negative family's dominant object at `n = 84..130` is a ridge that *grows* against the Parseval scale by `×1.06`
per level for forty levels with generic one-step Gauss sums, which the ridge mechanism (exact at its seeds) does
not explain; H survives it with `C = 8.7` at `ρ = 3^(-1/2)`, but the note should record the growth as OPEN rather
than attribute it to the coherent leading digits at `n = 112..128`.
