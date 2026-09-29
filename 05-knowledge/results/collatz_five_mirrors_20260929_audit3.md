# Independent audit (3) of section 2c of `collatz_five_mirrors_20260929.md`: the recursion to level 300, the excursion split, the random units

**Auditor:** independent session (Fable), 2026-09-29. **Audited state:** the note at commit `d61c39bc4` (line
numbers below refer to it): section 2c (lines 493–559), §0 item 1 (lines 85–100), the status rows 789–790 and
the "Next probes" sentence (lines 885–889); the scripts `collatz_five_mirrors_powers_of_two_20260929.py` (run at
300), `..._rate300_20260929.py`, `..._coherent_20260929.py`, `..._coherent120_20260929.py` and their outputs. The
untracked `collatz_five_mirrors_excursion_profile_20260929.py/.out` is not part of the commit and is used here
only as a cross-check of the band widths.

**Method.** `04-computation/experiments/collatz_five_mirrors_20260929_audit3.py` (61 s; output
`collatz_five_mirrors_20260929_audit3.out`; tags `[A1]`, `[C6]`, ... refer to its lines). Own numpy
implementation of the closed recursion to `h = 300` with `AMAX = 40` (the session's setting) and `AMAX = 60`
(certified), a reversed-summation run and a noise-injection run; own DP for `P_h` to 300; free, pinned and
fixed-`β` fits; an own DP for the excursion split, validated at `c = ∞` against the closed recursion; 100 random
units at `h = 30` and the exact coefficient distribution over all units at `h = 12, 14` from an own FFT law.

---

## (i) Claims checked

### (a) The rate to level 300

1. **The values to `h = 300` (`M(300) = 7.42·10^-11`, ratios, argmax).** VERIFIED. `[A1]`, `[A3]`: my
   `AMAX = 40` run agrees with the session's output at all 300 levels to `4.4·10^-7` (its seven printed digits);
   the argmax agrees except at `h = 5` (a `±2^s` tie inside a full period). `[A2]`: `AMAX = 60` agrees with
   `AMAX = 40` to `2.5·10^-12` relative at every level (argmax differences only at `h <= 9`, full-period ties).
2. **Is the level-300 value trustworthy under float64 and `a <= 40`?** YES, but not by the session's own
   setting. The dropped terms `sum_(a>40) 2^-a m_(n-1)(j-a)` involve exponents outside the computed window, for
   which only `|m| <= 1` is available a priori, so the rigorous truncation bound is `300·2^-40 = 2.7·10^-10 >
   M(300)`: `AMAX = 40` is not self-certifying beyond `h ≈ 250` `[A5]`. With `AMAX = 60` the truncation bound is
   `2.6·10^-16`, and the rounding bound is `50 eps sum_k maxwin|m_k| = 2.3·10^-14` (every summed term lies inside
   the window, so the per-level rounding is relative to the computed maximum, and errors add without
   amplification): `|error at 300| <= 2.4·10^-14`, i.e. `3.2·10^-4` relative, rigorous. The `AMAX = 40` values
   coincide with the certified ones to `2·10^-12`, so the session's numbers are right; the note should say the
   certification. Empirically the recursion averages errors down rather than accumulating them: a reversed
   summation order changes the values by `1.4·10^-14` and noise of `10^-12` injected at every level changes
   `M(300)` by `8·10^-12`, not `3·10^-10` `[A4]`. All conclusions of 2c survive unchanged.
3. **Doubling exponents `3.10, 5.27, 7.27, 9.19, 11.07, 12.97`; "slope `0.077` per level".** VERIFIED WITH
   CORRECTION. `[A6]` the exponents reproduce; the `0.077` is the two-point slope `(E(150→300) - E(50→100))/100`
   of the session's script; the least-squares slope over the six doublings is `0.0785` (endpoints `0.0790`),
   intercept `1.28`; `|log_2 3^(h*-1)| = 0.0793`; the slope corresponds to `r = 0.947`.
4. **Geometric-mean ratios of `M` (`0.9383, 0.9418, 0.9427, 0.9433`) and of `P_h` (`0.9375, 0.9404, 0.9413`).**
   VERIFIED `[A7]`, `[B1]`.
5. **Free fits over `100..300`: `M`: `β = 1.62, r = 0.9489` (residual `1.0%`); `P_h`: `β = 1.29, r = 0.9461`
   (`2.3%`); pinned to `0.9465`: `β = 1.14` (`5.9%`) and `1.38` (`2.5%`).** VERIFIED `[B3]` (`1.620, 0.94894,
   0.0100; 1.293, 0.94605, 0.0227; 1.144, 0.0586; 1.384, 0.0254`).
6. **"the exponential rate of `M` is `0.941–0.949`, bracketing `3^(h*-1)`" (lines 511–512; line 98).** VERIFIED
   WITH CORRECTION. The trade-off curve `[B3]` (`β` fixed, `r` fitted, max log-residual): `β = 1.00: 0.9458
   (8.1%)`, `1.25: 0.9471 (5.1%)`, `1.50: 0.9483 (2.1%)`, `1.75: 0.9496 (2.0%)`, `2.00: 0.9508 (4.4%)`; the pure
   geometric `0.9408` has residual `20%` `[B4]`. So the lower end `0.941` of the bracket is the pure-geometric
   misfit, and the fits with residual at most `2.5%` (the quality the pinned `P_h` fit achieves) give `r =
   0.948–0.950`; the doubling slope gives `0.947`; `3^(h*-1) = 0.9465` needs `β = 1.14` and fits six times worse
   than the free fit. For `P_h`, whose rate is `0.9465` by theorem, the free fit returns `0.9461` and pinning
   costs nothing (`2.3% → 2.5%`) `[B5]`: the same procedure applied to `M` returns `0.9489`. Honest reading: over
   `100..300` the exponential rate of `M` is `0.947–0.950` by every admissible estimator, `0.15–0.35%` above the
   no-descent rate, which sits at the lower edge; whether the gap is a prefactor effect or a genuine rate
   difference is undecidable from 300 levels. The same signal is the `M/P_h` drift (item 7).
7. **`M/P_h = 0.46 ± 0.02` up to `h = 180` (line 514), then `0.500, 0.519, 0.560, 0.579` at `200, 240, 280,
   300`; "`0.44–0.48` over `20..180`, `0.44–0.58` over `20..300`" (lines 91–93).** VERIFIED WITH CORRECTION.
   `[B2]`: range `0.443–0.481` over `20..180` (so "`± 0.02`" fails at `h = 41`, as audit 2 found; §0's
   "`0.44–0.48`" is right, line 514 is not), `0.443–0.584` over `20..300`; the values at `200..300` reproduce
   (`0.500, 0.519, 0.560, 0.579`). The rise by a factor `1.26` over `180 → 300` is what a rate difference of
   `0.2%` per level produces (item 6).
8. **"`s = h log_2 3 - 6.2 .. 7.1` throughout (`-6.49` at `h = 300`)" (lines 518–519).** VERIFIED WITH
   CORRECTION. `[A8]`: over all levels `20..300` the offset runs `-5.28 .. -7.22`; "`6.2..7.1`" is the range at
   the multiples of 25 that the rate script prints. The fact itself (a fixed fraction `2^s/3^h ≈ 2^-6` of the
   modulus) stands.
9. **"`M` and `P_h` have the same exponential order, with prefactors that differ by a slowly varying factor"
   (lines 515–517).** OBSERVED, correctly hedged as a description; but see item 6: the fits and the drift both
   point to a rate `0.002–0.003` above `P_h`'s, so "the same exponential order" is the conjecture, not a finding.

### (b) The excursion split

10. **The phase factorisation `e(2^s Y_h/3^h) = prod_j e((2^(s - A + P_(j-1)) mod 3^j)/3^j)`, including the
    modular reading for negative exponents.** VERIFIED. `Y_h = 2^-A sum_j 3^(h-j) 2^(P_(j-1))` with `2^-A` the
    inverse mod `3^h`; `3^(h-j) x / 3^h ≡ (x mod 3^j)/3^j mod 1`, and the inverse mod `3^h` reduces to the inverse
    mod `3^j`, so `2^(s-A+P)` with a negative exponent is the inverse power mod `3^j`, as both scripts compute.
    `[C1]`: my DP with no excursion constraint (`c = ∞`, `A` up to `2h + 8 sqrt(2h) + 40`, mass covered
    `1.000000`) reproduces the closed recursion to `7·10^-18` at `h = 20, s = 26` and `2·10^-19` at `h = 40,
    s = 57`; negative exponents `s - A + P` are exercised there.
11. **Completeness of the `A`-window `|A - s| <= 45 + c`.** VERIFIED. Band words have `h <= A = P_h < h log_2 3 +
    c <= s + 7.3 + c`, inside the window; `[C3]` at `h = 40, c = 6` the window and the full range `[h, lim_h]`
    contain the same 30 values of `A` and give identical sums and masses. The "full" coefficient is taken from
    the closed recursion, not from a windowed DP, so no truncation enters there.
12. **`c = 0` mass equals `P_h`.** VERIFIED `[C2]` (`2.986131·10^-3` both, difference `0`).
13. **The table: strict no-descent words carry `54%, 41%, 35%` (coherent fractions `0.25, 0.19, 0.15`); `c = 4`:
    `0.92, 0.93, 0.89`; `c = 6`: `1.008, 1.007, 0.997`, remainders `0.012, 0.108, 0.274`; `c = 16`: `0.001,
    0.013, 0.065`; `c = 24` at `120`: `0.002`; coherent fractions at `c = 6`: `0.020, 0.015, 0.013`; band masses
    `24, 30, 34 P_h`.** VERIFIED `[C4]`, `[C5]`: every quoted ratio reproduces to `5·10^-5` (masses `23.62,
    30.46, 33.62 P_h`).
14. **"The words with `E < c` reproduce the coefficient's magnitude at `c ≈ 6` at all three levels" (lines
    529–530); status row 790 "the band `E < 6` reproduces the magnitude".** VERIFIED WITH CORRECTION. `[C6]`:
    the magnitude ratio is within `4%` of `1` for every `c >= 5` at all three levels, but at `h = 120` it wanders
    (`1.04` at `c = 10`, `1.12` at `c = 12`, `1.05–1.06` at `c = 14–19`) while the vector remainder is `0.19,
    0.27, 0.29, 0.26, 0.19, 0.15, 0.13, 0.17` for `c = 5..12`: at `c = 6` the band vector is `27%` off and the
    magnitude agreement (`0.997`) is a coincidence of a rotated vector. The band that reproduces the coefficient as
    a vector to `10%` starts at `c = 5, 9, 15` (to `5%`: `c = 6, 10, 20`) at `h = 40, 80, 120`: the width grows
    roughly like `h/8`, which the note's next sentence ("the width of the band ... grows slowly with `h`")
    concedes but the headline "within about 6 bits" does not. The untracked profile output agrees (`10%` first
    reached at `m = 4` and `3` but not stably).
15. **"the strict no-descent words carry a third to a half of the coefficient" (line 888; lines 526–527).**
    VERIFIED WITH CORRECTION: `0.54, 0.41, 0.35` at `h = 40, 80, 120`, falling by about `0.1` per 40 levels; "a
    third to a half over `h = 40..120`, decreasing" is the fair statement (the trend is not mentioned).
16. **"so the descending words do not cancel; they add roughly in phase" (line 528).** VERIFIED: `|coh_0| +
    |rest| = 0.54 + 0.48 ≈ 1.02` of `|full|` at `h = 40` (`0.41 + 0.60`, `0.35 + 0.68`), and `arg(coh_0) -
    arg(full) = 0.21, 0.17` rad at `40, 80` (session's output).
17. **"the coefficient is a small residue of the band's mass, not a sum of aligned terms".** VERIFIED (coherent
    fraction `0.02, 0.015, 0.013` at `c = 6`).

### (c) The random units

18. **"twelve random units at `h = 30` give `|mu_hat| = 0.8·10^-8 .. 8·10^-8`, median `2.3·10^-8` against
    `0.84·3^-15 = 5.9·10^-8`" (lines 540–542).** VERIFIED WITH CORRECTION. The twelve values are genuine
    (uniformly random exponents `j_0`, i.e. uniformly random units, through the closed family); but the sentence
    compares a median with an rms. `[D1]` 100 own random units: rms `5.42·10^-8` against the level average
    `sqrt(0.70) 3^-15 = 5.83·10^-8` (mean `|mu_hat|^2 3^h = 0.60`, within the sampling error of a skewed
    variable), median `2.20·10^-8 = 0.41` rms, `60%` of the units below half the rms. `[D2]` the exact
    distribution over all `2·3^(h-1)` units at `h = 12, 14`: median `0.66, 0.62` rms, `35%, 38%` below half the
    rms (a Rayleigh law would give `0.83` and `22%`): the coefficient distribution over the units is skewed and
    gets more so with `h`. The twelve-unit rms `3.1·10^-8` is `1.9` below the level average, within the
    fluctuation of a heavy-tailed sample of twelve. "The typical coefficient is what square-root cancellation
    predicts" is right in the rms sense and should be stated that way.
19. **"`10^5` below the resonant window (`3.3·10^-3`)".** VERIFIED (`3.29·10^-3 / 2.3·10^-8 = 1.4·10^5`).

### (d) The targets and the wording

20. **(T1) "the sum ... over the words whose excursion exceeds `c` tends to `0` relative to the band sum as `c`
    grows, for every `h`" (lines 547–549).** VACUOUS as stated: for fixed `h` the set of words with excursion
    `> c` is empty once `c` exceeds the maximal excursion (`≈ h (2 - log_2 3)` for the words that matter), so the
    limit is trivial. The content is uniformity in `h`: `|full - coh_c| <= ε(c) |full|` with `ε(c) -> 0`
    uniformly in `h` — and that is exactly what the data put in doubt (item 14: at fixed `c = 6` the remainder
    grows `0.012, 0.108, 0.274`). Rewrite with "uniformly in `h`" (or with `c` allowed to grow like `o(h)`), and
    say that the data show the uniform version needs `c` growing at least like `h/8` on the observed range.
21. **(T2) "the sum over the words with excursion below `c` has modulus of the exponential order of `P_h` (its
    mass is `≍ P_h` by large deviations ...)" (lines 549–553).** Acceptable as OPEN. Two remarks: the band mass
    over strict mass grows on the observed range (`24, 30, 34` at `40, 80, 120`), so "`≍`" is a statement about
    exponential order only; and if (T1) holds only with `c` growing with `h`, (T2) with fixed `c` is not the
    right partner — the pair should be stated for a common `c = c(h)`.
22. **"the phase of a word at level `j` is `e((2^(s - T_j) mod 3^j)/3^j)`, `T_j` the suffix cost, which is near
    `1` iff the prefix sum `P_(j-1)` lies below `j log_2 3 + 6 - K`" (lines 544–547).** VERIFIED WITH
    CORRECTION: `s - T_j = s - A + P_(j-1)`, so the condition is `A - h log_2 3 + 6 <= P_(j-1) <= (A - h log_2 3)
    + j log_2 3 + 6 - K` (two-sided: below the lower end the exponent is negative and the phase is a "random"
    inverse power); the term `A - h log_2 3` is `O(c)` on the band but is not zero.
23. **Labels.** "VERIFIED numerics" (row 789, line 499): VERIFIED with the certification caveat of item 2;
    "OBSERVED law" (line 513): VERIFIED as a label, with the bracket corrected (item 6); "DIRECTION mechanism"
    (line 521–522): VERIFIED; "OPEN targets" (line 556): VERIFIED, with (T1) reworded (item 20). §0 item 1 (lines
    91–100): VERIFIED except the bracket (item 6). "`M(h) ≍ P_h` ... CONJECTURAL" (row 788): VERIFIED as a
    label; the data now lean against it (items 6, 7), which the row could say.
24. **"exponent window `12000` per level" (line 501).** VERIFIED (`40·300`, plus `640`).
25. **Next probes (lines 885–889): "section 2c refutes as a mechanism" the coherent no-descent family with fixed
    initial phases.** VERIFIED (coherent fractions `0.15–0.25`, phases spread).

---

## (ii) Corrections (old text -> new text; line numbers of commit `d61c39bc4`)

**C1 (lines 511–513; line 98; row 789).** "so the exponential rate of `M` is `0.941–0.949`, bracketing
`3^(h*-1)`, and the identification of the rate with the no-descent rate stays OBSERVED, now over three hundred
levels." -> "so the exponential rate of `M` is `0.947–0.950` by the admissible estimators (fits with residual at
most `2.5%` have `β = 1.5–1.75` and `r = 0.948–0.950`; the least-squares doubling slope `0.0785` gives `0.947`;
the pure geometric `0.941` has residual `20%` and the fit pinned to `3^(h*-1) = 0.9465` needs `β = 1.14` at six
times the free fit's residual): the no-descent rate sits at the lower edge of the bracket, `0.15–0.35%` below
the fitted rates, a gap that 300 levels cannot attribute to the prefactor or to the rate; the identification of
the rate with the no-descent rate stays OBSERVED only in this weaker sense." Line 98: "bracketed in
`0.941–0.949` by the prefactor degeneracy" -> "`0.947–0.950` by the admissible fits, the no-descent rate at its
lower edge". Row 789: "rate `0.941–0.949`" -> "rate `0.947–0.950` (`3^(h*-1)` at the lower edge)".

**C2 (lines 500–501).** "(`M(300) = 7.42·10^(-11)`; exponent window `12000` per level)" -> "(`M(300) =
7.42·10^(-11)`; exponent window `12000` per level; the truncation `a <= 40` is not self-certifying here, its a
priori bound `300·2^(-40) = 2.7·10^(-10)` exceeding `M(300)`; an `a <= 60` run certifies the values to
`3·10^(-4)` relative and agrees with the `a <= 40` run to `2·10^(-12)`)".

**C3 (lines 503–504).** "slope `0.077` per level against `|log_2 3^(h*-1)| = 0.079`" -> "least-squares slope
`0.0785` per level (two-point `0.077` between the doublings `50→100` and `150→300`) against `|log_2 3^(h*-1)|
= 0.079`, i.e. `r = 0.947`".

**C4 (line 514).** "The ratio `M/P_h` is `0.46 ± 0.02` up to `h = 180` and then drifts up" -> "The ratio
`M/P_h` stays in `0.44–0.48` up to `h = 180` and then drifts up".

**C5 (lines 515–517).** "`M` and `P_h` have the same exponential order, with prefactors that differ by a slowly
varying factor (`0.44–0.58` over `20 <= h <= 300`)" -> "over `20 <= h <= 300` the ratio moves within `0.44–0.58`;
the rise by a factor `1.26` over `180..300` is what a rate `0.2%` per level above the no-descent rate produces,
and whether `M` and `P_h` have the same exponential order (the conjecture `M ≍ P_h`) or `M` decays slightly
slower is not decided by these levels".

**C6 (lines 518–519).** "The resonant exponent is `s = h log_2 3 - 6.2 .. 7.1` throughout (`-6.49` at `h =
300`)." -> "The resonant exponent is `s = h log_2 3 - 5.3 .. 7.2` at every level `20..300` (`-6.2 .. -7.1` at
the multiples of 25; `-6.49` at `h = 300`)."

**C7 (lines 529–536; row 790).** "The words with `E < c` reproduce the coefficient's *magnitude* at `c ≈ 6` at
all three levels (`|coh_6|/|full| = 1.008, 1.007, 0.997`; at `c = 4`: `0.92, 0.93, 0.89`), while the vector
remainder `|full - coh_c|/|full|` decays with `c` more slowly as `h` grows (`c = 6`: `0.012, 0.108, 0.274`; `c =
16`: `0.001, 0.013, 0.065`; `c = 24` at `h = 120`: `0.002`): the deep-descending words' net contribution is a
phase rotation of decreasing size, and the width of the band that carries the coefficient grows slowly with
`h`." -> "The words with `E < c` reproduce the coefficient as a vector to `10%` from `c = 5, 9, 15` on (to `5%`
from `c = 6, 10, 20`) at `h = 40, 80, 120`: the band width grows roughly like `h/8`. The magnitude alone is
within `4%` for every `c >= 5` (`|coh_6|/|full| = 1.008, 1.007, 0.997`; at `c = 4`: `0.92, 0.93, 0.89`), but at
`h = 120` it wanders up to `1.12` (`c = 12`) while the vector remainder is `0.19–0.29` for `c = 5..9` (`c = 16`:
`0.065`; `c = 24`: `0.002`): the deep-descending words' net contribution is a rotation whose size decays with
`c` more slowly as `h` grows." Row 790: "the band `E < 6` reproduces the magnitude (`1.008, 1.007, 0.997`)" ->
"the band `E < c` reproduces the coefficient to `10%` from `c = 5, 9, 15` at `h = 40, 80, 120` (magnitude
within `4%` from `c = 5`)".

**C8 (lines 526–527; line 888).** "The strict no-descent words (`E < 0`) carry `54%, 41%, 35%` of
`|mu_hat_h(2^s)|` at `h = 40, 80, 120`" -> "... at `h = 40, 80, 120`, a share falling by about `0.1` per forty
levels"; line 888 "carry a third to a half of the coefficient" -> "carry a third to a half of the coefficient
over `h = 40..120`, a decreasing share".

**C9 (lines 539–542).** "The typical coefficient is what square-root cancellation predicts: twelve random units
at `h = 30` give `|mu_hat| = 0.8·10^(-8) .. 8·10^(-8)`, median `2.3·10^(-8)` against `0.84·3^(-15) =
5.9·10^(-8)`, `10^5` below the resonant window (`3.3·10^(-3)`)." -> "The typical coefficient is at the
square-root scale: random units at `h = 30` have rms `|mu_hat| ≈ 5·10^(-8)` (twelve units: `0.8·10^(-8) ..
8·10^(-8)`, rms `3·10^(-8)`; a hundred: rms `5.4·10^(-8)`) against the level average `0.84·3^(-15) =
5.9·10^(-8)`; the median is only `0.4` of the rms (`2.3·10^(-8)`) because the coefficient distribution over the
units is skewed (exact at `h = 14`: median `0.62` rms, `38%` of the units below half the rms), `10^5` below the
resonant window (`3.3·10^(-3)`)."

**C10 (lines 547–549, (T1)).** "the sum of `2^(-A) e(2^s Y_h/3^h)` over the words whose excursion exceeds `c`
tends to `0` relative to the band sum as `c` grows, for every `h`" -> "the sum of `2^(-A) e(2^s Y_h/3^h)` over
the words whose excursion exceeds `c` is at most `ε(c)` times the band sum, with `ε(c) -> 0` uniformly in `h`
(for fixed `h` the statement is empty; the data show the uniform version needs `c` of order `h/8` on `40..120`,
so the target may have to be stated with `c = c(h) = o(h)`)".

**C11 (lines 545–547).** "which is near `1` iff the prefix sum `P_(j-1)` lies below `j log_2 3 + 6 - K`" ->
"which is near `1` iff `A - h log_2 3 + 6 <= P_(j-1) <= A - h log_2 3 + j log_2 3 + 6 - K` (on the band `A -
h log_2 3` is between `-∞` and `c`; below the lower end the exponent is negative and the phase is a scrambled
inverse power)".

---

## (iii) Verdict

**SOUND WITH CORRECTIONS.**

Every number of section 2c reproduces with independent code: the recursion to `h = 300` (certified here with
`a <= 60` to `3·10^-4`; the session's `a <= 40` values agree to `2·10^-12` but are not self-certifying at
`h >= 250`), the doubling exponents, the geometric-mean ratios, the free and pinned fits, `P_h` to 300, the
excursion split at every listed `(h, c)` (the DP validated to `10^-18` at `c = ∞`, the phase factorisation and
its modular reading correct, the `A`-window complete), the band masses and the random units. The conclusions
survive the error analysis. The substantive corrections, in order:

1. **The rate bracket (item 6; C1).** "`0.941–0.949`, bracketing `3^(h*-1)`" presents the pure-geometric misfit
   (`20%` residual) as one end of the bracket; the admissible fits give `0.947–0.950`, with the no-descent rate
   `0.9465` at the lower edge (pinned fit six times worse than free). Together with the `M/P_h` rise from
   `0.46` to `0.58` over `180..300`, the data over 300 levels lean toward a rate slightly above `P_h`'s; the
   conjecture `M ≍ P_h` is not supported better by 300 levels than by 120, and the note should say so.
2. **The band width (items 14, 20; C7, C10).** "Reproduces the magnitude at `c ≈ 6`" is true but a coincidence
   of a rotated vector at `h = 120` (remainder `27%`); the band reproducing the coefficient to `10%` widens as
   `c = 5, 9, 15` at `h = 40, 80, 120`, roughly `h/8`. Consequently target (T1) as written ("for every `h`") is
   empty and must be a uniform statement, which the data show requires `c` growing with `h`.
3. **The certification of the level-300 numerics (item 2; C2).** The `a <= 40` truncation's a priori bound
   exceeds `M(300)`; the values are right (verified by `a <= 60`), and the note should record this.
4. **Smaller corrections (items 3, 7, 8, 15, 18, 22; C3–C6, C8, C9, C11):** the slope `0.077` is a two-point
   estimate (least squares `0.0785`); "`0.46 ± 0.02` up to `180`" contradicts §0's `0.44–0.48`; the exponent
   offset runs `5.3–7.2` over all levels; the strict share decreases; the random-unit comparison mixes a median
   with an rms (the distribution over units is skewed: median `0.4–0.6` rms); the phase-alignment condition is
   two-sided and carries the term `A - h log_2 3`.

Not checked: the untracked excursion-profile script beyond its band-width cross-check; the asymptotic prefactor
exponents (`3/2` expected for `P_h`; the fitted `1.3–1.6` are finite-range values); whether the argmax over all
units stays on the powers of two beyond `h = 18` (unchanged since audit 2).
