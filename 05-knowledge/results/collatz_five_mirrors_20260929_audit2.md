# Independent audit (2) of section 2b of `collatz_five_mirrors_20260929.md`: the powers of two to level 120 and the no-descent probability

**Auditor:** independent session (Fable), 2026-09-29. **Audited state:** the note at commit `7d59c6064` (line numbers
below refer to it): section 2b (lines 395–467), the sentences that cite it in the title, the header (lines 37–45),
§0 items 1 and 3 (lines 85–99, 119–123), §2 (lines 318–319, 330–344, 376–380), §6 (lines 640–652) and the status
rows (lines 684–686); the script `collatz_five_mirrors_powers_of_two_20260929.py` with its output, and the
untracked `collatz_five_mirrors_multiplier_families_20260929.py/.out` that §2b cites.

**Method.** `04-computation/experiments/collatz_five_mirrors_20260929_audit2.py` (62 s; output
`collatz_five_mirrors_20260929_audit2.out`; tags `[A2]`, `[C6]`, ... refer to its lines) implements the closed
recursion independently (numpy over an exponent window, phases from the exact integers `2^k mod 3^n`), checks it
against a brute-force word sum, reproduces the table, does a rigorous float64 error analysis plus an `AMAX = 50`
run and a 30-digit `mpmath` recomputation through the full dependency cone at `h = 30` and `60`, computes `P_h` by
an own dynamic programme (strict and non-strict inequality), fits the three decay models, and re-runs the
multiplier-family comparison to `h = 60`.

---

## (i) Claims checked

### (a) Closure and exactness

1. **The family `{2^j mod 3^n : j ∈ Z}` is closed under `t -> t 2^-a`.** PROVED (trivial: `2^j 2^-a = 2^(j-a)`).
   `[A1]`. Since `2` generates the units mod `3^n`, the family IS the unit group; the max over a window shorter
   than the period `2 3^(n-1)` is a lower bound for `M(n)`, over a full period it is `M(n)`.
2. **The recursion `m_n(j) = sum_a 2^-a e((2^(j-a) mod 3^n)/3^n) m_(n-1)(j-a)`, `m_0 = 1`, equals the word sum.**
   VERIFIED. `[A2]` against a brute-force sum over all `6^n` words (`n <= 4`, valuations `<= 6` on both sides,
   seven exponents): max deviation `3.2·10^-16`. `[A3]` `m_1(0) = (1/3)e(1/3) + (2/3)e(2/3)`, `|m_1| = 1/sqrt 3`.
3. **Inverse powers (`j < 0`) well defined.** VERIFIED. `[A4]` `2^-k 2^k = 1 mod 3^n` at `n = 5, 30, 120`.
4. **Truncation `a <= 40`, "error `<= 2^(-40)` per level".** VERIFIED and sharpened. The truncated step is a
   weighted average (weights `2^-a`, sum `1 - 2^-40`, unimodular phases), hence an `l^∞` contraction: errors add
   and never amplify, `|error at level n| <= n (2^-40 + 41 eps) = 1.1·10^-10` at `n = 120`, i.e. relative
   `2.7·10^-5` of the value `4.12·10^-6` `[C8]`. `AMAX = 40` against `50`: differences `<= 10^-14` `[C9]`.
5. **Exponent-range bookkeeping ("level `n` needs exponents down to `-40 (HMAX - n)`").** VERIFIED. `[A5]`:
   level `n` on `[-40(HMAX-n), JMAX]` uses level `n-1` on `[-40(HMAX-n+1), JMAX-1]`, inside its window (asserted in
   my code; the session's `prev.get(j - a, 0j)` default is never hit inside the reported window). "About `5000`
   exponents per level": `4800 + 2·120 + 20`.

### (b) Agreement with the FFT

6. **`max_j |m_h(j)|` reproduces `M(h)` to six digits for `h <= 18` with offsets `s - h = 1..5` for `9 <= h <=
   18`.** VERIFIED. `[B2]`, `[B3]`: max deviation `4.9·10^-7` (the FFT values are quoted to six decimals);
   offsets `+1, +2, +2, +2, +3, +3, +3, +4, +4, +5`.
7. **"for `h <= 8` the family covers all units and the maximum is `M(h)` by definition".** VERIFIED. `[B1]` the
   window (`4760` exponents at `h = 8`) exceeds the period `2 3^7 = 4374`; at `h = 9` (`13122`) it does not.
8. **"Beyond, `max_j |m_h(j)|` is a lower bound for `M(h)` (an equality if the maximum stays on the powers of
   two, true wherever checked)".** PROVED (subset maximum) / OBSERVED as hedged (item 26).

### (c) The values to `h = 120`, the fits, the numerics

9. **Table values `8.88e-3 .. 4.12e-6`, offsets `6 .. 64`, ratios `.905 .. .930`.** VERIFIED. `[C1]`, `[C2]`: my
   values are identical to the session's output to seven digits (`8.884593e-3`, `2.828142e-4`, `4.119742e-6`);
   the `2.3·10^-3` relative deviation from the table is its three-digit rounding; offsets and ratios identical.
10. **Local exponents `2.1, 2.3, 2.7, 3.1, 3.5, 4.4, 5.3, 6.1`.** VERIFIED. `[C4]` (`2.107, 2.308, 2.683, 3.102,
    3.539, 4.422, 5.272, 6.101`).
11. **"a geometric fit over `h = 40..120` has ratio `0.930` with maximal log-residual `0.12`, against `0.32` for
    the best shifted power law (exponent `6.9`, shift `20`)".** VERIFIED WITH CORRECTION. The geometric numbers
    reproduce (`0.93026`, `0.1208` `[C5]`). The shifted-power residual is an artefact of the grid bound on the
    shift: the session's grid stops at `c = 19.9` (residual `0.32`), mine at `c = 40` (`0.24`, exponent `8.4`,
    again at the boundary), and as `c -> ∞` a shifted power law `(1 + h/c)^-α` with `α ∝ c` tends to a geometric
    law, so "the best shifted power law" is not a competitor with a well-defined residual. The discriminating
    statistic is the doubling exponent `E(h) = log2(m(h)/m(2h))`: for `C h^-β r^h` it is `h |log2 r| + β`, linear
    in `h`; for a shifted power law it is bounded by `α` and saturates. Observed: linear with slope `0.083` per
    level and intercept `1.1` `[C4b]`; `|log2 0.9465| = 0.079`, `|log2 0.930| = 0.105`.
12. **"The decay is geometric".** VERIFIED WITH CORRECTION: geometric with a polynomial prefactor. `C h^-1.09
    (0.9437)^h` fits `40..120` to `0.9%` max log-residual `[C6]`, against `12%` for the pure geometric law
    (a systematic curvature: the pure ratio is `0.930` over `40..120` and `0.934` over `80..120`). The rate `r`
    of the prefactor fit (`0.944`) and the doubling-exponent slope (`0.083`, i.e. `r = 0.944`) both sit at the
    no-descent rate `0.9465`, not at the pure-geometric `0.930`.
13. **"local exponents ... grow without bound ... as `h |log_2 r|` for a geometric law" (lines 425–428).**
    VERIFIED WITH CORRECTION: as `h |log2 r| + β`, `β ≈ 1.1`; with `r = 0.930` the formula alone gives `1.0`
    at `h = 10` against the observed `2.1`. "Without bound" is the extrapolation of a linear trend seen over
    `h = 10..60`.
14. **"The audit's shifted power law ... fails beyond `h ≈ 30`" (lines 428–430).** VERIFIED. The `h <= 18` fit
    `C (h + 2.5)^-2.52` (`C = 22.3` from `M(17)`) gives `3.2e-3` at `h = 31` (measured `3.02e-3`), `1.75e-3` at
    `40` (`1.38e-3`), `6.6e-4` at `60` (`2.83e-4`): a factor `2` near `h = 55`, as audit 1 predicted (`h ≈ 31`
    for the separation of the two readings at that level of fit).
15. **"The resonant exponent grows linearly, `s/h = 1.52–1.53` at `h = 100–120`, approaching `log_2 3`" (lines
    430–431).** VERIFIED and sharpened: `s - h log2 3 = -5.5 .. -6.8` for every `20 <= h <= 120` `[C3]`, i.e.
    `s = h log2 3 - 6 ± 1` and `2^s/3^h = 0.009–0.019 ≈ 2^-6` `[C3b]`: the resonant frequency is a fixed fraction
    of the modulus (the character `e(t y/3^h)` with `t ≈ 3^h/64`). This amends my audit-1 correction C1, whose
    "`2^s/3^h` decreases geometrically" (now lines 318–319) holds only to `h ≈ 18`; beyond, the ratio stops near
    `2^-6` (item 37).
16. **Float64 rounding as the source of the decay.** EXCLUDED. Rigorous bound `1.1·10^-10` at `h = 120` `[C8]`;
    `AMAX = 40/50` agree to `10^-14` `[C9]`; a 30-digit `mpmath` recomputation of the argmax coefficient through
    its full dependency cone gives `|difference| = 9.7·10^-19` at `h = 30` and `3.1·10^-19` at `h = 60` `[C10]`.
    The decay from `0.58` to `4·10^-6` is real to at least ten digits.

### (d) The no-descent probability

17. **`P_h` values `1.93e-2 .. 9.31e-6`.** VERIFIED. `[D3]` own DP, three-digit rounding (`3.2·10^-3`).
18. **"`max_j |m_h(j)| / P_h = 0.46 ± 0.02` at every level `20 <= h <= 120`".** VERIFIED WITH CORRECTION. Range
    `0.443–0.481` `[D4]`; one level outside the band (`h = 41`: `0.481`); window means `0.459, 0.465, 0.456`
    over `20–40, 41–80, 81–120`; the last levels drift down (`0.449, 0.443` at `110, 120`) because the
    coefficient's prefactor decays slightly faster than `P_h`'s (local exponents `1.34` against `1.29` on
    `60 -> 120` `[D6]`). "`0.44–0.48`, drifting slowly" is the honest statement; "identity" is too strong.
19. **The rate `P_h^(1/h) -> e^(-I(log_2 3)) = 3^(h*-1) = 0.94650`, `h* = h(log_3 2) = 0.94996`.** The identity
    VERIFIED: with `p = log_3 2 = 1/log_2 3` and `h*` the binary entropy of `p` in bits, `(h* - 1) ln 3 = -ln p -
    ((1-p)/p) ln(1-p) - ln 3 = -I(1/p)` where `I(α) = α ln 2 + (α-1) ln(α-1) - α ln α` (numerically `0.054979`,
    `e^-I = 0.94650 = 3^(h*-1)` `[D2]`); `h* = 0.94996` is THM-4476's constant. The limit is a THEOREM, not an
    observation: `P_h <= P(S_h < ch) <= e^(-I(c) h)` (Chernoff, `c = log_2 3 < E a = 2`; holds at every `h`
    `[D7]`), and for the lower bound tilt the valuations to mean `c - ε`: `P_h >= e^(-h (I(c-ε) + |θ_(c-ε)| ε))
    P_θ(S_j < cj ∀ j <= h, S_h >= (c-2ε) h)`, the last probability bounded below by a constant (law of large
    numbers plus `P(a_1 = 1) = 1/2`), then `ε -> 0`. "The cheapest way to stay under the critical line is to walk
    along it" is the correct heuristic behind the tilted lower bound; adequate as a pointer, but the status should
    be PROVED (standard) rather than sit inside an OBSERVED paragraph.
20. **Finite-level behaviour of `P_h`.** `P_h^(1/h) = 0.848, 0.884, 0.908` at `h = 30, 60, 120` `[D5]`; the
    prefactor `P_h e^(I h)` decays like `h^-β` with `β = 1.27` fitted over `40..120` (locally `1.14 -> 1.29`
    `[D6]`; the expected asymptotic value is `3/2`, the local persistence exponent of the tilted mean-zero walk
    weighted by `e^(|θ*| X_h)`), so the level ratio approaches `0.9465` from below like `e^-I (1 - β/h)`.
21. **"the ratios `P_h/P_(h-1)` are `0.922, 0.931, 0.939, 0.948` at `h = 30, 40, 50, 60` and `0.949` at `120`,
    still rising toward `0.9465`" (lines 451–452).** REFUTED as a description. The single-level ratios oscillate
    with the integer part of `j log_2 3`: `0.903` at `70`, `0.910` at `80`, `0.919` at `90`, `0.944` at `100`,
    `0.946` at `110`, `0.949` at `120` `[D3]`; the quoted `0.948` and `0.949` already exceed `0.9465`. The
    geometric-mean ratio is `0.9325` over `61..120` and `0.9366` over `101..120` `[D5]`, rising toward `0.9465`
    from below as item 20 predicts.
22. **The strict inequality in `P_h`.** VERIFIED immaterial: `j log_2 3` is irrational for `j >= 1`, so `S_j < j
    log_2 3` iff `S_j <= floor(j log_2 3)`; strict and non-strict DPs agree exactly `[D1]`. The strict form is the
    natural one ("no descent" = `3^j > 2^(S_j)`, the pure-drift value stays above the start).
23. **The mechanism sentence (lines 452–458): "the coherent family is the no-descent set, whose members' phases
    `e(2^(s-S_j)/3^j)` are fixed roots of unity determined by the first few valuations (the constant) while the
    descending words cancel".** NOT ESTABLISHED (DIRECTION): no decomposition at any `h >= 20` supports it; the
    level-`h` coefficient involves the suffix sums (the last valuations; by the i.i.d. symmetry the suffix
    no-descent probability is the same `P_h`); for `s - S_j < j log_2 3` the integer `2^(s-S_j)` is small compared
    with `3^j`, so the phase is close to `1`, not a "fixed root of unity". Labelled "Reading" inside an OBSERVED
    paragraph; should carry DIRECTION explicitly.

### (e) Wording and status

24. **"PROVED for the closure and `M(h) >= max_j |m_h(j)|`" (lines 464–465).** VERIFIED (items 1, 8).
25. **"VERIFIED (exact recursion, truncation `2^(-40)`)" for the numbers (line 684).** VERIFIED as a label; the
    recursion is exact, the computation float64 with a rigorous error `1.1·10^-10` (item 16).
26. **The multiplier families (lines 433–441, 685): "at every level `h <= 60` the pure family `u = 1` gives the
    largest coefficient".** VERIFIED. `[E1]`, `[E2]`: pure family largest at every level `1..60`, the others at
    `0.33` of it for `h >= 40` (`u = 13`); ties at `h <= 9` (full-period coincidence). As a test of "the maximum
    stays on the powers of two" it samples seventeen windows of `~ 36 (60 - h) + 160` residues out of `2 3^(h-1)`
    units (every `u` is `2^k` for a 3-adic exponent `k`; the families are the pure family shifted by a 3-adic
    discrete log) `[E3]`; the note's OBSERVED label for the general statement is right.
27. **§0 item 1 (lines 95–97): "on the powers of two, which carry the maximum wherever that was checked, the
    sup-norm mixing rate of the 3-adic law is the no-descent rate (OBSERVED; CONJECTURAL for `M(h)`)".**
    OVERREACH in the main clause. What is observed is `max_j |m_h(j)| ≈ 0.44–0.48 P_h` to `h = 120`. "The
    sup-norm mixing rate of the 3-adic law" is (a) `M(h)`'s rate only if the maximum stays on the powers of two
    for all `h` (OBSERVED to `18` over all units, to `60` over a sample), and (b) "mixing rate" is a name for the
    decay of the maximal primitive Fourier coefficient, not the rate of the mixing estimate (2.3), which is an
    `ℓ^1` statement dominated by the bulk of square-root coefficients (S20 §2: `d(2, 18) = 0.74` against `M(3) =
    0.25`; S19: `d(m, 18) = 0.85 .. 0.35` for `m = 1..9`).
28. **Consequence (i) (lines 458–459): "on this family (2.3) holds with the geometric rate `3^(h*-1)` per level
    and cannot hold faster".** REFUTED as stated. (2.3) is not a statement on a family of frequencies. Through
    Proposition 2's inequality the family gives the LOWER bound `d(m, q) = ||mu_q - lift mu_m||_1 >= M(m+1) >=
    max_j |m_(m+1)(j)| ≈ 0.46 P_(m+1)` (PROVED / OBSERVED), so (2.3)'s `C_A m^-A` can never beat the no-descent
    rate — the "cannot hold faster" half is right. "Holds with the geometric rate" reverses Proposition 2: a
    geometric bound on `M(h)` is necessary for (2.3), not sufficient, and `d(m, q)` is far above `M(m+1)`.
29. **Consequence (ii) (lines 459–463).** Acceptable as CONJECTURAL; "are the same number" should read "the
    conjectured rate is `3^(h*-1)`, a function of `h*`".
30. **Consequence (iii) (line 463): "Mazur's `C* h^(-6409)` is a polynomial statement about a geometric
    quantity".** Conditional on the conjecture; "is" -> "would be".
31. **§6 (lines 643–647): "shows a geometric decay ... at the no-descent rate `3^(h*-1) = 0.9465` (OBSERVED). If
    the maximum stays on the powers of two (true wherever checked), (2.3) holds with a geometric constant".**
    First sentence VERIFIED WITH CORRECTION (the measured level ratio over `40..120` is `0.930`; the rate `0.9465`
    is the asymptotic rate supported by the doubling-exponent slope and the prefactor fit); second sentence
    REFUTED (item 28: the implication runs the other way).
32. **Title (line 1): "decaying at the no-descent rate `3^(h*-1)` on the powers of two to level 120".** The rate
    is inferred, the measured ratio `0.93`; "decaying like the no-descent probability (`0.46 P_h`) to level 120"
    is what was measured.
33. **§0 item 3 (lines 119–121): "is confirmed at scale by section 2b".** OVERREACH: section 2b shows a
    proportionality of two sequences, not that the no-descent words carry the coefficient; "consistent with".
34. **Header (lines 39–41): "the coefficient is `0.46` times the no-descent probability at every level
    `20..120`".** -> "`0.44–0.48` times" (item 18).
35. **§0 item 1 (line 94): "`P_h` decays at the rate `3^(h*-1) = 0.9465`".** True and PROVED (item 19); the
    finite-level ratio lies below it by the prefactor (item 20); no change needed beyond the status.
36. **Status row 686: "`M(h) ≍ P_h`, rate `3^(h*-1) = 0.9465` | CONJECTURAL".** VERIFIED as labelled.
37. **Audit-1 amendment (lines 318–319, my correction C1): "so `2^s/3^h` decreases geometrically (`0.12` at `h
    = 7`, `0.027` at `h = 14`, `0.022` at `h = 18`)".** VERIFIED WITH CORRECTION by section 2b's own data:
    the decrease stops near `2^-6` (`0.019` at `20`, `0.009–0.015` for `40..120` `[C3b]`); `s = h log_2 3 - 6 ± 1`.

---

## (ii) Corrections (old text -> new text; line numbers of commit `7d59c6064`)

**C1 (lines 458–459).** "Consequences: (i) on this family (2.3) holds with the geometric rate `3^(h*-1)` per
level and cannot hold faster;" -> "Consequences: (i) by Proposition 2's inequality the `ℓ^1` distances of (2.3)
satisfy `||mu_q - lift mu_m||_1 >= max_j |m_(m+1)(j)| ≈ 0.46 P_(m+1)` (PROVED / OBSERVED), so (2.3) cannot hold
with a rate faster than the no-descent rate; whether it holds with that rate is a statement about the bulk
coefficients, which dominate the distances (§2), and is not decided here;".

**C2 (lines 645–647).** "If the maximum stays on the powers of two (true wherever checked), (2.3) holds with a
geometric constant, which no published proof provides (Mazur's `C` is a three-fold tower, S19 audit), and the
sharp sup-norm form of the mixing estimate is `M(h) ≍ P_h`" -> "If the maximum stays on the powers of two (true
wherever checked), `M(h)` decays like `P_h`, which is the necessary condition for (2.3) that Proposition 2
extracts (a geometric bound on `M(h)` does not by itself give (2.3), whose `ℓ^1` distances are dominated by the
bulk coefficients), and the conjectured sharp form of the coefficient bound is `M(h) ≍ P_h`".

**C3 (lines 643–645).** "shows a geometric decay of `max_j |mu_hat_h(2^j)|` to `h = 120` at the no-descent rate
`3^(h*-1) = 0.9465` (OBSERVED)" -> "shows `max_j |mu_hat_h(2^j)|` decaying to `h = 120` like `C h^(-1.1) r^h`
with `r ≈ 0.944` (fit to `0.9%`; a pure geometric fit gives ratio `0.930` with `12%` residual), the doubling
exponents growing linearly with slope `0.083 ≈ |log_2 0.9465|`: geometric with a polynomial prefactor, at the
no-descent rate `3^(h*-1) = 0.9465` asymptotically (OBSERVED)".

**C4 (lines 422–428).** "The decay is geometric: a geometric fit over `h = 40..120` has ratio `0.930` with
maximal log-residual `0.12`, against `0.32` for the best shifted power law (which needs exponent `6.9` and shift
`20`), and the local exponents on successive doublings grow without bound, `2.1` (10→20), ..., `6.1` (60→120), as
`h |log_2 r|` for a geometric law." -> "The decay is geometric with a polynomial prefactor: `C h^(-1.1) r^h`, `r =
0.944`, fits `h = 40..120` to `0.9%` (a pure geometric fit has ratio `0.930` and residual `12%`; a shifted power
law is not a well-posed competitor over a bounded range, since with a large shift it tends to a geometric law),
and the local exponents on successive doublings, `2.1` (10→20), ..., `6.1` (60→120), grow linearly in `h` with
slope `0.083` per level and intercept `1.1`, as `h |log_2 r| + β` for `C h^(-β) r^h` (`|log_2 0.9465| = 0.079`);
a shifted power law would saturate them at its exponent."

**C5 (lines 443–445).** "`max_j |m_h(j)| / P_h = 0.46 ± 0.02` at every level `20 <= h <= 120` (`0.461` at
`20`, `0.463` at `50`, `0.462` at `100`, `0.443` at `120`)." -> "`max_j |m_h(j)| / P_h` lies in `0.44–0.48` at
every level `20 <= h <= 120` (window means `0.459, 0.465, 0.456` over `20–40, 41–80, 81–120`; `0.481` at
`h = 41`, `0.443` at `120`, the coefficient's prefactor decaying slightly faster than `P_h`'s)." Same change at
line 40 ("`0.46` times" -> "`0.44–0.48` times") and line 92 ("`0.46 ± 0.02`" -> "`0.44–0.48`").

**C6 (lines 445–452).** "The rate of `P_h` is the large-deviation rate of the valuation sum at the critical slope
`log_2 3`: `P_h^(1/h) -> e^(-I(log_2 3)) = 3^(h* - 1) = 0.94650` with `h* = h(log_3 2) = 0.94996` ... (standard:
the cheapest way to stay under the critical line is to walk along it; the ratios `P_h/P_(h-1)` are `0.922, 0.931,
0.939, 0.948` at `h = 30, 40, 50, 60` and `0.949` at `120`, still rising toward `0.9465`)." -> "The rate of `P_h`
is the large-deviation rate of the valuation sum at the critical slope `log_2 3`: `P_h^(1/h) -> e^(-I(log_2 3)) =
3^(h* - 1) = 0.94650` with `h* = h(log_3 2) = 0.94996` the binary entropy of `log_3 2` (PROVED, standard:
`P_h <= P(S_h < h log_2 3) <= e^(-I h)` by Chernoff, and the matching lower bound by tilting the valuations to
mean `log_2 3 - ε` — the cheapest way to stay under the critical line is to walk along it). The convergence is
slow because of a polynomial prefactor, `P_h e^(I h) ≈ h^(-1.3)` over `40..120` (expected `h^(-3/2)`): the
single-level ratios `P_h/P_(h-1)` oscillate with the integer part of `j log_2 3` (`0.903` at `70`, `0.949` at
`120`), and the geometric-mean ratio is `0.9325` over `61..120` and `0.9366` over `101..120`, approaching
`0.9465` from below."

**C7 (lines 452–458, "Reading").** Prefix "Reading (DIRECTION):" and replace "whose members' phases
`e(2^(s-S_j)/3^j)` are fixed roots of unity determined by the first few valuations (the constant)" -> "whose
members' phases `e(2^(s-T_j)/3^j)` (`T_j` the suffix sums, `2^(s-T_j)` small against `3^j` on the no-descent set)
stay near `1`, the last few valuations fixing the constant".

**C8 (line 463).** "(iii) Mazur's `C* h^(-6409)` is a polynomial statement about a geometric quantity." -> "(iii)
Mazur's `C* h^(-6409)` would then be a polynomial statement about a geometric quantity."

**C9 (lines 95–97).** "on the powers of two, which carry the maximum wherever that was checked, the sup-norm
mixing rate of the 3-adic law is the no-descent rate (OBSERVED; CONJECTURAL for `M(h)`)." -> "on the powers of
two, which carry the maximum wherever that was checked, the maximal Fourier coefficient decays like the no-descent
probability (OBSERVED to `h = 120`); that the decay rate of `M(h)` is the no-descent rate `3^(h*-1)` is
CONJECTURAL (it needs the maximum to stay on the powers of two), and it concerns the maximal primitive
coefficient, not the `ℓ^1` distances of (2.3)."

**C10 (title, line 1).** "decaying at the no-descent rate `3^(h*-1)` on the powers of two to level 120" ->
"decaying like the no-descent probability, `0.46 P_h`, on the powers of two to level 120".

**C11 (lines 119–121).** "is DIRECTION at these levels and is confirmed at scale by section 2b: to `h = 120` the
coefficient at the powers of two tracks the no-descent probability" -> "is DIRECTION at these levels and is
consistent with section 2b: to `h = 120` the coefficient at the powers of two is proportional to the no-descent
probability".

**C12 (lines 430–431).** "The resonant exponent grows linearly, `s/h = 1.52–1.53` at `h = 100–120`, approaching
`log_2 3 = 1.585`." -> "The resonant exponent is `s = h log_2 3 - 6 ± 1` for every `20 <= h <= 120`: the
resonant frequency is a fixed fraction `2^s/3^h ≈ 2^(-6)` (`0.009–0.019`) of the modulus, and `s/h` approaches
`log_2 3 = 1.585` (`1.53` at `h = 120`)."

**C13 (lines 318–319, the audit-1 correction).** "so `2^s/3^h` decreases geometrically (`0.12` at `h = 7`,
`0.027` at `h = 14`, `0.022` at `h = 18`)" -> "so `2^s/3^h` decreases (`0.12` at `h = 7`, `0.027` at `h = 14`,
`0.022` at `h = 18`) until it settles near `2^(-6)` (section 2b: `s = h log_2 3 - 6 ± 1` for `20 <= h <= 120`)".

**C14 (status row 684).** "`max_j |mu_hat_h(2^j)|` to `h = 120`: geometric, ratio `0.906 -> 0.930`, local
exponents `2.1 -> 6.1`; `= 0.46 P_h` at every level `20..120` | VERIFIED (exact recursion, truncation
`2^(-40)`); the identity OBSERVED" -> "`max_j |mu_hat_h(2^j)|` to `h = 120`: `C h^(-1.1) r^h`, `r ≈ 0.944`
(pure geometric ratio `0.906 -> 0.930`), doubling exponents `2.1 -> 6.1` with slope `0.083 ≈ |log_2 0.9465|`;
`= 0.44–0.48 P_h` over `20..120` | VERIFIED (exact recursion, truncation `2^(-40)`, float64 error `<= 1.1·10^(-10)`,
independently recomputed; 30-digit check at `h = 30, 60`); the proportionality OBSERVED".

**C15 (lines 89–92).** "decays geometrically (ratio `0.906 -> 0.930` over `h = 30..120`; local exponents `2.1,
2.3, 2.7, 3.1, 3.5, 4.4, 5.3, 6.1` on successive doublings, the signature of a geometric law)" -> "decays
geometrically up to a polynomial prefactor (`C h^(-1.1) r^h`, `r ≈ 0.944`; level ratio `0.906 -> 0.930` over `h =
30..120`; local exponents `2.1, ..., 6.1` on successive doublings, growing linearly with slope `0.083 ≈ |log_2
0.9465|`, the signature of a geometric law)".

---

## (iii) Verdict

**SOUND WITH CORRECTIONS.**

The closure is trivial and the recursion is exact (item 2); the numerical values, offsets, ratios and local
exponents to `h = 120`, the FFT agreement to `h = 18`, the no-descent probabilities and their ratio to the
coefficients all reproduce with an independent implementation (items 6–10, 17, 26); float64 rounding is excluded
rigorously and by a 30-digit recomputation (item 16); the identity `3^(h*-1) = e^(-I(log_2 3))` is right and
the rate statement for `P_h` is a theorem (item 19). The substantive corrections, in order of importance:

1. **The implication in consequence (i) and in §6 runs the wrong way** (items 28, 31; C1, C2). The family gives a
   LOWER bound on the `ℓ^1` distances of (2.3) (they cannot decay faster than the no-descent rate); it does not
   make (2.3) "hold with a geometric constant", and "(2.3) holds on this family" is not a statement about (2.3).
   Proposition 2 extracts a necessary condition from (2.3); section 2b measures that condition.
2. **The rate is asymptotic, the measured decay has a polynomial prefactor** (items 11–13, 21; C3, C4, C6, C14,
   C15). Over `40..120` the level ratio is `0.930` and a pure geometric law misfits by `12%`; `C h^(-1.1) r^h`
   fits to `0.9%` with `r ≈ 0.944`, and the doubling exponents grow linearly with slope `0.083`, against
   `|log_2 0.9465| = 0.079` — this, not the fit residual of a shifted power law (a grid artefact), is the
   evidence for the no-descent rate. The `P_h` ratios sentence is wrong: single-level ratios oscillate and two of
   the quoted ones exceed `0.9465`; the geometric-mean ratio approaches it from below.
3. **"`0.46 ± 0.02` at every level" fails** (item 18; C5): the range is `0.443–0.481`, one level outside the
   band, with a slow downward drift; "`0.44–0.48`, an approximate proportionality" is the observation.
4. **Status and wording** (items 19, 23, 27, 30, 32, 33; C6–C11): the `P_h` rate is PROVED (standard LD), not
   OBSERVED; the mechanism sentence is DIRECTION; "the sup-norm mixing rate of the 3-adic law is the no-descent
   rate", "confirmed at scale", the title's "decaying at the no-descent rate" and consequence (iii)'s "is" state
   more than was measured.
5. **The resonant frequency is a fixed fraction of the modulus** (item 15, 37; C12, C13): `s = h log_2 3 - 6 ±
   1`, `2^s/3^h ≈ 2^(-6)` for `20 <= h <= 120`; this also amends the audit-1 correction "`2^s/3^h` decreases
   geometrically", which holds only to `h ≈ 18`.

Not checked: the completeness of the `3x + k` material and the rest of the note (unchanged since audit 1); the
asymptotic value `3/2` of the prefactor exponent of `P_h` (expected, not proved here); the identity `M(h) ≍ P_h`
beyond the seventeen sampled windows (CONJECTURAL, as labelled).
