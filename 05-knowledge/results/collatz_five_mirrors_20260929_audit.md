# Independent audit of `collatz_five_mirrors_20260929.md` (S20: five mirrors, Fourier profile `M(h)`, Propositions 1–6)

**Auditor:** independent session (Fable), 2026-09-29. **Audited state:** the note as of commit `b37902bfa`
(545 lines; line numbers below refer to it) with the scripts `collatz_five_mirrors_20260929.py`,
`..._fourier_deep_20260929.py`, `..._costsplit_20260929.py` and their `.out` files, which were the assignment.
The authoring session committed four times while this audit ran (`3d4239156`, `ddad1c079`, `0d018330f`,
`b37902bfa`: the `h = 14` decomposition, the profile to level 18, reversal on the cycles of `3x + k`, and a
Proposition 6). That material is covered only by the quick checks of section J below; everything else refers to
the current text. Context used: the S19 note and its audit, `mazur.txt` (the paper's extracted text), the four
extracted paper texts (`paper2..5.txt`), and the arXiv HTML of Viaclovsky 2609.33785 (fetched).

**Method.** `04-computation/experiments/collatz_five_mirrors_20260929_audit.py` (2 minutes, < 2 GB; output
`collatz_five_mirrors_20260929_audit.out`; tags `[A1]`, `[H5]`, ... below refer to its lines) recomputes everything
with its own code and imports nothing from the session's scripts: the law by the FORWARD recursion
`z = 2^-a (3y + 1) mod 3^n` (exact rationals to level 6 using the period `L = 2 3^(n-1)` of `2^-a`; float64 to
level 14), its FFT profile, Proposition 2's inequalities, the Gauss sums with an exhaustive test of the stated
inequality, the same-length collisions with valuation caps 24 and 70 at every depth `<= 5`, the reciprocity
identity symbolically (sympy), the cycles and the `3x + 139` orbit, the cost decomposition at `h = 10` by an own
joint recursion, the cycle-sum identity of Proposition 6, and the decay models on `h <= 18`. Levels 15–18 of the
profile are taken from the session's `fourier_deep` / `fourier_deep18` outputs (8 and 20 GB; not rerun), which
agree with my level-14 values to `4·10^-7`.

---

## (i) Claims checked

### A. Proposition 1 and the consistency of the laws

1. **Consistency: `Y_n mod 3^h` has the law `mu_h`.** VERIFIED. `[A3]` exact for all `1 <= m < n <= 6`; `[A4]`
   float64 to `3.5·10^-14` for all `m < n <= 14`; `[A5]` `Y_n mod 3^h = F_h(a_(n-h+1), ..., a_n)` on 2000 random
   word pairs (`3 Y_k mod 3^h` depends only on `Y_k mod 3^(h-1)`), so the law of `Y_n mod 3^h` is that of `Y_h`
   for i.i.d. valuations. The citation "S19 §7 and the audit's item 13" is correct.
2. **Proposition 1 (lines 248–254).** VERIFIED. The proof is the two-line reduction; `[B2]` `max |M_n(h) -
   M_14(h)| = 1.1·10^-16` over `n = 10, 12`, same argmax up to sign.
3. **"the script confirms the same profile at `n = 4, 6, ..., 14, 17`" (line 253–254).** VERIFIED from the
   `.out` files and `[B4]` (`4.3·10^-7` against the deep output to `h = 14`); the level-18 output agrees too.
4. **`M(1) = 1/sqrt 3` exactly (line 279).** VERIFIED. `[B3]`; `|(1/3) e(1/3) + (2/3) e(2/3)|^2 = 5/9 + (4/9)
   cos(2π/3) = 1/3`.
5. **"the maximum is always at a power of two" (lines 279–280).** VERIFIED WITH CORRECTION: up to sign.
   `|mu_hat(-u)| = |mu_hat(u)|`; my FFT argmax at `h = 9` is `18659 = -2^10 mod 3^9` `[B1]` and the session's
   deep18 output has `19 = -2^3`, `65 = -2^4`, `1528787 = -2^16` at `h = 3, 4, 13`. Write `u = ±2^s`.
6. **"`s = h + 3` to `h + 5`" (line 66; line 280 "for `h >= 7`"; header line 30 and status line 530 "`2^(h+3)`,
   `2^(h+4)`").** REFUTED as a description of the table. From the argmax row and `[B6b]`: `s - h = 0` for
   `h <= 6`, then `1, 1, 1` (`h = 7–9`), `2, 2, 2` (`10–12`), `3, 3, 3` (`13–15`), `4, 4` (`16–17`), `5` (`18`):
   a step of one every two or three levels; "`h + 3` to `h + 5`" holds only for `h = 13..18`. Hence
   `2^s/3^h` is not "`≈ 0.03–0.05`" (line 281) but decreases geometrically: `0.117, 0.078, 0.052, 0.069, 0.046,
   0.031, 0.041, 0.027, 0.018, 0.024, 0.016` for `h = 7..17` `[B6]` (`0.022` at `h = 18`).

### B. Proposition 2 and the quotation of (2.3)

7. **(2.3) quoted correctly (lines 256–258).** VERIFIED. Mazur (2.2)–(2.3), `mazur.txt` 120–144: `rho_q = (2/3)
   3^q mu_q`, `||rho_q - rho_m ∘ π_(q,m)||_q <= (2/3) C_A m^-A` for `1 <= m <= q`, with the paper's own gloss
   "`C_A` bounds the `ℓ^1` distance between the probability law `mu_q` and the uniform lift of `mu_m` to `G_q`;
   the factor 2/3 comes from (2.2)". With `||.||_q` the mean over `G_q`,
   `mean_(G_q) |rho_q - rho_m ∘ π| = (2/3) ||mu_q - lift mu_m||_1` exactly `[C3]` (`0.367274 = 0.367274` at
   `q = 8, m = 3`), so the S20 form is the same statement; no normalisation factor reaches the conclusion.
8. **The lift has no primitive Fourier coefficients (lines 259–262).** VERIFIED. `[C1]` max `|lift^hat| = 0.0`
   at conductor `> 3^m` for `m = 1..7` at `N = 14` (the lift is constant on the fibres of reduction mod `3^m`;
   the note's `sum_(k=0)^2 e(uk/3) = 0` is the case `m = h - 1`).
9. **`|mu_hat_h(u)| <= ||mu_h - lift mu_(h-1)||_1` for primitive `u`.** VERIFIED. `[C2]` ratio at most `0.79`
   over `h <= 12`; `[C1]` `max_(cond > 3^m) |mu_hat_N| = M(m+1)` exactly, as Proposition 1 predicts.
10. **Conclusion `M(h) <= C_A (h-1)^-A`, `h >= 2`.** VERIFIED (`q = h`, `m = h - 1 >= 1`).
11. **"`d(m, N) >= M(m+1)` ... far from tight (`0.74` against `0.25` at `m = 2`)" (lines 310–311).** VERIFIED.
    S19's `N = 18` distance `0.739` against `M(3) = 0.2522`; at `N = 14` `[C1]`: `0.721` against `0.252`
    (ratio `2.9`, rising to `6.5` at `m = 7`).
12. **"Mazur's Lemma 8.1 is the case `A = 6`"; the quote "primitive Fourier coefficient bound `C* h^(-6409)` at
    level `h`" (lines 263–266; header lines 19–21).** VERIFIED WITH CORRECTION. `mazur.txt` 793–796: Lemma
    8.1 is `||rho_q - rho_m ∘ π||_q <= 2C/(3 m^6)`; the quoted phrase is in the paragraph "Quantitative
    dependence of the mixing estimate" following the lemma in §8.1 (lines 804–811), in full "the resulting
    primitive Fourier coefficient bound is `C* h^-6409` at level `h` for characters not factoring through level
    `h - 1`", i.e. a bound on `M(h)`. §2 attributes it correctly to "his Section 8.1"; the header attributes it
    to "Lemma 8.1". "`C*` a three-fold exponential tower" is the S19 audit's item 80 (not re-derived here).
13. **"Tao's Proposition 1.14 in Mazur's restatement", "from memory".** NOT CHECKED against Tao's text;
    consistent with the S19 note and audit (item 88 there).

### C. The profile, the masses, Parseval

14. **`M(h)`, `h <= 18` (lines 64–65, 273).** VERIFIED to `h = 14` (`[B1]`, `4.3·10^-7`); `h = 15..18` NOT
    RECOMPUTED (the session's FFTs at 8 and 20 GB; consistent with every fit below).
15. **Fourier mass per level `0.667, 0.476, 0.462, 0.464 .. 0.472`; typical `|mu_hat|^2 3^h = 0.70`.** VERIFIED.
    `[B1]` `0.66667, 0.47619, 0.46157, 0.46421, ..., 0.47025` (`h <= 14`); typical `0.692 → 0.705`.
16. **Parseval: "the second moment grows by `0.31` per level (S19), so the Fourier mass per conductor level is
    `(3/2)·0.31 = 0.466`, the table's constant" (lines 302–305; line 80).** VERIFIED WITH CORRECTION.
    `sum_t |mu_hat_n|^2 = 3^n sum mu_n^2 = (3/2) E_units[rho_n^2]` holds to `10^-5` at every level `[B8]`, and
    with Proposition 1 it is an exact level-by-level identity `mass(n) = (3/2)(E_n - E_(n-1))` (`0.4762 =
    0.4762`, ..., `0.4703 = 0.4703`). `E_n = 2.667, 3.911, 5.162` at `n = 6, 10, 14` as in S19; the increment
    is `0.308 → 0.3135`, slowly increasing, so the mass is `0.462 → 0.472`, not a constant.
17. **`|mu_hat_12(2^s)|`, `s = 0..40`.** VERIFIED. `[B5]` identical to the main output.
18. **"the generic size `3·10^(-4)` (`= 0.84 · 3^(-17/2)`)" (line 287).** REFUTED as a number. `0.84 · 3^(-8.5)
    = 7.4·10^-5` `[B7]`; the deep output prints `0.0000–0.0005` outside `s = 13..27`; `3·10^-4` would be
    `3.4 · 3^(-17/2)`.
19. **The bump at `s = 21` (`0.0125`), `0.0111` at `20`, `0.0106` at `23`; at level 18 centred at `s = 23`
    (`0.0112`) (lines 285–288).** VERIFIED from the two deep outputs (`0.0124` at `s = 22`, level 17: the
    centre is between `21` and `22`).
20. **Cost decomposition at `h = 10`, `t = 2^12` (lines 313–319).** VERIFIED. `[G1–G4]` own joint recursion:
    `|mu_hat| = 0.03828`, generic `t = 7`: `0.00261`; `A = 14: 0.0168`, `A = 16: 0.0144`, `A = 13: 0.0095`,
    `A = 15: 0.0072`, `A = 17: 0.0044`; `sum |.| = 0.0658`, coherent `0.581`; `a = 1..4`: `0.0225, 0.0120,
    0.0051, 0.0019`; the masses `0.0436`, `0.0764` are the exact `C(A-1, 9) 2^-A` `[G5]`; `|.|`-weighted
    `A/h = 1.493` `[G3]`. The `h = 14` numbers (lines 319–324) are read from `costsplit14` only (NOT RECOMPUTED,
    5 GB).
21. **The rate function (lines 326–329): `I(α) = α H(1/α) - α ln 2`, `I(1.48) = -0.093`, mass `e^(-0.093 h)`,
    ratio `0.91`.** VERIFIED as a Cramér rate: `α H(1/α) - α ln 2 = -(α ln 2 + (α-1) ln(α-1) - α ln α)`, and
    `I(1.48) = 0.0933`, `I(1.35) = 0.1632` `[G6]` (the note's sign convention is the exponent of the mass). At
    finite `h` the exact `-ln P(A = 1.35h)/h` is `0.31, 0.29, 0.27, 0.22` at `h = 10, 14, 20, 40` `[G7]`
    (subexponential prefactors), so the "`≈`" is a rate, not a probability.

### D. The interpretive claims on the decay

22. **Local exponents `1.76, 1.86, 1.99, 2.08, 2.10` (lines 70, 290–291).** VERIFIED. `[H3]` `1.756, 1.861,
    1.988, 2.079, 2.102`.
23. **"A fixed power law is excluded by the steepening" (lines 69–70, 291–292).** VERIFIED WITH CORRECTION.
    Only a pure power law `C h^-α` is excluded (best fit `α = 2.06` on `7..18`, max residual `4.5%` `[H4]`).
    A shifted power law `C (h + 2.5)^-2.52` fits `h = 7..18` to `1.8%` max / `1.0%` rms `[H5]` — better than
    the geometric `C r^h`, `r = 0.862`, on `11..18` (`4.3%` / `2.7%` `[H6]`) — and reproduces the steepening
    exactly (local exponents `1.86, 1.94, 2.01, 2.06, 2.10`) and the ratios at `14..18` (`0.854, 0.862, 0.869,
    0.876, 0.881`); so does `C h^-1.72 (0.971)^h` (`2.1%` `[H7]`). A rising local exponent is what any
    polynomial law with a shift produces; the sentence is literally true and does not support its use.
24. **"the ratio ... rises from `0.65` to about `0.87–0.89` and levels off over `h = 14..18`" (lines 67–69);
    "are flat over `14..18` at `0.87 ± 0.02`" (lines 293–294).** REFUTED as a description. Window means
    `0.781, 0.823, 0.868, 0.873` over `h = 6–9, 10–13, 14–17, 14–18`; least-squares slope `+0.009` per level
    over `11–18` `[H2]`; `1 - ratio` runs `0.203 → 0.106` over `10..18` `[H10]`; the ratio at `h = 18`
    (`0.894`) is the largest so far. No plateau; the shifted power law predicts exactly this slow rise.
25. **"a geometric decay at rate about `0.87–0.89` per level fits levels `11–18` and would satisfy (2.3) with
    room" (lines 71–73); "compatible with Proposition 2 (super-polynomial) with room" (lines 296–297); §6 "if
    the ratio stays below `1`, the decay is geometric and (2.3) holds with a usable constant" (lines 497–499).**
    REFUTED as inferences. The data to `h = 18` are better fitted by `C (h + 2.5)^-2.5`, which violates the
    conclusion of Proposition 2 for every `A >= 3`; the two readings separate by a factor `2` only at `h ≈ 31`
    `[H9]` (`3^31` classes: out of reach of the FFT). The new level 18 is an out-of-sample test that goes the
    wrong way for the note: fitted on `7..17` / `11..17`, the shifted power law predicted `M(18) = 0.01096`
    (ratio `0.880`) and the geometric law `0.01047` (ratio `0.857`); measured `0.01119` (ratio `0.894`) `[H5b]`.
    The §6 conditional is empty as written (a ratio below `1` tending to `1` is polynomial decay; it must read
    "below a fixed `r < 1`"). The note's own line 299–300 ("OBSERVED, eighteen levels; the asymptotic regime is
    not established, and a ratio creeping to `1` is not excluded") is the honest sentence and contradicts
    lines 67–73, 293–297 and 497–499.
26. **"`3.75 h^(-2)` fits `7 <= h <= 14` to `4%` and undershoots by `2–4%` at `15–18`" (lines 292–293).**
    VERIFIED WITH CORRECTION. `[H8]` residuals `+0.9% .. -3.8%` on `7–14`; at `15–18` the fit lies `+2.3%,
    +1.7%, +3.7%, +3.5%` ABOVE `M(h)`: `M` undershoots the fit, not the fit `M`.
27. **§0.3 (lines 88–95, no label): "the low-cost words ... are the obstruction to fast sup-norm mixing", "the
    resonance is a fixed large-deviation family ... the observed `0.87` is that rate times the coherence loss";
    §2 lines 329–335 "the observed sup-norm ratio `0.87` is this rate times a coherence loss of about `0.96`
    per level ... whose rate (`≈ 0.87` per level) is what a sharp version of (2.3) would have to compute".**
    NOT ESTABLISHED: heuristic stated as mechanism. The decompositions at `h = 10, 14` are FINITE-EXACT; the
    "coherence loss `0.96`" is defined as `0.87/0.91` and explains nothing; that the sup-norm ratio equals a
    large-deviation rate times anything is not derived. Label DIRECTION (and note that under the polynomial
    reading of item 23 the "rate" is not a constant at all).
28. **"its density is in no `L^p` for `p >= 2` if this continues (S19 §7, OBSERVED)", "`1/f`" (lines 305–307).**
    VERIFIED as labelled.

### E. Proposition 3 (Gauss sums)

29. **`G_j(t) = c_j sum_(r=1)^(L_j) 2^-r e(t 2^-r/3^j)`, `c_j = 1/(1 - 2^-L_j)`, as the one-step factor (lines
    386–389).** VERIFIED (group `a` by its class mod `L_j`; the frequency recursion itself holds to `2.7·10^-14`
    for eight `t` including non-units, all `a` summed `[A6]`).
30. **(i) `2^-r mod 3^j = (1 + m_r 3^j)/2^r`, `m_r = -3^-j mod 2^r` (lines 391–394).** VERIFIED: `1 + m_r 3^j
    ≡ 0 mod 2^r`, the quotient is an integer in `(0, 3^j)` congruent to `2^-r`; numerically to `10^-10` on ten
    random units per `j <= 12` `[D1]`.
31. **The table (lines 397–401): `sup |G_j| = 0.577, ..., 0.9997` at `t = 1, 8, 8, 16, ..., 18659, ...`; mean
    `0.5430`; `2.9%` above `0.9`.** VERIFIED `[D1]` (`0.57735, 0.58175, 0.78935, 0.88682, 0.94362, 0.97140,
    0.98589, 0.99371, 0.99674, 0.99852, 0.99928, 0.99968`; mean `0.5430` from `j = 7`; share `0.0286`). The
    argmax is `±2^(j+1) mod 3^j` in every case (my scan returns the conjugate representatives `1, 1, 19, 16,
    ...` of the note's `1, 8, 8, 16, ...`); §0.4's "at `t = 2^(j+1)`" is right up to sign.
32. **(ii) "`|G_j(2^s)| >= 1 - 2π 2^(-K) - 2^(-s)` when `2^s <= 3^j/2^K`" (lines 394–396).** REFUTED as
    stated. Exhaustive test over `j <= 12`, `s <= 40`, all admissible `K`: 178 violations, at `s = 1..7`
    (`62, 46, 32, 20, 12, 5, 1`), the worst `|G_9(2)| = 0.207` against the bound `0.499` (`K = 13`) `[D2]`;
    `|G_j(2)| = 0.21–0.38` and `|G_j(4)| = 0.60–0.67` for `j = 4..12` `[D3]`. The head `sum_(r<=s)` has modulus
    at least `(1 - 2^-s) cos(2π 2^-K)` and the tail `sum_(r>s)` at most `2^-s`, so the provable bound is
    `1 - 2π 2^-K - 2^(1-s)` (zero violations). The conclusion `sup_t |G_j(t)| -> 1` survives (let `s, K -> ∞`
    with `2^(s+K) <= 3^j`); `c_j >= 1` only helps and the truncation at `r = 60` is immaterial.
33. **"Mixing is a multi-step (renewal) phenomenon, as in Tao's proof" (lines 101–103).** The "no uniform
    one-step gap" is PROVED (item 32); "as in Tao's proof" NOT CHECKED (consistent with Mazur's
    "Fourier–renewal proof", `mazur.txt` 804).

### F. Proposition 4 (same-length spread) and the collision table

34. **The inequality chain (lines 356–362).** VERIFIED: `Y(w) ≡ Y(w') mod 3^n` gives `3^n | C_w 2^A' - C_w'
    2^A`; the difference is nonzero (item 36); `0 < |.| < 2^(A+A') 3^d/2` from `C_w <= 2^A (3^d - 1)/2`; hence
    `3^(n-d) < 2^(A+A'-1)`.
35. **`C_w` odd; `C_w <= 2^A (3^d - 1)/2`; `Y(w) = C_w 2^-A`.** VERIFIED `[E4]` (3000 random words; the sharp
    bound is `C_w <= 2^(A-d)(3^d - 2^d)`).
36. **"`(d, C_w)` determines `w` (`a_1 = v_2(C_w - 3^(d-1))`, then recurse)" (lines 359–360).** REFUTED as
    stated; proof repairable. `C_w` does not involve `a_d`: `C_(1,2) = C_(1,3) = C_(1,9) = 5` `[E5]`. The
    recursion recovers `a_1, ..., a_(d-1)` and `a_d = A - (a_1 + ... + a_(d-1))` `[E6]`, so `(d, C_w, A)`
    determines `w`; in the proof `A = A'` is already established (both carries odd), so the argument stands.
37. **"the corresponding depth-`d` predecessors of any two integers in one class mod `3^n` are distinct nodes"
    (lines 352–354).** VERIFIED but vacuous as phrased (predecessors of distinct targets are distinct integers
    regardless); the content is that their residues `Y(w) mod 3^n` are distinct for low-cost words.
38. **The table (lines 364–372; valuations `<= 24`; first colliding depth `2, 2, 3, 3, 3, 3, 4, 4, 4`; minimal
    `A + A' = 29, 31, 32, 32, 56, 56, 58, 61, 63`).** VERIFIED as defined `[E1]` (identical, same classes), and
    the values are the true minima at those depths: unchanged with valuations `<= 70` `[E2]`, and any pair with
    a valuation `> 70` has `A + A' >= 70 + 3d` `[E3]`.
39. **"The actual first collisions in the tree of `1` need cost sums `29..63` at depths `2..4`" (lines
    107–108) and the table's definition "minimal `A + A'` among the colliding pairs at the first colliding
    depth" (lines 365–366).** REFUTED as a statement about the tree; the stated definition produces a cap
    artefact. `[E2]`: with valuations `<= 70` the first colliding depth is `2` for `n = 6..9` (sums `49, 73,
    99, 99`) and `3` for `n = 10..12` (`70, 80, 127`; the last two are upper bounds), and the minimal colliding
    sum at depth `d` DECREASES with `d` for `n >= 8` (`n = 8`: `56, 43, 44` at `d = 3, 4, 5`; `n = 10`: `58,
    50`; `n = 12`: `63, 59` at `d = 4, 5`). The quantity to tabulate is "minimal colliding `A + A'` at depth
    `d`" for each `d`.
40. **"three to five times the bound" (line 108), "loose by a factor `3–5`" (line 374).** REFUTED as numbers.
    The table's own ratios are `7.0, 5.4, 5.6, 4.4, 6.3, 5.3, 5.5, 5.0, 4.6` (`29/4.17` ... `63/13.7`); over all
    depths `2..5` and `n <= 12` `[E2]` the factor is never below `4.4` (`32/7.34` at `n = 7`, `d = 3`) and
    reaches `9–13` where the bound is small (`n = 5..7`, `d = 4, 5`; at `n = 4`, `d >= 4` the bound is `<= 1`).
41. **"the mixing of section 2 happens only once the exponentially many words of typical cost exhaust the
    `2·3^(n-1)` units (depth `≈ n log 3/log(4/3) = 3.8 n`)" (lines 378–380).** REFUTED. No derivation is
    given; `n log 3/log(4/3)` is the depth at which `(3/4)^d = 3^-n` (a size statement about typical
    predecessors), not a residue statement. The entropy count (`4^d` effective words of length `d`) gives
    `d ≈ n log 3/log 4 = 0.79 n`, and directly `[E7]`: the depth-`d` layers (valuations `<= 70`) hit every unit
    class mod `3^n` from `d = 2` (`n = 4`), `d = 3` (`n = 5, 6, 7`), `d = 4` (`n = 8`) on, i.e. from `d ≈ n/2`,
    against "`3.8 n`" (`15, 19, 23, 27, 31`).

### G. Proposition 5 (carry reciprocity) and the reversed cycles

42. **`C_(rev w)(u, v) = u^(d-1) v^A C'_w(1/u, 1/v)` (lines 413–418).** VERIFIED as a polynomial identity
    (sympy, 120 random words `[F1]`); the index substitution `i = d + 1 - j` is correct.
43. **`w = (1,2)`: `C_w = 5`, `C'_w = 14`, `C_rev = 7 = 3·8·C'_w(1/3, 1/2)` (lines 420–421).** VERIFIED `[F2]`.
44. **The seven-cycle: orbit, word `(1,1,1,2,1,1,4)`, `C_w = 2363`, `y_0 = -17`; reversed `(4,1,1,2,1,1,1)`,
    `C_rev = 13801`, `y_0' = -13801/139`, not a rotation, not an integer; `{-5,-7}` a rotation with `y_0' =
    -7`; `-1` fixed (lines 431–435).** VERIFIED `[F3]` (orbit `-17 → -25 → -37 → -55 → -41 → -61 → -91 → -17`;
    `139 = 3^7 - 2^11`, `gcd(13801, 139) = 1`).
45. **"a rational cycle of `3x+1`... i.e. of the map `x -> (3x + 139)/2^v` on integers" (lines 433–435).**
    VERIFIED `[F4]`: `-13801 → -2579 → -3799 → -5629 → -4187 → -6211 → -9247 → -13801`, word `(4,1,1,2,1,1,1)`.
46. **"the 2-adic class of the source, `x ≡ -C_w 3^(-d) mod 2^(A+1)`" (line 424; line 115–116 without
    modulus).** REFUTED as stated. `3^d x = 2^A y - C_w` with `y` odd gives `x ≡ -C_w 3^-d mod 2^A` (2000
    sources `[F5]`) and `x ≡ (2^A - C_w) 3^-d mod 2^(A+1)`, never `-C_w 3^-d mod 2^(A+1)` (the S19 audit's item 1
    has the right class).
47. **"Reversal is an involution on rational cycles that preserves `(k, A)`" (lines 117–118).** VERIFIED
    (trivial).
48. **"Reversal therefore does not preserve integrality ... one more way the negative cycles are special objects
    of the 3-adic law (S19, Theorem B)" (lines 435–438).** VERIFIED as a finite fact; Theorem B correctly
    cited; the reading is DIRECTION.

### H. The five mirrors (section 1)

49. **Viaclovsky summary (lines 143–153).** VERIFIED against the arXiv HTML: `y^2 = x^3 + t x + 1` (1.1);
    fibres `III*, I_1, I_1, I_1` (abstract); Mordell–Weil infinite cyclic generated by `P = (0, 1)` (Lemma
    2.1); `M = O_S(P - O)`, `M^×` its complement of the zero section; `Φ = λ Φ_0`, `λ` small, lifting
    translation by `-2P` (2.8); Mumford fillings (Prop. 3.1); the four-torsion section `a = ℓ/4` and `Y_∞ =
    A/<ρ~>` (3.12), one multiple fibre of multiplicity four; `T_* T_1 T_2 T_3 = I`, `T_*^4 = I`, `(T_i - I)^2 =
    0` (4.12); "Primitivity forces `d^2 = 1`" (Prop. 4.2). The wording "an order-four free affine action" NOT
    CHECKED literally.
50. **Viaclovsky mirror (lines 153–165).** Fair as labelled ("shape only; nothing to compute"). `S_a` is a
    contraction of `Z_3` with ratio `|3|_3 = 1/3` and the 3-adic law is the stationary law of the IFS with
    weights `2^-a` (the S19 reference should be §2/§7 rather than "§4", minor). "The inverse histories are the
    quotient by the lift" and "Terras primitivity is the primitivity that pins `d`" are analogies without
    content.
51. **Narode (lines 167–170): `x^4 + 3x^2 + 1` (`g = t^2 + 1`, `f(1) = f(-1) = 5`, `Δ = 400`, monogenic)
    refutes the necessity of "`f(1) f(-1)` squarefree"; Theorem 1.4 as stated.** VERIFIED (`paper2.txt`,
    Example 1.3, Theorem 1.4).
52. **Narode mirror (lines 172–178).** Fair: one exact identity (Proposition 5), one finite fact, and the
    honest "no degree-halving analogue was found".
53. **Merca (lines 180–186): `a(n)` the signed count (sign `(-1)^(r_1 + r_2)`) of `n = (2r_1+1)x + (2r_2+1)y`,
    `x < y`; the double Lambert series; Theorem 1.1 `a(n) = (σ(n) - λ(n))/4 - σ(n/2)/2 + σ(n/4)`; `a(n) >= 0`;
    Theorem 1.3 `A = P·S`, `S = sum_(r>=1) (-1)^(T_(r-1)) T_(⌈r/2⌉) q^(T_(r+1))`; Theorem 1.6 `a(12n + 11) ≡ 0
    mod 3`.** VERIFIED (`paper3.txt`). "theta series over triangular numbers with signs `(-1)^(T_(r-1))`" omits
    the coefficients `T_(⌈r/2⌉)` (Merca: "theta-type correction series"), minor; "Jacobi's two- and four-square
    counts" is the proof's route, correct.
54. **"the analogue of Merca's divisor recursion" (line 197).** REFUTED as a reference: Merca has a divisor-sum
    FORMULA (Theorem 1.1) and a linear RECURRENCE from the `pod` theta identity (Corollary 1.5, `sum_r
    (-1)^(T_r) a(n - T_r) = ...`); there is no "divisor recursion".
55. **Merca mirror, "Verdict: the frame is right" (lines 188–199).** OVERREACH. `mu_n(y) = sum_w 2^(-A(w))
    [Y_n(w) ≡ y]` is positive by definition, so Merca's phenomenon (positivity invisible in a signed
    definition) has no counterpart; the Fourier expansion is a signed representation, not a definition; and
    the open `liminf H_n(1) > 0` is an asymptotic statement about a positive sequence. Honest verdict: remark.
56. **Dalfó–Fiol–Reyes (lines 201–207): steps `(+a,+b), (+a,-b), (-a,+b)`, not `(-a,-b)`; sector-constrained
    non-symmetric distance; `3ℓ + 1` vertices at distance `ℓ`; `N(k) = (3k^2 + 5k + 2)/2`; not attained for
    `k > 1` because the optimal tiles do not tessellate; `u_i a + v_i b ≡ 0 mod N`, `N = det`.** VERIFIED
    (`paper4.txt`, §1–2, (1)).
57. **DFR mirror: "the frequency-side recursion ... with the sector restriction that the parity of `a` is fixed
    by the class of the target mod `3`" (lines 209–214).** REFUTED as stated. The frequency recursion sums over
    ALL `a >= 1` (`[A6]` verifies it so); the parity restriction `a ≡ ε(z) mod 2` belongs to the space-side
    parents' recursion of the law (S19 §5). The "sector" of the analogy sits on the wrong side of the
    transform. The rest (exponential ball growth `2^(A-1)`; "the Moore bound is Proposition 4"; "it produced
    the Fourier profile") is fair, with items 39–41 applying to "far from attained" and "no tiling obstruction
    at all" (a metaphor, unlabelled).
58. **Lyu (lines 221–228): `τ(H)`, `p(H)`, `p <= τ <= r_tr(p)`; Theorem 1.1 (connected, girth `> g`, `ω_ao =
    ω_ro = 3`, `τ = r_tr(p)`); Corollary 1.2 (a forbidden `F` with a cycle in `U(F)` gives no polynomial bound,
    so only forests can); Theorem 1.3(2): `F_alt`-free gives `τ <= 3p - 2`, `= p` when triangle-free.** VERIFIED
    (`paper5.txt`).
59. **Lyu mirror (lines 228–237).** Fair as "remark". The atlas rows it points to are NOT CHECKED here; "words
    with all `a_j <= 2` and mean below `log_2 3` exist" is trivially true of words.

### I. Status hygiene and cross-references

60. **Labels PROVED for Propositions 1–5 (lines 24–27).** VERIFIED for 1, 2, 4, 5 (with the proof repairs of
    items 36 and 46, which do not touch the statements); Proposition 3(ii) is PROVED only after the constant is
    corrected (item 32); its stated inequality is false.
61. **"FINITE-EXACT (the Fourier profile of the law to level 17 in float64 ...)" (lines 27–28), "Computed to
    `h = 18` (FINITE-EXACT, float64)" (line 63), "FINITE-EXACT, float64 FFT" (line 268), status line 530.**
    VERIFIED WITH CORRECTION: under the S19 audit's convention (its C25) a float64 computation is VERIFIED, not
    FINITE-EXACT; here it is VERIFIED against exact rationals to level 6 and an independent forward recursion
    to level 14 (`4·10^-7`). The header still says "level 17" and the status table "`h <= 17`", "checked at
    levels 4–17" (lines 528, 530) after the level-18 extension.
62. **OBSERVED stated as established.** Lines 67–73 ("levels off", "would satisfy (2.3) with room"); 88–95
    (the unlabelled mechanism sentence); 293–297 ("flat", "compatible ... with room"); 329–335 ("coherence loss
    of about `0.96`", "whose rate (`≈ 0.87` per level)"); 107–108 and 374–380 ("the actual first collisions",
    "three to five times", "no tiling obstruction at all", "depth `≈ 3.8 n`"); 497–499 (the empty conditional).
    Items 24, 25, 27, 39–41.
63. **Numbers trace to output files.** VERIFIED except: "`3·10^(-4)`" (no source; wrong, item 18), "three to
    five" (wrong, item 40), "`3.8 n`" (no source, item 41), "`0.96`" (defined as a quotient, item 27); the
    rate `e^(-0.093 h)` and "`0.74` against `0.25`" are computed correctly (items 21, 11).
64. **S19 cross-references.** VERIFIED: Theorem A (harmonic mass), Corollary A2 (`H_(n+1) = H^(1) + 2H^(2)`),
    Theorem B (cycle resonances), the consistency (S19 §7; audit item 13), the second-moment slope `0.31` (S19
    table `2.667, 3.911, 5.162, 6.419`; increments `0.310–0.314`), the `ℓ^1` distances `0.854 .. 0.354` for
    `m = 1..9` at `N = 18`, "`C` a three-fold tower" and "`m` of the order of `2^(2^(2^8697))`" (S19 audit item
    80). "S19 §4" for the IFS is loose (item 50). The barrier atlas §§0, 7, 8 and THM-4514/THM-4476 as cited
    NOT CHECKED.
65. **Reproduction lines (§7).** VERIFIED for the main and costsplit scripts (items 3, 17, 20, 31, 38, 43–45
    rerun their content independently); the deep scripts at 17/18 and the two new scripts NOT RERUN.

### J. Material added during the audit (commits `3d4239156` .. `b37902bfa`), quick checks only

66. **Level 18: `M(18) = 0.0112`, ratio `0.894`, argmax `2^23`, bump centred at `s = 23`.** VERIFIED as read
    from `fourier_deep18_20260929.out` (`0.011187`, `0.8942`, `8388608 = 2^23`); NOT RECOMPUTED. Its effect
    on the note's reading is item 25: it is the largest ratio so far and lands on the shifted power law's
    prediction, not the geometric one.
67. **The `h = 14` decomposition (lines 319–330; `costsplit14` output).** VERIFIED as read (`A = 21: 0.0083`,
    coherent `0.560`, `|.|`-weighted `A/h = 1.481`, `a = 1: 0.0110`); NOT RECOMPUTED; the interpretation is
    item 27.
68. **Reversal on the integer cycles of `3x + k` (lines 440–465; `reversal` output).** Spot-checked by hand:
    `3x+13`: `227 → 347 → 527 → 797 → 601 → 227` (word `(1,1,1,2,3)`) and `259 → 395 → 599 → 905 → 341 → 259`
    (word `(1,1,1,3,2)`, a rotation of `(3,2,1,1,1) = rev(1,1,1,2,3)`); `3x+37`: `23 → 53 → 49 → 23` (`(1,2,3)`)
    and `29 → 31 → 65 → 29` (`(2,1,3)`, a rotation of `(3,2,1)`): the two named pairs are genuine reversal
    pairs. "Every cycle of `3x+1` is self-dual except the seven-cycle": VERIFIED for the four known cycles
    (`(2)`, `(1)`, `(1,2)` are rotation-symmetric). The completeness of the survey (all cycles reached from
    `|x| <= 5·10^4`, cap `10^15`) and the other `k` NOT CHECKED. "`3(kx) + k = k(3x + 1)`" is right.
69. **Proposition 6 (lines 467–486): `S_(rev w) = S_w` for the cycle-sum polynomial `S_w = sum_j C_(rot_j
    w)`; equal traces `2499`, `125`, `-327`.** VERIFIED. The proof's mechanism is right (the coefficient of
    `u^(d-i)` in `S_w` is `sum_j v^(b_(j,i-1))` over the cyclic blocks of length `i - 1`, and reversal maps
    cyclic blocks of `w` to cyclic blocks of `rev w` with the same sums); `[F6]` the identity on 300 random
    words with rational `(u, v)`, and the seven-cycle trace `-327 = -327` from the two words; `227 + 347 + 527 +
    797 + 601 = 259 + 395 + 599 + 905 + 341 = 2499` and `23 + 53 + 49 = 29 + 31 + 65 = 125` by hand. "The sum of
    squares is not invariant" NOT CHECKED. The closing sentence "a strong constraint on any candidate family"
    is DIRECTION.

---

## (ii) Corrections (old text -> new text; line numbers of commit `b37902bfa`)

**C1 (lines 66–67; 280–281; header 29–30; status 530; the wiring's "`2^(h+3)`").** "the maximum sits at `u =
2^s` with `s = h + 3` to `h + 5` (`2^17` at `h = 14`, `2^21` at `h = 17`, `2^23` at `h = 18`)" -> "the maximum
sits at `u = ±2^s` with `s - h = 0` for `h <= 6` and `s - h = 1, 2, 3, 4, 5` on `h = 7–9, 10–12, 13–15, 16–17,
18` (a step every two or three levels, OBSERVED; `2^17` at `h = 14`, `2^21` at `h = 17`, `2^23` at `h = 18`)".
Line 280–281: "The maximum is always at a power of two, `u = 2^s` with `s = h + 3` to `h + 5` for `h >= 7`, i.e.
`t/3^h = 2^s/3^h ≈ 0.03–0.05`" -> "The maximum is always at a power of two up to sign (`|mu_hat(-u)| =
|mu_hat(u)|`), `u = ±2^s` with `s - h` stepping up by one every two or three levels (`s - h = 1` at `h = 7`, `5`
at `h = 18`), so `2^s/3^h` decreases geometrically, `0.12` at `h = 7`, `0.027` at `h = 14`, `0.022` at `h =
18`". Header: "the maxima sit at the powers of two `2^(h+3)`, `2^(h+4)`" -> "the maxima sit at the powers of
two `±2^s`, `s - h` growing slowly". Status line 530: "argmax at `2^(h+3)`, `2^(h+4)`" -> "argmax at `±2^s`,
`s - h = 0 .. 5`".

**C2 (lines 67–73).** "the ratio `M(h)/M(h-1)` rises from `0.65` to about `0.87–0.89` and levels off over `h =
14..18` (`0.867, 0.851, 0.885, 0.868, 0.894`). A fixed power law is excluded by the steepening (local exponents
`1.76, 1.86, 1.99, 2.08, 2.10` on the doublings `5→10, 6→12, 7→14, 8→16, 9→18`); a geometric decay at rate about
`0.87–0.89` per level fits levels `11–18` and would satisfy (2.3) with room." -> "the ratio `M(h)/M(h-1)` rises
from `0.65` to `0.89` with no plateau (window means `0.78, 0.82, 0.87` over `h = 6–9, 10–13, 14–18`; slope
`+0.009` per level over `11–18`; `0.894` at `h = 18` is the largest yet). The steepening local exponents (`1.76,
1.86, 1.99, 2.08, 2.10` on the doublings `5→10, ..., 9→18`) exclude a pure power law `C h^(-α)` but not a
shifted one: `C (h + 2.5)^(-2.5)` fits `h = 7–18` to `1.8%` (a geometric `C r^h`, `r = 0.86`, fits `11–18` to
`4.3%`), reproduces both the steepening and the ratios, and predicted the level-18 value from levels `<= 17`
(`0.0110` against the geometric `0.0105`; measured `0.0112`); the two readings separate by a factor `2` only
near `h = 31`. Eighteen levels therefore do not decide between a geometric decay (which would give (2.3) with a
usable constant) and a polynomial one (which would contradict Proposition 2); OBSERVED."

**C3 (lines 291–297).** "steepening, so no fixed power law; `3.75 h^(-2)` fits `7 <= h <= 14` to `4%` and
undershoots by `2–4%` at `15–18`. The ratios rise from `0.65` to `0.87` and are flat over `14..18` at `0.87 ±
0.02` (`0.894` at `h = 18`; they fluctuate with the `s = h+3`/`h+4`/`h+5` alternation of the argmax). Reading: a
geometric decay at rate about `0.87–0.89` per level, which is compatible with Proposition 2 (super-polynomial)
with room, and with Tao's proposition, whose constants say nothing below `m` of the order of `2^(2^(2^8697))`."
-> "steepening, so no pure power law `C h^(-α)`; a shifted power law `C (h + 2.5)^(-2.52)` fits `7 <= h <= 18`
to `1.8%` with local exponents `1.86, 1.94, 2.01, 2.06, 2.10` and ratios `0.854, 0.862, 0.869, 0.876, 0.881` at
`14..18`; `3.75 h^(-2)` fits `7 <= h <= 14` to `4%` and `M(h)` lies `2–4%` below it at `15–18`. The ratios rise
from `0.65` to `0.89` and are still rising over `14..18` (mean `0.873` against `0.823` over `10..13`; the `±0.02`
fluctuation follows the steps of `s - h`). Reading: either a geometric decay near `0.87` per level or a
polynomial decay of exponent about `2.5` with a shift; the data to `h = 18` cannot separate them (factor `2` at
`h ≈ 31`), and only the first is compatible with Proposition 2. Tao's proposition says nothing below `m` of the
order of `2^(2^(2^8697))`." (Lines 299–300 then stand as the conclusion.)

**C4 (lines 497–499).** "The measured `M(h)` to `h = 17` decays with a ratio rising to `0.87`; if the ratio
stays below `1`, the decay is geometric and (2.3) holds with a usable constant, which no published proof
provides" -> "The measured `M(h)` to `h = 18` decays with a ratio rising to `0.89`; if the ratio stays below a
fixed `r < 1`, the decay is geometric and (2.3) holds with a usable constant, which no published proof provides;
if it keeps rising to `1` as the equally good shifted power law `(h + 2.5)^(-2.5)` predicts, (2.3) fails".

**C5 (line 396).** "so `|G_j(2^s)| >= 1 - 2π 2^(-K) - 2^(-s)`" -> "so `|G_j(2^s)| >= (1 - 2^(-s)) cos(2π
2^(-K)) - 2^(-s) >= 1 - 2π 2^(-K) - 2^(1-s)` (head of modulus at least `(1 - 2^(-s)) cos(2π 2^(-K))`, tail at
most `2^(-s)`; with `2^(-s)` in place of `2^(1-s)` the bound fails, e.g. `|G_9(2)| = 0.207`)".

**C6 (lines 359–360).** "and `(d, C_w)` determines `w` (`a_1 = v_2(C_w - 3^(d-1))`, then recurse on `(C_w -
3^(d-1))/2^(a_1)`)" -> "and `(d, C_w, A)` determines `w` (`C_w` does not involve `a_d`: `a_1 = v_2(C_w -
3^(d-1))`, recurse on `(C_w - 3^(d-1))/2^(a_1)` for `a_2, ..., a_(d-1)`, then `a_d = A - a_1 - ... - a_(d-1)`)".

**C7 (line 424; lines 115–116).** "`x ≡ -C_w 3^(-d) mod 2^(A+1)`" -> "`x ≡ -C_w 3^(-d) mod 2^A` (the class mod
`2^(A+1)` is `(2^A - C_w) 3^(-d)`, `y` being odd)"; line 115–116 "the 2-adic source `x ≡ -C_w 3^(-d)`" -> "the
2-adic source `x ≡ -C_w 3^(-d) mod 2^A`".

**C8 (lines 107–111).** "The actual first collisions in the tree of `1` need cost sums `29..63` at depths
`2..4` for `n = 4..12`, three to five times the bound: the "tessellation obstruction" of the paper (optimal
tiles that fail to tile) has no counterpart; the Syracuse residues stay injective far beyond the Moore regime."
-> "In the tree of `1` with valuations `<= 24`, the first same-length collisions mod `3^n` appear at depths
`2..4` with minimal cost sums `29..63` for `n = 4..12`, `4.4–7` times the bound (with valuations `<= 70` there
are collisions at depth `2` already for `n <= 9` and at depth `3` for `n <= 12`, at larger sums; and the minimal
colliding sum at depth `d` decreases with `d`: `n = 8`: `56, 43, 44` at `d = 3, 4, 5`; `n = 12`: `63, 59` at `d
= 4, 5`); at every depth `d <= 5` and `n <= 12` the minimal colliding `A + A'` exceeds the bound by a factor of
at least `4.4`: the Syracuse residues stay injective well beyond the Moore regime, and no counterpart of the
paper's tessellation obstruction was found (an analogy, not a theorem)."

**C9 (lines 364–380).** Redefine the table as "minimal colliding `A + A'` at depth `d` (valuations `<= 70`;
exact for every entry `<= 70 + 3d`)" with rows `d = 2..5` from `[E2]` (`n = 4`: `29, 22, 26, 29`; `5`: `31, 32,
29, 29`; `6`: `49, 32, 29, 33`; `7`: `73, 32, 36, 40`; `8`: `99*, 56, 43, 44`; `9`: `99*, 56, 55, 44`; `10`: `-,
70, 58, 50`; `11`: `-, 80*, 61, 52`; `12`: `-, 127*, 63, 59`; `*` = upper bound) and the bounds `(n - d) log_2 3
+ 1`; line 374 "The bound is loose by a factor `3–5`" -> "The bound is loose by a factor of at least `4.4`
(`4.4–7` on the tabulated first depths)"; lines 378–380 "and the mixing of section 2 happens only once the
exponentially many words of typical cost exhaust the `2·3^(n-1)` units (depth `≈ n log 3/log(4/3) = 3.8 n`)" ->
delete, or "(the depth-`d` layers with valuations `<= 70` cover every unit class mod `3^n` from `d ≈ n/2` on for
`n <= 8`; the entropy count `4^d ≈ 3^n` gives `d ≈ 0.79 n`)".

**C10 (line 287).** "and the generic size `3·10^(-4)` (`= 0.84 · 3^(-17/2)`) outside `s = 13..27`" -> "and the
generic size `7·10^(-5)` (`= 0.84 · 3^(-17/2)`) outside `s = 13..27`".

**C11 (lines 303–305; line 80).** "the second moment grows by `0.31` per level (S19), so the Fourier mass per
conductor level is `(3/2)·0.31 = 0.466`, the table's constant" -> "with Proposition 1 the mass at level `h` is
exactly `(3/2)(E_units[rho_h^2] - E_units[rho_(h-1)^2])`; the second-moment increment is `0.308 → 0.314` over
`h = 3..17` (S19: `0.31`), so the mass is `0.462 → 0.472`, slowly increasing".

**C12 (lines 88–95; lines 329–335).** Line 88–89 "the low-cost words, exponentially rare but phase-coherent,
are the obstruction to fast sup-norm mixing." -> "... (FINITE-EXACT at `h = 10` and `14`); the reading that the
low-cost words are the obstruction to fast sup-norm mixing is DIRECTION."; lines 93–95 "the resonance is a fixed
large-deviation family, whose mass rate at `A/h = 1.48` is `e^(-0.093 h)` (ratio `0.91` per level); the
observed `0.87` is that rate times the coherence loss." -> "the resonant cost ratio is stable at two levels;
the cost family `A ≈ 1.48 h` has large-deviation mass `e^(-0.093 h)` (ratio `0.91` per level); whether the
sup-norm ratio is this rate times a coherence factor is DIRECTION (the quotient `0.87/0.91 = 0.96` is a
definition, not a mechanism), and under a polynomial reading of the profile there is no constant rate at all."
Lines 329–335: "the observed sup-norm ratio `0.87` is this rate times a coherence loss of about `0.96` per
level. This is the mechanism behind the slow sup-norm decay ... whose rate (`≈ 0.87` per level) is what a sharp
version of (2.3) would have to compute." -> "the observed sup-norm ratio `0.87` would be this rate times a
coherence factor of about `0.96` per level if the decay is geometric (DIRECTION). This is the candidate
mechanism ... whose contribution is what a sharp version of (2.3) would have to bound."

**C13 (lines 211–214).** "is a twisted geometric walk on the cyclic unit group `<2> ≅ Z/(2·3^(n-1))` — a
weighted circulant with generator `2^(-1)` — with the sector restriction that the parity of `a` is fixed by
the class of the target mod `3`; the ball growth" -> "is a twisted geometric walk on the cyclic unit group `<2>
≅ Z/(2·3^(n-1))` — a weighted circulant with generator `2^(-1)`, all steps `a >= 1` allowed (the parity
restriction `a ≡ ε(z) mod 2` lives on the space side, in the parents' recursion of S19 §5); the ball growth".

**C14 (lines 196–199).** "is the analogue of Merca's divisor recursion. Verdict: the frame is right; a closed
formula for `H_n(1)` is not in reach (it would be a formula for the tree of `1`)." -> "is the analogue of
Merca's linear recurrence (his Corollary 1.5). Verdict: remark; `mu_n` is positive by definition (a sum of the
weights `2^(-A(w))`), so Merca's positivity phenomenon has no counterpart, and the open `liminf H_n(1) > 0` is
an asymptotic question about a positive sequence; a closed formula for `H_n(1)` is not in reach (it would be a
formula for the tree of `1`)."

**C15 (lines 19–21).** "Mazur 2026 (Lemma 8.1: "primitive Fourier coefficient bound `C* h^(-6409)` at level
`h`")" -> "Mazur 2026 (Lemma 8.1, and the paragraph following it in §8.1: "primitive Fourier coefficient bound
`C* h^(-6409)` at level `h` for characters not factoring through level `h - 1`")".

**C16 (lines 27–28, 63, 268, 528, 530).** "FINITE-EXACT (the Fourier profile of the law to level 17 in float64;
...)" -> "VERIFIED (the Fourier profile of the law to level 18 in float64, agreeing with the rational law to
level 6 and with an independent forward recursion to level 14 to `4·10^(-7)`; ...)"; "Computed to `h = 18`
(FINITE-EXACT, float64)" -> "(VERIFIED, float64)"; line 268 "FINITE-EXACT, float64 FFT" -> "VERIFIED, float64
FFT"; line 528 "checked at levels 4–17" -> "4–18"; line 530 "`M(h)`, `h <= 17` ... FINITE-EXACT (float64 FFT)"
-> "`M(h)`, `h <= 18` ... VERIFIED (float64 FFT)". FINITE-EXACT stays for the collision table, the cycles, the
Gauss-sum table, Proposition 6's checks.

**C17 (lines 398–399).** "at `t = 1, 8, 8, 16, 32, 64, 128, 256, 18659, 2048, 4096, 8192`" -> "at `t = ±2^(j+1)
mod 3^j` (`1, 8, 8, 16, 32, 64, 128, 256, 18659 = -2^10, 2048, 4096, 8192`)".

**C18 (line 157).** "(S19 §4)" -> "(S19 §2, §7)" (minor).

---

## (iii) Verdict

**SOUND WITH CORRECTIONS.**

Propositions 1, 2, 4 and 5 are correct as stated and their proofs hold, with two repairs that do not touch the
statements: Proposition 4's "(d, C_w) determines w" must read "(d, C_w, A)" (item 36), and the 2-adic source
class is mod `2^A`, not `2^(A+1)` (item 46). Proposition 3(ii) is true, but the inequality written as its proof
is false (178 violations, `|G_9(2)| = 0.207` against `0.499`); the tail costs `2^(1-s)` (item 32). Proposition
6, added during the audit, is right (item 69). Every number recomputed reproduces: the law (exact to level 6,
float to 14), the profile `M(h)` and its argmax up to sign, the masses, the Parseval identity (sharpened to an
exact level-by-level identity), the Gauss-sum table, the collision table as defined, the cost decomposition, the
cycles and the `3x + 139` orbit, the traces. The (2.3) and Lemma 8.1 quotations and the S19 cross-references
are faithful; the five paper summaries are accurate against their texts.

Substantive corrections, in order of importance:

1. **The decay reading is not supported by the data (items 23–25; C2–C4).** A shifted power law
   `C (h + 2.5)^(-2.5)` fits `h = 7..18` better (`1.8%`) than the geometric law (`4.3%`), reproduces the
   "steepening" local exponents and the ratios `0.85–0.89`, predicted the new level 18 (`0.0110` against
   `0.0105` geometric; measured `0.0112`), and would contradict Proposition 2's conclusion; the readings separate
   only near `h = 31`. "Levels off"/"flat" is false (the ratio rises `0.78 → 0.82 → 0.87` over successive
   windows and `0.894` at `h = 18` is the largest yet); "would satisfy (2.3) with room" and "if the ratio stays
   below `1`, the decay is geometric" must go. Lines 299–300 already say the honest thing and should govern
   §0.1, §2 and §6.
2. **The argmax law is misdescribed (item 6; C1).** `s - h = 0, 1, 2, 3, 4, 5` in steps of two or three levels,
   not "`h + 3` to `h + 5`" (true only for `h >= 13`); `2^s/3^h` decays geometrically instead of sitting near
   `0.03–0.05`. The wiring propagates "`2^(h+3)`".
3. **The collision table measures its valuation cap (items 39–41; C8, C9).** With valuations `<= 70` the first
   colliding depth drops to `2` (`n <= 9`) or `3` (`n <= 12`) and the minimal colliding sum decreases with
   depth; the looseness factor is at least `4.4` (`4.4–7` on the tabulated rows), not `3–5`; "depth `≈ 3.8 n`"
   is unsupported (full coverage of the unit classes by `d ≈ n/2` for `n <= 8`).
4. **Proposition 3(ii)'s inequality (item 32; C5)** and **Proposition 4's recovery step (item 36; C6)** are
   wrong as written; both proofs are repaired by one line.
5. **The 2-adic source class (item 46; C7):** mod `2^A`.
6. **The resonance "mechanism" (item 27; C12)** is a heuristic stated as fact, and its "coherence loss `0.96`"
   is a quotient, not a measurement.
7. **Mirrors (items 54–57; C13, C14):** the DFR "sector restriction" is on the wrong side of the Fourier
   transform; the Merca "frame is right" is an overreach (`mu_n` is positive by definition) and "divisor
   recursion" is not in Merca.
8. **Slips and labels (items 12, 16, 18, 26, 61; C10, C11, C15–C18):** `3·10^(-4)` is `7·10^(-5)`; the mass per
   level is not a constant; "undershoots" is inverted; the Lemma 8.1 attribution; float64 is VERIFIED, not
   FINITE-EXACT; header and status table still say "level 17".

Not checked (outside reach): the literal text of Tao's Proposition 1.14; levels 15–18 of the profile and the
`h = 14` decomposition (the session's FFTs, consistent with all fits); the completeness of the `3x + k` cycle
survey and its `k` other than 1, 13, 37; "the sum of squares is not invariant"; the barrier atlas rows and
THM-4514/THM-4476 as cited for the Lyu mirror; the wording "order-four free affine action" for Viaclovsky's
logarithmic transform.
