# Audit of `collatz_three_mirrors_20260929.md` (S23), sections 0–5 and 7

**Auditor:** independent session, 2026-09-29. Own code only:
`04-computation/experiments/collatz_three_mirrors_20260929_audit.py` (Python with three embedded C programs
compiled at run time into a temp directory outside the repo) → `05-knowledge/results/collatz_three_mirrors_20260929_audit.out`
(236 s on the 64 GB machine; `python3 04-computation/experiments/collatz_three_mirrors_20260929_audit.py 16 19 30`).
The repo's scripts were read but not run or imported. Method, briefly: the law `mu_n` on `Z/3^n` was recomputed for
`n <= 16` by a forward DP that embeds level `n-1` into level `n` (the map `y -> 3y+1` from `Z/3^(n-1)` into `Z/3^n` is
injective, so the valuation mixture is a gather with no duplicate accumulation; truncation `a <= 60`; the projection
`mu_n -> mu_(n-1)` was checked to `2e-16` at every level and `H_2(1) = 8/7`, `H_2(-1) = 22/7` exactly). From the law,
`S_n(psi_j) = E[psi_j(Y_n)]` was computed directly in discrete-log coordinates and `m_n(k) = mu_hat_n(2^k)` by an FFT
over `Z/3^n`; Proposition 2 was implemented three ways (sequential Python filter, an exact circular-convolution closed
form, a sequential C program from level 0 streaming level 19 with a 100-step warm-up); the census was rerun in Python
integers (box A) and in C with residues mod `3^30` (box B); the basins by a bitmask sieve in C on `[1, 2^30]` that
records every target on an orbit (so nestedness is tested, not assumed), cross-checked by a pure-Python sieve on
`[1, 2^22]`. Web check: the Krasikov–Lagarias theorem statement was read from arXiv math/0205002 (Theorem 6.1).

**Overall verdict: SOUND WITH CORRECTIONS.** Every computed number in sections 2–4 that I recomputed is reproduced
(the identities to `1e-9`, the maxima to all seven printed digits at `n = 1..18`, the census histograms and deepest
lines line by line, the basin densities to all six printed decimals at every dyadic range, the level-19 maximum). The
two theorems are correct as stated. The corrections are textual: one wrong constant in a displayed formula
(section 2, "partial sums"), one formula off by a factor 2 (`psi_(±2)(Y_n)`), one misreading of a normalisation
(the "imprimitive characters" explanation of the Fourier mass), one general rule that is off by one at three of the
four levels it is meant to summarise (`s = floor(h log_2 3) - 6`), one wrong historical sentence (the
Krasikov–Lagarias threshold "astronomically beyond" the checked range), one overstated precision ("five digits" for the
parity identity), "geometric to 1% through depth 12" (true through depth 10), and two forward references to an audit
(section 6) that is still a placeholder. No claim typed PROVED fails; the NUMEROLOGY items are typed as such; the
DIRECTION readings are heuristics, one of which (section 0, item 3) is presented without its type.

---

## 1. Numbered claims and verdicts

### Theorem 1 and the spectrum (section 2)

1. **Theorem 1(i)**, `sum_(k mod L_n) m_n(k) conj(psi_j(2^k)) = tau(conj psi_j) S_n(psi_j)` for primitive `psi_j`. **HOLDS.**
   Own derivation: substituting `u = 2^k` turns the left side into `sum_y mu_n(y) sum_u conj psi_j(u) e(uy/3^n)`, and for a
   unit `y` the substitution `u = v y^(-1)` gives `psi_j(y) tau(conj psi_j)` — this holds for *every* character of
   `(Z/3^n)^×`, primitive or not (the note's restriction to primitive characters is harmless; it matters only for
   `|tau| = 3^(n/2)`). Numerically (own `S` from the law in discrete-log coordinates, own `F = FFT(m)`, own
   `tau = FFT(omega)`): `max_j |F - tau S| / 3^(n/2) <= 2.5e-16` at every level `n = 1..16`, and `|tau|/3^(n/2) = 1`
   to twelve digits for every primitive `j`. The claim "for imprimitive characters the Gauss sum vanishes" (script
   comment) is correct for `n >= 2` (`max|tau| <= 1.1e-9`, i.e. FFT rounding of size `3^n · 1e-16`; `max|F| <= 4e-16`),
   and at `n = 1` the trivial character has `tau = -1`, `F = -1`, `S = 1`, so (i) holds there too. `S_n(psi_j)` itself
   does not vanish at imprimitive `j`: it equals the lower-level moment, `max_j |S_n(psi_(3j')) - S_(n-1)(psi_j')| <= 3.3e-16`
   at every level (this is the "primitive iff `3 ∤ j`" statement, verified). The note handles the imprimitive
   characters correctly by summing over the primitive ones only in (ii).
2. **Theorem 1(ii)**, `sum_(psi prim) S_n(psi) = L_n mu_n(1) - (L_n/3) mu_(n-1)(1) = rho_n(1) - rho_(n-1)(1)` and the
   parity version, `n >= 2`; the `n = 1` case `-1/3`, `+1/3`; the constant `L_n = (2/3) 3^n`. **HOLDS.** Own proof:
   `sum_prim psi(y) = L_n 1[y ≡ 1 (3^n)] - L_(n-1) 1[y ≡ 1 (3^(n-1))]` with `L_(n-1) = L_n/3` for `n >= 2` (the
   characters induced from level `n-1` number `L_(n-1)`); consistency gives `mu_n(y ≡ 1 mod 3^(n-1)) = mu_(n-1)(1)`.
   At `n = 1` the induced characters number `1`, not `L_1/3 = 2/3`, which is why the general formula gives `0` there
   while the true sum is `-1/3`; the note states the `n = 1` case separately and its formula `rho_n(±1) = 1 ∓ 1/3 +
   sum_(m=2)^n ...` is right (`rho_1(1) = 2/3`, `rho_1(-1) = 4/3`). Numerically both sides agree to the nine printed
   decimals at every level `2..16` (e.g. `n = 16`: `-0.015796933` and `+213.433932660`, both ways); the accumulated
   sums give `rho_n(1) = 0.761905, 0.661379, 0.618418, 0.642868, 0.636947, 0.573573, 0.516384, 0.464503, 0.424971,
   0.394279, 0.358243, 0.333431, 0.315492, 0.305906, 0.290110` for `n = 2..16`, each equal to S19's independent
   `mazur_harmonic_mass_deep18` value at all five printed digits, and `rho_16(-1) = 640.2286` against S19's `640.2`.
   The corollary `E[chi_(-3)(Y_n)] = -1/3` holds (`S_1(chi_(-3)) = -1/3` exactly, `mu_1(2) = 2/3`).
3. **"Both identities reproduce S19's independently computed values to five digits at every level 10..16"** (section 0
   and the status line). **HOLDS WITH CORRECTION.** True for `rho_n(1)`. For the parity identity S19 prints `rho_n(-1)`
   to four significant digits (`56.12, 84.23, 126.4, 189.6, 284.5, 426.8, 640.2`), so its differences are known to
   3–4 digits only (`28.11, 42.2, 63.2, 94.9, 142.3, 213.4`), which is what the note's own section 2 lists. Also the
   listed `rho_12(1) = 0.35825` is a rounding slip: the accumulated sum is `0.358243` (S19: `0.35824`).
4. **The rms `0.8452, 0.8321, ..., 0.8408` (`n = 2..16`), `rms^2 = 0.70 =` the typical `|mu_hat|^2 3^n`.** **HOLDS.**
   All fifteen values reproduced to four decimals. The identity behind it: `sum_prim |S|^2 3^n = sum_j |F_j|^2 = L_n
   sum_k |m_n(k)|^2` (Parseval on the cycle, `F` vanishing at imprimitive `j`), so `rms^2 = (3/2) sum_(units) |mu_hat_n|^2
   = 3^n/L_n ×` S20's "mass at level `n`" (`0.4714` at `n = 16`; `1.5 × 0.4714 = 0.7070` ✓).
5. **"Fourier mass per level 0.709 at n = 19 (S20: 0.462 → 0.472 was the mass on the primitive characters; the
   all-units figure here includes the imprimitive ones and converges to 0.71)."** **HOLDS WITH CORRECTION.** The
   number is right (`3^19 sum_k |m_19(k)|^2 / L_19 = 0.7089`; `sum_(units) |mu_hat_19|^2 = 0.4726`), the explanation is
   wrong: both figures are sums over the *same* set, the units of `Z/3^n` (the family `m_n(k)` covers the units and
   nothing else; the imprimitive additive frequencies `3 | t` are not on the cycle at all). The factor `1.5 = 3^n/L_n`
   is the normalisation of the fullperiod script's "mass" (a mean over the `L_n` units times `3^n`), not extra
   characters. `0.7084/0.4722 = 1.500` at `n = 18`.
6. **`P(z > 1, 2, 3) = 0.27, 0.13, 0.074` against `0.37, 0.14, 0.05`; largest `|S| 3^(n/2) = 1.00, 1.41, ..., 9.11`
   (`n = 2..16`), "roughly `n/2`".** **HOLDS.** Reproduced: `0.266, 0.127, 0.074` at `n = 16`; maxima `1.000, 1.411,
   1.586, 2.022, 2.206, 2.823, 2.945, 3.604, 4.085, 4.347, 5.751, 6.035, 7.409, 7.786, 9.111`. In rms units the maximum
   is `10.8` at `n = 16` (`max/rms = 1.18, 1.70, 1.90, 2.42, 2.64, 3.38, 3.52, 4.31, 4.88, 5.19, 6.86, 7.19, 8.82, 9.27, 10.84`).
7. **"The maximum of `L_n` independent Rayleigh variables would be `0.84 sqrt(ln L_n) = 3.5` at `n = 16`."** **HOLDS**
   (as an order of magnitude; the count is loose). With `|S|^2/mean ~ Exp(1)` the expected maximum of `N` samples is
   `rms sqrt(ln N + gamma)`: `0.84 sqrt(ln L_16) = 3.48` with `N = L_n`, `3.44` with the correct count `N = (2/3) L_n`
   of primitive characters, `3.43` with `N = (1/3) L_n` independent conjugate pairs (`|S_j| = |S_(-j)|` since `mu_n` is
   real). The observed `9.11` is 2.6 times any of these; the conclusion (non-Rayleigh) stands.
8. **"The top eight at `n = 16` are all `8.9–9.1`, a cluster, not an outlier."** **HOLDS WITH CORRECTION.** The top eight
   entries are four conjugate pairs (`9.11` at `±1516613`, `8.90` at `±1398515`, `8.86` at `±505708`, `8.86` at
   `±335633`); the statement is about the top four moduli, each carried by a pair `psi_(±j)`.
9. **Which characters carry the maximum: `psi_(±2)` at `n = 3, 4, 6`, `psi_(±8)` at `n = 5`, `psi_(±4)` next at `n = 4`,
   `psi_(±80)` at `7`, `psi_(±278) = psi_(∓2^12)` at `8` (`278 = L_8 - 4096 = 2(3^7 - 2^11)`), then `±521, ±1193, ±1127,
   ±4219, ±40162, ±71333, ±154112 = ±2^9·301, ±1516613`.** **HOLDS.** All reproduced from the law (my argmax list:
   `2, 2, 8, 2, 80, 278, 521, 1193, 1127, 4219, 40162, 71333, 154112, 1516613` for `n = 3..16`; second at `n = 4` is
   `±4` with `1.42`). `4374 - 4096 = 278 = 2·139` ✓.
10. **The formula `psi_(±2)(Y_n) = e(∓a_n/(2·3^(n-1))) e(±log_4(1 + 3z)/3^(n-1))`.** **HOLDS WITH CORRECTION** (factor
    2). `psi_2(y) = e(2 log_2 y / L_n) = e(log_2 y / 3^(n-1))`, and `log_2(Y_n) = -a_n + 2 log_4(1 + 3z)`, so
    `psi_(±2)(Y_n) = e(∓a_n/3^(n-1)) e(±2 log_4(1 + 3z)/3^(n-1))`. The displayed formula (with `2·3^(n-1)` and a single
    `log_4`) is `psi_(±1)(Y_n)`. The expression for `z = Y_(n-1)` as a sum over suffix costs is right.
11. **Parity split: `mean_odd S = -mean_even S` to three digits (`∓1.116·10^(-5)` at `n = 16`), explained by
    Theorem 1(ii) with `sum_even - sum_odd = rho_n(-1) - rho_(n-1)(-1) ≈ 0.325 (3/2)^n` and `sum_even + sum_odd =
    O(0.05)`; the per-character margin `≈ 0.5 (1/2)^n`, below the Parseval scale by `0.866^n`.** **HOLDS.** Reproduced:
    `-1.1157e-5 / +1.1155e-5` at `n = 16`, sum of the means `-1.7e-9`; `#odd = #even = L_n/3`; `0.325 (3/2)^n / L_n =
    0.4875 (1/2)^n`; `(1/2)/3^(-1/2) = 0.866`.
12. **"Theorem C's question `liminf H_n(1) > 0` is the question whether the partial sums `1 + sum_(m<=n) T_m`, `T_m =
    (3/2) sum_(prim mod 3^m) S_m(psi)`, stay away from zero; the `T_m` are `-1/2, +0.143, -0.151, ...`."** **FAILS as
    written (constant), the list of `T_m` HOLDS.** With `T_1 = -1/2` the partial sum `1 + T_1 + T_2 = 0.643`, but
    `H_2(1) = 8/7 = 1.143`. Since the increment identity fails at `m = 1` (item 2), the correct statement is
    `H_n(1) = 1 + sum_(2<=m<=n) T_m` (equivalently `3/2 + sum_(1<=m<=n) T_m`). The sixteen `T_m` values are reproduced
    to the printed digits.
13. **"each a sum of `L_m` moments of size `0.84·3^(-m/2)` whose random-phase size would be `0.69`."** **HOLDS WITH
    CORRECTION.** The number of primitive characters is `(2/3) L_m`, giving a random-phase size `sqrt((2/3) L_m) · 0.84 ·
    3^(-m/2) = 0.56` (the note's `0.69 = sqrt(L_m) 0.84 3^(-m/2)` uses all `L_m` characters). Observed `|sum| = 0.0158`
    at `m = 16`; the point (near-complete cancellation) is unchanged.

### Proposition 2 and the level-19 extension (section 2)

14. **Proposition 2: the recursion `m_n(k) = sum_a 2^(-a) omega_n(k-a) m_(n-1)(k-a)` on `Z/L_n` is the one-pole filter
    `m_n(k) = (g(k-1) + m_n(k-1))/2`, its periodic solution is unique, and the transient is `2^(-steps)`.** **HOLDS.**
    Own derivation: conditioning on the last valuation gives `mu_hat_n(t) = sum_a 2^(-a) e((t 2^(-a) mod 3^n)/3^n)
    mu_hat_(n-1)(t 2^(-a) mod 3^(n-1))` (the inverse of 2 mod `3^n` reduces to the inverse mod `3^(n-1)`), and
    `m_n(k-1) = 2[m_n(k) - g(k-1)/2]`; the homogeneous solutions `c 2^k` are not `L_n`-periodic, so the periodic
    solution is unique, and from any state the error halves per step. Numerically, one recursion step from the law's
    `m_(n-1)` reproduces the law's `m_n` to `max|diff| <= 2.5e-16` at every level `1..16` (closed form) and `1..10`
    (sequential filter); the sequential C program run from level 0 agrees elementwise with the law-derived `m_12` and
    `m_16` to `8.7e-17` and `7.5e-17` (no valuation truncation, no accumulated error).
15. **The maxima agree with the S20 FFT at every level `1..18` to seven digits (`0.5773503, 0.3779236, ..., 0.0144095
    (16), 0.0125107 (17), 0.0111873 (18)`).** **HOLDS.** All eighteen values reproduced to all seven digits by the C
    recursion (difference `0.0` at the printed precision) and to `5e-8` by the law's FFT (`n <= 16`); the argmax is
    `±2^s` with S20's `s` at every level (`2^s` is one of `2^k`, `2^(k + L/2)` for the argmax `k`: True at all 18).
16. **`M(19) = 0.00982` at `k = 24` (`= 19 log_2 3 - 6.1`), mirror `-2^24` second, `k = 25, 23` next; "the maximum
    over all `774,840,978` units is on `±2^s` with `s = floor(19 log_2 3) - 6`"; Fourier mass `0.709`.** **HOLDS**
    (independently recomputed: the shipped `fullperiod_max` output ends at `n = 18` and contains no level-19 line, so
    at audit time this claim rested on no output in the repo). Own C recursion, level 19 streamed: `M(19) = 0.0098157`
    at `k = 24` and its mirror `387420513 = 24 + L_19/2`; then `±2^25` (`0.0097694`), then `2^23`'s mirror (`0.0087743`);
    `24 - 19 log_2 3 = -6.11`; mass `0.7089`; `L_19 = 774,840,978` ✓.
17. **The general rule in section 0: "the maximum over all units is still at `±2^s`, `s = floor(h log_2 3) - 6`."**
    **HOLDS WITH CORRECTION.** `floor(h log_2 3) - 6 = 19, 20, 22` at `h = 16, 17, 18`, but the maxima sit at `2^20,
    2^21, 2^23`; the formula holds at `h = 19` only (offsets `-5.36, -5.94, -5.53, -6.11` at `h = 16..19`). S20's rule
    `s = h log_2 3 - 6 ± 1` is the correct general statement; the note's own Proposition 2 paragraph uses it correctly.

### The census and the critical-rate identity (section 3)

18. **The null model: exactly one sign works, `P(depth >= 1) = 1`, `P(depth >= d) = 3^(-(d-1))`, Poisson counts with
    mean `N 3^(-(d-1))`.** **HOLDS, with a remark that changes the reading of the histogram.** `u 2^Q` is a unit mod 3,
    so exactly one of `u 2^Q ∓ 1` is divisible by 3 (never both): depth `>= 1` is a tautology, and the conditional law
    is `3^(-(d-1))` under uniformity on the class. But for fixed `u` the map `Q -> u 2^Q mod 3^d` is a bijection onto
    the units on every window of `L_d = 2·3^(d-1)` consecutive `Q`, so whenever `W >= L_d` each multiplier contributes
    `2 floor(W/L_d) + O(1)` lines of depth `>= d` *deterministically*: box B's ratios `1.000` at `d <= 8` (`L_8 = 4374 <
    10^4`) are equidistribution of the powers of 2, not evidence for a Poisson model (own measurement: the per-`u`
    count of depth-`>= 6` lines has mean `41.15` and variance `0.13`, where Poisson would give `41`; at `d = 9`,
    `L_9 = 13122 > W`, mean `1.52`, variance `0.25`, still sub-Poisson since each `u` contributes at most one hit per
    target residue). The Poisson comparison is informative for `d >= 10` in box B and `d >= 8` in box A. Likewise
    "lines of depth `>= 5` with `Q <= 60` exist at every `Q`" is forced: every `Q` has exactly `4` or `5` such
    multipliers `u <= 1000` (two per block of `486` consecutive odd `u` prime to 3).
19. **Box B numbers: `N = 16,670,000`; ratios `1.000, ..., 0.951, 0.956, 1.100, 1.291, 1.291, 1.291, 3.87` (`d = 1..17`);
    cumulative `276/282.3 (0.65)`, `97/94.1 (0.40)`, `37/31.4 (0.18)`, `14/10.5 (0.17)`, `5/3.5 (0.27)`, `2/1.16
    (0.32)`, `1/0.39 (0.32)`, none deeper; deepest `1187·2^5031 ≡ 1 mod 3^17`, `2441·2^8384 ≡ -1 mod 3^16`, depth 15 at
    `(3025, 846, -)`, `(1685, 8663, -)`, `(55, 423, +)`, the nine depth-14 lines.** **HOLDS.** Own C census with
    residues mod `3^30` (cap 30 instead of 20): the histogram is identical to the repo's line by line (`11113333,
    3704445, 1234814, 411603, 137201, 45735, 15230, 5102, 1690, 571, 179, 60, 23, 9, 3, 1, 1`, nothing at `d >= 18`),
    the Poisson tails are `0.654, 0.396, 0.178, 0.171, 0.272, 0.324, 0.321` for `d >= 11..17`, and the deepest lines
    are the same pairs with the same signs; all depth-14/15/16/17 lines were also verified exactly in Python integers
    (`v_3(1187·2^5031 - 1) = 17`, `v_3(2441·2^8384 + 1) = 16`, etc.). The expected maximum depth `log_3 N + O(1) =
    15.1` against the observed `17` is consistent.
20. **"the depth histogram is geometric to 1% through depth 12".** **HOLDS WITH CORRECTION.** Through depth 10 (ratios
    within `1.1%`); at `d = 11, 12` the ratios are `0.951, 0.956` (`5%` off, within the Poisson fluctuations `±7%`,
    `±13%`), and `d <= 8` is deterministic (item 18).
21. **Box A: `3` lines of depth `>= 14` against `0.42` (`P ≈ 0.009`), `4` of depth `>= 13` against `1.25` (`0.04`);
    the wide box dissolves the excess.** **HOLDS.** Own Python census: identical histogram (`444000, 148000, 49334,
    16444, 5485, 1815, 614, 204, 74, 22, 2, 2, 1, 2, 1`), `P = 0.0089` and `0.0386`; the deep lines `(55, 423, +)`,
    `(917, 1557, -)`, `(521, 1365, -)`, `(907, 455, +)`, `(301, 124, -)`, `(41, 271, -)`, `(497, 41, -)`, `(413, 1806, +)`
    and the window list (`(497, 41)` depth 11; `(205, 5)`, `(133, 17)`, `(647, 27)`, `(901, 30)`, `(997, 35)`, `(379, 45)`,
    `(491, 58)` depth 8) reproduced.
22. **Exact checks `v_3(55·2^423 + 1) = 15`, `v_3(13·2^154 - 1) = 7`, `v_3(2^486 - 1) = 6`, `v_3(2^480 - 1) = 2`.**
    **HOLDS.** Verified in Python integers; the last two also by LTE (`v_3(2^k - 1) = 1 + v_3(k/2)` for even `k`:
    `1 + v_3(243) = 6`, `1 + v_3(240) = 2`).
23. **The identity `e^(-I) 3^(θ*/ln 2) = log_2 3 - 1`.** **HOLDS** (numerically to `1e-16`; algebraically `3^(θ*/ln 2)
    = e^(θ* m)` and `e^(-I) = e^(-θ* m)(m-1)` by the definition `I := θ* m - ln(m-1)`). It is the definition of `I`
    rewritten and holds for any `θ*`, so it carries no information about the law; the typing "PROVED (a one-line
    identity)" is accurate but the word "exactly" in section 0 suggests a coincidence where there is none. The reading
    (the Chernoff amplitude law `M(n) u^(-0.438)` extrapolated to `u ≈ 3^n` gives `(m-1)^n`, above the Parseval scale
    `3^(-n/2)`, so the exponent must steepen; H's `1.3%` margin `= 1 - 3^(-1/2)/(m-1)`) is a heuristic: it uses
    `M(n) ≈ e^(-nI)`, which S20/S21 type CONJECTURAL, and the families `u 2^j` overlap heavily. Section 3 types it
    DIRECTION; section 0 item 3 states it ("so Parseval forces ...") without a type — see the corrections.

### The saddle chain (section 4.1)

24. **`T^(2j)(1 + 4^j t) = 1 + 3^j t`; the chain `1 + 3^j 4^(6-j)` = `4097, 3073, 2305, 1729, 1297, 973, 730`;
    `1729 = 1 + 12^3 = 1 + 3^3 4^3`; `1 + 12^m` the midpoint of the chain of `1 + 4^(2m)` (`13, 145, 1729, 20737, 248833,
    2985985`); `17 = 1 + 4^2` heads `17, 13, 10`.** **HOLDS.** `x ≡ 1 mod 4` gives `T(x) = 2 + 3·2^(2j-1) t` (even for
    `j >= 1`) and `T^2(x) = 1 + 3·4^(j-1) t`; checked for `j <= 7`, `t < 200`, the chain step for `m <= 11`, the
    `T^2`-iterates of `4097`, and `T^(2m)(1 + 4^(2m)) = 1 + 12^m` for `m <= 8` by iteration. (That `1 + 12^m` is the
    `j = m` point of the chain of `1 + 4^(2m)` is the identity `1 + 3^m 4^m = 1 + 12^m`, i.e. true by definition once
    Proposition 4 is in hand.)
25. **Factorisations and the numerology: `4097 = 17·241`, `3073 = 7·439`, `2305 = 5·461`, `1729 = 7·13·19` (all
    `≡ 1 mod 6`, `= 13·Φ_6(12)`), `1297 = 6^4 + 1` prime, `973 = 7·139`, `139 = 3^7 - 2^11`, `4·3^5 + 1 = 7(3^7 - 2^11)`,
    `730 = 3^6 + 1`, `1729 = 3^6 + 10^3`.** **HOLDS** (all verified with sympy; `Φ_6(12) = 133 = 7·19`). Typed
    NUMEROLOGY where the note says so; correct.
26. **The `-17` cycle `-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34`, seven odd steps, clock `139`.**
    **HOLDS.** The `T`-orbit returns to `-17` after 11 steps (7 odd, 4 even); `2^11 - 3^7 = -139`; the cycle constant
    is `c = -17 (2^11 - 3^7) = 2363 = 17·139`.
27. **The orbits of `1729` and `27` merge at `137` (positions 10 and 12): `1729 -> 2594 -> 1297 -> 1946 -> 973 -> 1460 ->
    730 -> 365 -> 548 -> 274 -> 137`.** **HOLDS** (orbit of `1729`: 68 `T`-steps, max `4616`; of `27`: 70 steps; first
    common element `137` at positions 10 and 12; `137 = 1 + 8·17`).
28. **"The chain is the `j = 6` row of Theorem 2.3 of the extended-Collatz note (`730, 973, 1297, 1729, 2305, 3073, 1024`
    ...), read forward."** **HOLDS** (faithful quote of that table; its last entry `1024 = (3073 - 1)/3` is the
    backward-greedy `C`-preimage, the "growth block", not the chain point `4097 = 1 + 4·1024`).

### The basins and Krasikov–Lagarias (sections 4.2, 4.3)

29. **The density table on `[1, 2^30]`: `0.003339, 0.003524, 0.004163, 0.004395, 0.004881, 0.005795, 0.014459, 0 (26
    numbers), 0.299245`; stability from `2^24` on (`0.004388, ..., 0.004395`); `|B(1729) ∩ [1, 2^30]| = 4,719,191`; the
    two codes agree at `2^21`.** **HOLDS.** Own bitmask sieve in C (records every target on the orbit): every density
    at every dyadic range `2^20..2^30` equals the repo's to all six decimals, `|B(1729)| = 4719191`, `B(27) = {27·2^k :
    k <= 25}` (26 numbers, listed), and the pure-Python sieve on `[1, 2^22]` gives the same nine densities as the C code
    and the repo at `2^22`.
30. **Nestedness `B(4097) ⊂ ... ⊂ B(730)`; the increments are the side entries `(2^k a - 1)/3`, `k >= 4`: `0.00023` at
    `1729` (through `9221, 36885, ...`), `0.0087` at `730` (through `3893 = 1 + 4·973`).** **HOLDS.** Tested rather than
    assumed: zero violations of `B(chain_j) ⊂ B(chain_(j+1))`, of `B(chain) ⊂ B(137)` and of `B(27) ⊂ B(137)` on
    `[1, 2^30]`; `|B(1729) \ B(2305)| = 249739` (density `0.000233`), `|B(730) \ B(973)| = 9303565` (`0.008665`);
    `(2^4·1729 - 1)/3 = 9221`, `(2^6·1729 - 1)/3 = 36885`, `(2^4·730 - 1)/3 = 3893 = 1 + 4·973`.
31. **"`B(27)` is the doubling ray of `27` (`27 ≡ 0 mod 3` has no odd preimage): the famous starting value has an
    empty tree above it."** **HOLDS** (`27·2^k` is never `(3m+1)/2`; the sieve finds exactly the 26 numbers `27·2^k
    <= 2^30`).
32. **`dens B(137) = 0.299245`; the path "`137 -> 103 -> 155 -> 233 -> 350 -> ... -> 577 -> ... -> 5`"; `B(5) = 0.938`
    (mac-mini's `e_2`).** **HOLDS WITH CORRECTION.** Density reproduced. The path skips `206 = T(137)`: `137 -> 206 ->
    103 -> 155 -> 233 -> 350 -> 175 -> ...`. `e_2 = 0.93796` is what the necklace note reports.
33. **Krasikov–Lagarias as quoted ("for every `a ≢ 0 mod 3` there is `X_0(a)` with `|B(a) ∩ [1, X]| >= X^0.84` for `X >=
    X_0(a)`"), and `B(a)` as the right set.** **HOLDS.** The theorem (arXiv math/0205002, Theorem 6.1; Acta Arith. 109
    (2003)) reads: for each positive `a ≢ 0 (mod 3)`, `π_a(x) := |{1 <= n <= x : some T^(j)(n) = a}|` satisfies
    `π_a(x) >= x^0.84` for all sufficiently large `x >= x_0(a)`, with `T(n) = n/2`, `(3n+1)/2` — the same map as the
    note's sieve, and `B(a) ∩ [1, X]` is exactly the set counted (`a` itself included). For the odd roots the count is
    the same under the Collatz map `C`; for the even chain point `730` the `C`-basin additionally contains the doubling
    ray of `243 = 3^5` (`C(243) = 730`, skipped by `T`), a density-zero difference (23 numbers below `2^30`).
34. **`(2^30)^0.84 = 3.85·10^7`; the count is `12.2%` of the bound, rising by `2^0.16 = 1.117` per doubling (`0.046` at
    `2^21`, `0.098` at `2^28`, `0.122` at `2^30`); `X_0(1729) > 2^30` (PROVED by the count).** **HOLDS.** `X^0.84 =
    3.8544·10^7`, ratio `0.1224`; the ratios `0.0458, ..., 0.0980, 0.1095, 0.1224` at `2^21..2^30` reproduced. Since
    the inequality fails at `X = 2^30`, every threshold `x_0(1729)` for which it holds for *all* `x >= x_0` exceeds
    `2^30`: correct, and a direct consequence of the theorem's asymptotic form (the typing "PROVED (by the count)" is
    accurate; it is a finite-exact count plus the theorem's quantifier).
35. **"the inequality first holds near `X = 2^30 · 1.117^(-log(0.122)/log(1.117)) ≈ 2^49.3 ≈ 7·10^14` (OBSERVED
    extrapolation)".** **HOLDS WITH CORRECTION** (arithmetic). `-ln(0.1224)/ln(1.1173) = 18.94` doublings, i.e.
    `2^48.9 ≈ 5.4·10^14`; equivalently `X* = 0.004395^(-1/0.16) = 5.39·10^14 = 2^48.94`. "Near `2^49`" is fine; the
    decimals `2^49.3 ≈ 7·10^14` are not.
36. **"The theorem is asymptotic; for the root `1729` its threshold is astronomically beyond any range in which Collatz
    has been checked exhaustively at the time of the theorem (`2^68`, Barina, is beyond it — but the theorem's own
    constants are not explicit)."** **FAILS.** The extrapolated threshold `2^48.9 ≈ 5·10^14` lies *inside* the range
    verified exhaustively before the theorem (Oliveira e Silva 1999: all `n < 3·2^53 ≈ 2.7·10^16 = 2^54.6`) and far
    below Barina's `2^68 ≈ 3·10^20`; the sentence contradicts its own parenthesis. The valid point is only that the
    theorem's `x_0(a)` is ineffective and that, for this root, the bound is far from the truth (a density `0.0044`
    against a sublinear `X^0.84`).

### Attribution to the preprints, the parallel note, and statuses (sections 0, 1, 5, 7)

37. **The description of the three preprints** (Chocian, arXiv:2607.23177 / 2607.27503 / 2608.08724; `b_(chi,j) = f
    B_(1, chi ω^(-j)) mod p`; `mu_p = (1 + z zeta_p)/(1 + conj(z) zeta_p)`, `z = -zeta_3^2`; `sum_m P_m(X) Y^m = -log(1 -
    X(1 - e^(-Y)))`; `P_(p-j)(h) - P_(p-j)(1-h) = -(2h-1)(j-1)! b_j`; the twelve lines at `p = 67, 103, 139, 199, 241,
    271, 331, 337 (twice), 409, 421, 457`, seven at regular primes; the Gauss-sum identity `sum_t conj(chi)(t) P_m(h_t)
    = tau(conj chi) B_(m,chi)/(m·m!)`, `h_t = zeta_f^t/(zeta_f^t - 1)`; conductor 5, eleven lines; `27,508` zero lines
    over `55,121` pairs; `chi^2 = 0.05, 3.93`; digits uniform; conjugate quartic characters never vanishing at a common
    index; no association with classical irregularity; eight non-simple zeros, one of depth three at `(19, 37, 16)`).
    **HOLDS** against the text files `p1.txt`, `p2.txt`, `p3.txt` (Theorem 4.1's table, the abstracts). Two small
    omissions: paper 2 states its identity for *odd* `m`; the twelve-line catalogue is stated for `p < 500`, `p ≡ 1
    mod 6` (the note says "below 500"). Dates and identifiers match.
38. **Claims about the mac-mini necklace note**: THM-4515/4516/4517 and HYP-9165 exist there; the `(4,3,17)` statement
    is recorded UNVERIFIED there; that note names `(3,5,7)` as the smallest open Beal signature and lists `(3,4,n)`
    solved only for `n = 4, 5` from arXiv:2412.11933v2; `e_2 = 0.938`; the clock-shadow observation; the addendum on
    arXiv:2609.26996. **HOLDS** (all present in
    `collatz_necklace_20260929_fair_splits_power_clocks_basins.md`, sections 2.5, 3; the addendum is in its 2.5).
39. **"mac-mini's clock shadow argument ... is confirmed by the audit of section 6" (section 0, item 5) and "is
    confirmed by the audit below" (section 5); "adds an independent audit of its three theorems (section 6)" (header);
    "Audit of mac-mini's THM-4515–4517: section 6" (status line).** **FAILS as the text stands:** section 6 is
    `AUDIT_PLACEHOLDER`. These sentences assert a result that the note does not contain. (The clock-shadow statement
    itself was not audited here; it is out of scope.)
40. **Status typing (section 7 table and the status line).** **HOLDS WITH CORRECTIONS** listed above: the PROVED items
    are proved (Theorem 1, Proposition 2, the chain identities, the tautological rate identity, `X_0(1729) > 2^30`);
    VERIFIED/FINITE-EXACT items are reproduced; the level-19 line was VERIFIED by this audit rather than by a shipped
    output; the DIRECTION reading of section 3 is presented untyped in section 0 item 3; "five digits" applies to
    `rho_n(1)` only; and the claim "the primitive maxima over all units to level 18 against the S20 FFT to seven
    digits, extended to level 19" is now supported (item 16).

Not audited (out of scope or not decidable here): section 6 (placeholder), the history remark "the first
implementation ran two passes ... all maxima 12% low", the truth of the `(4,3,17)` statement (the note's UNVERIFIED
typing is appropriate), and the AUTHOR-CLAIMED arXiv:2609.26996. Remark on obligation (b): level 20 does not need
`37 GB` — streaming the last level as done here needs only `m_19` in memory (`12.4 GB` complex128, `6.2 GB`
complex64) plus residues generated by doubling mod `3^20` (no int64 overflow, unlike the block multiplication).

---

## 2. Textual corrections (exact phrase → replacement)

1. Section 0, item 1 and the status line: "Both identities reproduce S19's independently computed values to five
   digits at every level `10..16`" → "The seed-1 identity reproduces S19's independently computed `rho_n(1)` to five
   digits at every level `2..16`; the parity identity reproduces S19's differences to the 3–4 digits S19 prints".
2. Section 2, "The identities, checked": "`rho_n(1) = 0.42497, 0.39428, 0.35825, 0.33343, ...`" → "`... 0.35824 ...`".
3. Section 2, Proposition 2 paragraph: "Fourier mass per level `0.709` at `n = 19` (S20: `0.462 -> 0.472` was the mass
   on the *primitive* characters; the all-units figure here includes the imprimitive ones and converges to `0.71`)" →
   "Fourier mass per level `0.709` at `n = 19` (this is `3^n/L_n = 3/2` times S20's `0.462 -> 0.472`, the same sum
   `sum_(units) |mu_hat_n|^2` normalised per unit; it equals the rms² of the character spectrum)".
4. Section 0, item 2: "the maximum over all units is still at `±2^s`, `s = floor(h log_2 3) - 6`" → "the maximum over
   all units is still at `±2^s`, `s = h log_2 3 - 6 ± 1` (`s = 20, 21, 23, 24` at `h = 16..19`; `floor(h log_2 3) - 6`
   is exact at `h = 19` only)".
5. Section 2, spectrum paragraph: "the top eight at `n = 16` are all `8.9–9.1`, a cluster, not an outlier" → "the top
   four conjugate pairs `psi_(±j)` at `n = 16` are all `8.9–9.1`, a cluster, not an outlier".
6. Section 2, spectrum paragraph: "`psi_(±2)(Y_n) = e(∓a_n/(2·3^(n-1))) e(±log_4(1 + 3z)/3^(n-1))`" →
   "`psi_(±2)(Y_n) = e(∓a_n/3^(n-1)) e(±2 log_4(1 + 3z)/3^(n-1))`".
7. Section 2, "What this changes": "whether the partial sums `1 + sum_(m<=n) T_m`, `T_m = (3/2) sum_(prim mod 3^m)
   S_m(psi)`, stay away from zero" → "whether the partial sums `H_n(1) = 1 + sum_(2<=m<=n) T_m`, `T_m = (3/2)
   sum_(prim mod 3^m) S_m(psi)`, stay away from zero (the `m = 1` term `T_1 = -1/2` is already inside the `1 = H_1(1)`)".
8. Same paragraph: "each a sum of `L_m` moments of size `0.84·3^(-m/2)` whose random-phase size would be `0.69`" →
   "each a sum of `(2/3) L_m` moments of size `0.84·3^(-m/2)` whose random-phase size would be `0.56` (`0.84` for `T_m`)".
9. Section 0, item 3: "the depth histogram is geometric to `1%` through depth 12 and the tail is Poisson-consistent at
   every depth" → "the depth histogram is geometric to `1%` through depth 10 (forced by the equidistribution of the
   powers of two for `d <= 8`, where `L_d <= W`) and the tail is Poisson-consistent at every depth `>= 10`".
10. Section 0, item 3: "The critical rate of Theorem C is the Chernoff amplitude law of the multiplier families extended
    to `u ≈ 3^n`: `e^(-I) 3^(θ*/ln 2) = log_2 3 - 1` exactly, so Parseval forces the multiplier exponent to steepen
    beyond `-0.438` for large `u`, and H's margin ..." → "The identity `e^(-I) 3^(θ*/ln 2) = log_2 3 - 1` (the
    definition of `I` rewritten) reads the critical rate of Theorem C as the Chernoff amplitude law of the multiplier
    families extended to `u ≈ 3^n` under `M(n) ≈ e^(-nI)` (CONJECTURAL); on that reading Parseval would force the
    multiplier exponent to steepen beyond `-0.438` for large `u`, and H's margin ... (DIRECTION)".
11. Section 3, census paragraph: add after "Box B's depth histogram against `N (2/3) 3^(-(d-1))`: ratios `1.000, ...`":
    "(the ratios at `d <= 8` are forced: `L_d <= W`, so each multiplier contributes `2W/L_d + O(1)` lines
    deterministically; the Poisson test begins at `d >= 10`)".
12. Section 4.2: "the path `137 -> 103 -> 155 -> 233 -> 350 -> ...`" → "the path `137 -> 206 -> 103 -> 155 -> 233 ->
    350 -> ...`".
13. Section 4.3: "near `X = 2^30 · 1.117^(-log(0.122)/log(1.117)) ≈ 2^49.3 ≈ 7·10^14`" → "near `X = 2^30 ·
    1.117^(-log(0.122)/log(1.117)) = 2^(30 + 18.9) ≈ 2^48.9 ≈ 5·10^14` (equivalently `0.0044^(-1/0.16)`)"; and in
    section 0 item 4 and section 7 "`X_0(1729) ≈ 2^49`" may stay.
14. Section 4.3: "The theorem is asymptotic; for the root `1729` its threshold is astronomically beyond any range in
    which Collatz has been checked exhaustively at the time of the theorem (`2^68`, Barina, is beyond it — but the
    theorem's own constants are not explicit)." → "The theorem is asymptotic with an ineffective `x_0(a)`; for the
    root `1729` the extrapolated threshold `≈ 2^49` lies inside the range checked exhaustively before the theorem
    (`3·2^53`, Oliveira e Silva 1999, CITED) and far below Barina's `2^68`, so for this root the sublinear bound is
    weaker than the truth throughout the verified range."
15. Header ("adds an independent audit of its three theorems (section 6)"), status line ("Audit of mac-mini's
    THM-4515–4517: section 6"), section 0 item 5 ("is confirmed by the audit of section 6") and section 5 ("is
    confirmed by the audit below"): until section 6 is written, replace by "(audit of THM-4515–4517: pending, section
    6)" and drop "confirmed".
16. Section 2, Proposition 2 paragraph: "at level 19 `M(19) = 0.00982` at `k = 24`" — keep, and add "(fullperiod
    output; the level-19 line of that output was not present at audit time and is independently reproduced by the
    audit: `0.0098157`, mass `0.7089`)".
17. Section 1, second paper: "`sum_t conj(chi)(t) P_m(h_t) = tau(conj chi) B_(m,chi)/(m·m!)`" → "`... (odd m)`".
18. Section 7 table, row "the multiplicative spectrum ...": "VERIFIED (`n <= 16`)" — keep; row "Proposition 2 ...
    `M(n)` agrees with the S20 FFT at `n = 1..18` to seven digits" — keep and add "level 19 VERIFIED by the audit".

---

## 3. Decisive evidence, in one table

| item | own result | note |
|---|---|---|
| Thm 1(i), `max_j |F - tau S| / 3^(n/2)`, `n <= 16` | `<= 2.5e-16`; `|tau| = 3^(n/2)` to 12 digits | PROVED ✓ |
| Thm 1(ii), both identities, `n = 2..16` | equal to 9 decimals; `rho_n(1)` = S19 to 5 digits at all `n` | PROVED, VERIFIED ✓ (parity: 3–4 digits vs S19) |
| Prop 2, one step from the law's `m_(n-1)` | `<= 2.5e-16` (closed form, `n <= 16`; sequential, `n <= 10`) | PROVED ✓ |
| C recursion from level 0 vs law, elementwise | `8.7e-17` (`n = 12`), `7.5e-17` (`n = 16`) | exact ✓ |
| `M(n)` vs S20 FFT, `n = 1..18` | identical at 7 digits; argmax `±2^s` with S20's `s` | ✓ |
| `M(19)` | `0.0098157` at `±2^24`, then `±2^25` (`0.0097694`), `∓2^23` (`0.0087743`); offset `-6.11`; mass `0.7089` | `0.00982`, `k = 25, 23` next ✓ |
| rms, `P(z > 1,2,3)`, maxima, argmax characters, parity means, `n <= 16` | all reproduced | ✓ |
| census box B (cap 30) | histogram, cumulatives, Poisson tails, deepest lines identical | ✓ (d ≤ 8 deterministic) |
| census box A | identical; `P = 0.0089`, `0.0386` | ✓ |
| seeds `55·2^423+1`, `13·2^154-1`, `2^486-1`, `2^480-1` | `v_3 = 15, 7, 6, 2` | ✓ |
| `e^(-I) 3^(θ*/ln 2) - (m-1)` | `-1.1e-16`; margin `1.30%` | ✓ (tautological) |
| chain, factorisations, merge at 137 (10, 12), `-17` cycle | all exact | ✓ |
| basins `[1, 2^30]`, nine densities, all dyadic ranges | identical to 6 decimals; `|B(1729)| = 4719191`; `B(27)`: 26 numbers; nestedness violations `0` | ✓ |
| Krasikov–Lagarias | Theorem 6.1 of math/0205002 as quoted; ratio `0.1224`; threshold `2^48.9` | quote ✓; `2^49.3` → `2^48.9`; "astronomically beyond" ✗ |
