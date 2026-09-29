# Independent audit of `mazur_positive_density_20260928.md` (Mazur digest, Theorems A–C, spike profile, seed-1 test)

**Auditor:** independent session (Fable), 2026-09-28. **Audited state:** the note at commit
`0a7d62a99` (working copy identical to HEAD; 41755 bytes, mtime 22:15:16) and the scripts
`mazur_positive_density_20260928.py`, `mazur_harmonic_mass_20260928.py`,
`mazur_harmonic_mass_deep_20260928.py`, `mazur_seed1_test_20260928.py`,
`mazur_seed1_deep_20260928.py`, `kaprekar_cuboid_20260928.py` with their `.out` files, against the
paper's extracted text (`mazur.txt`, 1035 lines). The note was rewritten twice by the authoring
session while this audit ran (the level-38 and level-48 extensions of the seed-1 sequence, commits
`df907c0f0` and `0a7d62a99`); everything below refers to the final text, and the depth-decomposition
section (levels 19–48, `mazur_seed1_deep*_20260928.out`) is audited as well.

**Method.** `04-computation/experiments/mazur_positive_density_20260928_audit.py` (27 s, no
imports from the session's scripts; output `mazur_positive_density_20260928_audit.out`) recomputes
the 3-adic Syracuse law exactly to level 7 by the forward recursion (integer numerators) and, as a
cross-check, by the parents recursion; computes it in float64 to level 14 by the vectorised parents
recursion (a different algorithm from the session's discrete-log FFT); checks Terras classes by
brute-force scans, the affine identity, the converse and the sign statement; enumerates tree layers
with truncated valuations for eight seeds; verifies the cycle words, the resonance lower bounds, the
path-sum profile `D` by dynamic programming on the closure of `-1`, the harmonic-sum inequality,
`Pi(x)` to `2·10^5`, the identity `1/m = 3^n 2^-A Pi(m)` exactly, the empirical density
`54.22%`, the global statistics to level 14, the Kaprekar cycles and Euler bricks, Lemma 7.3's
arithmetic, and the Section 8.1 constants in iterated logarithms. Bracketed tags `[A1]`, `[C9]`,
... refer to lines of the `.out` file.

---

## (i) Claims checked

### A. Theorem A (harmonic mass of a tree layer)

1. **Terras class (proof step (i)).** VERIFIED. `[A1]` full scan for 165 words (`n <= 4`): the odd `x` with word `w` in `(-2^(A+3), 2^(A+3))` are exactly the 8 elements of one class mod `2^(A+1)`, namely `(2^A - C_w) 3^-n`; `[A1b]` 900 random words of length 5–8. The induction step of the note is right: `x = c + 2^(a_1+1) t` gives `x_1 = c' + 6t`, and `t mod 2^B -> x_1 mod 2^(B+1)` is a bijection onto the odd residues.
2. **Affine identity `3^n x + C_w = 2^A S^n(x)`, `C_w` odd.** VERIFIED. `[A2]` 3000 random odd `x` of both signs, `n <= 12`; `C_w = 3^(n-1) + (even terms)`.
3. **`S^n(x) = Y_n(w) mod 3^n`, `Y_n(w) = C_w 2^-A` = the recursion `Y_(k+1) = 2^-a(3Y_k + 1)`.** VERIFIED. `[A2]`; `C_w 2^-A = sum_j 3^(n-j) 2^-(a_j+...+a_n)`, Tao's `F_n`.
4. **Converse (iii).** VERIFIED. `[A3]` every word `n <= 4`, `a <= 4`, every odd `y` in `(-4·3^n, 4·3^n)` with `3 !| y`: 1360 cases where `y = Y_n(w)` give an odd integer `x` with word `w`, `S^n(x) = y`, sign of `y`; 58960 cases where `y != Y_n(w)` give a non-integer.
5. **Sign statement (v).** VERIFIED. `[A3]`, `[A4]`; `3x + 1 <= -2` for `x <= -1`.
6. **The identity for `y = 1, 5, 7, 11, -1, -5, -7, -17`, `n <= 5`.** VERIFIED. `[A4]` truncated layer sums (`amax = 8, 12, 16, 20`) increase to the exact `3^n mu_n(y mod 3^n)` with gaps `10^-2 -> 10^-6`.
7. **`3 | y` gives `0` on both sides.** VERIFIED (`S` never outputs a multiple of 3; `Y_n = 2^-a_n mod 3`).
8. **`H_1 = 1`, `H_2 = 8/7`, `H_3 = 1376/1387`, `H_4..6 = 0.927627, 0.964301, 0.955420`.** VERIFIED. `[A5]` exact law (also `H_7 = 0.860359`).
9. **Corollary A1 (means `1`, `2/3`, `4/3`, zero on nonunits).** VERIFIED. `[A7]` at `n = 1, 6, 10, 14`.
10. **Corollary A2 (`H_(n+1) = H^(1) + 2H^(2)`), depth-1 split `16/21, 4/21, 1/21`.** VERIFIED. `[A5]`.
11. **First-arrival masses `1/4, 0.393, 0.135, 0.184, 0.269, 0.232` and `H_n = sum (3/4)^(n-j) H_j^first`.** VERIFIED. `[A5]`, `[C5]`: `H_n^first = H_n - (3/4) H_(n-1) = 0.25, 0.392857, 0.134926, 0.183575, 0.268581, 0.232...`.
12. **"The recursion reproduces each next first-arrival mass exactly."** VERIFIED (children of first-arrival nodes are first arrivals; `1` is only its own child).
13. **Consistency (`Y_n mod 3^m` has the law of `Y_m`).** VERIFIED. `[A6]` exact to level 7; proof: `Y_n mod 3^m = F_m(a_(n-m+1), ..., a_n)` with i.i.d. valuations.
14. **Forward and parents recursions agree.** VERIFIED. `[A6]` exact to level 5; `[D]` pointwise; float law vs exact `4·10^-17`.
15. **§0.1 "his weighted inverse histories ... are exactly the harmonic mass of the Syracuse tree of `M`".** VERIFIED WITH CORRECTION. Mazur's `omega(w)` is a summand of `H_d(M)` (as §4 "Reading" says); his sums `Z_N(M)` run over the restricted family (central first-crossing words, unit intermediate endpoints, terminal restriction; paper §3, §5), not over the full layer.
16. **Paper §2: "these transfers agree with the actual inverse-orbit sums".** VERIFIED (session's P3; my item 4).

### B. Theorem B (cycle resonances)

17. **Cycle data: `-1` (`k = 1, A = 1`), `{-5, -7}` (`k = 2, A = 3`, word `(1, 2)`), `-17` (`k = 7, A = 11`, word `(1,1,1,2,1,1,4)`), `y_0 = C_w/(2^A - 3^k)`.** VERIFIED. `[B]`.
18. **`mu_(km)(y_0) >= 2^(-Am)`, `rho >= (2/3)(3^k/2^A)^m`.** VERIFIED. Proof valid (`y_0 ∈ T_(km)(y_0)` with word `w^m`); numerically at all levels `<= 14`.
19. **`Y_n(1^n) = (3/2)^n - 1 = -1 mod 3^n`.** VERIFIED. `[B]` exact, `n <= 14`.
20. **Negative cycles have `3^k > 2^A`.** VERIFIED (`y_0 (2^A - 3^k) = C_w > 0`).
21. **`s = 0.3691, 0.0536, 0.0086`.** VERIFIED WITH CORRECTION. `1 - (11/7) log_3 2 = 0.008539`, i.e. `0.0085`, not `0.0086`.
22. **"The limit density has singularities of order `|y - y_0|_3^(-s)`" (§0.2, §5).** VERIFIED WITH CORRECTION. Only the lower bound is proved (ball averages `>= (2/3) 3^(sn)`); the limit's existence is CONDITIONAL; write "of order at least".
23. **`rho_n(-1)/(3/2)^n -> 0.9748`.** VERIFIED WITH CORRECTION. The ratio is increasing (`0.9712` at 8, `0.97456` at 14, `0.9748` at 18) with halving increments; `0.9748` is the level-18 value, the limit is `>= 0.9748` (about `0.975`).
24. **"the `{-5, -7}` cycle grows like `(9/8)^(n/2)`, the seven-cycle of `-17` like `(2187/2048)^(n/7)`".** VERIFIED WITH CORRECTION. Only lower bounds are proved. Observed: `rho(-5), rho(-7)` grow by `1.16–1.20` per two levels at `n <= 18` (bound `1.125`); `rho_n(-17)` is non-monotone through level 18 (`1.79, 1.34, 1.46, 1.87, 2.38, 2.24, 2.11, 2.47, 1.98` for `n = 10..18`): no growth law is visible for `-17`.
25. **`sup rho_n >= (2/3)(3/2)^n`, no level-uniform bound.** VERIFIED.
26. **The class of `-1` is the maximal atom at every level.** VERIFIED to level 14 `[B]`; levels 15–18 not scanned (the note now says so).
27. **`rho_18(-5) = 6.3`, `rho_18(-7) = 7.6`, `-17` at `2.0–2.5`.** VERIFIED against `mazur_harmonic_mass_deep18_20260928.out` (my level-14 values `4.677, 5.558, 2.383` coincide with the session's).

### D. The spike profile

28. **Recursion `rho_(n+1)(z) = 3 sum_(a = eps(z) mod 2) 2^-a rho_n((2^a z - 1)/3 mod 3^n)`.** VERIFIED. Derived from the definition (`3Y_n + 1 = 2^a z mod 3^(n+1)` forces the parity of `a` and `Y_n = (2^a z - 1)/3 mod 3^n`); `[D]` exact pointwise check; my entire float law is built on it.
29. **`D(-1/2) = 1/2`, `D(-1/4) = D(1/8) = D(11/16) = D(49/32) = 3/4`, `D(-1/8) = D(1/16) = D(11/32) = 3/8`, `D(-1/16) = 3/16` as path sums.** VERIFIED. `[D]` dynamic programming on the closure (denominators `<= 2^9`), `D(-1) = 1`, `D = 0` off the closure.
30. **`x_j = 3^(j+1)/2^(j+2) - 1` is the all-ones forward orbit of `-1/4`, `x_j = -1 mod 3^(j+1)`, `D(x_j) = 3/4`.** VERIFIED. `[D]` (`S_1(x_j) = x_(j+1)`; path sums give `3/4` for `j <= 7`). The induction for all `j` uses "`D = 0` off the closure", a heuristic, as the note's OBSERVED label allows.
31. **Numerical agreement of the ratios with `D`.** VERIFIED. `[D]` level 14: the nine named points within `3·10^-4` (`0.5000, 0.7499, 0.7499, 0.7498, 0.7497, 0.3750, 0.3750, 0.3749, 0.1875`); over all 33 closure points with denominator `<= 2^6` the maximal deviation is `0.022, 0.011, 0.006, 0.003` at levels `8, 10, 12, 14` (at `x = 11/64`): convergence at roughly a factor 2 per two levels, consistent with four decimals at 17–18.
32. **"The second tier of atoms at level `n` is the set of classes `x_j mod 3^n`, `j < n`, each at `3/4`; `n = 8`: `-1: 24.9`, then five classes at `18.7`".** VERIFIED WITH CORRECTION. `[D]` at levels 8, 10, 12 the classes with `rho >= 0.7 rho(-1)` are exactly the `n - 1` shadows `x_j`, `j <= n - 2` (`x_(n-1) = -1 mod 3^n` is the top atom itself), at ratios `0.746–0.752`; at `n = 8` there are seven such classes (the note's "five" is the truncated top-6 list).
33. **Lower bound `rho_n(x) >= (2/3)(3/2)^n 2^(depth - A_0)` on the closure (PROVED).** VERIFIED: the word `(1^(n-depth), path)` lands in `x mod 3^n` (each `S_a` multiplies a discrepancy `3^(n-depth) u` by `3/2^a`).
34. **"`D = 0` off the closure."** NOT PROVED (heuristic; correctly inside the OBSERVED label).

### C. Theorem C (the seed-1 test)

35. **Statement.** VERIFIED as correctly stated (hypothesis: positive lower natural density for some `C`; conclusion: Cesàro means of `n^(1/6) H_n` bounded below, hence `limsup > 0`).
36. **(a) Partial summation.** VERIFIED: `N(t) >= c t` for `t >= X_0` gives `sum_(x∈F, x<X) 1/x >= ∫_(X_0)^X N(t)/t^2 dt >= c ln(X/X_0)`.
37. **(b) Odd parts.** VERIFIED: `tau(2^k m) = k + tau(m)`; `n(m) <= tau(m)/2 <= tau(m)` (`[C7]` odd `m <= 10^5`); `sum_k 2^-k = 2`; `m = 1` gives the additive `2`.
38. **(c) `1/m = 3^n 2^-A Pi(m)`.** VERIFIED. `[C3]` exact for `m = 3, 9, 27, 97, 871, 993`; telescoping `2^(a_j) x_j = 3 x_(j-1)(1 + 1/(3x_(j-1)))`.
39. **Distinctness of `x_0, ..., x_(n-1)` (all `> 1`).** VERIFIED. `[C4]` odd `m <= 20001`; a repeat forces periodicity, hence an earlier `1`.
40. **`sum_(k<=n) 1/(2k-1) <= 1 + (1/2) ln n` for all `n >= 1`.** VERIFIED. `[C1]` to `n = 10^6` (slack `0` at `n = 1`, `0.01824` at `10^6`, limit `1 - ln 2 - gamma/2`); proof: the increment `1/(2n+1) - (1/2) ln(1 + 1/n)` is `<= 0` because `ln(1 + 1/n) = 2 artanh(1/(2n+1)) >= 2/(2n+1)`, so the difference is non-increasing from its value `0` at `n = 1`.
41. **`Pi <= e^(1/3) n^(1/6)`.** VERIFIED. `[C2]` ratio `<= 0.7643`.
42. **`G_n <= e^(1/3) n^(1/6) H_n^first <= e^(1/3) n^(1/6) H_n`.** VERIFIED. `[C5]` `n <= 6`.
43. **(d) Cesàro/limsup.** VERIFIED: `sum_(n<=N) n^(1/6) H_n >= (c/(2e^(1/3)))(N/C - ln X_0) - e^(-1/3)`; `limsup a_n >= liminf` of the Cesàro means.
44. **`max Pi = 1.2531` at `x = 993`, mean `1.167`, `Pi/(e^(1/3) n^(1/6)) <= 0.764`.** VERIFIED. `[C2]` (`1.25314` at `993`, mean `1.1666` over `x <= 2·10^5`, `0.7643`).
45. **`Pi(9) = 1.2486` (`sum 1/x_j = 0.68`), `Pi(27) = 1.1989`.** VERIFIED against `mazur_seed1_test_20260928.out` (identity checked exactly at `m = 9, 27`).
46. **"`H_n = o(n^(-1/6))` would refute it, for every `C`."** VERIFIED.
47. **No claim that Mazur's theorem is false or that `H_n -> 0`.** VERIFIED: verdict OPEN; "decay is excluded through `n = 48`" is a statement about lower bounds `>= 0.36` at `n <= 48`, correct.
48. **`n^(1/6) H_18 = 0.665`.** VERIFIED (`0.6645`).
49. **"Every ratio `H_(n+1)/H_n` from `n = 7` to `n = 21` is below `1`."** VERIFIED (in fact from `n = 5`: `H_6/H_5 = 0.991`).
50. **"falling about `3%` per level from level 8 to level 22".** REFUTED as a number. `[C10]` average decline `5.1%` per level (`6.2%` over 8–18, `2.6%` over 18–22).
51. **Parents form at `z = 1` and the fixed-point extrapolation `H_∞ ≈ 0.35`.** VERIFIED. `[F]` `rho_14(1) = 3 sum 4^-j rho_13(R_j)` to `10^-6`; extrapolation `0.416` (level 14), `0.379` (17), `0.355` (18): drifting, as the note says.
52. **Data row `n <= 18`.** VERIFIED to level 14 (`rho_n(1)` identical to five decimals `[F]`); levels 15–18 not recomputed here (the session's FFT, consistent with its own parents-form assertion at every level).

### I. The depth decomposition (levels 19–48)

53. **Identity `H_(m+d)(1) = sum_(y ∈ T_d(1)) 3^d 2^(-A_d(y)) H_m(y)`.** VERIFIED (Theorem A at the depth-`d` ancestor; the fibres `{x : S^m(x) = y}` are disjoint; weights multiply).
54. **The computed values are lower bounds (up to float error).** VERIFIED (all terms nonnegative).
55. **"lower bounds whose loss is at most the pruned weight"; "pruning loss `<= 3·10^(-2)`".** REFUTED as a bound. The loss is `sum_(pruned z at depth d') W(z) H_(m+d-d')(z)`, and `H` is not bounded by `1` (it reaches `2·10^3` on the class of `-1` at level 18); the validation shows loss `≈` pruned weight (ratios `1.00` at `d = 10, 15, 18`, `0.6` at `d = 2`), an empirical statement, exactly the "heuristic correction" the session's own script prints.
56. **"minimum `0.370` at `n = 22`", "bottoms at `0.370` near `n = 22`", "the minimum is `H_22 = 0.3697`", "rises to `0.544` at `n = 43`", "`0.513` at `n = 48`", "`n^(1/6) H_n ≈ 1.0` for `n >= 41`".** REFUTED (the minimum) / VERIFIED (the rest). `[C9]`, `[C10]`: in both runs the minimum of the computed sequence is `0.3634` at `n = 27` (`0.36336`, corrected `0.36338`); `n = 26, 27, 28, 29` all lie below `H_22 = 0.36971`; the sequence rises monotonically only from `n = 28`. Maximum `0.54371` at `n = 43`, `H_48 = 0.51337` (corrected `0.542`), `n^(1/6) H_n = 0.98–1.02` for `n >= 41`: correct.
57. **Validation deficits equal to pruned weights; the two runs agree at `n = 38` up to their pruned weights.** VERIFIED from the outputs (`0.44510 + 0.0103 ≈ 0.45187 + 0.0036`).
58. **"Two routes (`(d, 18)` and `(d+1, 17)`) agree to `10^(-3)`."** VERIFIED WITH CORRECTION: to `2.4·10^-3` (depth-20 run, `n = 36`) and `3.6·10^-3` (depth-30 run, `n = 44`: `0.54034` vs `0.53670`), the second route carrying one more pruned level.
59. **Class shares: "`f^(0)` sits at `0.34–0.35` from depth 6 on"; "from depth 22 on all three shares are `1/3 ± 0.02`".** VERIFIED WITH CORRECTION: `0.325–0.354` from depth 6; `0.355` at depth 30 (`± 0.022`).
60. **`rho_2 = (16, 32, 22, 8, 4, 44)/21` on `1, 2, 4, 5, 7, 8 mod 9`; the poor classes `5, 7 mod 9` are those whose dominant child is a leaf.** VERIFIED. `[A8]` exact.
61. **`R_7 = 5461 = 7 mod 9`, its `a = 2` child `7281` is a leaf, `rho_18(R_7) = 0.047`; `rho_18(5) = 0.28`.** VERIFIED (arithmetic; values from the level-18 output).

### E. Global statistics of the law

62. **`rho_n(1)` for `n <= 10` (`0.66667, 0.76190, 0.66138, 0.61842, 0.64287, 0.63695, 0.57357, 0.51638, 0.46450, 0.42497`).** VERIFIED. `[F]` identical to five decimals.
63. **`E[rho^2]` over units `2.667` (`n = 6`), `3.911` (10), `5.162` (14).** VERIFIED (`2.6669, 3.9105, 5.1618`).
64. **Entropy deficit `0.935`, `1.062`, `1.131` at `n = 6, 10, 14`.** VERIFIED (`0.9351, 1.0615, 1.1309`).
65. **Medians `0.587, 0.533, 0.511`; shares below `0.1`: `0.043, 0.058, 0.066`; `||f_n - f_(n-1)||_1 = 0.263, 0.178, 0.130`.** VERIFIED.
66. **"`E[rho_n^2]` grows by `0.31` per level" (`0.31 n + 0.8`).** VERIFIED (increments `0.3103 -> 0.3135`, slowly increasing).
67. **"the entropy deficit converges (increments ... ratio `0.9`), so the information dimension is `1` and `E_units[rho log rho] -> 0.84`".** VERIFIED WITH CORRECTION. OBSERVED only: the increment ratios rise with `n` (`0.835` at `n = 8`, `0.879` at 14, `≈ 0.9` at 18), so the data do not separate convergence from a slow (logarithmic) divergence. `E_units[rho ln rho] = deficit - ln(3/2)` (exact identity, `[F]`) is `0.767` at `n = 18`; `0.84` is the extrapolated limit under geometric convergence of the deficit to `≈ 1.25`.
68. **"information dimension `1`".** VERIFIED as OBSERVED (needs only a sublinear deficit).
69. **"decays like `n^(-0.93)`".** VERIFIED WITH CORRECTION. Range-dependent: my fits `-0.83` (6–14), `-0.95` (10–14), `-1.01` (12–14); the session's values give `-1.04` (10–18), `-1.16` (14–18). "About `1/n`, steepening" is what the data show.
70. **"tail `P(rho > t) ~ t^(-2)`".** VERIFIED WITH CORRECTION: heuristic; at level 14 the local exponents run `1.3 -> 2.35` across `t = 1..32` and `t^2 P(rho > t)` peaks near `t = 8`; consistent, not established.
71. **Fine-scale distances at `N = 18` (`0.854 ... 0.354`, decreasing in `m`).** VERIFIED against the output; my `N = 14` values `0.840 ... 0.307` have the same shape.
72. **§7: density "infinite on the forward rational closure of every negative cycle, not in `L^2`, with finite entropy".** VERIFIED WITH CORRECTION. The local averages `f_n(x) -> ∞` on the closure points are PROVED (items 18, 33); `f_∞ ∉ L^2` follows under (2.3) from Jensen (`∫ f_n^2 <= ∫ f_∞^2`) plus the OBSERVED unbounded growth of `E[rho^2]`; "finite entropy" is the OBSERVED (extrapolated) convergence of item 67. Label the two as OBSERVED inside the CONDITIONAL paragraph.

### F. The paper digest (sections 1–3)

73. **Theorem 1.1 as stated (count of `1 <= n < X`, natural log, `c`, `X_0` of Section 8).** VERIFIED (`mazur.txt` 39–46).
74. **Corollary 1.2 (`262/25 = 10.48`).** VERIFIED (line 48).
75. **Author and version: "M. Mazur, ..., v2, September 2026".** REFUTED. The paper is by **Lech Mazur** ("L. Mazur"), "Date: September 6, 2026; version 2.1" (lines 3, 55); the session's own script header has it right ("L. Mazur ... v2.1, 2026-09-06"). The wrong initial is repeated in the atlas §8, the synthesis and the ledger.
76. **Lemma 2.1 (`<|T_w g|>_(t+d) = 2^(-A) <|g|>_t`, with and without absolute values; single-point fibres; `3^d` cancels).** VERIFIED (lines 158–164).
77. **(2.3) as quoted (`<= (2/3) C_A m^-A`, `1 <= m <= q`; "finite-distribution form of Tao's Prop. 1.14"; exponents 6 and 2; Lemma 8.1 gives `C` for exponent 6).** VERIFIED (lines 127–144, 793–798).
78. **(3.9): "the number of central histories landing in a class is bounded by a polynomial in the generation".** VERIFIED WITH CORRECTION. The bound is `(2W* + 1)(K* + 4n + 1)(2^(b_0+1) + 16^(b_0)(9/16)^(b_n))` with `W* <= 2n b_0^(3/5) G^n`, `G = 2013/2000`, i.e. `(n+1)^3 G^n` in (3.10): polynomial times a slowly growing exponential (lines 311–330).
79. **Section 8.1 formulas: `A* = 6409`, `E* = 2170`, `L* = 280`, `D_exp = 2^(8192 A* 2^(3E*))`, `C* = (32 A* D*)^A*`, `N = 20000(ceil(log_2 F) + 64)`.** VERIFIED as quoted (lines 768–792, 824–829). `F = 2467 b 16^b (C+1)` versus the script's `2^467 b^(16b) (C+1)`: NOT CHECKED (the extraction "2467b16b" is ambiguous; immaterial for the magnitudes).
80. **Magnitudes: "`D_exp ≈ 2^(2^6536)`", "`log_2 log_2 C* ≈ 6548`", "`log_2 N ≈ 6563`", "explicit constant `2^(2^6536)`-sized", "says nothing below `m` of the order of `2^(2^6536)`", "`c^(-1)` beyond `2^(2^(2^6535))`", "`X_0 ... likewise`".** REFUTED except for `D_exp` and the lower bound on `c^(-1)`. `[H]`: `D* = max(D_1, D_2, D_3)` is dominated by `D_sc >= P*^10`, where `P* = g^(R*-1)(T*)` is a cubic map composed `R* - 1 ≈ 2^8696` times (`T* = 10 A* 2^(3E*)`, `R* = 2^E*(T* + 24586) + 1`), so `log_2 log_2 log_2 P* ≈ 8696.6`. Hence `log_2 log_2 log_2 C* ≈ 8697` (not `log_2 log_2 C* ≈ 6548`), `log_2 log_2 N ≈ 8697` (not `log_2 N ≈ 6563`), the mixing coefficient `C` is a three-fold tower `2^2^2^8697` and Lemma 8.1 is vacuous below `m` of that order (not `2^(2^6536)`); `c^(-1) = 256 M^2 m/3` with `log_2 M ≈ 6q ≈ 960·280·1.01^N` is a four-fold tower `2^2^2^2^8697` (so "beyond `2^(2^(2^6535))`" is true but understates by a level), and `X_0 = 32(2^B M + 1)` with `B ≈ 42 beta(J)`, `J ≈ 32000 beta(N)`, is a five-fold tower, one level above `c^(-1)`, not "likewise". The session's `part8` assumed `D* = D_exp`.
81. **"an explicit certificate, not a claim of numerically usable constants".** VERIFIED (line 1010).
82. **"`523/50 = 10.46` is `3/log(4/3) = 10.428` with rounding losses".** VERIFIED. `[H]` (7.9): `kappa(L + h) + 1 = 7.25012 <= (523/50) L = 7.25032` with `kappa = 34881/10000`, `h <= ln 3 + 1/12288` (margin `2·10^-4`); replacing `kappa` by `1/ln(4/3)` gives exactly `3/ln(4/3) = 10.4282`.
83. **Lemma 4.1 / Proposition 4.2 as summarised (seeds `R_j` permute `Z/3^q`; a fixed odd `M >= 16 b_0`, `3 !| M`, `S(M) = 1`).** VERIFIED (lines 366–377; LTE proof `v_3(R_j - R_k) = v_3(j - k)`, `[H]`).
84. **"coverage (Proposition 6.3): the source charge `x omega(w) <= M`".** VERIFIED WITH CORRECTION. The source charge is Lemma 6.1 (`3^d x < 2^A M`, hence `x omega(w) <= M`), coverage of all large scales is Lemma 6.2, and Proposition 6.3 is the count (lines 545–582).
85. **Section 9 quotations; axioms `propext`, `Classical.choice`, `Quot.sound`.** VERIFIED (lines 1013–1022).
86. **Platform status (386 files, 61,169 lines, toolchain, mathlib hash, commit `830b9d3f`, the Lean statement text, the companion page and its three theorem names).** NOT CHECKED (no platform access; not in the PDF).
87. **Krasikov–Lagarias `x^0.84` (from memory).** VERIFIED (line 50). Omission: the paper's introduction cites Mazur's own certified `pi_a(X) >= X^0.90` (ref. [3], July 2026) and `>= c_a X^0.901`; the note's "the Krasikov–Lagarias `x^0.84` row would be replaced by positive density" skips the already-existing `x^0.90` step (CITED via this paper).
88. **Tao's Proposition 1.14 is literally (2.3).** NOT CHECKED against Tao's text (not available here). Consistent with my recollection: Tao's `Osc_(m,n)(X)` is the `ell^1` distance between the law of `X` on `Z/3^n` and the uniform lift of its projection to `Z/3^m`, and Proposition 1.14 bounds it by `m^-A`; with consistency (item 13) this is (2.3) up to the factor `2/3`.
89. **"(2.3) as stated is the `L^1`-Cauchy property of the Haar densities, hence absolute continuity of the 3-adic law, with `||f_∞ - f_m||_1 <= C_A m^-A`."** VERIFIED. `||rho_q - rho_m∘pi||_q = (2/3)||mu_q - lift(mu_m)||_1 = (2/3)||f_q - f_m||_(L^1(Haar))`; `f_m = E[f_n | level m]` (item 13); `sup_(q>=m) ||f_q - f_m||_1 <= C_A m^-A` is Cauchy in the complete space `L^1`; the limit of `mu_n(B)` on every ball is `∫_B f_∞`. CONDITIONAL on (2.3), as labelled.
90. **P1–P8 table rows.** VERIFIED that the `.out` supports P1–P5, P7 (and P8's formulas). P6 is vacuous as run: the seed `R_40`, depth 3, valuations `<= 5` yield 5 histories in 5 distinct `(D, A)` classes (`[H]`), so the spacing and equal-endpoint assertions never compare two endpoints (both statements are trivialities anyway). P2's assertion `x.denominator != 1 or x <= 0 or any(w % 2 == 0 for w in [1])` has dead branches (harmless).
91. **"`54.22%` of `n < 10^6` satisfy `tau(n) <= 10.46 ln n`".** VERIFIED (`[C6]` `0.5422`).

### G. Atlas typing (§7 of the note; atlas §0 rule)

92. **DRIFT O.** VERIFIED as defensible: the conclusion (positive density reaching `1` in log time) is conjecturally false for `5x+1`; the proof uses the drift through `A ≈ 2d`.
93. **STICKY O.** VERIFIED as defensible, with the same hedge as the atlas's Tao row: the aliquot analogue is *conjecturally* false (parity lock: `s(2^a m)` with `m` a non-square is even, so reaching `1` needs a visit to the density-zero set `2^a·square`; unprovable at present), and the feature named (size-free exact inverse-orbit weights) is real (item 16).
94. **SHEET B, DEFECT B, INTEGRAL B, DIM B.** VERIFIED as defensible under the conclusion rule (`3x-1` conjecturally has positive-density log-time basins; null sets and rational cycles do not touch a density statement; no descent certificate). The companion theorem invoked for SHEET is CITED, NOT CHECKED.
95. **UNIFORM B "(presumably runs unchanged for `3x+k`; not checked)".** Asserted as a presumption; acceptable, but the atlas §8 row should carry the label CONJECTURAL, since the explicit constants and the seeds are map-specific.
96. **"The pattern of atlas §7 holds: a residue-averaging mechanism."** VERIFIED (atlas §7: "every mechanism that overcomes STICKY is a residue-averaging mechanism, blind to SHEET and DIM").
97. **§7(iii): "the S18 observation that `80%` of visited values above the start range are `2 mod 3` is the class mean `4/3` of Corollary A1 seen from the orbit side".** VERIFIED WITH CORRECTION. The class mean `4/3` says `P(Y_n = 2 mod 3) = 2/3` (`67%`), the unconditional share; `80%` is a conditional statistic (values above the start range were reached by small valuations, and `a = 1` gives residue `2`); the identification is loose.

### H. Status hygiene

98. **Header: "FINITE-EXACT (the 3-adic Syracuse law to level 18, rational to level 6)".** VERIFIED WITH CORRECTION: levels 7–18 are float64 FFT, verified against exact rationals only to level 6 (the status table says "18 (float)" correctly); the header should not call them FINITE-EXACT.
99. **§0.2: "The relative profile of the `-1` spike ... is exact:"** OBSERVED stated as fact (item 31, 34); §5 labels it correctly.
100. **§0.2: "grows like", "singularities of order", "`-> 0.9748`".** Lower bounds and finite-level values stated as asymptotics (items 22–24).
101. **§0.4 / §5: "converges", "`-> 0.84`", "`n^(-0.93)`", "`t^(-2)`".** Extrapolations stated as facts (items 67, 69, 70).
102. **§6 and status table: "loss is at most the pruned weight", "pruning loss `<= 3·10^(-2)`", "FINITE-EXACT evaluation".** A heuristic stated as a bound; a float lower-bound computation labelled FINITE-EXACT (item 55).
103. **§0.3, §6, title, status table: the minimum at `n = 22`.** Arithmetic/reading slip (item 56).
104. **§0.3: "`3%` per level".** Wrong number (item 50).
105. **§5: `s = 0.0086`.** Rounding slip (item 21).
106. **Numbers without a backing file.** The platform figures (item 86) only; every other number in the note traces to an `.out` file or to the paper.
107. **Kaprekar cycles `d = 2..6`, `9 | image`, digit-multiset invariant; Euler bricks `(44,117,240)`, `(240,252,275)` below 300, no integral space diagonal.** VERIFIED (`[G]`).
108. **§8 (inverter): `Phi(a) = -sum 2^(a_1+...+a_(j-1))/3^j`; periodic words give odd-denominator rationals; a rational `p/q` has the word of `p` under `x -> (3x+q)/2^v`; THM-4476 corollaries (2), (3), (5), (8) as quoted.** VERIFIED (2-adic limit of `(2^A x_n - C_w)/3^n`; `S(p/q) = (3p+q)/(q 2^v)`; the numbered corollaries of `THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md` say what the note attributes to them).

---

## (ii) Corrections (old text -> new text), current note at `0a7d62a99`

**C1 (title, line 1).** "`H_n(1)` bottoms at `0.370` near `n = 22`, rises to `0.54` by `n = 42` and stays near `0.52` through `n = 48`" -> "`H_n(1)` falls to `0.363` at `n = 27` (local minimum `0.370` at `n = 22`), rises to `0.54` by `n = 42` and stays near `0.52` through `n = 48` (lower bounds)".

**C2 (§0.3, lines 84–85).** "falling about `3%` per level from level 8 to level 22, where it bottoms at `0.370`." -> "falling about `5%` per level on average from level 8 to level 22 (`0.370`), with its minimum `0.363` at `n = 27`."

**C3 (§0.3, lines 87–88).** "loss `<= 3·10^(-2)`, measured against the exact `H_18` at every split" -> "loss about the pruned weight, `3·10^(-2)` at `n = 48`, an empirical calibration against the exact `H_18` (the pruned subtrees are not bounded by their weight; their mass is measured, not bounded)".

**C4 (§6, lines 419–421).** "lower bounds whose loss is at most the pruned weight" -> "lower bounds whose loss is empirically about the pruned weight (validation at `m + d = 18`: deficit/pruned weight `= 1.00` at `d = 10, 15, 18`); no bound on the loss is proved, since `H_18(z)` on the pruned nodes is not bounded by `1`".

**C5 (§6, line 426).** "the minimum is `H_22 = 0.3697`; from `n = 22` to `n = 43` the ratios are mostly above `1`" -> "`H_22 = 0.3697` is a local minimum; the minimum of the computed sequence is `H_27 = 0.3634` (`n = 26..29` all lie below `H_22`); from `n = 28` to `n = 43` every ratio is above `1`".

**C6 (§6, lines 473–474).** "The sequence has a minimum `0.370` at `n = 22`, rises to `0.544` at `n = 43`" -> "The sequence has its minimum `0.363` at `n = 27` (after a local minimum `0.370` at `n = 22`), rises to `0.544` at `n = 43`".

**C7 (status table, line 631).** "minimum `0.370` at `n = 22`, rising to `0.544` at `n = 43`, `0.513` at `n = 48` (pruning loss `<= 3·10^(-2)`)" -> "minimum `0.363` at `n = 27` (local minimum `0.370` at `n = 22`), rising to `0.544` at `n = 43`, `0.513` at `n = 48` (lower bounds; pruning loss about `3·10^(-2)` at `n = 48`, calibrated, not bounded)". The same three corrections (minimum at `n = 27`; "about", not "`<=`"; "3%" -> "5%") apply to the atlas §8 (lines 223–226), `PROBLEM-LEDGER.md` line 597, the synthesis lines 1539–1541 and the two commit messages.

**C8 (§6, line 450, "FINITE-EXACT evaluation").** "PROVED identity, FINITE-EXACT evaluation" -> "PROVED identity; float64 lower bounds with an empirical loss estimate".

**C9 (§2.7, lines 205–208).** "`D_exp = 2^(8192 A* 2^(3 E*))` (about `2^(2^6536)`), `C* = (32 A* D*)^A*` (`log_2 log_2 C* ≈ 6548`), `N = 20000(ceil(log_2 F) + 64)` with `F = 2467 b 16^b (C+1)` (`log_2 N ≈ 6563`), and `c^(-1)` beyond `2^(2^(2^6535))`." -> "`D_exp = 2^(8192 A* 2^(3 E*))` (about `2^(2^6536)`) but `D* = max(D_1, D_2, D_3)` is dominated by `D_sc >= P*^10` with `P* = g^(R*-1)(T*)`, a cubic map composed `R* ≈ 2^8696` times, so `log_2 log_2 log_2 D* ≈ 8697`; `C* = (32 A* D*)^A*` (`log_2 log_2 log_2 C* ≈ 8697`), `N = 20000(ceil(log_2 F) + 64)` (`log_2 log_2 N ≈ 8697`); `c^(-1)` is a four-fold tower `2^(2^(2^(2^8697)))` and `X_0` a five-fold one."

**C10 (§3, P8 row, line 224).** "`D_exp ≈ 2^(2^6536)`, `log_2 log_2 C* ≈ 6548`, `log_2 N ≈ 6563`, `c^(-1) > 2^(2^(2^6535))`" -> "`D_exp ≈ 2^(2^6536)`; `D* ≈ D_sc`, `log_2^(3) C* ≈ 8697`, `log_2^(2) N ≈ 8697`, `log_2^(4) c^(-1) ≈ 8697`, `log_2^(5) X_0 ≈ 8697` (the script's `part8` took `D* = D_exp`; corrected in the audit script, section H)".

**C11 (§0.4, line 106).** "whose explicit constant here is `2^(2^6536)`-sized" -> "whose explicit constant here is a three-fold tower `2^(2^(2^8697))`".

**C12 (§5, lines 375–376).** "says nothing below `m` of the order of `2^(2^6536)`" -> "says nothing below `m` of the order of `2^(2^(2^8697))`". Propagate to the atlas §8 line 204 ("an explicit coefficient of size `2^(2^6536)`" -> "of size `2^(2^(2^8697))`") and line 195 / ledger line 589 (`c^(-1) > 2^(2^(2^6535))` is true; better `c^(-1) ≈ 2^(2^(2^(2^8697)))`).

**C13 (§1, lines 12–13).** "M. Mazur, *Explicit Positive-Density Collatz Convergence in Logarithmic Time*, v2, September 2026" -> "L. Mazur (Lech Mazur), *Explicit Positive-Density Collatz Convergence in Logarithmic Time*, version 2.1, September 6, 2026". Same in atlas §8 line 188, synthesis line 1506, ledger line 587 ("v2").

**C14 (§0.2, lines 70–72).** "the `{-5, -7}` cycle grows like `(9/8)^(n/2)`, the seven-cycle of `-17` like `(2187/2048)^(n/7)`. The limit density has singularities of order `|y - y_0|_3^(-s)`" -> "the `{-5, -7}` cycle grows at least like `(9/8)^(n/2)` (observed factor `1.16–1.20` per two levels through `n = 18`), the seven-cycle of `-17` at least like `(2187/2048)^(n/7)` (no growth visible through `n = 18`). Where the limit density exists (§7) it has singularities of order at least `|y - y_0|_3^(-s)`".

**C15 (§0.2, lines 74–75).** "The relative profile of the `-1` spike on its forward rational closure is exact:" -> "The relative profile of the `-1` spike on its forward rational closure is, to four decimals at levels 17 and 18 and as a proved lower bound:".

**C16 (§0.2 line 69 and §5 line 322).** "`rho_n(-1) = 0.9748 (3/2)^n`" / "`rho_n(-1)/(3/2)^n -> 0.9748`" -> "`rho_n(-1) = 0.9748 (3/2)^n` at `n = 18` (the ratio increases with `n`; limit about `0.975`)".

**C17 (§5, line 311).** "`s = 0.3691` at `-1`, `0.0536` at `-5` and `-7`, `0.0086` at the seven points" -> "... `0.0085` at the seven points" (also §0.2 line 73: "`0.3691, 0.0536, 0.0086`" -> "`0.3691, 0.0536, 0.0085`").

**C18 (§0.4, lines 99–101, and §5 lines 364–366).** "Its entropy deficit `n ln 3 - H(mu_n)` converges (..., increments shrinking geometrically; information dimension `1`, `E[rho log rho] -> 0.84` on units)" -> "Its entropy deficit `n ln 3 - H(mu_n)` appears to converge (..., increments shrinking, but their ratios rise from `0.84` at `n = 8` to `0.9` at `n = 18`, so a slow divergence is not excluded); `E_units[rho ln rho] = deficit - ln(3/2) = 0.767` at `n = 18`, about `0.84` if the deficit converges to about `1.25`".

**C19 (§5, lines 366–367).** "The `L^1` distance between consecutive Haar densities decays like `n^(-0.93)`." -> "The `L^1` distance between consecutive Haar densities decays roughly like `1/n` (fitted exponent `-0.93` over all levels, `-1.04` over `10..18`, `-1.16` over `14..18`)."

**C20 (§5, line 363, "tail").** "its density has a tail `P(rho > t) ~ t^(-2)`" -> "its density is consistent with a tail `P(rho > t) ~ t^(-2)` (heuristic; local exponents `1.3–2.3` at level 14)".

**C21 (§5, lines 340–345, second tier).** "the set of classes `3^(j+1)/2^(j+2) - 1 mod 3^n`, `j < n`" -> "`j <= n - 2`"; "(`n = 8`: `-1: 24.9`, then five classes at `18.7`, all of the form ...)" -> "(`n = 8`: `-1: 24.9`, then the seven shadow classes `j = 0..6` at `18.6–18.7`, of which the top-6 list shows five)".

**C22 (§2.4, lines 192–193).** "so the number of central histories landing in a class is bounded by a polynomial in the generation" -> "so the number of central histories landing in a class is bounded by `(n+1)^3 G^n`, `G = 2013/2000` ((3.10)), a polynomial times a slowly growing exponential".

**C23 (§2.6, lines 199–201).** "then **coverage (Proposition 6.3):** the source charge `x omega(w) <= M`" -> "then the **source charge (Lemma 6.1)** `x omega(w) <= M`, **coverage of all large scales (Lemma 6.2)** and the **count (Proposition 6.3)**".

**C24 (§0.1, lines 55–57).** "his weighted inverse histories from a seed `M` with weights `omega(w) = 3^d 2^(-A(w))` are exactly the harmonic mass of the Syracuse tree of `M`" -> "his weighted inverse histories from a seed `M`, weights `omega(w) = 3^d 2^(-A(w))`, are summands of the harmonic mass of the Syracuse tree of `M` (he sums them over his restricted family of central first-crossing histories)".

**C25 (header, lines 31–32).** "FINITE-EXACT (the 3-adic Syracuse law to level 18, rational to level 6; ...)" -> "FINITE-EXACT (the 3-adic Syracuse law, rational, to level 6; ...) + VERIFIED (the law in float64 to level 18, agreeing with the rational law to `10^(-12)` where both exist; ...)".

**C26 (§6, line 466, "Two routes ... agree to `10^(-3)`").** -> "agree to `2.5·10^(-3)` (depth-20 run) and `4·10^(-3)` (depth-30 run), the second route carrying one more pruned level".

**C27 (§6, line 432).** "the leaf share `f^{(0)}` sits at `0.34–0.35` from depth 6 on" -> "`0.33–0.35` from depth 6 on"; "from depth 22 on all three shares are `1/3 ± 0.02`" -> "`1/3 ± 0.022`".

**C28 (§7(iii), lines 543–545).** "the S18 observation that `80%` of visited values above the start range are `2 mod 3` is the class mean `4/3` of Corollary A1 seen from the orbit side" -> "the class mean `4/3` of Corollary A1 is `P(Y_n = 2 mod 3) = 2/3` seen from the orbit side (two thirds of the arriving valuations are odd); the S18 figure of `80%` above the start range is the conditional enrichment of `a = 1` arrivals, not the class mean".

**C29 (§7 "For the repo" (i), line 529).** "The Krasikov–Lagarias `x^0.84` row (predecessors of `1`) would be replaced by positive density" -> "The Krasikov–Lagarias `x^0.84` row (predecessors of `1`), already superseded by Mazur's certified `pi_a(X) >= X^0.90` (his ref. [3], July 2026, cited in this paper's introduction), would be replaced by positive density".

**C30 (§7, lines 498–500, absolute continuity paragraph).** "Its density is unbounded (Theorem B), infinite on the forward rational closure of every negative cycle, not in `L^2` (the second moment grows linearly), with finite entropy `int f log f`." -> "Its local averages are unbounded (Theorem B) and tend to infinity on the forward rational closure of every negative cycle (PROVED); it is not in `L^2` if the second moment keeps growing (OBSERVED through level 18; Jensen), and has finite entropy if the deficit converges (OBSERVED)."

**C31 (§3, P6 row, line 222).** "| P6 | endpoint spread (3.9) ... | enumeration | holds |" -> "| P6 | endpoint spread (3.9) ... | enumeration (vacuous as run: 5 histories in 5 `(D, A)` classes, no two endpoints ever compared; both statements are trivial) | holds |".

**C32 (atlas §8 typing).** "UNIFORM **B** (presumably runs unchanged for `3x+k`, `3 ∤ k`; not checked)" -> "UNIFORM **B, CONJECTURAL** (the mechanism should transfer to `3x+k` with seeds `(4^j - k)/3` or `(2·4^j - k)/3` and GGM's mixing; the explicit constants and seeds are map-specific; not checked)".

---

## (iii) Verdict

**SOUND WITH CORRECTIONS.**

Theorems A, B (lower-bound form) and C are correct as stated and their proofs hold step by step
(items 1–14, 17–20, 35–43); the spike-profile recursion and path sums are right and the numerics
reproduce to two–four decimals at levels 8–14 (items 28–33); every number tested (`rho_n(1)` to
level 10 and 14, `8/7`, `1376/1387`, class means, second moments, entropy deficits, medians,
`1.2531` at `993`, `54.22%`, Kaprekar cycles, Euler bricks) reproduces (items 44, 62–66, 91, 107);
the digest of Theorem 1.1, Corollary 1.2, Lemma 2.1, (2.3), Lemma 4.1/Prop. 4.2, Lemma 6.1 and
Section 9 is faithful (items 73–77, 81–85); the absolute-continuity derivation is valid conditional
on (2.3) (item 89); the typing is defensible under the atlas rule (items 92–96). Nothing in the note
asserts that Mazur's theorem is false or that `H_n -> 0`, and the sequence data are honestly labelled
OPEN.

Substantive corrections, in order of importance:

1. **The paper's constants are understated by an exponential level** (item 80; corrections C9–C12): the session's `part8` set `D* = D_exp`, but `D* = D_sc` contains `P* = g^(R*-1)(T*)`, a cubic map composed `2^8696` times. `C` is `2^2^2^8697`, not `2^2^6536`; `log_2 log_2 N ≈ 8697`, not `log_2 N ≈ 6563`; `c^(-1)` is a four-fold tower and `X_0` a five-fold one. The qualitative conclusion (no numerical range) is unchanged; the numbers in §0.4, §2.7, §3, §5, the atlas §8 and the ledger are wrong.
2. **The minimum of the seed-1 sequence is at `n = 27` (`0.363`), not `n = 22` (`0.370`)** (item 56; C1, C2, C5–C7): four levels (`26–29`) lie below `H_22` in the session's own output; the monotone rise starts at `n = 28`. Title, §0, §6, the status table, the atlas, the ledger, the synthesis and two commit messages carry the wrong location.
3. **"Loss at most the pruned weight" is a calibration, not a bound** (item 55; C3, C4, C7, C8): the pruned subtrees' masses are unbounded a priori (`H_18` reaches `2·10^3` on resonant residues); the deficits equal the pruned weights empirically at `m + d = 18`. The computed values remain rigorous lower bounds, so "no decay through `n = 48`" stands; the "FINITE-EXACT evaluation" label does not.
4. **Author and version** (item 75; C13): Lech Mazur, version 2.1 of September 6, 2026, not "M. Mazur, v2".
5. **Asymptotics stated for what are lower bounds or finite-level values** (items 22–24, 31, 67, 69, 70, 72; C14–C20, C30): "grows like", "singularities of order", "`-> 0.9748`", "is exact", "converges", "`-> 0.84`", "`n^(-0.93)`", "`t^(-2)`", "finite entropy". Each is OBSERVED or a proved one-sided bound; the phrasing in §0 and parts of §5 and §7 presents them as established.
6. **Smaller misquotations and slips** (items 21, 32, 50, 58, 59, 78, 84, 87, 97; C17, C21–C23, C26–C29): `s = 0.0085`; the second tier is `j <= n - 2` with seven classes at `n = 8`; "3% per level" is 5%; the two routes agree to `2.5–4·10^(-3)`; (3.9) is `(n+1)^3 G^n`, not a polynomial; the source charge is Lemma 6.1; Mazur's `x^0.90` already supersedes the Krasikov–Lagarias row; the `80%` figure is not the class mean.
7. **Labels** (items 90, 95, 98; C25, C31, C32): the level-18 law is float (VERIFIED against exact to level 6), not FINITE-EXACT; P6 is vacuous as run; UNIFORM B is conjectural.

Not checked (outside reach): the platform's Lean statement, file counts, toolchain and the companion
theorem (item 86); the literal text of Tao's Proposition 1.14 (item 88, consistent with recollection);
the `F = 2467 b 16^b (C+1)` reading (item 79, immaterial); levels 15–18 of the law (only the
session's FFT; its internal assertions and the parents-form identity hold there).
