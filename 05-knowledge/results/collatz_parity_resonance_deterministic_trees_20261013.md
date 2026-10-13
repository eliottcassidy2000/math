# Twisted transfer law: the step-count characters of Collatz orbits hear the near-returns of log2 3 by congruence class, and the deterministic trees are silent

2026-10-13, opus-2026-10-13-S26 (worktree `codex/session-papers-collatz-hn-jc-bsd`, continuing S25). Owner prompt: *think deterministic trees, square triangular numbers and prior work adjacent to the golden ratio and go down many tangents of seemingly unrelated work and merge ideas into your continued exploration towards proofs; work next steps as they emerge; prioritize critical mathematics.*

**Status.** PROVED: the exact class-mean identities of the inverse tree on `Z_2 x Z_3`; the closed-form characteristic polynomial of the transfer operator on every space of locally constant functions of `(n mod 2^j, n mod 3^k)`, hence the spectral silence of every finite deterministic truncation; the resonance classification of the scalar law for every character `exp(2 pi i (alpha tau + beta omega)/m)` with its small-defect root laws; the Cesaro-shadow formula (all in THM-4629). FINITE-EXACT: `tau`, `omega` to `2^30`; sympy verification of the closed form in 30 cases; near-root scans at levels 9, 27, 81. NUMERICAL: every amplitude, frequency, decay and phase below. CITED: HYP-9230 (mac-mini), the Codex note of the same day (precedence for the `(-1)^tau` case and the qualitative mean-zero theorem), THM-4627, THM-4581, THM-4528, THM-3335. OPEN: the Haar hypothesis at infinite depth; the amplitudes. Audit: section 11.

**Precedence.** While this session ran, Codex pushed [collatz_signed_resonance_pair_mass_20261013b.md](collatz_signed_resonance_pair_mass_20261013b.md) (commit `ad10738e20`): the signed transfer law for `(-1)^tau`, HYP-9230's root equation on the half-integer lattice, the root table `2.5, 14.5, 26.5, 79.5, 173.5, 332.5`, the rule `Q odd, P even`, a census to `2^34` at 1024 and 4096 bins per octave, and a PROVED qualitative theorem (the ordinary mean of the ROOT-supported parity is zero, from THM-4581's almost-sure coalescence). That case is theirs by first push; this note cites it and builds on it. New here: the law for every character of `(tau, omega)` with its congruence classification (the odd-step count and the Collatz-map count hear the two other parity classes of near-returns; the `(Z/3)^2` characters refine by residue); the exact `Z_2 x Z_3` hierarchy and the theorem that it is spectrally silent at every finite depth, with an explicit characteristic polynomial; the small-defect expansion of the root; 29 lines of three parities with fitted decays against the exact roots; the class-resolved test; the Cesaro shadow, which turns the lines into the ordinary mean and predicts the Codex census value at `2^34` from data at `2^30`.

## 1. Setting and the exact hierarchy

`T(n) = n/2` or `(3n+1)/2`; `tau(n)` the `T`-steps to `1`; `omega(n)` the odd steps among them; `tau + omega` the steps of the Collatz map `C`. `theta = log2 3 = 1.5849625`, `u = log2 n`. A character `chi(n) = zeta^(alpha tau(n) + beta omega(n))`, `zeta = e^(2 pi i/m)`, satisfies `chi(n) = zeta^alpha chi(Tn)` for even `n` and `chi(n) = zeta^(alpha+beta) chi(Tn)` for odd `n`. The three parities are `(alpha, beta, m) = (1,0,2)`, `(0,1,2)`, `(1,1,2)`; the trivial character is HYP-9230's class density.

**Exact identities (THM-4629, statement 1).** Halving maps the integers `n = r (mod 2^j)`, `n = rho (mod 3^k)` of a log-bin at `u` onto the class `(r/2 mod 2^(j-1), rho/2 mod 3^k)` of the bin at `u - 1`; the odd step maps them onto the class `((3r+1)/2 mod 2^(j-1), c~(rho) mod 3^(k+1))` of the bin at `u + theta - 1`, where `c~(rho) = (3 rho' + 1)/2` with `rho'` the odd representative of `rho` modulo `2 3^k`. So the class mean at `(j,k)` is expressed through class means at `(j-1, k)` and `(j-1, k+1)`: every step consumes a binary digit, every odd step produces a ternary digit, and the exact hierarchy closes at no finite level. (This is the classical fact that `T` is 2-adically expanding and 3-adically contracting, in the form relevant to means.) The hypothesis that closes it at `(j,k)`, that the mean function depends on at most `j` binary and `k` ternary digits of `n`, is what every finite deterministic truncation of the inverse tree assumes; call it `H_(j,k)`. `H_(0,0)` is the scalar law
`m(u) = (1/2) zeta^alpha m(u-1) + (1/2) zeta^(alpha+beta) m(u + theta - 1)`.

## 2. The scalar law: resonances by congruence class and the root

A mode `e^(su)`, `s = 2 pi i f + w`, solves the scalar law iff `zeta^alpha e^(-s) + zeta^(alpha+beta) e^(s(theta-1)) = 2`. At `w = 0` both terms are unimodular, so a root near the axis needs both near `1`: `f = alpha/m (mod 1)` and `f theta = -beta/m (mod 1)`. Writing `f = Q/m`: `Q = alpha (mod m)` and `Q theta = P + m theta_f` with `P = -beta (mod m)`, `theta_f` small. The equation at such `f` is `e^(-w) + e^(2 pi i theta_f) e^(w(theta-1)) = 2`, HYP-9230's root equation with `theta_f` in place of the defect `theta` of the unsigned law. The general rule:

> the character `(alpha, beta)` modulo `m` resonates at `f = Q/m` exactly for the near-returns `2^P ~ 3^Q` with `Q = alpha` and `P = -beta` modulo `m`.

For `m = 2`: `(-1)^tau` hears `(Q odd, P even)`; `(-1)^omega` hears `(even, odd)`; `(-1)^(tau+omega)` hears `(odd, odd)`; the trivial character hears `(even, even)`, which is every integer `f` with small `||f theta||` (HYP-9230). The near-returns with `|Q theta - P| < 0.06`, `Q <= 730`, sorted by class:

| class | near-returns `(Q, P)` | temperament reading |
|---|---|---|
| `(-1)^tau`: odd, even | (5,8) (29,46) (53,84) (159,252) (241,382) (347,550) (453,718) (559,886) (665,1054) | 53-TET (Mercator), 665-TET |
| `(-1)^omega`: even, odd | (12,19) (36,57) (70,111) (94,149) (118,187) (200,317) (224,355) (306,485) (388,615) (412,653) (494,783) (518,821) (730,1157) | 12-TET, 94-, 306-TET |
| `(-1)^(tau+omega)`: odd, odd | (7,11) (17,27) (41,65) (65,103) (147,233) (253,401) (359,569) | 41-TET, 7-TET |
| trivial: even, even | (24,38) (82,130) (106,168) (212,336) (424,672) (506,802) | HYP-9230's 12, 41, 53, 106, 212, 253 |

**The root.** Continued from `w(0) = 0`,
`w = 2 pi i theta_f/(2 - theta) - 4 pi^2 theta_f^2/(2 - theta)^3 + O(theta_f^3)`
(derivation in THM-4629, proof of 3). The linear term is the half-drift `1/(2 - theta) = 2.4094` of HYP-9258 (the Codex note made the same identification); the quadratic term is the decay: per octave `exp(-552.2 theta_f^2)`. Its range (mpmath against the exact root):

| `theta_f` | exact `Re w` | `-4 pi^2 theta_f^2/(2-theta)^3` | exact drift | `theta_f/(2-theta)` |
|---|---|---|---|---|
| 0.0005 | -0.000138 | -0.000138 | +0.001204 | +0.001205 |
| 0.001 | -0.000550 | -0.000552 | +0.002404 | +0.002409 |
| 0.002 | -0.002173 | -0.002209 | +0.004779 | +0.004819 |
| 0.005 | -0.012580 | -0.013805 | +0.011489 | +0.012047 |
| 0.01 | -0.041139 | -0.055220 | +0.020711 | +0.024094 |
| 0.02 | -0.107510 | -0.220881 | +0.033092 | +0.048188 |
| 0.05 | -0.274797 | -1.380504 | +0.049988 | +0.120471 |

The factor `2` at `theta_f = 0.02` is the "factor-2 discrepancy" that HYP-9230's audit A2 resolved by using the exact root; here it is the third-order term. Exact roots at the predicted lines (`collatz_parity_three_characters_20261013.out`, Part 1; the Codex table agrees to all printed digits for `tau`):

| `(-1)^tau` | `f_eff` | decay | `(-1)^omega` | `f_eff` | decay | `(-1)^(tau+omega)` | `f_eff` | decay |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.4485 | 0.4864 | 1 | 1.0567 | 0.6587 | 0.5 | 0.5390 | 0.4323 |
| 2.5 | 2.4549 | 0.8088 | 6 | 6.0203 | 0.9611 | 1.5 | 1.4420 | 0.5859 |
| 9.5 | 9.5520 | 0.7353 | 18 | 18.0404 | 0.8475 | 3.5 | 3.5491 | 0.7694 |
| 14.5 | 14.4689 | 0.9097 | 47 | 46.9850 | 0.9787 | 8.5 | 8.4606 | 0.8551 |
| 26.5 | 26.5036 | 0.9988 | 153 | 152.9982 | 0.9997 | 20.5 | 20.4822 | 0.9701 |
| 332.5 | 332.5001 | 1.0000 | | | | 73.5 | 73.4880 | 0.9863 |
| | | | | | | 179.5 | 179.5019 | 0.9997 |

The three parities pick out the three classical temperaments: the odd-step parity the 12-tone comma (`3^12/2^19`), the Collatz-map parity the 41-tone comma (`3^41/2^65`), the `T`-step parity the 53-tone comma of Mercator (`3^53/2^84`) and the 665-tone comma. The golden ratio enters through THM-4528: the parity words of the Collatz map `C` are golden-mean words, and the `C`-step count `tau + omega` is the third character, the one that hears the `(odd, odd)` commas.

## 3. Data: the three parities to 2^30

`tau` and `omega` were computed exactly to `2^30` (int16, 10 s with numba; maximum `tau` 616). Bin means of `(-1)^tau`, `(-1)^omega`, `(-1)^(tau+omega)` at 1/64 octave over `2^16..2^30` and at 1/768 octave over `2^20..2^30`; per-octave least-squares amplitude at a frequency (cosine, sine, constant); fitted decay = exponential of the slope of log amplitude against octave.

**Twenty-nine lines.** The table is THM-4629, statement 4 (`collatz_parity_highlines_decay_20261013.out`). Summary: for the 24 lines of amplitude above `0.02` the fitted decay per octave differs from the exact root by at most `0.006` (median `0.001`); `(-1)^omega` at `6.0203`: `0.9611` fitted against `0.9611` predicted; `(-1)^(tau+omega)` at `8.4606`: `0.8551` against `0.8551`; `(-1)^tau` at `26.5036`: `0.9997` against `0.9988`. The four lines of amplitude below `0.02` (`332.5`, `18`, `206`, `3.5`) are at the fit's noise floor; `206` (`(412,653)`, defect `0.0046`, predicted decay `0.9972`) is present at amplitude `0.008` in all ten octaves, a small initial amplitude, not an absent line.

**Seen directly.** `(-1)^tau` bin means at 1/212 octave over two Mercator periods above `2^29`: `+0.41, +0.29, +0.34, -0.02, -0.25, -0.24, -0.33, +0.03, +0.29, +0.24, +0.20, -0.05, -0.43, -0.20, -0.20, +0.07`. The parity of the total stopping time is not equidistributed on any scale computation reaches; it swings by `0.8` every `1/26.5` octave.

**The 1/64-octave scan** (`collatz_parity_spectrum_hires_20261013.out`, `collatz_parity_three_characters_20261013.out`) sees the same lines attenuated by `sinc(f/64)` (`0.262` at `26.5` for a true `0.354`) and, for `(-1)^omega`, the track `0.5456, 0.5060, 0.4676, 0.4325, 0.3994, 0.3691, 0.3410, 0.3278` at `2^16, ..., 2^29` against the prediction `0.5451, 0.5055, 0.4674, 0.4328, 0.3995, 0.3689, 0.3408, 0.3277` from the root `0.9611` and the amplitude at `2^16` (three decimals over thirteen octaves). The joint damped-sinusoid fits of S25's low mode give `f = 2.4550`, decay `0.823` (root: `2.4549`, `0.8088`), `14.470`, `0.905` (`14.4689`, `0.9097`) and `26.5025`, `1.000` (`26.5036`, `0.9988`). The side peaks at `28.0`, `25.0` (and `4.5`, `7.5` beside `6`, `45.5` beside `47`) are the `+-1.5`-cycle sidelobes of the one-octave rectangular window (`sinc(1.5) = 0.21`; `0.354 x 0.21 = 0.075` as observed). A transient peak at the Nyquist frequency `32` of the 64-per-octave analysis of `(-1)^omega` and `(-1)^(tau+omega)` (`0.22`, `0.30` at `2^16`, decaying to `0.09`, `0.13`) has no near-return of its class and is absent at 768 per octave: an alias, recorded, not interpreted.

## 4. Data: the characters of (Z/3)^2

`chi_(alpha,beta)(n) = exp(2 pi i (alpha tau + beta omega)/3)`, top octave, 1/768 octave, complex least squares (`collatz_parity_mod3_characters_20261013.out`; the script reports the larger of the amplitudes at `+f` and `-f`, so it conflates a character with its conjugate, and the rule is read for the pair: `Q = +-alpha`, `P = -+beta`):

| pair | peaks `f (Q,P) amplitude` |
|---|---|
| `omega mod 3` | 4.00 (12,19) 0.290; 151 (453,718) 0.090; 49 (147,233) 0.076; 8 (24,38) 0.056; 102 (306,485) 0.045 |
| `tau mod 3` | 17.75 (53,84) 0.121; 35.25 (106,168) 0.076; 5.5 (17,27) 0.069; 70.75 (212,336) 0.054; 164.75 (494,783) 0.048; 23.25 (70,111) 0.042; 141.25 (424,672) 0.041; 341.25 (1024,1623) 0.031 |
| `(1,1)/(2,2)` | 31.25 (94,149) 0.067; 21.75 (65,103) 0.066; 168.75 (506,802) 0.044; 243.25 (730,1157) 0.038; 74.75 (224,355) 0.037; 137.25 (412,653) 0.037; 9.75 (29,46) 0.034; 84.25 (253,401) 0.032 |
| `(1,2)/(2,1)` | 13.75 (41,65) 0.144; 27.25 (82,130) 0.109; 39.5 (118,187) 0.084; 80.25 (241,382) 0.063; 66.75 (200,317) 0.048; 119.75 (359,569) 0.033 |

All 26 peaks are near-returns in the predicted class (checked by `Q = 3f`, `P = round(Q theta)`, residues); the other entries of the top-eight lists are sidelobes (`2.5`, `5.5` of `4.00`; `25.75`, `12.25` of `27.25`, `13.75`) and the window's `0.5` trend. The grid of the scan is `0.25`, so `35.25` reads `35.33 + drift`, etc.

## 5. Deterministic trees are silent

**Theorem (THM-4629, statement 2).** On `W_(j,k)`, the locally constant functions of `(n mod 2^j, n mod 3^k)`, which are exact invariant subspaces of the transfer operator, the characteristic polynomial is
`det(I - M(s)) = [1 - (1/2) zeta^alpha e^(-s) - (1/2) zeta^(alpha+beta) e^(s(theta-1))] prod_(i=1)^(k) [1 - ((1/2) zeta^alpha e^(-s))^(2 3^(i-1))]`
for every `j`. The first factor is the scalar law; the others vanish only at `Re s = -ln 2`, i.e. modes decaying by half per octave, at the frequencies `alpha/m - l/(2 3^(i-1))`.

*Proof idea.* Filter by 3-adic depth: the odd branch is zero on every graded piece (its image class modulo `3^i` depends only on `rho mod 3^(i-1)`), leaving the even branch, which is `zeta^alpha e^(-s)` times (2-adic lumping) tensor (the permutation `rho -> rho/2` on the characters of exact order `3^i`). The 2-adic lumping has, by Sylvester's determinant identity applied level by level, the single nonzero eigenvalue `1/2` at every depth `j`; the permutation is one cycle of length `2 3^(i-1)` because `2` is a primitive root modulo `3^i`. At depth `0` the same Sylvester reduction collapses the 2-adic levels to the two classes even/odd, whose matrix has eigenvalues `a' + b'` and `0`. Full proof in THM-4629.

*Checks.* sympy at `(j,k) = (0,1), (1,0), (2,0), (3,0), (1,1), (2,1), (1,2), (3,1), (2,2), (0,2)` for all three parities: the closed form matches in all 30 cases (`collatz_parity_mixed_levels_20261013.out`). The explicit factorisations at `(0,1)`, `(0,2)`: `+-(a/2 - 1)(a/2 + 1)(scalar)` and `(a/2 - 1)^2 (a/2 + 1)^2 (scalar)(a^2/4 - a/2 + 1)(a^2/4 + a/2 + 1)` for `(-1)^tau`, with the signs permuted for the other two. Numeric near-root scans at levels `9, 27` (all three) and `81` (`tau`): every near-root off the scalar list sits at `Re s = -0.69 = -ln 2` (`collatz_parity_3adic_hierarchy_20261013.out`).

*Data.* The nine classes `n mod 9` on the top octave carry every line with the global amplitude: maximal relative spread `0.1%-0.7%` and phase spread `<= 0.006` rad for the 20 lines of amplitude `>= 0.03`; the `26.5036` line is `0.3536..0.3540` in all nine classes; the aggregates mod 3 agree to `0.0001` (`collatz_parity_classes_mod9_20261013.out`).

**Wall.** The theorem is a no-go for a natural proof strategy: refining the inverse tree to any finite 2-adic and 3-adic depth, the deterministic structure that Terras's periodicity and the 3-adic classes provide, cannot change the slow spectrum, so it cannot prove the modal expansion for the integers either; the expansion is exactly the Haar hypothesis at infinite depth (that the mean functions are locally constant in the 6-adic digits at the scales considered). The theorem says nothing about mean functions depending on infinitely many digits: that is where the open content lives, and the data (section 4 of the Codex note, section 5 here) say that at the computed scales the dependence is below `0.5%`.

## 6. The Cesaro shadow: the ordinary mean is the lines, reduced by the Cesaro factor

For a line `a(x) cos(2 pi f x + phi)` in `x = log2 n`, the prefix mean `(1/N) sum_(n<=N)` is `a(log2 N) cos(2 pi f log2 N + phi - arctan om)/sqrt(1 + om^2)` with `om = 2 pi f/ln 2`, up to the decay correction and `O(1/N)` (THM-4629, statement 6). For the Mercator line, `1/sqrt(1 + om^2) = 0.00416`: predicted Cesaro amplitude `0.354 x 0.00416 = 0.00147`.

`collatz_parity_cesaro_shadow_20261013.out`: the nine `(-1)^tau` lines fitted jointly on `[2^29, 2^30)` (residual rms `0.062` against a signal rms `0.296`; constant `+0.00014`) and extrapolated with the predicted decays give, for the exact prefix means at 256 points per octave over `2^24..2^30`: exact rms `0.001301`, maximum `0.00358`; model rms `0.001298`; rms of the difference `0.000177`. The exact prefix mean oscillates at `26.50` cycles per octave with amplitude `0.001484` (`14.5`: `0.00048`). At the powers of two the exact means alternate, `+0.00154, -0.00203, +0.00225, -0.00238, +0.00233, -0.00223, +0.00208` (`2^24..2^30`), model `+0.00161, -0.00174, +0.00227, -0.00209, +0.00238, -0.00202, +0.00218`. The model at `2^34`: `+0.00159` (envelope `0.00225`, range `-0.00197..+0.00226` over `2^33.9..2^34.1`); the Codex census gives the exact value `+0.001452` at `2^34` (prefix sum `24,953,018`). The extrapolation over five octaves from data to `2^30` lands within the model's own rms error.

Two readings. (i) The alternation of the dyadic-block means of `(-1)^tau` found in S25 (HYP-9258 (iii), "forced by `tau(2n) = tau(n) + 1`") is the half-integer lattice: every line of `(-1)^tau` has half-integer `f`, so every one changes sign from `2^k` to `2^(k+1)`; the lines of `(-1)^omega` have integer `f` and do not alternate (`omega(2n) = omega(n)`). (ii) The Codex theorem that the ordinary mean is `o(1)` is approached, in the model, at the rate of the Mercator line, `0.9988^(log2 N) = N^(-0.00179)`: the prefix mean falls below `0.001` only beyond `2^350`. A limit theorem about the parity of `tau` is invisible to every census.

## 7. Tangents merged and typed

- **Exponential Pell against Pell (square-triangular numbers).** The resonances are the near-solutions of `X^2 - 3 Y^2 ~ 0` with `X = 2^(P/2)`, `Y = 3^((Q-1)/2)`, an exponential Diophantine approximation with no exact solutions (Mihailescu), and each near-solution gives a damped line with decay quadratic in the defect (section 2). The true Pell equation `x^2 - 2y^2 = 1` behind the square triangular numbers `1, 36, 1225, 41616, ...` (THM-3335; section 4 of the Codex companion note) has exact solutions with the exact log-period `log2 (3 + 2 sqrt 2)^2 = 5.0857` octaves and no damping. The convergents of `sqrt 2` and of `log2 3` share the denominators `5, 12` and the semiconvergent `29` (`(29,46) = (12,19) + (17,27)`); this is a coincidence of continued-fraction prefixes (`[1; 2, 2, 2, ...]` against `[1; 1, 1, 2, 2, 3, 1, 5, ...]`) and is typed NUMEROLOGY.
- **Golden control against Collatz (reading, ANALOGY).** The Codex companion proves universal entry into a finite box for a map `G` on `Z[phi]` whose two branches both contract in the real embedding (factors `1/2` and `phi/2`). In the language of this note, a one-sided walk in `log` size has no near-returns, so no resonances, and a contracting norm exists; Collatz's odd branch expands by `3/2`, the walk is two-sided with step ratio `log2(3/2)` irrational, the near-returns `2^P ~ 3^Q` exist, and the resonance lines are the spectral signature of the missing contracting norm. The golden ratio's own appearance in Collatz is THM-4528's: the `C`-parity words are golden-mean words, and the `C`-step parity is the character of the `(odd, odd)` commas.
- **Deterministic trees.** Terras's periodicity makes every 2-adic depth an exact identity, and the 3-adic classes are produced by the odd step; both are spectrally silent (section 5). The only non-silent deterministic input is the pair of step sizes `(1, log2(3/2))`.
- **Beatty words.** THM-4514's complementary Beatty sequences describe tight spine blocks of individual orbits; the resonance law is about ensemble means and does not use them. No link beyond the common irrational `theta` is claimed.
- **Baker-type bounds** (CITED, not used): lower bounds for `|Q theta - P|` bound the defect of every line from below and hence its decay from above; they are far too weak (`Q^(-13)` scale) to say anything at `Q <= 10^3`.

## 8. What is proved about the integers, and what is not

| statement | type |
|---|---|
| exact class identities on `Z_2 x Z_3`; the hierarchy closes at no finite level | PROVED (THM-4629, 1) |
| closed-form characteristic polynomial on every `W_(j,k)`; finite deterministic truncations have the scalar slow spectrum | PROVED (THM-4629, 2); FINITE-EXACT in 30 symbolic cases |
| resonance classification by `(Q mod m, P mod m)`; root = HYP-9230's with `theta_f`; `w = 2 pi i theta_f/(2-theta) - 4 pi^2 theta_f^2/(2-theta)^3 + ...` | PROVED (model) |
| Cesaro-shadow formula | PROVED (model) |
| 29 lines: frequencies to `0.1%`, decays to `0.005` per octave; 26 `(Z/3)^2` peaks; class silence mod 9; prefix means to `14%` rms; `2^34` prediction | NUMERICAL (data exact to `2^30`) |
| ordinary mean of the ROOT-supported `(-1)^tau` is zero | PROVED, Codex note (from THM-4581) |
| the modal expansion for the integers; the amplitudes | OPEN (Haar hypothesis at infinite depth) |

## 9. Next steps

1. **The amplitude problem.** The amplitudes are the projection of the small-`n` data onto the bounded modes of an advance-delay equation (THM-4629, scope). Target: a spectral projection formula from the exact means at `2^8..2^14` (where the law is already accurate to the bin resolution) to the amplitudes at `2^29`; a correct prediction of `0.354` for the Mercator line would make the amplitudes a theorem of the model.
2. **A rate for the signed residue error.** The Codex note's identity (7) bounds the error of the scalar law by the surviving pair mass; the error is a sum over the class means modulo `3` at the larger scale, and section 5 says those are the global mean plus modes decaying by half per octave in the model. The open step is the same Haar hypothesis.
3. **`m >= 4` characters and the full near-return census.** The rule predicts every line of every character; a census of `(Z/m)^2` for `m = 4, 5, 6` would test the rule at residues where no temperament reading exists.
4. **Golden analogue.** Compute the parity spectra of THM-4528's base-`phi` map and of the Codex control map `G`: the law predicts no lines for `G` (one-sided walk) and lines at the near-returns of `log_phi 2` for the base-`phi` Collatz.

## 10. Files

Scripts and outputs in `04-computation/experiments/`: `collatz_parity_roots_20261013.py`, `collatz_parity_mod3_roots_20261013.py`, `collatz_parity_spectrum_hires_20261013.py`, `collatz_parity_three_characters_20261013.py`, `collatz_parity_3adic_hierarchy_20261013.py`, `collatz_parity_mod3_characters_20261013.py`, `collatz_parity_highlines_decay_20261013.py`, `collatz_parity_classes_mod9_20261013.py`, `collatz_parity_cesaro_shadow_20261013.py`, `collatz_parity_mixed_levels_20261013.py` (each with `.out`; bin arrays as `.npy`). Run with `python3 -P <script>`; `tau`, `omega` to `2^30` need 4 GB.

## 11. Audit

(filled after the independent adversarial audit)
