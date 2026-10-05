# Five Fourier experiments on the Syracuse law (2026-10-04): the cold-frequency rate is universal across fixed units AND across multipliers `q` (it is the incoherent rate `1/sqrt3` of the 2-adic recursion, not the Parseval rate `q^(-1/2)`), so H1 is a statement about the binary digits of `3^(-n)`; the Fourier mass is uniform in the real frequency coordinate (the "hot dyadic frequencies" reading is refuted); the density has a `t^(-2)` tail and `E[rho^2]` grows `0.31` per level (P7 supported), which is exactly THM-4263's uniform-integrability condition for density transport

**Session:** opus, `collatz-synthesis-20261004` (opus-2026-10-04-S1, cycle 2 of the
owner's experiment loop: "test the most relevant hypotheses repeatedly, inching
toward a proof; when stuck, search niche past concepts; expand a web of
connections; reframe"), 2026-10-04.
**Inherits (cited):** the coalescence note
([`collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md`](collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md),
predictions P6, P7), [HYP-9166](../hypotheses/HYP-9166-h1-frequency-one-coefficient-decays-below-critical.md)
(H1, mac-mini; the closed window recursion `f_n(k) = sum_a 2^-a omega_n(k-a) f_(n-1)(k-a)` and its rescaled form,
`collatz_h1_20260929_mu1_rescaled.py`, whose phase device `(2^-M mod 3^n)/3^n = frac(R/2^M) + 2^-M 3^-n`,
`R = -3^-n mod 2^(M+1)`, is reused with a general unit and a general multiplier),
[THM-4519](../../01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md),
[THM-4520](../../01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md),
the S19 Mazur digest ([`mazur_positive_density_20260928.md`](mazur_positive_density_20260928.md): the reference
density `rho_n = (2/3) 3^n mu_n`, its atoms and `E[rho_n^2]`), the S20 five-mirrors note (the resonant maxima at
`+-2^s`), the S23 three-mirrors note (the full-period recursion), and
[THM-4263](../../01-canon/theorems/THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport.md)
(finite-factor density transport iff uniformly integrable fibre weights, condition (14)).

**Status: FINITE-EXACT / VERIFIED (every number below is from the scripts named in section 6, each with a
direct-enumeration control at small levels); PROVED (two identities: the phase sequence of the recursion is the
2-adic digit string of `-u q^-n`; the Cauchy-Schwarz form of the incoherent rate); OBSERVED (the rates);
REFUTED (the real-coordinate concentration hypothesis of E3); SUPPORTED (P7, uniform integrability);
DIRECTION (the digit reframing of H1). H1 remains OPEN; nothing here is a Collatz step.**

---

## 0. What is new, in one screen

1. **E1 (P6, VERIFIED to level 800).** For the fixed units `u = 1, 5, 7, 11, 13, 17` the coefficients
   `|mu_hat_n(u)|` of the 3-adic Syracuse law decay at one rate: least-squares rates over `n = 200..800` are
   `0.5680, 0.5695, 0.5687, 0.5709, 0.5703, 0.5704` (block rates `0.553-0.580`), all at or below the Parseval
   rate `3^(-1/2) = 0.5774` and far below the critical `0.5850` of H1. The ridge events are frequency-specific
   (`u = 1`: the known wave at `n = 127..135` reaching `10.4` Parseval units; `u = 13`: `6.9` at `n = 34..36`;
   `u = 11`: `2.6` at `n = 24`). Truncation control `A = 40` vs `60`: relative difference `< 3 10^-9` to `n = 400`.
2. **E4 (the drift control; the main finding).** The same recursion for the `q`-adic law of `qx+1` gives, for
   `q = 3, 5, 7, 11` and units `u = 1, 5` (`u = 2` for `q = 5`), rates over `n = 200..600` of `0.5721, 0.5698;
   0.5735, 0.5738; 0.5716, 0.5718; 0.5720, 0.5734` -- **the same `0.57` for every multiplier**, while the
   Parseval rates are `0.5774, 0.4472, 0.3780, 0.3015`. So the cold-frequency rate is not `q`'s Parseval scale;
   it is the incoherent rate of the 2-adic recursion, `(sum_a 4^-a)^(1/2) = 1/sqrt3 = 0.57735`, which for `q = 3`
   happens to coincide with the Parseval rate. For `q = 5, 7, 11` the fixed-frequency coefficients are
   astronomically ABOVE the Parseval scale (`q^(n/2)|mu_hat_n(1)| = 10^58, 10^99, 10^152` at `n = 600`): the
   `q`-adic law of a drifting map equidistributes at fixed frequencies only at the 2-adic rate `0.57`.
3. **The identity behind it (PROVED, one line).** `omega_n(-a) = e((u 2^-a mod q^n)/q^n) = e(0.b_a b_(a-1) ... b_1
   + u 2^-a q^-n)`, where `b_i` is the `i`-th binary digit of the 2-adic number `-u q^-n`. So `mu_hat_n(u) =
   f_n(0)` is a functional of the binary digit strings of `q^-1, q^-2, ..., q^-n`, consecutive strings being
   related by the 2-adic multiplication by `q^-1` (for `q = 3`: the carry automaton of `x -> 3x` run backwards).
   **H1 is therefore a statement about the binary digits of the powers of `3^-1` in `Z_2`** -- the 2-versus-3
   transversality that the thread has repeatedly identified as the missing mechanism, in Fourier form; and the
   `q`-universality of the rate says the mechanism uses nothing about `3` beyond oddness.
4. **E5 (the random-digit control; section 4).** Replacing the digit strings by i.i.d. bits, or by the digits of
   a random odd multiplier, gives the rates recorded in section 4: this decides whether mac-mini's measured deficit
   (`0.569` against `0.5774`, i.e. `3^(n/2)|mu_hat_n(1)|` falling from `9 10^-2` to `3 10^-15` over `n =
   200..2500`) is arithmetic or a property of the recursion.
5. **E3 (REFUTED hypothesis).** The Fourier mass at level `n <= 13`, in the real coordinate `theta = t/3^n`, is
   uniform to `+-5%` in 64 bins; the mass within a quarter-spacing of the dyadic rationals `j/2^m` equals its
   Haar expectation to `1%` for `m = 3..8`; the fixed units `|t| <= 100, 1000, 10000` carry their Haar share
   (ratios `0.94-1.46`). The resonant maxima `+-2^s` are isolated spikes (`|mu_hat| = 0.022-0.026` at
   `n = 12, 13`, `30x` the rms) carrying `0.1%` of the mass each. So "the law has structure at the real scale
   `1/64`" is false; the powers of two are special as 2-adically simple frequencies, not as real positions.
6. **E2 (P7 supported; THM-4263's condition named).** `E[rho_n^2]` grows by `0.3077 ... 0.3135` per level
   (`n = 2..14`, slowly rising), `E[rho_n^2.5]` grows geometrically (`x1.17` per level), `E[rho_n^1.5]` has
   increments decaying like `n^-1.1` (`0.0282` at `n = 14`), `E[rho_n^1.9]` increments decay like `0.98^n`; the
   tails `P(rho_n > t)` are stable in `n` at fixed `t` with local exponents `1.85, 2.05, 2.35, 2.89` at
   `t = 4, 8, 16, 32` (level 14), i.e. `P(rho > t) ~ c t^-2` with `c ~ 0.5` in the bulk. The maximal atom is
   the class of `-1` at every level (`284.5 = 0.9748 (3/2)^14`). THM-4263's uniform-integrability quantity
   `E[rho_n 1(rho_n > M)]` is `0.284, 0.145, 0.065, 0.024, 0.006, 0.001` at `M = 4, ..., 128` (level 14), still
   rising in `n` at fixed `M` (the law has not converged at these `M` by level 14) but of the size `2c/M` the
   `t^-2` tail predicts. **Reading:** THM-4263 says density-one statements transport from uniform 3-adic targets
   to uniform 2-adic sources iff (14) holds; (14) is the `L^1`-uniform tail of the Syracuse density, which the
   `t^-2` law satisfies with room (`t^-1` would not). The Mazur digest's absolute continuity (conditional on
   (2.3)) is the stronger property; (14) is the one that transport needs.

---

## 1. E1: the cold-frequency rate across fixed units (`q = 3`)

Script `collatz_fixed_frequency_rates_20261004.py` (`N = 800`, `A = 40`; output `.out`). The recursion
`f_n(k) = sum_(a=1)^A 2^-a omega^(u)_n(k-a) f_(n-1)(k-a)` on the window `[-(N-n)A, 0]`, rescaled by `3^(n/2)`,
with `omega^(u)_n(j) = e((u 2^j mod 3^n)/3^n)`; `mu_hat_n(u) = f_n(0)`. Control: at `n = 3, 5, 7` the recursion
agrees with a direct enumeration of the law (valuations `<= 30`) to ten digits for `u = 1, 5, 7`.

| `u` | 200-block medians of `3^(n/2)|mu_hat_n(u)|` (`n = 1..800`) | rate `200..400` | `400..600` | `600..800` | `200..800` | local maxima `>= 1` (`n`, value) |
|---|---|---|---|---|---|---|
| 1 | `1.5e-1, 7.4e-3, 9.5e-4, 4.4e-6` | 0.5741 | 0.5802 | 0.5534 | **0.5680** | (127, 3.88), (129, 7.22), (131, 10.38), (135, 2.62) |
| 5 | `1.0e-1, 3.5e-3, 3.5e-4, 2.0e-5` | 0.5629 | 0.5698 | 0.5663 | **0.5695** | (68, 2.69) |
| 7 | `4.4e-2, 5.3e-3, 1.9e-4, 9.4e-6` | 0.5681 | 0.5795 | 0.5754 | **0.5687** | (11, 1.80), (49, 1.27) |
| 11 | `1.4e-1, 2.6e-3, 1.3e-4, 3.5e-5` | 0.5652 | 0.5685 | 0.5567 | **0.5709** | (24, 2.57), (32, 2.06), (180, 1.57) |
| 13 | `1.0e-1, 4.6e-3, 1.6e-4, 4.0e-5` | 0.5666 | 0.5748 | 0.5653 | **0.5703** | (34, 6.89), (36, 6.98), (72, 1.62) |
| 17 | `1.2e-1, 3.7e-3, 3.0e-4, 3.1e-5` | 0.5627 | 0.5761 | 0.5693 | **0.5704** | (27, 2.08), (92, 1.28) |

The `u = 1` row reproduces mac-mini's `n = 127..135` wave exactly (`3.882, 7.221, 10.377`). The medians of the
ratios `|mu_hat_n(u)|/|mu_hat_n(1)|` over 200-blocks are `O(1)` (`0.17` to `8.8`): the same exponential rate with
independent fluctuations. **P6 VERIFIED to level 800 for six units**; the universal constant itself remains
OBSERVED (`0.569 +- 0.003`).

## 2. E4: the drift control (`q = 3, 5, 7, 11`)

Script `collatz_fixed_frequency_drift_control_20261004.py` (`N = 600`, `A = 40`). The `q`-adic law of `qx+1`:
`Y_n = 2^-a (q Y_(n-1) + 1) mod q^n`, the same recursion with `omega^(u,q)_n(j) = e((u 2^j mod q^n)/q^n)`,
rescaled by `q^(n/2)`. Control: direct enumeration at `n = 4, 5` for `q = 3, 5, 7`, ten digits.

| `q` | Parseval `q^(-1/2)` | `u` | rate `100..300` | `300..600` | `200..600` | ratio to Parseval | `q^(n/2)|mu_hat_600(u)|` |
|---|---|---|---|---|---|---|---|
| 3 | 0.5774 | 1 | 0.5656 | 0.5720 | 0.5721 | 0.991 | `~10^-3` |
| 3 | 0.5774 | 5 | 0.5677 | 0.5713 | 0.5698 | 0.987 | `~10^-4` |
| 5 | 0.4472 | 1 | 0.5768 | 0.5765 | 0.5735 | 1.282 | `2 10^58` |
| 5 | 0.4472 | 2 | 0.5764 | 0.5765 | 0.5738 | 1.283 | `2 10^58` |
| 7 | 0.3780 | 1 | 0.5715 | 0.5743 | 0.5716 | 1.512 | `3 10^99` |
| 7 | 0.3780 | 5 | 0.5712 | 0.5665 | 0.5718 | 1.513 | `8 10^97` |
| 11 | 0.3015 | 1 | 0.5685 | 0.5739 | 0.5720 | 1.897 | `5 10^152` |
| 11 | 0.3015 | 5 | 0.5694 | 0.5772 | 0.5734 | 1.902 | `3 10^152` |

**Reading.** The rate is a property of the 2-adic side (the geometric weights `2^-a` and the phases, which are
the binary digits of `-u q^-n`), not of `q`. Under independent random phases the recursion gives
`E|f_n(0)|^2 = sum_a 4^-a E|f_(n-1)(-a)|^2`, i.e. the rate `(1/3)^(1/2) = 0.57735` -- the Cauchy-Schwarz bound
`|f_n(0)| <= 3^(-1/2) (sum_a |f_(n-1)(-a)|^2)^(1/2)` is the incoherent model. For `q = 3` this incoherent rate
equals the Parseval rate `3^(-1/2)` of the law, which is why H1 (`rho < 0.585`) has only the `1.3%` margin and why
the question looked `3`-specific. It is not: the same `0.57` governs `5x+1`, whose law is nowhere near
Parseval-flat at fixed frequencies.

**Consequence for H1.** H1 asks `|mu_hat_n(1)| <= C rho^n` with `rho < log_2 3 - 1 = 0.585`. The data (here to
`800`, mac-mini to `2500`) say the true rate is `1/sqrt3` or slightly below, for every `q`. A proof of
"rate `<= 1/sqrt3 + eps`" from the structure of the recursion alone would prove H1 with `eps < 0.0076`. What such
a proof needs is a decorrelation statement for the phase sequence `a -> e(0.b_a ... b_1)` built from the binary
digits of `3^-n`, uniformly in `n` -- a 2-versus-3 digit statement. This is the Fourier form of the
transversality mechanism named by the procgen synthesis (HYP-9127's cube-swap value, THM-4469's Mahler bridge)
and by the S13 directions note; it is now one explicit sequence.

## 3. E3: where the Fourier mass sits in the real coordinate (REFUTED hypothesis)

Script `collatz_fourier_mass_real_coordinate_20261004.py` (`n <= 13`, full FFT of the dense law). The
hypothesis tested: the Parseval mass at level `n` is carried by frequencies `t` whose real position
`theta = t/3^n` is near a dyadic rational of small denominator (the resonant maxima sit at `2^s/3^n ~ 2^-6`).
Result: the 64-bin mass distribution over `theta` is flat (top bins `0.0163-0.0169` against uniform `0.0156`);
the mass within a quarter spacing of `{j/2^m}` is `0.496-0.506` for `m = 3..8` against the uniform `0.500`; the
fixed units `|t| <= 10^2, 10^3, 10^4` carry `0.94-1.46` of their Haar share at `n = 11..13`. The top coefficients
at `n = 13` are `+-2^16 (0.0221), +-2^15 (0.0218), +-2^17 (0.0188), +-2^14 (0.0185)` (S20's resonant family, with
`theta = 0.041, 0.021, 0.082, 0.010`) and the lifts of the level-12 maxima (`t = 465905 = -2^16 mod 3^12`, a
different unit at level 13, `0.0157`): isolated spikes of `0.1%` of the mass each. **REFUTED:** the law has no
structure at the real scale `1/64`; the powers of two are special as 2-adically simple frequencies (the closed
family of the recursion), not as real positions. The cold fixed units of E1 are at their Haar share at level 13
(`(0.57/0.5774)^13 = 0.84`) and only fall below it at large `n` (`0.1` by `n = 180`).

## 4. E5: the random-digit control

Script `collatz_fixed_frequency_random_digits_20261004.py` (`N = 1500`, `A = 40`, three seeds). Modes: the real
digits of `-3^-n` and `-5^-n`; i.i.d. random bits at every level; the digits of `-Q^-n` for a random odd 40-bit
`Q`. Values are `3^(n/2)|f_n(0)|` (rescaled by the incoherent rate).

**Results (`N = 1500`, `A = 40`; least-squares rates and their ratio to `1/sqrt3`).**

| digit source | rate `300..1500` | ratio | rate `750..1500` | ratio |
|---|---|---|---|---|
| real `-3^-n` (the Collatz case) | **0.5689** | 0.9854 | **0.5695** | 0.9864 |
| real `-5^-n` | 0.5736 | 0.9934 | 0.5745 | 0.9950 |
| i.i.d. bits, seed 0 / 1 / 2 | 0.5741 / 0.5720 / 0.5723 | 0.994 / 0.991 / 0.991 | 0.5726 / 0.5733 / 0.5731 | 0.992 / 0.993 / 0.993 |
| random odd 40-bit multiplier `Q`, seed 0 / 1 / 2 | 0.5736 / 0.5721 / 0.5727 | 0.993 / 0.991 / 0.992 | 0.5755 / 0.5739 / 0.5725 | 0.997 / 0.994 / 0.992 |

**Reading.** (i) Even with i.i.d. digits the recursion decays at `0.573 +- 0.001`, `0.8%` below the
incoherent rate `1/sqrt3`: the window entries `f_(n-1)(-a)` are correlated (they share their past), so the
Cauchy-Schwarz/incoherent model is an upper heuristic, not the rate. (ii) The real digits of `3^-n` give a
further `0.6%` deficit (`0.569`, mac-mini's value to `n = 2500`), below every random model in both windows,
while the real digits of `5^-n` sit at the random level. So the Collatz case shows an additional arithmetic
cancellation of about `0.6%` per level, specific to `q = 3` among the cases run (OBSERVED; three seeds per
model; the wider control of section 4a settles the significance). (iii) Whatever its origin, for H1 only the
inequality matters: every model, random or arithmetic, sits at `0.569-0.576 < 0.585`; H1's margin against the
critical rate is `1.6-2.8%`, and what a proof must show is that the Collatz digits are no worse than random
in this one functional.

### 4b. The mean-square rate is exactly incoherent in the i.i.d. model (PROVED); the rest is a Jensen gap

**Proposition (PROVED, one line).** In the i.i.d.-digit model (fresh uniform bits at every level) the window
phases `omega(-d) = e(theta_d)`, `theta_d = 0.b_d b_(d-1) ... b_1`, are pairwise uncorrelated: for `d < e` the top
bit `b_e` enters `theta_e - theta_d` with coefficient `1/2` and nowhere else, so `E e(theta_d - theta_e)` contains the
factor `(1 + e(1/2))/2 = 0`. Hence `E|f_n(k)|^2 = sum_a 4^-a E|f_(n-1)(k-a)|^2` and, from `f_0 = 1`, **`E|f_n(k)|^2 =
3^-n` exactly** for every window position `k` and every `n`: the root-mean-square rate of the i.i.d. model is
exactly `1/sqrt3`. (Monte Carlo, `collatz_fixed_frequency_iid_second_moment_20261004.py`, 200 seeds, `n <= 100`:
`E[3^n |f_n(0)|^2] = 0.90, 0.79, 1.00, 2.78, 1.11, 0.53` with heavy-tailed errors `0.1-2`; consistent.)

**The Jensen gap.** The typical value `exp(E log 3^(n/2)|f_n(0)|)` of the same model decays (`0.61, 0.50, 0.39,
0.28, 0.26, 0.18` at `n = 10..100`), so the almost-sure rate of one realisation is below the rms rate: `0.573` at
`n <= 1500` (E5) against `0.5774`. This is the ordinary gap between the mean-square and the Lyapunov growth of a
product of fluctuating factors. H1 concerns one realisation (the real digits), so its natural rate is the typical
one, `0.569`; the mean-square rate `1/sqrt3` is an upper envelope that the sup over `n` can touch only at the
ridges.

**E5d (window energy along the real digits; `collatz_fixed_frequency_window_energy_20261004.py`, `N = 600`,
`W = 100`).** The per-level ratio of the window energy `e_n = sum_(k=-W)^0 |f_n(k)|^2` in rescaled units
(incoherent `= 1`) has mean `0.994, 1.019, 0.997` for the real digits of `3^-n, 5^-n, 7^-n` and `1.021, 0.994,
1.004` for three i.i.d. seeds (medians `0.96-0.98`, heavy-tailed): **in mean square the real digits are
incoherent within the noise, like random digits**; the fitted energy rates `0.5723, 0.5746, 0.5718` (real) and
`0.5727, 0.5716, 0.5729` (i.i.d.) do not separate `q = 3` at this length, and the typical rates at `N = 600`
(`0.5706, 0.5741, 0.5710` real; `0.5715, 0.5724, 0.5727` i.i.d.) are within `0.002` of each other. The `0.6%`
deficit of section 4 is therefore a long-run effect of the fluctuation structure (the ridges), not of the second
moment; the wider control of section 4a (eight seeds, `q = 3, 7, 11, 13, 19`, `N = 1500`) decides whether it is
real and whether it is a `3 mod 8` phenomenon (`11 = 3 mod 8` against `13 = 5 mod 8`).

**What this makes of H1 (DIRECTION, typed).** H1 (`|mu_hat_n(1)| <= C rho^n`, `rho < 0.585`) decomposes into
(a) *mean-square incoherence of the real digit string:* the window energy decays at `(1/3)^n` up to
`exp(o(n))`, i.e. the covariance of the twisted vector `D_n f_(n-1)` with the geometric symbol (THM-4521(3) on the
window) has `o(n)` partial sums -- an autocorrelation statement about the binary digits of `3^-n`; and
(b) *sub-exponential ridges:* the excursions of `3^(n/2)|f_n(0)|` above its typical size are polynomial in `n`
-- the ridge inventory of S22/S23 (depths `v_3(u 2^Q -+ 1)`, Poisson), a statement about the 3-adic digits
of powers of `2`. (a) holds in mean square in the random model exactly and for the real digits within noise to
`n = 600`; (b) holds to `n = 2500` (largest excursion `10.4` at `n = 131`, nothing above `0.08` past `200`).
The two mirror endgames of the procgen synthesis (low binary digits of `3^A u` for Q1, low ternary digits of
`2^K w` for Q2) are exactly (a) and (b).

## 5. E2: `L^p` moments, tails, and THM-4263's condition (14)

Script `collatz_syracuse_law_lp_moments_20261004.py` (`n <= 14`, `A = 40`; the dense law by the exact recursion;
controls: `rho_n(-1) = 1.333, 2.095, 3.198, 4.846, 7.319, 11.018` and `E[rho_6^2] = 2.6669`, `rho_6(1) = 0.637`
reproduce S19).

| `n` | `E[rho^1.5]` | `E[rho^1.9]` | `E[rho^2]` | increment of `E[rho^2]` | `E[rho^2.5]` | `rho_n(-1)` | `rho_n(1)` |
|---|---|---|---|---|---|---|---|
| 6 | 1.4694 | 2.3329 | 2.6669 | +0.3108 | 5.693 | 11.02 | 0.637 |
| 8 | 1.5722 | 2.7780 | 3.2878 | +0.3106 | 8.732 | 24.89 | 0.516 |
| 10 | 1.6560 | 3.2036 | 3.9105 | +0.3116 | 12.737 | 56.12 | 0.425 |
| 12 | 1.7258 | 3.6125 | 4.5353 | +0.3126 | 17.936 | 126.40 | 0.358 |
| 14 | 1.7846 | 4.0062 | 5.1618 | +0.3135 | 24.610 | 284.50 | 0.315 |

Tails `P(rho_14 > t)` for `t = 2, 4, 8, 16, 32, 64, 128, 256`: `0.110, 0.0359, 0.0100, 0.00241, 4.7e-4, 6.4e-5,
4.7e-6, 3.1e-7`; stable in `n` at `t <= 16` from level 9 on; local exponents `1.61, 1.85, 2.05, 2.35, 2.89, 3.76,
3.91` (the steepening beyond `t = 32` is the finite-level cutoff at the maximal atom `0.975 (3/2)^n`).
`E[rho_n 1(rho_n > M)]` at level 14: `0.284, 0.145, 0.065, 0.024, 0.0064, 0.0010` for `M = 4, 8, 16, 32, 64, 128`
(levels 9..14 at `M = 32`: `0.0028, 0.011, 0.018, 0.015, 0.023, 0.024`, still rising).

**Verdict on P7.** SUPPORTED: `t^-2` tail in the bulk, linear second moment, geometric `p = 2.5`; the convergence
of `E[rho^p]` for `p` just below `2` is too slow to confirm at level 14 (increments `0.98^n` at `p = 1.9`), and the
UI tails at fixed `M` have not saturated. A level-18 run (the S19 FFT machinery) would settle `M <= 64`.

---

## 6. Reproduction

```bash
python 04-computation/experiments/collatz_fixed_frequency_rates_20261004.py 800 40          # E1, ~40 s
python 04-computation/experiments/collatz_fixed_frequency_drift_control_20261004.py 600 40  # E4, ~20 s
python 04-computation/experiments/collatz_fourier_mass_real_coordinate_20261004.py 13 40    # E3, ~30 s
python 04-computation/experiments/collatz_syracuse_law_lp_moments_20261004.py 14 40         # E2, ~10 s
python 04-computation/experiments/collatz_fixed_frequency_random_digits_20261004.py 1500 40 3   # E5
```
Outputs beside the scripts (`.out`; the `.log` files are the raw runs).

## 7. Verdicts

| claim | status |
|---|---|
| P6: one rate for six fixed units, `0.569 +- 0.003`, below Parseval | VERIFIED to `n = 800` |
| the rate is `q`-independent (`q = 3, 5, 7, 11`) | VERIFIED to `n = 600` (OBSERVED constant `~0.57`) |
| the phase sequence is the digit string of `-u q^-n`; incoherent rate `1/sqrt3` by Cauchy-Schwarz | PROVED (one line each) |
| H1 reframed as a binary-digit statement about `3^-n` | DIRECTION (exact reformulation of what a proof must control) |
| real-coordinate concentration of the Fourier mass | REFUTED (`n <= 13`) |
| P7 (`L^p` iff `p < 2`); THM-4263 (14) for the Syracuse density | SUPPORTED, not confirmed at level 14 |
| H1, Collatz | OPEN |
