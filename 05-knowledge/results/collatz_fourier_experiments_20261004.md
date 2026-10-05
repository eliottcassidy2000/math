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

### 4a. The wide control (eight seeds per random model; `q = 3, 7, 11, 13, 19`; `N = 1500`)

`collatz_fixed_frequency_random_digits_20261004_wide.log` (same script, `8` seeds, `q` list):

| digit source | rate `300..1500` | rate `750..1500` |
|---|---|---|
| real `-3^-n` | **0.5689** | **0.5695** |
| real `-7^-n` | 0.5727 | 0.5726 |
| real `-11^-n` (`11 = 3 mod 8`) | 0.5722 | 0.5719 |
| real `-13^-n` (`13 = 5 mod 8`) | 0.5724 | 0.5726 |
| real `-19^-n` | 0.5711 | 0.5709 |
| i.i.d. bits, 8 seeds | `0.5741, 0.5720, 0.5723, 0.5725, 0.5726, 0.5738, 0.5716, 0.5719`: mean **0.5726**, sd `0.0009` | mean `0.5728`, sd `0.0015` (one seed at `0.5697`) |
| random odd multiplier, 8 seeds | `0.5736, 0.5721, 0.5727, 0.5729, 0.5733, 0.5731, 0.5720, 0.5727`: mean **0.5728**, sd `0.0005` | mean `0.5732`, sd `0.0011` |

**Verdicts.** (i) The random models agree with each other: `0.5727 +- 0.001`, i.e. `0.8%` below `1/sqrt3`,
the Jensen gap of section 4b. (ii) `q = 7, 11, 13` sit at the random level; **the "`3 mod 8`" sub-hypothesis
(`11` like `3`, `13` like `5`) is REFUTED.** (iii) `q = 3` is `0.0037` below the i.i.d. mean (`4` standard
deviations of the seed scatter; `7` of the random-multiplier scatter) on `300..1500`, and `0.0033` below on
`750..1500` (`2.2` and `3.4` sd); the six fixed units `u = 1, 5, 7, 11, 13, 17` of `q = 3` rerun to `N = 1500`
(`collatz_fixed_frequency_rates_20261004.log`; the `N = 800` run kept as `_N800`) give `0.5691, 0.5699, 0.5692,
0.5712, 0.5709, 0.5697` on `200..1500`: mean **`0.5700 +- 0.0009`** against the twelve i.i.d. runs' `0.5723 +-
0.0008` -- a difference of `0.0023`, five standard errors: **the digits of `3^-n` produce a typical decay about
`0.4%` per level faster than random digits, for every unit tried, and `3` is the only multiplier among `3, 5, 7,
11, 13` to do so** (`19` is marginal, `0.5711`; the `q = 5, 7` unit families are in section 4c). OBSERVED; mechanism OPEN (candidates: the fluctuation structure -- section 4b and
E5e -- or an arithmetic anticorrelation specific to the carry automaton of `x -> 3x`). For H1 only the sign
matters: the Collatz digits are at least as cancelling as random ones in this functional.

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

### 4c. Where the `q = 3` excess does NOT come from (E5e, E6), and the level-to-level phase law (PROVED)

**E5e (fluctuations; `collatz_fixed_frequency_fluctuations_20261004.py`, `N = 1500`, window `300..1500`).**
Detrended log-variance of `3^(n/2)|f_n(0)|`: real `3^-n` `6.19`, `5^-n` `5.98`, `7^-n` `3.52`, `11^-n` `3.80`;
four i.i.d. seeds `6.86, 4.65, 2.77, 4.20`. Largest upward excursion: `+9.1` (`q = 3`, `n = 889`), `+9.7` (i.i.d.
seed 0, `n = 740`); counts of `n` with `v_n >= 3 x` trend: `389, 414, 274, 336` (real) against `349, 348, 326,
365` (i.i.d.). **The Collatz digits are not more fluctuating than random ones**; the Jensen gap does not single
out `q = 3`. The twelve i.i.d. runs now pooled (E5b's eight and E5e's four) give `0.5723 +- 0.0008` on
`300..1500`.

**E6 (within-level digit statistics; `collatz_digit_phase_autocorrelation_20261004.py`, `n = 60`, `M = 2^20`
digits).** For the real strings of `-3^-60, -5^-60, -7^-60, -11^-60`: digit balance `0.4991-0.5000`, phase
autocorrelations `|C_s| = |(1/M) sum_d e(theta_(d+s) - theta_d)|` at every lag `s <= 40` of size `10^-3 =
M^(-1/2)` (rms `9.6-10.0 x 10^-4`, exactly the i.i.d. level `9.5-9.9 x 10^-4`), bit correlations at lags `1..8`
of size `10^-3`. **Along one level the real digits are indistinguishable from random.** (The digits of `q^-n mod
2^M` are the digits of the power `q^(2^(M-2) - n) mod 2^M`, so this is the Dupuy-Weirich-type averaged
equidistribution of the low binary digits of powers of `q`, seen on one exponent.)

**The level-to-level law (PROVED, one line).** With `R_n = -u q^-n mod 2^d` and `theta_(n,d) = R_n/2^d`,
`q R_(n+1) = R_n mod 2^d`, so `q theta_(n+1,d) = theta_(n,d) mod 1`: **`e(theta_(n+1,d))` is a `q`-th root of
`e(theta_(n,d))`**, the branch being fixed by `R_n mod q`. For Collatz the phase at depth `d` of level `n+1` is
a cube root of the phase at depth `d` of level `n`. The random models differ exactly here: i.i.d. digits have
unrelated levels, the random-multiplier model relates levels by a `Q`-th root with `Q ~ 2^40` (effectively
unrelated), and the `x q^-n` family (`u` varying, E1/E5f) keeps the `q`-th-root law with a random start. The
`q = 3` excess therefore lives in the cross-level coupling of the cube-root law with the geometric-weight
recursion (section 4a: six units of `q = 3`, all below the random models; `q = 5, 7, 11, 13` single units at the
random level; the `q = 5, 7` unit families are run in `collatz_fixed_frequency_unit_families_20261004.py`).
Mechanism OPEN; the cleanest next statement would be the typical rate of the `q`-th-root random model as a
function of `q`.

**The adjacent-depth identity (PROVED, one line) and a failed naive prediction.** Along one level the phases
obey the odometer law `theta_(d-1) = 2 theta_d (mod 1)`, so for `q = 2^j +- 1` the `q`-th-root law reads
`omega_n(d) = omega_(n+1)(d)^q = omega_(n+1)(d-j) omega_(n+1)(d)^(+-1)`: the level-`n` phase at depth `d` is a
product of two level-`(n+1)` phases at depths `d` and `d-j`. **`q = 3` (`j = 1`) is the only root order for
which the two depths are adjacent**, i.e. for which the cross-level law pairs the two terms of the recursion
with the largest weights (`2^-1` and `2^-2`); `q = 5` pairs depths two apart, `q = 7, 9` three apart, `q = 15, 17`
four apart, and `q = 11, 13` are sums of three phases. This is a precise sense in which the Collatz tower has
the strongest cross-level coupling, and it would predict an excess decreasing like `2^-j` (`q = 5`: half of
`q = 3`'s). The data refute the naive form: `q = 5` sits `0.0008` ABOVE the i.i.d. level (six units,
`0.5731 +- 0.0007`), not `0.001` below. So the identity locates where the coupling is strongest but does not by
itself give the sign or size of the effect; typed DIRECTION.

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

### 4d. Per-level increments, and the root-order reading (E5h; DIRECTION)

From the saved E5e arrays (`n = 300..1500`): the per-level log-increments of `3^(n/2)|f_n(0)|` have mean
`-0.0133` (`q = 3`; rate `0.5697`), `-0.0087` (`5`), `-0.0059` (`7`), `-0.0046` (`11`), `-0.0066..-0.0100`
(four i.i.d. seeds), with variances `0.81, 0.68, 0.60, 0.88` (real) and `0.68-0.74` (i.i.d.), skew `~0`, and
heavy-tailed squared ratios (means `4-17`). The variance does not order the rates (`q = 11` has the largest
variance and the highest rate), so the `q = 3` excess is not a per-level fluctuation effect either.

**Root-order reading.** The random models differ only in how consecutive levels are coupled: unrelated
(i.i.d. digits, `0.5723 +- 0.0008`), coupled by a `Q`-th root with `Q ~ 2^40` (`0.5728 +- 0.0005`), coupled by a
cube root with a random start (the `u`-family of `q = 3`, `0.5700 +- 0.0009`). The `q = 5, 7` families and the
composite and prime-power root orders (section 4f) all sit at `0.5725 +- 0.0007`: the typical rate does not
depend on the root order at all, and `q = 3` is the single exception.
ANALOGY (typed, lost coordinate named): the tower of levels related by `q`-th roots is a Cartier-type tower
(the `q`-th-root-of-Frobenius structure of THM-4210's Rule-30 carriers), and THM-4210's lesson -- bounded
truncations of a Frobenius/Cartier carrier do not decide an all-scale statement -- is the shape of H1's
difficulty; nothing transfers as a method (the Collatz tower is on phases, not on power series).

### 4f. The root order does not matter either: one rate for every `q` (E5g, E5i, E5j) -- and a retracted reading

`collatz_fixed_frequency_unit_families_20261004.py` (six units each, `N = 1500`): `q = 5`: `0.5731 +- 0.0007`;
`q = 7`: `0.5723 +- 0.0009` -- the i.i.d. level (`0.5723 +- 0.0008`), with `q = 3` (`0.5700 +- 0.0009`) the only
prime clearly below and `q = 5` marginally above.

**Retraction (same session, self-caught).** The first version of this section read the runs of
`collatz_fixed_frequency_root_order_20261004.py` and `..._root_order2_...py` as giving `0.3305` for `q = 9, 15,
21, 25, 33, 45, 49` and pre-registered a "digit-counting law" `rate = c^(Omega(q))`. Those two scripts' reporting
line applied the `1/sqrt3` rescaling a second time (the vector is stored rescaled by `3^(n/2)`, the fit already
removes it); the pre-registered control `q = 5` in the second script came out `0.3309` against its known `0.5731`
(E4, E5g) and exposed the error. Dividing the printed values by `3^(-1/2)` (archived logs
`*_buggy_reporting.log`; the corrected scripts rerun to the same values, outputs `.out`):

| `q` | units | corrected rates | mean |
|---|---|---|---|
| 9 | 1, 5, 7, 11 | `0.5731, 0.5724, 0.5723, 0.5723` | `0.5725 +- 0.0004` |
| 15 | 1, 7, 11, 13 | `0.5737, 0.5724, 0.5705, 0.5730` | `0.5724 +- 0.0013` |
| 21 | 1, 5, 7, 11 | `0.5712, 0.5742, 0.5728, 0.5731` | `0.5728 +- 0.0012` |
| 27 | 1, 5, 7, 11 | `0.5743, 0.5714, 0.5728, 0.5733` | `0.5730 +- 0.0012` |
| 25 | 1, 2, 7 | `0.5721, 0.5723, 0.5711` | `0.5718 +- 0.0007` |
| 49 | 1, 2, 11 | `0.5731, 0.5731, 0.5728` | `0.5730 +- 0.0002` |
| 45 | 1, 2, 7 | `0.5726, 0.5726, 0.5737` | `0.5730 +- 0.0006` |
| 33 | 1, 2, 7 | `0.5712, 0.5712, 0.5724` | `0.5716 +- 0.0007` |
| 5 (control) | 1, 2, 7 | `0.5733, 0.5733, 0.5728` | `0.5731 +- 0.0003` |
| 3 (control, buggy run) | 1, 5, 7, 11 | `0.5690, 0.5700, 0.5693, 0.5714` | `0.5699 +- 0.0011` |

-- every `q >= 5` at the universal `0.5725 +- 0.0007`, and the `q = 3` control reproduces the anomaly
(`0.5699`), so the buggy run was internally consistent and only its last reporting line was wrong.
So the typical cold rate is the same constant for every odd `q`, prime, prime power or composite, and the only
exception remains `q = 3` at `0.5700`. The CRT remark stands as an identity (the frequency-`1` character of
`Z/15^n` is a product of a `3`-adic and a `5`-adic character), but the product's typical rate is not the product
of the typical rates: a `q`-adic digit advanced per level costs the same `0.5725` whether it is one digit of one
prime, two digits of one prime, or one digit each of two primes. The lesson is logged in MISTAKES (2026-10-04,
opus Fourier note): a control value that contradicts an earlier measurement must be checked before any
pre-registered prediction is read as confirmed -- here `q = 25` and `q = 49` "confirmed" the law for an hour
before the control line was read.

### 4g. The anomaly localized: it is the cube-root coupling of consecutive levels (E5k, hybrid towers)

`collatz_fixed_frequency_hybrid_20261004.py` (`N = 1500`, `A = 40`, rates over `300..1500`; reference: real
`0.5700 +- 0.0009`, i.i.d. `0.5723 +- 0.0008`):

| tower | `u = 1, 5, 7` (or three seeds) | mean |
|---|---|---|
| A: real digits of `-u 3^-n` for `n > 20`, i.i.d. digits for `n <= 20` | `0.5689, 0.5705, 0.5693` | **`0.5696`** |
| A: real for `n > 100`, i.i.d. for `n <= 100` | `0.5701, 0.5698, 0.5692` | **`0.5697`** |
| B: i.i.d. for `n > 20`, real for `n <= 20` | `0.5721, 0.5734, 0.5728` | `0.5727` |
| B: i.i.d. for `n > 100`, real for `n <= 100` | `0.5726, 0.5730, 0.5727` | `0.5727` |
| C: a genuine string of `-u_n 3^-n` at every level, with a fresh random unit `u_n` per level | `0.5725, 0.5735, 0.5715` | `0.5725` |

**Reading.** (A) The excess is a steady-state property of the high levels (it survives replacing the first 20
or 100 levels, where the phase strings are periodic inside the window, by random digits). (B) It disappears
when the high levels are random. (C) **It disappears when every level is a genuine `3^-n` digit string but the
unit changes from level to level** -- i.e. when consecutive levels are no longer related by the cube-root law
`e(theta_(n+1,d))^3 = e(theta_(n,d))` of section 4c. So the Collatz-specific `0.4%` per level extra cancellation
is produced by the exact cross-level relation between the digits of `u 3^-n` and `u 3^-(n+1)` (the carry
automaton of `x -> 3x` linking consecutive levels), not by the digit strings themselves (E6) and not by their
fluctuations (E5e). The `q`-th-root couplings for `q >= 5` produce no excess (section 4f). Why the cube root in
particular interacts with the geometric weights `2^-a` while higher roots do not is OPEN; the ensemble
experiment `collatz_fixed_frequency_ensemble_rms_20261004.py` (sixteen units, `q = 3` and `5`) tests whether the
excess is a mean-square anticorrelation or a Jensen effect (results in 4h).

### 4h. Mean square or Jensen? The sixteen-unit ensembles (E5l) and their calibration

`collatz_fixed_frequency_ensemble_rms_20261004.py` (`N = 1000`, `A = 40`, window `300..1000`, sixteen units
prime to `3q`): the ensemble mean square `E_u[q^n |mu_hat_n(u)|^2]` and the ensemble mean log.

| `q` | ensemble rms rate | ensemble typical rate | per-unit typical | Jensen gap per level |
|---|---|---|---|---|
| 3 | **`0.5703`** | `0.5692` | `0.5692 +- 0.0012` | `+0.0019` |
| 5 | `0.5735` | `0.5728` | `0.5728 +- 0.0013` | `+0.0013` |
| i.i.d. digits, 16 seeds (true rms rate exactly `0.5774`, E5c) | `0.5732` | `0.5722` | `0.5722 +- 0.0011` | `+0.0017` |
| random odd multiplier, 16 seeds | `0.5730` | `0.5726` | `0.5726 +- 0.0011` | `+0.0006` |

**Calibrated reading: INCONCLUSIVE.** The sixteen-member "ensemble rms rate" of the i.i.d. model is `0.5732`
although its true mean-square rate is exactly `0.5774` (E5c): with ensemble mean squares this heavy-tailed
(100-block means for `q = 3`: `9e-3, 5e-5, 5e-5, 2e-2, 8e-9, 3e-7, 3e-9`) the estimator inherits the Jensen bias of
the members, about `-0.004`. Against that baseline `q = 3` sits `0.003` lower in the ensemble statistic, exactly
as in the typical rate (`0.5692` against `0.5722`); so the sixteen-unit ensemble cannot separate a true
mean-square anticorrelation from a typical-rate effect, and the `1.2%` reading of the first version of this
section is withdrawn. What is exact instead (PROVED, Parseval): over the FULL frequency ensemble `t mod 3^n`,
`E_t |mu_hat_n(t)|^2 = sum_y mu_n(y)^2 = (3/2) 3^-n E_units[rho_n^2]`, the collision probability -- rate exactly
`1/sqrt3` with the prefactor `(3/2) E[rho_n^2] ~ 0.47 n` (E2), against the prefactor `1` of the i.i.d. model. The
Collatz ensemble therefore has the same exponential mean-square rate as the random model and a heavier tail
(the linearly growing collision prefactor is the forward closure of the `-1` spike), which by Jensen is
consistent with a lower typical rate; but the fixed small units are a null subfamily of that ensemble, so this
is a heuristic for the `0.4%`, not a derivation. The excess itself is robust (E5f: five standard errors;
E5k: localized to the cube-root coupling); its mechanism stays OPEN.

**E5m (per-level coherence ratio; inconclusive).** `collatz_fixed_frequency_coherence_20261004.py` measures
`kappa_n = |f_(n+1)(0)|^2 / sum_a 4^-a |f_n(-a)|^2` along one run (mean `1` under random phases whatever the
input). Means over `n = 300..1500`: real `3^-n` (`u = 1, 5, 7`) `0.985, 0.961, 0.971`; real `5^-n` `0.985`; real
`7^-n` `1.179`; i.i.d. `1.026, 1.012, 1.016`; medians `0.89-0.97` everywhere except `7^-n` (`1.15`); typical values
`exp(E log kappa) = 0.73-0.79` (`0.98` for `7^-n`). The ratio is a heavy-tailed random variable whose mean over
`1200` levels scatters by `+-0.1` (the `7^-n` value), so a `0.4%` per-level effect is not resolvable this way;
the `q = 3` means below `1` are suggestive only. The decisive test is the two-hundred-member ensemble
(`collatz_fixed_frequency_large_ensemble_20261004.py`, i.i.d. calibration at the same size; results in 4i).

### 4e. The web around the cold rate (connections found by the niche search; all typed)

| repo thread | the object there | the map to the cold-frequency problem | preserved | lost / sidecar | type |
|---|---|---|---|---|---|
| [THM-4263](../../01-canon/theorems/THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport.md) finite-factor density transport | condition (14): uniformly integrable fibre weights | the Syracuse reference density `rho_n` is the fibre weight of the factor `Z/3^n -> words`; (14) is `E[rho_n 1(rho_n > M)] -> 0` uniformly | density-one transport target -> source | pointwise statements (Collatz needs them); the tail law itself | EXACT reading; E2 SUPPORTS (14) via the `t^-2` tail |
| [Q1 mirror](collatz_procgen_20260922_q1_mirror.md) section 5 | the endgame reads the low binary digits of `3^A u` along Beatty exponents | the window phases ARE the low binary digits of `3^-n = 3^(2^(M-2)-n) mod 2^M` | the digit object (powers of `3` in base `2`) | the Beatty exponent selection (Q1 needs specific `A`); here all `n` enter with geometric weights | EXACT identification of the object; the two statements differ in the quantifier |
| Dupuy-Weirich (CITED there) | low `q`-adic digits of `p^n` equidistributed on average over `n` | E6: equidistribution of the digits of one `3^-n` to the noise floor, `M = 2^20` | averaged balance | the pointwise/along-the-recursion statement (a) needs more than balance: the weighted cross-level sums | CITED + FINITE-EXACT consistency |
| [THM-3848](../../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md) / the Mahler `3/2` frontier | `frac(xi (3/2)^n)`, the safe-prefix tree, loneliness `2/5` of the mixed-power speed row | the digits of `3^n mod 2^d` are `frac(3^n/2^d) 2^d`: the phase tower is Mahler's object read at depth `d` | the `3^n mod 2^d` arithmetic | Mahler fixes `xi` and varies `n`; here `n` and `d` both vary with geometric weights; no Z-number enters | ANALOGY (shared object, different predicate) |
| [THM-4210](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md) Rule 30 Cartier tree | even/odd Cartier lift, Frobenius/Cartier carrier, all-scale admissibility | the level-to-level `q`-th-root law is a Cartier-type tower on phases | "bounded truncations do not decide an all-scale statement" | the power-series structure; nothing transfers as a method | ANALOGY |
| [THM-485](../../01-canon/theorems/THM-485-two-temperatures-viswanath.md) Viswanath's constant | Lyapunov exponent of random Fibonacci products via the Stern-Brocot stationary measure; the golden-mean shift = Zeckendorf | the typical cold rate (`0.5723` i.i.d., `0.5700` Collatz) is a Lyapunov exponent of a random product of THM-4520's operators; the Collatz parity language is the golden-mean shift (THM-4528) | "a typical rate below the mean-square rate, computed from an invariant measure" | the product here acts on an infinite window, not on `R^2`; no stationary measure identified | DIRECTION: a Viswanath-type exact value for the i.i.d. rate would turn (a) of HYP-9176 into a computable constant |

## 5b. What the day established about the cold rate (consolidated; the owner's experiment loop, cycles 2-3)

| statement | status | evidence |
|---|---|---|
| the frequency-`u` coefficient of the `q`-adic Syracuse law, `u` fixed, decays at one typical rate `0.5725 +- 0.0007` for every odd `q >= 5`, prime, prime power or composite, and for every unit tried | OBSERVED (`N = 1200-1500`; `q = 5, 7, 9, 11, 13, 15, 19, 21, 25, 27, 33, 45, 49`; 6 units at `5, 7`; 3-4 units elsewhere) | E4, E5g, E5i/E5j (corrected) |
| the same rate for i.i.d. digit strings (`0.5723 +- 0.0008`, twelve seeds) and for random odd multipliers (`0.5728 +- 0.0005`, eight) | OBSERVED | E5, E5b |
| Collatz, `q = 3`: `0.5700 +- 0.0009` over six units, five standard errors below -- the only exception | OBSERVED; mechanism OPEN | E1, E5f |
| the rate is not `q`'s Parseval scale `q^(-1/2)` (for `q >= 5` the fixed units are astronomically above it) | OBSERVED | E4 |
| the phases are the binary digits of `-u q^-n`; consecutive levels are `q`-th roots | PROVED (identities) | section 0.3, 4c |
| in the i.i.d. model `E|f_n(k)|^2 = 3^-n` exactly; the typical rate is below the rms rate by a Jensen gap | PROVED + OBSERVED | 4b, E5c |
| the Collatz digits are neither correlated within a level nor more fluctuating; the `q = 3` excess lives in the cross-level coupling | OBSERVED (null results) | E6, E5e, E5h |
| the mean-square (window-energy) rate of the real digits is incoherent within noise | OBSERVED (`n <= 600`) | E5d |
| "the cold rate counts adic digits, `c^(Omega(q))`" | REFUTED (reporting bug, caught by the `q = 5` control; MISTAKES) | 4f |
| "the Fourier mass concentrates near dyadic real frequencies" | REFUTED | E3 |
| H1 = (a) digit incoherence in mean square + (b) polynomial ridges | HYP-9176 (CONJECTURED) | 4b |

**Where this leaves H1 (HYP-9166).** The rate H1 needs, `rho < 0.585`, is beaten by every model and every
`q` by `1.6-2.8%`; the one Collatz-specific fact is a `0.4%` extra cancellation, in the right direction. A proof
would have to show (a) for the digits of `3^-n` and (b) for the ridge inventory; neither is attempted here.

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
