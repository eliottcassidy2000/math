---
id: HYP-9176
title: "The H1 decomposition. (a) Digit incoherence: the fixed-frequency window energy e_n = sum_{k=-W}^0 |mu_hat_n(2^k)|^2 of the 3-adic Syracuse law decays at (1/3)^n exp(o(n)), i.e. the covariance of the twisted vector with the geometric symbol has o(n) partial sums -- an autocorrelation statement about the binary digits of 3^-n (CONJECTURED; holds exactly in mean in the i.i.d.-digit model, PROVED, and within noise for the real digits to n = 600); (b) ridge control: the excursions of 3^(n/2)|mu_hat_n(1)| above its typical size are polynomial in n (CONJECTURED; the ridge inventory of S22/S23, Poisson in the 3-adic zero lines u 2^Q = -+1 mod 3^k; largest 10.4 at n = 131, nothing above 0.08 past n = 200 to 2500). (a) and (b) together give H1 with rho = 1/sqrt3 + eps < 0.585. Universality: the typical rate of |mu_hat_n(u)| is the same for every fixed unit u and every odd multiplier q (0.5727 +- 0.001 for random digits; q = 5, 7, 11, 13 at that level), and the Collatz digits of 3^-n are about 0.4% per level MORE cancelling (0.5691-0.5712 over six units at N = 1500, mean 0.5700 +- 0.0009; OBSERVED, mechanism OPEN)"
status: >
  OPEN (CONJECTURAL). Evidence: the five Fourier experiments of
  collatz_fourier_experiments_20261004.md (E1 six units to n = 800; E4 four
  multipliers to 600; E5/E5b sixteen random-digit runs to 1500; E5c the exact
  second moment 3^-n of the i.i.d. model, PROVED one line; E5d window energies
  to 600; E5e/E6: the Collatz digits are neither more fluctuating nor
  correlated within a level; E5f: six units of q = 3 to n = 1500 give
  0.5700 +- 0.0009 against twelve i.i.d. runs at 0.5723 +- 0.0008, five
  standard errors; the level-to-level law e(theta_(n+1,d))^q = e(theta_(n,d))
  is PROVED and locates the q = 3 excess in the cross-level coupling; E5k
  hybrid towers CONFIRM the localization: real digits at the high levels with
  random digits below keep the excess (0.5696), random digits above lose it
  (0.5727), and genuine 3^-n strings at every level but with a fresh unit per
  level -- breaking the cube-root relation between consecutive levels -- lose
  it too (0.5725); the root orders q = 5..49, prime, prime power or composite,
  all sit at 0.5725 +- 0.0007, so the cube root alone interacts with the
  geometric weights; E5l sixteen-unit ensembles are INCONCLUSIVE on whether
  the excess is a mean-square anticorrelation (the estimator's Jensen bias,
  calibrated on the i.i.d. model, equals the effect); exact: over the full
  frequency ensemble the mean square is the collision probability
  (3/2) 3^-n E[rho_n^2], rate 1/sqrt3 with a linearly growing prefactor
  (Parseval, PROVED); and for the RANDOM-START tower (uniform 2-adic start,
  all levels by the Pascal identity theta_{N-m,d} = sum_i C(m,i) theta_{N,d-i},
  PROVED) the second moment is EXACTLY 3^-n at every level for every odd q
  (the path value Phi_q(a) = sum_i q^(n-i) 2^-D_i mod 1 is injective on
  paths, PROVED), while the fourth moment exceeds the i.i.d. model's for every
  q (weighted additive energy; x3 at n = 4): same mean square as random
  digits, exponentially heavier tails = the ridges; with INTEGER starts prime
  to q (the Collatz family proper) the scaled fourth moment 9^n E|f|^4 at
  n = 10..40 is 2.5-9 for q = 3 and 1.3-6 for q = 5 against 5-19 for i.i.d.
  digits -- no heavier tail (E5q, 2000 starts; the first version sampled
  units not prime to q and read a 50-100x inflation that was projectivity,
  audit C3, MISTAKES) -- while two hundred integer units prime to 3 give the
  typical rate 0.5694 +- 0.0002 against 0.5721 +- 0.0001 for two hundred
  i.i.d. seeds (E5n): the excess is real; E5r: the UNIFORM-start cube-root
  tower (uniform deepest string, levels by the exact Pascal law, no
  arithmetic input, second moment exactly 3^-n) has the same typical rate
  0.5698 +- 0.0009 while the uniform-start fifth-root tower has 0.5724 +-
  0.0005 -- so the excess is a Lyapunov-exponent property of the step-1
  Pascal coupling itself (mechanism localized to an exact random model; its
  value not derived, OPEN); E7: the
  frequency-one coefficient to n = 5000 decays at 0.5698 with excursions at
  most ~90x its trend and sup_(n>=200) |mu_hat_n(1)|/0.585^n = 0.0046,
  consistent with (b)) and
  mac-mini's |mu_hat_n(1)| to 2500 (HYP-9166). What it adds to
  HYP-9166: H1's rate is not 3's Parseval scale but the incoherent rate of the
  2-adic window recursion (identical for 5x+1, whose law is nowhere near
  Parseval-flat); the margin 0.577 vs 0.585 is a near-coincidence at q = 3;
  the two sub-statements are the two mirror digit endgames of the procgen
  synthesis (low binary digits of 3^A u, low ternary digits of 2^K w; Q1 and
  Q2), now as one explicit sequence each. Cheapest tests: (a) the window
  covariance partial sums along the real recursion to n = 2500 (bounded?
  sqrt n? -c n?); (b) the ridge maxima of 3^(n/2)|mu_hat_n(1)| to n = 10^4
  against the inventory's polynomial law; (c) the q = 3 excess: pooled rate
  over 12 units to n = 1500 against 16 random runs; (d) a proof of (a) for the
  i.i.d. model beyond the mean (the typical rate 0.5727 as a Lyapunov exponent
  of the random product of THM-4520's operators).
source: opus-2026-10-04-S1 (worktree collatz-synthesis-20261004), the Fourier experiments note (cycle 2 of the owner's experiment loop)
related:
  - 05-knowledge/hypotheses/HYP-9166-h1-frequency-one-coefficient-decays-below-critical.md
  - 01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md
  - 01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md
  - 01-canon/theorems/THM-4521-collatz-base-cycles-cohere-to-every-level-half-turn-intertwining-and-exact-energy-transfer.md (the energy identity, item 3)
  - 05-knowledge/results/collatz_fourier_experiments_20261004.md
  - 05-knowledge/results/collatz_procgen_20260922_q1_mirror.md (section 5: the low binary digits of 3^A u)
scripts: 04-computation/experiments/collatz_fixed_frequency_{rates,drift_control,random_digits,iid_second_moment,window_energy,fluctuations}_20261004.py
---

# HYP-9176 -- the H1 decomposition: digit incoherence plus ridge control

**The object.** `f_n(k) = mu_hat_n(2^k)` for `k <= 0` obeys
`f_n(k) = sum_(a>=1) 2^-a omega_n(k-a) f_(n-1)(k-a)` with
`omega_n(-d) = e(0.b_d b_(d-1) ... b_1 + 2^-d 3^-n)`, `b_i` the binary digits of
the 2-adic number `-3^-n` (PROVED identity, Fourier note section 0.3). H1 is
`|f_n(0)| <= C rho^n`, `rho < log_2 3 - 1 = 0.585`.

**(a) Mean-square incoherence.** In the i.i.d.-digit model the phases are
pairwise uncorrelated (the top bit `b_e` enters `theta_e - theta_d` with
coefficient `1/2`), so `E|f_n(k)|^2 = 3^-n` exactly. For the real digits the
window energy's per-level ratio is `1/3` times `(1 + cov_n)` with `cov_n` the
covariance of the twisted vector's spectrum with the geometric symbol
(THM-4521(3) on the window); measured mean ratios `0.994 (q = 3), 1.019 (5),
0.997 (7)` in rescaled units to `n = 600`. Conjecture: `sum_(m<=n) cov_m = o(n)`.

**(b) Ridges.** The typical value of `3^(n/2)|f_n(0)|` decays like
`(0.5727/0.5774)^n` in the random model and `(0.569/0.5774)^n` for the real
digits (the Jensen gap between the mean-square and the almost-sure rate);
the sup over `n` is attained at ridge arrivals (`n = 131`, `10.4`; HYP-9166's
inventory by depth `v_3(u 2^Q -+ 1)`), which the three-mirrors census finds
Poisson. Conjecture: the ridge heights are `n^O(1)`.

**Why (a) gives H1, and what (b) adds (logic corrected after the audit, C5).** (a) alone implies H1:
`|f_n(0)|^2 <= e_n = 3^-n exp(o(n))` gives `|mu_hat_n(1)| <= C_eps (1/sqrt3 + eps)^n` for every `eps > 0`, and
`1/sqrt3 + eps < 0.585` for `eps < 0.0076`. (b) is not needed for the rate: it describes where the sup over `n`
of `|mu_hat_n(1)|/rho^n` is attained (the ridges), hence the size of the constant `C`, and it is what makes the
typical rate sit below `1/sqrt3`. H1 does not imply (a). The two statements are the mean-square and the sup
forms of one incoherence.

**Non-consequences.** H1 itself does not prove Collatz (THM-4519: it gives
Theorem C's resonant rate `3^(h*-1)`); the universality across `q` says the
cold-frequency decay is not where the drift enters.

**Update 2026-10-05 (opus, `pascal-lyapunov-20261005`; [note](../results/collatz_pascal_tower_lyapunov_20261005.md)).** The uniform-start tower is a twisted Pascal triangle and `f_N(0)` is the Fourier transform at a Haar frequency of the law of the random backward Syracuse iterate (PROVED); its Lyapunov exponent is `-0.56234 +- 0.00006` (rate `0.56988 +- 0.00004`, levels `1000..6000`), against `-0.5574 +- 0.0001` (rate `0.5727`) for the step-2 tower, i.i.d. digits and random multipliers: the excess of (b) is `0.0049` per level, flat from level `1000` to `6000`. The exact Collatz coefficients `mu_hat_n(u)` of 200 integer units `u <= 600` prime to `3` (window kernel, `A = 40`, levels `300..1500`) have typical rate `0.56981 +- 0.00011` (`0.56989 +- 0.00015` on `600..1500`) against `0.56939 +- 0.00011` (`0.56946 +- 0.00015`) for the uniform-start tower computed the same way: the integer-start family and the uniform-start tower agree to `0.0004`, both at the step-1 value, so the excess uses nothing of the digits of `3^-n`. It is carried by valuations `>= 3`, grows smoothly as the coin favours long halving runs, is specific to `q = +-3` (`-3`: `0.0138`), and is invisible to every pair statistic (the mean square over every dyadic shell of the frequency is the same for every `q`, PROVED): a real-contraction mechanism was conjectured and REFUTED in the session. For (a) nothing changes: the mean square is exactly `3^-N`; for (b) the typical rate is now known to four decimals and its mechanism is a four-point (not pair) property of the path phases.
