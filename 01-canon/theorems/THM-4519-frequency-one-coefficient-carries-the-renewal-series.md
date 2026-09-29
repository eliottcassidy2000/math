---
id: THM-4519
title: "The frequency-one coefficient carries the renewal series of the resonant Syracuse Fourier coefficient: exactly c_J = 2^-s mu_hat_J(1) W_{h,J}(s) with W the unweighted top-part phase sum over compositions (|W| <= C(s,h-J)), so Theorem C of the five-mirrors note holds under the single-sequence hypothesis H1: |mu_hat_n(1)| <= C rho^n, rho < log_2 3 - 1; the level phases omega_n(k) = e(2^k/3^n) on Z/(2 3^(n-1)) have a flat Gauss-sum spectrum on the primitive characters; FINITE-EXACT: the identity to h = 2000, |mu_hat_n(1)| <= 0.03 (0.58)^n for 200 <= n <= 1200 with sup_n |mu_hat_n(1)|/0.585^n = 1.87 at n = 131, the resonant coefficient 0.32 h^-1 e^-hI on 200 <= h <= 600 bending toward h^-3/2, the growing remnant explained as a near-neutral multiplicative random walk of ridge transport (energy-weighted slope 4/3) exiting through frequency one, and the ridge inventory by depth v_3(u 2^Q -+ 1) matching the residue-class count (2/3) W/(2 3^(k-1)) per multiplier"
status: >
  PROVED: Lemma R'' (exact identity), Theorem C' (conditional on H1), the
  Gauss-sum spectrum of the level phases (standard, CITED), the energy-weighted
  transport slope 4/3 under random phases, the residue-class arithmetic of the
  seeds (2 a primitive root mod 3^k; LTE for u = 1). FINITE-EXACT: J-series and
  identity check (relative 10^-13..4 10^-14) at h = 40..2000 with valuation
  truncation A = 40..190; mu_hat_n(1) to n = 1200; the negative family to n = 420
  (m <= 80) and n = 140 (m <= 160); seed inventory u <= 200, Q <= 700. Truncation
  control: mu_hat_n(1) from A = 60 and A = 120 agree to 6e-14 relative for all
  n <= 600 (A = 100 vs 120: identical), so the truncated words cancel like the
  rest; the rescaled recursion v_n = 3^(n/2) f_n (no underflow) gives mu_hat_n(1) to
  n = 2500: block rates 0.5721, 0.5668, 0.5703, 0.5683 over 200..600, .., 1800..2500,
  3^(n/2)|mu_hat_n(1)| median 9e-2 -> 3e-15 across the ten 250-blocks, no value above
  0.078 past n = 200: exponentially below the Parseval rate.
  EMPIRICAL: the prefactor law R(h) = |mu_hat_h(2^s*)| h^(3/2) e^(hI) = 4.53,
  5.62, 6.44, 7.84, 8.79, 9.87, 10.42, 10.72 at h = 200, 300, 400, 600, 800, 1200,
  1600, 2000 (R/sqrt h = 0.320, 0.324, 0.322, 0.320, 0.311, 0.285, 0.261, 0.240;
  local prefactor exponent 0.97 -> 1.37 drifting toward 3/2); H1 constants; the transport
  statistics (per-level log-gain std 0.35, geometric mean 0.984 over the wave).
  OPEN: H1. NOT independently audited. Consequences for S22: H -> H1 (the ridges
  of the negative family matter only at exponent 0); the gamma_J limits exist
  for fixed J (gamma_5 = 2.05 e^(-0.29 i), gamma_7 = 2.60 e^(1.06 i), ...) but
  grow like J 0.987^J, so the constant of (T2) is a quantity of h >> 10^4 and
  M(h) is 0.32 h^-1 e^-hI in the accessible range (M/P_h ~ sqrt h, S21's rise);
  the growing remnant is a +1.1 sigma excursion of a driftless transport whose
  only one-step-coherent level is n = 128 (3^128 = 1 mod 2^9), exiting at n =
  131 as the sup of the H1 ratio; the small-exponent negative family is 100x
  below the Parseval scale for 140 <= n <= 420; arrivals of seeds farther than
  Q ~ 200 are dead on arrival.
source: collatz-necklace-20260929 session (mac-mini), part 2, 2026-09-29; owner directive: prove H in any form with rate below 0.585, explain the growing remnant, the J-term limits at larger h, the ridge inventory by v_3(u 2^Q -+ 1)
depends_on:
  - 05-knowledge/results/collatz_five_mirrors_20260929.md (section 2d: Lemma R', the mass law, Theorem C; the closed frequency recursion)
related:
  - 05-knowledge/results/collatz_h1_20260929_frequency_one_jseries_remnant_ridges.md (the full note)
  - 05-knowledge/results/mazur_positive_density_20260928.md (Theorem A: the harmonic mass; Proposition 2 of the five-mirrors note)
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md (the exponent h*)
note: 05-knowledge/results/collatz_h1_20260929_frequency_one_jseries_remnant_ridges.md
scripts:
  - 04-computation/experiments/collatz_h1_20260929_jseries.py
  - 04-computation/experiments/collatz_h1_20260929_mu1_track.py
  - 04-computation/experiments/collatz_h1_20260929_remnant.py
  - 04-computation/experiments/collatz_h1_20260929_arrivals.py
  - 04-computation/experiments/collatz_h1_20260929_ridge_inventory.py
outputs: 05-knowledge/results/collatz_h1_20260929_{jseries_h40_validation, jseries_h40..h2000, gamma_table_h200_1600, mu1_track_n1200, mu1_rescaled_n2500, mu1_trunc_check, remnant_n140, arrivals_n420, ridge_inventory_u200_Q700}.out (scripts collatz_h1_20260929_{mu1_trunc_check, mu1_rescaled, gamma_table}.py)
script_sha256: 718d812623df56e77860262a4ba4755adb5111301014a6bffe2094f8129bd4d8 (jseries), e28d8c07a0cd7e3bdba447d0526da1d3d0a41494ea4e92073282fd2e996a8c08 (mu1_track), 44fae1b567bf02b19dbc9ff56ec9a1c4cdc8fe6ddd1d6fcf0db22c648c705651 (remnant), f864993cb5209556056084d5d527314368b4a31d78a2f4364fdd8ae3c308921c (arrivals), 9db9060775ec5a3c622c71923ddd843206390d5042ff2e6acf8fce1cedb76fff (ridge_inventory)
hash_basis: raw LF bytes
audit: NOT independently audited; internal controls: the h = 40 tables reproduce S22 to four digits, the identity sum_J c_J = f_h(s) holds to 4e-14 at every h, the seed counts match the residue-class prediction, and the n = 128 coherent step is the audit's |G_128| = 0.967.
---

# THM-4519 -- the frequency-one coefficient carries the renewal series

**PROVED (identity, Theorem C', structure) + FINITE-EXACT + EMPIRICAL; not independently audited.** Full note: [collatz_h1_20260929_frequency_one_jseries_remnant_ridges](../../05-knowledge/results/collatz_h1_20260929_frequency_one_jseries_remnant_ridges.md).

## 1. Lemma R'' and Theorem C'

With S22's coordinates (`kappa_j = s - T_j`, `J(w) = #{j : kappa_j < 0}`, `c_J` the contribution of the words with `J(w) = J`):

`c_J = 2^-s mu_hat_J(1) W_{h,J}(s)`, `W_{h,J}(s) = sum_{u in Z_{>=1}^{h-J}, T(u) <= s} prod_{j=J+1}^{h} omega_j(s - T_j(u))`, `|W_{h,J}| <= binom(s, h-J)`,

because the crossing step `a_J = kappa + b` and the bottom word combine into `2^-kappa sum_b 2^-b omega_J(-b) mu_hat_{J-1}(2^-b) = 2^-kappa mu_hat_J(1)` (the frequency recursion at `t = 1`). Hence `|c_J| <= mass_h(J) |mu_hat_J(1)|`, and **Theorem C'**: if `|mu_hat_n(1)| <= C rho^n` with `rho < log_2 3 - 1` (H1), then `|mu_hat_h(2^s)| <= e^{-hI} e^{theta* delta} [1 + C (m-1)^{-1}/(1 - rho/(m-1))]` for all `h` and `0 <= s = hm - delta`. H1 is weaker than S22's H (`|mu_hat_J(1)| <= Ñ_{J-1}`): the negative family enters only through its value at the exponent `0`.

**Gauss-sum structure.** On `Z/L_n`, `L_n = 2 3^{n-1}`, `hat omega_n(xi) = tau(bar psi_xi)` with `psi_xi(2^k) = e(k xi/L_n)`; `|hat omega_n(xi)| = 3^{n/2}` for `3 ∤ xi`, `0` for `3 | xi` (`n >= 2`). The frequency recursion is `f_n = G * (omega_n lift f_{n-1})`: a fixed geometric convolution composed with a spectrally flat multiplier, and H1 concerns the entry `f_n(0)`.

## 2. The data

* H1 (`n <= 1200`, and to `n = 2500` by the rescaled recursion): `sup_n |mu_hat_n(1)|/0.585^n = 1.87`, `/0.58^n = 5.7`, `/3^{-n/2} = 10.3`, all at `n = 131`; over `n >= 200` the three suprema are `< 0.02`, `0.028`, `0.076`. Least-squares rates `0.5685` (`20..1200`), `0.5576` (`900..1200`); block rates `0.567–0.572` on `200..2500`, i.e. `|mu_hat_n(1)|` decays at about `0.569` per level, below the Parseval rate, with `3^{n/2}|mu_hat_n(1)| ~ 3 10^-15` by `n = 2500`. The normalised energy of the coefficients at `2^-m`, `m <= 80`, is `0.0–0.5` for all `140 <= n <= 420` against a Parseval value `56`: the small-exponent negative family empties out.
* The resonant coefficient at `s* = round(h log_2 3) - 6`: `R(h) = |mu_hat_h(2^{s*})| h^{3/2} e^{hI} = 3.16, 3.71, 3.97, 4.53, 5.62, 6.44, 7.84, 8.79, 9.87, 10.42, 10.72` at `h = 40, 80, 120, 200, 300, 400, 600, 800, 1200, 1600, 2000`; `R/sqrt h = 0.320 ± 0.004` on `200..600`, then `0.311, 0.285, 0.261, 0.240`. So `M(h) = 0.32 h^{-1} e^{-hI}` on `200 <= h <= 600` and the prefactor exponent drifts `1.0 -> 1.37` toward the fixed-`J` value `3/2`; `M/P_h ≍ sqrt h` in this range is S21's observed rise. Each `gamma_J = lim c_J/(h^{-3/2} e^{-hI})` exists (moduli and phases fixed to `1%` between `h = 600` and `800` for `J <= 10`), the terms behave like `J (0.987)^J` (mass `(m-1)^{-J}` against `|mu_hat_J(1)| ≍ 3^{-J/2}` against top coherence `≈ 3J/h`), and the Gaussian `J`-window of width `0.6 sqrt h` truncates the series far below its peak at `J ≈ 77` for every `h <= 10^4`.
* The growing remnant: the ridge slope is the energy-weighted mean step `4/3` of the recursion (random-phase transport with step law `3 4^{-a}`); along the argmax path the incoherent term is `0.26–0.29` of the previous peak energy (amplitude `0.9` per level), the cross term is positive at `43` of `51` levels, the per-level gain has geometric mean `0.984` and log-standard-deviation `0.35`; the packet energy grows `57 -> 483` over `n = 87..128` (true amplification), the only one-step-coherent level being `n = 128` (`3^128 = 1 mod 2^9`, gain `1.435`), and the packet exits through the exponent `0` (`3^{n/2}|mu_hat_n(1)| = 3.9, 7.2, 10.4` at `n = 127, 129, 131`). The growth (`+2.6` in the log over `44` levels) and the dominant ridge's decay (`-4.4` over `270` levels) are `1.1` and `0.8` standard deviations of one near-neutral multiplicative random walk.
* The inventory: depth `k = v_3(u 2^Q -+ 1)`; for fixed `(u, sign)` the seeds of depth `>= k` form one residue class mod `2 3^{k-1}`; counts over `134` pairs and `Q <= 700`: `389, 131, 50, 11, 0, 0, 0, 1, 0, 0, 1` for `k = 5..15` against the predicted `386.6, 128.9, 43.0, 14.3, 4.8, 1.6, 0.5, 0.18, 0.06, 0.02, 0.007`; S22's dominant ridge `(55, +, 423)` is the depth-`15` seed, `(41, -, 271)` the depth-`12` one, and the growing remnant's seed `(13, -, 154)` has depth `7`. Predicted arrivals at frequency one at `n ≈ 190, 215, 253, 263, 315, 331, 361, 377, 440, 536, 564` give at most `0.28` (`n = 189`): seeds farther than `Q ≈ 200` arrive dead.

## 3. Boundary

H1 is OPEN; no unconditional bound with rate below `1` on `|mu_hat_n(1)|` is proved. The constant of (T2) is not computable in the accessible range. The transport gain statistics are measured, not derived. Collatz is OPEN.
