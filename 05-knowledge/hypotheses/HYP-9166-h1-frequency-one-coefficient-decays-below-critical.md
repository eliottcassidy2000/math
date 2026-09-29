---
id: HYP-9166
title: "H1: the frequency-one Fourier coefficient of the 3-adic Syracuse law decays at a rate below log_2 3 - 1: |mu_hat_n(1)| = |E e(Y_n/3^n)| <= C rho^n with rho < 0.58496 — the single-sequence hypothesis under which Theorem C of the five-mirrors note holds (THM-4519), with the data |mu_hat_n(1)| <= 5.7 (0.58)^n for all n <= 1200 and <= 0.03 (0.58)^n for 200 <= n <= 1200"
status: >
  OPEN (CONJECTURAL). Implied by S22's hypothesis H (the weighted negative-power
  norm) and implying Theorem C's sharp rate 3^(h*-1) for the resonant window of
  the powers of two (THM-4519). FINITE-EXACT support to n = 1200: the sup of
  |mu_hat_n(1)|/rho^n is 1.87, 5.7, 10.3, 17.7 for rho = 0.585, 0.58, 0.5774,
  0.575, all attained at n = 131 (the arrival of the Q = 154 ridge at frequency
  one); over n >= 200 the same suprema are < 0.02, 0.028, 0.076, 0.67; the
  least-squares rate over 900..1200 is 0.5576. The stronger form rho = 3^(-1/2)
  (square-root cancellation) fails only at n = 127..135 and holds with constant
  10.3. Structure: mu_hat_n(1) = f_n(0) for f_n = G * (omega_n lift f_(n-1)), a
  geometric convolution composed with a Gauss-sum multiplier (flat spectrum on
  the primitive characters); the only mechanism found for a large value is the
  arrival of a ridge seed 2^-Q = -+u mod 3^k, which is enumerable (depth
  v_3(u 2^Q -+ 1), residue classes mod 2 3^(k-1)) and dead on arrival for
  Q >~ 200 (EMPIRICAL). Cheapest tests: mu_hat_n(1) to n = 4000 (the runs to
  2000 are in progress at close-out); the arrivals of the depth-7/8 seeds at
  Q = 729 (n ~ 564: 0.01 observed), 1458 (n ~ 1125: 0.00); a proof of
  |mu_hat_n(1)| <= C rho^n for any rho < 1 would already be new.
source: collatz-necklace-20260929 session (mac-mini), part 2, 2026-09-29
related:
  - 01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md
  - 05-knowledge/results/collatz_h1_20260929_frequency_one_jseries_remnant_ridges.md
  - 05-knowledge/results/collatz_five_mirrors_20260929.md (section 2d, hypothesis H and Theorem C)
---

# HYP-9166 -- H1: the frequency-one coefficient decays below the critical rate

`Y_n` is Tao's Syracuse random variable mod `3^n` (`Y_0 = 0`, `Y_{k+1} = 2^-A (3 Y_k + 1)`, `A` geometric of mean `2`). The hypothesis is `|E e(Y_n/3^n)| <= C rho^n` with `rho < log_2 3 - 1`. By THM-4519 it is exactly what the renewal series of the resonant coefficient needs: `c_J = 2^-s mu_hat_J(1) W_{h,J}(s)` with `|W| <= binom(s, h-J)`. The critical rate `0.585` is `1.3%` above the Parseval rate `3^{-1/2}`; the data sit at or below the Parseval rate except during the `n = 131` arrival. A proof would be a square-root-cancellation statement for one explicit character sum over the 2-adic digits of `3^{-n}` and is not attempted here.
