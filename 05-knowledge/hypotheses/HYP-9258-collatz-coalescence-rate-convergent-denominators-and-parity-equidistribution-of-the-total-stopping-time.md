---
id: HYP-9258
title: "Collatz coalescence along lags: (i) the density of n <= N with tau(n) != tau(n+1) decays like a power of 1/log N (THM-4581's T^(-1/2) with T of order log N) and the mean run of consecutive integers with equal total stopping time grows accordingly; (ii) the nonzero differences tau(n+h) - tau(n) are small integer combinations 19a + 8b of the two near-return lengths 8 (3^5/2^8) and 19 (3^12/2^19) of log_2 3, the 19-direction favoured over the 8-direction (19, 38, 27, 46, 8, 57, 65, 84, 11, 30 the ten commonest at every lag h <= 4; 84% of the differences below the Frobenius number 125 lie in the semigroup a, b >= 0); (iii) the parity of tau is log-periodically biased: the dyadic-block means of (-1)^tau alternate in sign (forced by tau(2n) = tau(n) + 1) with a non-decaying amplitude 0.007 over 2^25..2^29 and a beat node near 2^21.5, the 1/16-octave means have rms 0.044 at 2^28 (down from 0.082 at 2^18) with one dominant mode at 2.455 cycles per octave, near half the drift 4.82 of tau per octave (or the doubling mode 5/2; unresolved at resolution 0.09), so E[(-1)^tau] over [1, N] has not converged at the scales computed even though the two-point correlations of (-1)^tau tend to +1 at every fixed lag (THM-4627)"
status: >
  OPEN. PROVED context: THM-4627 (break density -> 0; correlations -> 1); the identity tau(2n) = tau(n) + 1, which forces the even
  half of each dyadic block to carry the negated parity mean of the previous block (E_j = -B_(j-1), verified exactly). NUMERICAL:
  break density 0.555, 0.537, 0.522, 0.510, 0.500, 0.490, 0.481 at N = 2^16, 2^18, ..., 2^28; parity flips exactly half the breaks
  (0.5003); mean run length 1.80 -> 2.08; longest run 120 below 2^28; the fitted constant in c (ln N)^(-1/2) drifts from 1.85
  (2^16) to 2.12 (2^28), so (i) is stated as a power of 1/log N with the exponent unknown. FINITE-EXACT: the difference censuses
  to 2^27 at lags 1..4 (commonest |d|: 19, 38, 27, 46, 8, 57, 65, 84, 11, 30, 76, 16, 73, 54; the share of |d| <= 125 inside the
  semigroup <8,19>_{>=0} is 0.8416 at every lag; the commonest values outside it, 11, 30, 3, 49, 14, 22, 68, are 19a - 8b with
  a <= 4, b <= 3). NUMERICAL (parity): dyadic-block means B_j of (-1)^tau alternate in sign with |B_j| = 0.0056, 0.0065, 0.0070,
  0.0070, 0.0068 at j = 25..29 after a beat node near j = 21.5 (|B_j| = 0.020 at j = 16); the odd part O_j alternates with the same
  amplitude; 1/16-octave means have rms 0.082, 0.072, 0.065, 0.059, 0.055, 0.051, 0.049, 0.047, 0.045, 0.044, 0.044 per octave from
  2^18 to 2^29; least-squares spectrum in log_2 n: dominant mode f = 2.455 cycles per octave (amplitude 0.042) with side structure
  at 2.54 and 2.375, next 1.52; the drift of tau is 4.82 steps per octave (measured 4.76-4.83 from 2^19 to 2^29), so f is at half
  the drift to within the resolution 0.09; the variance of tau is 1106 (2^18) to 1867 (2^29), so no Gaussian fluctuation model
  can produce a parity bias of this size, so it comes from the tree relations (inference, not a proof of a mechanism). Cumulative ordinary means
  oscillate between -0.0024 and +0.0023 at 2^27..2^29 without converging; the logarithmic mean +0.033, +0.028, +0.023 at 2^20,
  2^24, 2^29 is the fixed initial segment diluted as 1/log N. HEURISTIC: the mechanism of (ii), a segment of q steps with p odd
  steps, p/q a convergent of log_2 3, returns to nearly the same size (3^5/2^8 = 0.949, 3^12/2^19 = 1.014) and can land on the
  other orbit at time offset q; and of (iii), the parity of round(4.82 x) at x = log_2 n, which the tree relations keep coherent
  across a bin although tau's fluctuations are large.
session: opus-2026-10-13-S25
source: 05-knowledge/results/five_manuscripts_chowla_collatz_residue_bsd_20261013.md (sections 3.3-3.4)
scripts:
  - 04-computation/experiments/collatz_tau_runs_20261013.py (+ .out)
  - 04-computation/experiments/collatz_merge_lags_20261013.py (+ .out)
  - 04-computation/experiments/collatz_tau_diffs_20261013.py, collatz_tau_oscillation_20261013.py, collatz_tau_spectrum_20261013.py, collatz_tau_drift_20261013.py (+ .out)
related:
  - THM-4627 (the density-one coalescence and the correlation limits)
  - THM-4581 (P(no merge by T) between c T^(-1/2) and C T^(-1/2) (log T)^2, both bounds at sketch level)
  - HYP-9217 (debt recurrence law, exponent 1/2)
  - HYP-9230 (class-density resonance spectrum; the natural home of a log-periodic mechanism, not yet connected)
  - Garner, On heights in the Collatz 3n+1 problem, Discrete Math. 55 (1985) 57-64 (consecutive integers with equal total stopping time; CITED via the search of 2026-10-13); LaDue, Clusters of integers with equal total stopping times in the 3x+1 problem, Fibonacci Quart. 56 (2018), arXiv:1709.02979 (merge criterion for consecutive even-odd pairs; the identity sigma_inf(8i+4) = sigma_inf(8i+5), which is the depth-3 mass 1/8 of THM-4627); Elia and Tucker, Consecutive integers and the Collatz conjecture, arXiv:1511.09141 (counterexamples to Garner's conjectural family); Lagarias's survey attributes the clustering to coalescence of trajectories; no source found gives the density of n with tau(n) = tau(n+1), nor a parity bias
---

# HYP-9258 — rates and residues of Collatz coalescence along lags

## Statements

(i) **Rate.** `beta(N) := (1/N) #{n <= N : tau(n) != tau(n+1)}` tends to `0` (THM-4627) and `beta(N) = N^(o(1))`-slowly: conjecturally `beta(N) ~ c (log N)^(-gamma)` with `gamma` between `1/2` (THM-4581's merge law with `T` of order `log N`) and the observed slower decay; the mean length of a maximal run of consecutive integers with equal `tau` is `1/beta(N)`.
(ii) **Residues.** Among `n` with `tau(n+h) != tau(n)`, the difference is `19a + 8b` with small integers `a, b` (both signs), the `19`-direction favoured over the `8`-direction and the frequency falling with the size of the coefficients: the ten commonest absolute values at every lag `h <= 4` are `19, 38, 27, 46, 8, 57, 65, 84, 11, 30`, i.e. `(a,b) = (1,0), (2,0), (1,1), (2,1), (0,1), (3,0), (3,1), (4,1), (1,-1), (2,-1)`; below the Frobenius number `125` of the semigroup `<8, 19>`, `84.2%` of the differences lie in the semigroup `a, b >= 0`, and the commonest exceptions `11, 30, 3, 49, 14, 22, 68` are `19a - 8b` with `a <= 4`, `b <= 3`. The general law should be that unequal-time coincidences occur at time offsets that are signed combinations of the convergent and semiconvergent denominators of `log_2 3`, the shortest first.
(iii) **Parity bias.** `E[(-1)^tau]` over `[1, N]` has not converged at the scales computed (`|E| <= 0.0024`, oscillating): the parity of `tau` carries a log-periodic bias in `x = log_2 n` with one dominant mode at `f = 2.455` cycles per octave, near half the drift `4.82` of `tau` per octave (the alternative, the doubling mode at `5/2`, is not excluded at the resolution `0.09`), and a slowly decaying envelope (rms `0.082 -> 0.044` per octave from `2^18` to `2^29`); the dyadic-block means alternate in sign, as the identity `tau(2n) = tau(n) + 1` forces for the even half, with an amplitude `0.007` that does not decay from `2^25` to `2^29` and a beat node near `2^21.5`. Conjecture: the envelope decays to zero (so the logarithmic mean of `(-1)^tau` is `0` and the ordinary mean is `o(1)` only along subsequences), the dominant frequency is exactly half the drift, and the beat is the interference of the drift mode `s/2` with the doubling mode at the half-integer `5/2`.

## Cheapest decisive tests

- (i): compute the un-merged mass of the `(y, y+1)` pair chain at depths `K = 26, ..., 40` exactly (dyadic rationals) and fit the tail; compare with `beta` at `N = 2^K`.
- (ii): census of `tau(n+h) - tau(n)` for `h = 2, 3, 4` and of the time offsets of first unequal-time coincidences; test whether `65` (`3^41/2^65`) and `84` appear next.
- (iii): 1/16-octave means of `(-1)^tau` to `2^34` (int16 storage, 16 GB, or streamed); fit the envelope exponent and resolve the two candidate frequencies `2.41` (half drift) and `2.5` (doubling) by the beat period; connect to HYP-9230's resonance modes.

## Why

THM-4627 makes `(-1)^tau` locally constant along lags on a density-one set; the three statements quantify the complement: how fast the exceptional density decays, what the exceptional differences are, and how the global average, which is not locally random, oscillates. The parity bias is the first one-point statistic in this thread that is visibly deterministic at the scale of `2^29`: a Gaussian model of `tau` with variance `1867` would make `E[(-1)^tau]` smaller than `10^(-1000)`.
