---
id: HYP-9288
title: "Energy cone of the signed ROOT transport: the energy of clock_duality_oct15's odd transport of the Collatz ROOT current never exceeds its depth-0 value by more than a logarithmic factor (pointwise: sup_d E_d(N) <= C N log N; averaged over octaves: sup_d T_d(X) <= X^2 (a + b log X), strong form T_d(X) <= T_0(X) = (51/7 + o(1)) X^2); in the termination regime this is a mean-square maximal inequality for the stopping-time shell partial sums of the current; in the ideal regime it is dissipation at every depth (the ideal energy is non-increasing in d)"
status: >
  OPEN (CONJECTURE). Calibrated: the averaged form, even with a logarithmic budget, gives a budget of the same type for
  the census energy (THM-4658), hence the Mercator pole of HYP-9283(B) and limsup log B(X)/log X >= 0.99821 (beyond
  Krasikov-Lagarias). Evidence (FINITE-EXACT, THM-4667 statement 7): strong form (B') holds for every octave
  X = 2^2..2^22 (exception X = 2); the pointwise ratio sup_d E_d(N)/N is at most 10.06 for every N < 2^23 (at N = 7910);
  (C) holds for d <= 6. REVISED 2026-10-17 before any test of the revised forms: the first draft asked for an absolute
  pointwise constant (8.33 from 320 sampled N). The complete scan gives 10.06, and the random-merge model predicts
  log log growth of sup_N R(N)^2/N, so the pointwise form now carries a logarithmic factor.
session: opus-2026-10-17-S41
source: 05-knowledge/results/root_transport_energy_20261017.md
scripts:
  - 04-computation/experiments/root_transport_incremental_trace_20261017.py (every N < 2^23), root_transport_energy_profile_20261017.py, root_transport_octave_trace_20261017.py
  - 04-computation/experiments/root_transport_depth_ladder_20261017.py, root_transport_sibling_anatomy_20261017.py
related: [THM-4667, THM-4658, HYP-9283, HYP-9284, HYP-9285, THM-4659, THM-4660]
---

# HYP-9288 — energy cone of the signed ROOT transport

## Statement

The notation is that of THM-4667:
- `c_0` are the initial seed coefficients of the ROOT current `R(N) = H(4N) - H(N)`, and `T` is the odd transport;
- `E_d(N) = ||T^d c_0||^2` is the energy, and `T_d(X) = Σ_(X <= N < 2X) E_d(N)` is the octave trace;
- `Vroot_d(N) = c_d(1)` is the ROOT coefficient, the current of the sources with at most `d` odd transports;
- `E_d^ideal` is the ideal energy ladder.

The conjecture has three forms.

- **(A, pointwise)** `sup_d E_d(N) <= C N log(N + 2)` for all `N`. This gives `|R(N)| <= sqrt(C N log N)`.
- **(B, averaged)** `sup_d T_d(X) <= X^2 (a + b log_2 X)` for all dyadic `X`, with `b` of order at most `10^3`.
  - **Strong form (B').** `T_d(X) <= T_0(X)` for all `X >= 4` and all `d`; here `T_0(X) = (51/7 + o(1)) X^2`.
  - **Shell form.** By THM-4667 statement 9, (B) contains the mean-square maximal inequality
    `sup_d Σ_N Vroot_d(N)^2 <= X^2 (a + b log_2 X)` for the stopping-time shell partial sums of the current.
- **(C, ideal dissipation)** `E_d^ideal` is non-increasing in `d`. Weaker: the ideal cross terms are summable. Equivalently,
  the ideal sibling covariance (THM-4667 statement 8, in the ideal limit) has non-negative, or summable, aggregate.

## Why these forms

**The model.** In the random-merge model (independent fair relative signs at every merge) the energy is a martingale.
Doob's inequality then gives (A) and (B) in probability.

**Why the logarithmic factors.** Absolute constants are not expected. `R(N)^2/N` has chi-square-like fluctuations with
about one effective value per octave, so `sup_N R(N)^2/N` grows like `log log N`. The logarithmic factors cost nothing
for the calibration: HYP-9283(B) certifies the pole under logarithmic budgets.

**The actual transport.**
- It is dissipative in the ideal regime: negative cross terms at every computed depth, carried by sibling
  anti-correlation.
- For large currents its late energy is the squared ROOT coefficient.
- So (B) is a maximal inequality for the shell partial sums, a census-type statement whose endpoint is the current-energy
  budget.

## Evidence

| octave `X` | `2^9` | `2^10` | `2^11` | `2^12` | `2^13` | `2^14` | `2^15` | `2^16` | `2^17` | `2^18` | `2^19` | `2^20` | `2^21` | `2^22` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `D_X/X^2` (= final `T`) | 1.02 | 1.06 | 1.95 | 3.98 | 2.95 | 0.66 | 2.04 | 0.78 | 0.99 | 1.60 | 0.84 | 0.21 | 0.23 | 4.54 |
| `max_N sup_d E_d/N` | 6.41 | 4.87 | 5.82 | 10.06 | 9.76 | 5.48 | 6.25 | 4.89 | 5.74 | 4.86 | 4.86 | 4.86 | 4.86 | 6.27 |

- In every octave the trace maximum is `T_0/X^2 ≈ 7.2857`, at `d = 0`.
- The largest final value relative to `T_0` is 0.62, at `X = 2^22`.
- **The endpoint `d = ∞` of (B') holds further out.** S40's exact current gives `D_X <= 4.54 X^2` for every octave
  `2^8..2^27`; the maximum is at `2^22`, the same octave as here.

## Tests

* **Larger N.** Complete incremental scans beyond `2^23`, whose memory is about 10 nodes per source. Octave `2^22` came
  closest to (B'), so `2^23..2^26` are the natural next octaves.
* **Ideal dissipation for all `d`.** Prove (C). `d = 7` is being computed. The cross-term ratio is about 0.6, which
  suggests a spectral gap of the two-path correlation operator: compare S20's THM-4565 and the band-sliced Syracuse
  characteristic function of THM-4667 statement 10.
* **A Lyapunov weight on the cone of actual vectors.** The ideal transport is not a contraction even on balanced vectors
  (norm 1.28-1.35), and pure power weights fail. A weight must use the cone generated from `c_0`.
