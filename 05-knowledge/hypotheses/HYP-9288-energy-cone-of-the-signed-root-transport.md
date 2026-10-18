---
id: HYP-9288
title: "Energy cone of the signed ROOT transport: the energy of clock_duality_oct15's odd transport of the Collatz ROOT current never exceeds its depth-0 value by more than a logarithmic factor (pointwise: sup_d E_d(N) <= C N log N; averaged over octaves: sup_d T_d(X) <= X^2 (a + b log_2(X/2^37)), strong form T_d(X) <= T_0(X) = (51/7 + o(1)) X^2); this implies a uniform mean-square bound for the stopping-time shell partial sums of the current; in the ideal regime it is dissipation at every depth (the ideal energy is non-increasing in d)"
status: >
  OPEN (CONJECTURE). Calibrated: the averaged form, with a budget inside HYP-9283(B)'s certified region, gives a budget
  of the same type for the census energy (THM-4658). Hence the Mercator pole of HYP-9283(B) and
  limsup log B(X)/log X >= 0.99821 along a sequence of scales (beyond Krasikov-Lagarias). In the form
  a + b log_2 X only b <= 45.7 (2^37) or b <= 338 (2^40) is certified (S41 audit a7). Evidence (FINITE-EXACT, THM-4667 statement 7): strong form (B') holds for every octave
  X = 2^2..2^23 (exception X = 2); the pointwise ratio sup_d E_d(N)/N is at most 10.06 for every N < 2^24 (at N = 7910);
  (C) holds for d <= 12, with cross-term ratios rising to 0.83-0.90 (slower than geometric). REVISED 2026-10-17 before any test of the revised forms: the first draft asked for an absolute
  pointwise constant (8.33 from 320 sampled N). The complete scan gives 10.06, and the random-merge model predicts
  log log growth of sup_N R(N)^2/N, so the pointwise form now carries a logarithmic factor.
  Scope repaired by clock_duality_oct15: ROOT projection retains the unrooted complement; the liminf current-budget
  implication is unconditional, but a shell-only bound is weaker than total energy. Ideal dissipation includes both
  INTRA and INTER; ordinary Syracuse comparisons retain boundary kernels.
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

- **(A, pointwise)** `sup_d E_d(N) <= C N log(N + 2)` for all `N`. This gives `|R(N)| <= sqrt(C N log(N+2))` without assuming universal convergence.
- **(B, averaged)** `sup_d T_d(X) <= X^2 (a + b log_2(X/2^37))` for all dyadic `X >= 2^37`, with `(a, b)` inside
  HYP-9283(B)'s certified region. For example a constant `a <= 1806.6`, or `(2, <= 1194.9)`.
  - **Strong form (B').** `T_d(X) <= T_0(X)` for all `X >= 4` and all `d`; here `T_0(X) = (51/7 + o(1)) X^2`.
  - **Uniform shell bound.** By THM-4667 statement 9, (B) implies `sup_d Σ_N Vroot_d(N)^2 <= X^2 (a + b log_2(X/2^37))`
    for the stopping-time shell partial sums.
  - **Shell maximal inequality.** The stronger `Σ_N sup_d Vroot_d(N)^2 <= ...`, with the maximum inside the sum, is not
    implied by (B). It is observed in every octave to `2^22`.
- **(C, ideal dissipation)** `E_d^ideal` is non-increasing in d, equivalently `INTRA_d+INTER_d<=0` at every ideal depth.
  A weaker sufficient target for bounded ideal energy is summability of the positive parts of these total increments.
  Sibling covariance alone is not equivalent: the between-block term INTER must also be paid.

By repaired THM-4667 statements4--5, the ROOT coefficient eventually equals R(N) for every finite initial support,
even if unrooted sources exist. Hence `D_X<=liminf_d T_d(X)`. This is the unconditional direction needed for (B).
The shell form is necessary for (B), but is not proved sufficient for it.

## Why these forms

**The model.** In the random-merge model (independent fair relative signs at every merge) the energy is a martingale.
Doob gives `Pr(sup_d E_d>=lambda E_0)<=1/lambda`. This is a tail inequality for that model, not a deterministic
proof of (A) or (B).

**Why the logarithmic factors (HEURISTIC model).** Absolute constants are not expected in this model. Its analogue of `R(N)^2/N` has chi-square-like fluctuations with
about one effective value per octave, suggesting log log growth of its running maximum. No such law is proved for the actual ROOT current. The logarithmic factors cost nothing
for the calibration, provided the budget is written as in HYP-9283(B), `a + b log_2(X/2^37)`.

**The actual transport.**
- It is dissipative in the ideal regime: negative cross terms at every computed depth, carried by sibling
  anti-correlation.
- The tested large-current rows are ROOT-dominated late; this is finite evidence.
- (B) implies the uniform shell bound and the current-energy budget. Its converse additionally needs uniform control
  of nonterminal coefficient energy. The full ideal bound remains open. Matching is exact except at ROOT and `O(1)` small
  seeds for `d <= 10` (THM-4667 statement 7).

## Evidence

| octave `X` | `2^9` | `2^10` | `2^11` | `2^12` | `2^13` | `2^14` | `2^15` | `2^16` | `2^17` | `2^18` | `2^19` | `2^20` | `2^21` | `2^22` | `2^23` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `D_X/X^2` (= final `T`) | 1.02 | 1.06 | 1.95 | 3.98 | 2.95 | 0.66 | 2.04 | 0.78 | 0.99 | 1.60 | 0.84 | 0.21 | 0.23 | 4.54 | 1.51 |
| `max_N sup_d E_d/N` | 6.41 | 4.87 | 5.82 | 10.06 | 9.76 | 5.48 | 6.25 | 4.89 | 5.74 | 4.86 | 4.86 | 4.86 | 4.86 | 6.27 | 4.86 |

- In every octave the trace maximum is `T_0/X^2 ≈ 7.2857`, at `d = 0`.
- The largest final value relative to `T_0` is 0.62, at `X = 2^22`.
- **The shell form holds with the maximum inside the sum.** `Σ_N max_d Vroot_d(N)^2 <= 5.20 X^2 < T_0` for every
  octave `2^9..2^22`. The ratio to `D_X` ranges over 1.07-5.18; it exceeds the martingale constant 4 only where
  `D_X` is tiny, i.e. where the final cancellation is strong.
- **The endpoint `d = ∞` of (B') holds further out.** S40's exact current gives `D_X <= 4.54 X^2` for every octave
  `2^8..2^27`; the maximum is at `2^22`, the same octave as here.

## Tests

* **Larger N.** Complete incremental scans beyond `2^23`, whose memory is about 10 nodes per source. Octave `2^22` came
  closest to (B'), so `2^23..2^26` are the natural next octaves.
* **Ideal dissipation for all `d`.** Prove (C). It is verified for `d <= 12`: every total increment INTRA + INTER is
  negative.
  - The cross-term ratios rise, from 0.51 at `d = 4` to 0.83-0.90 at `d = 10..12`, so the decay is slower than
    geometric. A first reading of a spectral gap from `d <= 6` was withdrawn.
  - The natural objects are the joint band/frequency law of the primitive affine characteristic function (repaired
    THM-4667 statement 10). This preserves the A, C boundary operators and the band of clock_duality's equation (5).
    S20's two-path correlations (THM-4565) are the other natural object.
* **A Lyapunov weight on the cone of actual vectors.** The ideal transport is not a contraction even on balanced vectors
  (norm 1.28-1.35), and pure power weights fail. A weight must use the cone generated from `c_0`.
