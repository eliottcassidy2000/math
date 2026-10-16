---
id: HYP-9283
title: "Constant dyadic energy budgets fail eventually, logarithmic budgets suffice: for the Collatz mod-3 discrepancy H, the normalized octave profile H(N)/sqrt(N) behaves like a short-memory stationary process in log N (the exact two-octave memory of THM-4658 is almost cancelled), so every constant budget (E_X <= 2X^2, the clock-mesh coarse target C <= (4/9)X^2, |H| <= 1.5 sqrt(x)) is expected to fail at some large X; meanwhile the nonzero Mercator residue is interval-certified under budgets E_X <= X^2 (a + b log_2(X/2^34)) up to (64, 16), and any such budget forces limsup log B(X)/log X >= 0.99821 > 0.84"
status: >
  Part (A) HEURISTIC with NUMERICAL support: a Gaussian AR surrogate fitted to the exact census 2^14..2^34 gives
  P(some dyadic E_X > 2X^2 with X in [2^34, 2^K]) = 5.7%, 16%, 31%, 57%, 88% for K = 40, 50, 64, 100, 200.
  Part (B) CONDITIONAL, rigorously interval-checked: the log-budget implication extends Theorem B3 of the amplitude note.
  Part (C) PROVED implication (from B3's continuation argument): any such budget implies the population bound.
  The energy budgets themselves are OPEN.
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_octave_memory_census_20261015.py
  - 04-computation/experiments/collatz_octave_memory_analysis_20261015.py (output collatz_octave_memory_analysis_20261015.json/.out)
  - 04-computation/experiments/collatz_pole_certificate_logenergy_20261015.py (output collatz_pole_certificate_logenergy_20261015.json)
related: [THM-4658, HYP-9282, collatz_amplitude_boundary_20261013c.md, collatz_clock_mesh_20261015a.md, collatz_terminal_mesh_20261015a.md]
---

# HYP-9283: constant energy budgets fail, logarithmic budgets suffice

**Setting.** `H(N) = Σ_(n<=N) w_n t(n)` and `E_X = Σ_(X<=n<2X) H(n)²`. Theorem B3 of the amplitude note certifies a nonzero simple pole of `F(s) = Σ t(n) n^(−s)` at the Mercator zero `s_0 ≈ 0.99821 + 240.248i` of `D(s)`, under either of two hypotheses:
- (7) `|H(x)| <= sqrt(x) log_2(x)/4`;
- (9a) `E_X <= 2X²` for every dyadic `X >= 2^34`.

The clock-mesh note proves the fine-scale part of (9a) and reduces it to a coarse-mean budget `C(P_b) <= (4/9)X²`.

## (A) Constant budgets are expected to fail eventually (HEURISTIC)

From the exact census (`2^34`; THM-4658, HYP-9282):
- **Gaussian-like marginal.** The normalized profile `g(u) = H(2^u)/sqrt(2^(u+1))` has a Gaussian-like marginal (kurtosis 3.20), with mean `+0.30` (the bias of HYP-9282) and standard deviation 0.28.
- **Short memory.** Its autocorrelation is 0.69 at 1/4 octave, 0.39 at 1/2, 0.06 at 1, and 0.01 at 2 octaves.
- **The memory is cancelled.** THM-4658 gives `H(4N) = H(N) + Δ(N)`, yet the regression of `g(u+2)` on `g(u)` has slope 0.006. The level-2 increment `Δ` almost exactly cancels the inherited value, with octave correlations `κ = Corr(H(N), Δ(N))` between −0.97 and +0.91.
- **Spectrum.** 92% of the spectral power of `g` lies below one cycle per octave. There is no line; the density observable `S(N)/N` shows the Mercator line at 26.5 cycles per octave as a positive control.
- **Energy split.** Window energies `E/X²` over 1,217 windows have mean 0.51, 99% quantile 1.77 and maximum 1.85. The coarse mean carries 84% of the energy and refinement 16%.

A stationary non-degenerate Gaussian process with decaying correlations has an infinite supremum. So in this model every constant budget is eventually exceeded. A Gaussian AR(8) surrogate at quarter-octave resolution, fitted on `2^14..2^34` with the bias included and refinement energies resampled, gives, from `2^34` to `2^K`:

| `K` | 40 | 50 | 64 | 100 | 200 |
|---|---:|---:|---:|---:|---:|
| P(some dyadic window has `E_X > 2X²`) | 5.7% | 16.4% | 30.7% | 56.6% | 87.7% |
| P(some window at any phase has `E > 2X²`) | 9.4% | 27.6% | 46.6% | 75.9% | 97.4% |
| median of the maximum window energy | 1.11 | 1.63 | 1.95 | 2.34 | 2.77 |

**Consequence.** Hypothesis (9a), the stronger pointwise `|H| <= 1.5 sqrt(x)`, and the clock-mesh/terminal-mesh target `C <= (4/9)X²` on all dyadic `X >= 2^34` are each expected to fail at some large `X`, although probably not soon. A proof aiming at a constant budget is aiming at a statement that is likely false. The pointwise hypothesis (7) has a logarithmic factor and is consistent with the model.

## (B) Logarithmic budgets suffice (CONDITIONAL, interval-certified)

The finite numerator of B3 is unchanged. The future tail is recomputed under

`(9b)  E_X <= X² (a + b j)` for every `X = 2^(34+j)`, `j >= 0`,

by Cauchy–Schwarz on each dyadic block. The `j`-sum is exact in Arb for `j < 400` and bounded rigorously beyond. The normalized residue `2|r| 2^(29(σ_0 − 1))` lies in:

| `(a, b)` | `(2, 0)` | `(4, 1)` | `(8, 2)` | `(16, 4)` | `(32, 8)` | `(64, 16)` |
|---|---|---|---|---|---|---|
| interval | [0.3215, 0.3885] | [0.2964, 0.4137] | [0.2722, 0.4379] | [0.2379, 0.4721] | [0.1895, 0.5205] | [0.1210, 0.5890] |

All are nonzero. `(2, 0)` reproduces the published constant-budget interval. In the clock-mesh formulation, a coarse-mean budget `C(P_b) <= X²(a' + b' j)` with the proved fine allowance gives `E_X/X² <= 4(a' + b' j) + 4/27 + O(X^(−1/2))`. For example `(a', b') = (15, 4)` lies inside the certified `(64, 16)` row.

**Reading.** The amplitude conclusion does not need a constant budget. A budget about 30 times larger than 2, growing linearly in `log X`, still certifies the pole. In the model of (A), such a budget fails with negligible probability.

## (C) Calibration (PROVED implication)

Assume any budget in (B), or (7). Then:
- `E(s)` continues holomorphically to `Re s > 1/2`;
- `F` has a nonzero pole at `s_0`;
- `S(X) = Σ_(n<=X) t(n)` is not `O(X^β)` for any `β < Re s_0 = 0.99821`.

Since `|S(X)| <= B(X) = #{n <= X : n reaches 1}`, this gives `limsup log B(X)/log X >= 0.99821`. The best unconditional lower bound for every large `X` is `B(X) >= X^0.84` (Krasikov–Lagarias 2003). So proving any of these budgets, even with the logarithmic slack, would beat the known population exponent along a sequence of scales. That calibrates the remaining coarse-mean obligation as a major open problem, not a technical gap.

## Tests

- Census octaves `2^35..2^40`. The model predicts a constant-budget exceedance there with probability about 6% (dyadic windows) and a median maximum window energy of about 1.1.
- Any proof strategy for the coarse means should be checked against the near-cancellation of the two-octave memory (`κ` up to `±0.97`) and against the persistent biases of HYP-9282.
