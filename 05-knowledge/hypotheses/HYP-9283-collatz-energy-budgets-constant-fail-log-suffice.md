---
id: HYP-9283
title: "Energy budgets for the Collatz mod-3 discrepancy H: for dyadic X >= 2^37, any constant budget E_X <= a X^2 with a <= 1806.6, or logarithmic budget E_X <= X^2 (a + b log_2(X/2^37)) inside the certified region (e.g. a = 2, b <= 1194.9), certifies the nonzero Mercator residue of F(s) = sum t(n) n^(-s) and forces limsup log B(X)/log X >= 0.99821 > 0.84 (Krasikov-Lagarias); the stationary-noise forecast that constant budgets fail eventually was refuted out of sample"
status: >
  UPDATE 2026-10-16 (opus-2026-10-16-S38): new exact values e_X = 0.139 at 2^40 and 1.142 at 2^41 (single run, integrity
  checked; replication pre-registered as R5): the S38 pre-registered Q3 (e_X <= 2) PASSED and Q4 (e_X <= 1) FAILED.
  THM-4660 shows that the much weaker hypothesis (W), E_X = O(X^(2+eps)), already gives the continuation and the
  Riesz-mean explicit formula. (D) The certificate at the 1054/665 zero rho* = 0.99999921 + 3014.02i (Codex's 332.5
  line) is not certified at N = 2^42 (interval [-0.00304, 0.0230] under (2, 0)); its radius scales like N^(-1/2), and
  N = 2^43 is pre-registered (R4, expected margin about 8%). A certified nonzero residue there would give
  limsup log B(X)/log X >= 0.99999921 under the budget beyond 2^43.
  (A) REFUTED as a forecast (2026-10-15, same session): the Gaussian AR surrogate fitted on 2^14..2^34 predicted a
  median maximum window energy 0.77 for window starts in (2^33, 2^36]; the exact census to 2^37 gives 0.117
  (predictive probability of a maximum this small 0.1%); the claim "constant budgets are expected to fail" is
  withdrawn and the question is OPEN in both directions. (B) CONDITIONAL, rigorously interval-checked (Arb) at
  N = 2^34, 2^37 and 2^40; the intervals are nested (a code check implied by the data in between, not a test of the budget; second referee
  E1); certified constant budget up to
  14348 X^2 for dyadic X >= 2^40.
  (C) PROVED implication. FINITE-EXACT: E_X/X^2 <= 1.769 for every dyadic 4 <= X <= 2^39 (maximum at 2^24; e = 3.25 at X = 2);
  E_X/X^2 <= 0.096 for X = 2^33..2^37, then 0.931 at 2^38 and 0.078 at 2^39. The budgets themselves are OPEN. AUDITED (S37 referee; I8, I9, I10 applied).
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_octave_memory_analysis_20261015.py (in-sample surrogate)
  - 04-computation/experiments/collatz_census_extension_analysis_20261015.py (out-of-sample check of the surrogate)
  - 04-computation/experiments/collatz_pole_certificate_logenergy_20261015.py (N = 2^34 intervals)
  - 04-computation/experiments/collatz_pole_certificate_general_N_20261015.py, collatz_boundary_moments_2p37_20261015.py (N = 2^37 intervals)
  - 04-computation/experiments/collatz_pole_certificate_margins_20261015.py (largest certified budgets)
related: [THM-4658, HYP-9282, HYP-9284, collatz_amplitude_boundary_20261013c.md, collatz_clock_mesh_20261015a.md, collatz_terminal_mesh_20261015a.md, collatz_annular_overlap_20261015c.md]
---

# HYP-9283: energy budgets - what certifies the pole, and what the data say

**Setting.** `H(N) = Σ_(n<=N) w_n t(n)` and `E_X = Σ_(X<=n<2X) H(n)²`. Theorem B3 of the amplitude note (renewal_shapes_oct11) certifies a nonzero simple pole of `F(s) = Σ t(n) n^(−s)` at the Mercator zero `s_0 ≈ 0.99821 + 240.248i` of `D(s) = 1 + 2^(−s) + (3/2)^s/3`, under either of two hypotheses:
- (7) `|H(x)| <= sqrt(x) log_2(x)/4` for `x >= 2^34`;
- (9a) `E_X <= 2X²` for every dyadic `X >= 2^34`.

The clock-mesh note (clock_duality_oct15) proves the fine-scale part of (9a) and reduces it to a coarse-mean budget `C(P_b) <= (4/9)X²`.

## (A) The forecast that constant budgets fail (REFUTED out of sample)

**What was claimed (in-sample, 2^14..2^34).** The normalized profile `g(u) = H(2^u)/sqrt(2^(u+1))` looked like a short-memory stationary process: mean +0.30, sd 0.28, point autocorrelation 0.06 at one octave. A Gaussian AR(8) surrogate at quarter-octave resolution, with the bias included and refinement energies resampled, then gave P(some dyadic `E_X > 2X²` before `2^K`) = 5.7%, 16%, 31%, 57%, 88% for `K` = 40, 50, 64, 100, 200. Since a stationary non-degenerate Gaussian process has infinite supremum, every constant budget would eventually fail in that model.

**Fragility (referee, I8).** The probabilities depended on the persistent bias (1.0% instead of 32% before 2^64 with the bias removed), on the innovation scale (8%/58% at ×0.8/×1.2), and on the fit range (40% fitted on `u >= 20`); 77 quarter-octave samples fixed 9 parameters, and the kurtosis 3.20 has standard error about 0.95.

**What happened out of sample (FINITE-EXACT census to 2^37).**
- Exact dyadic energies `E_X/X²`: 0.095, 0.060, 0.080, 0.037 for `X` = 2^33, 2^34, 2^35, 2^36, the four smallest of the whole record from 2^20 (maximum 1.769 at 2^24).
- Window energies for starts in `(2^33, 2^36]`: mean 0.072 (in-sample mean 0.51), maximum 0.117. The refitted surrogate predicts a median maximum of 0.77 (10%-90%: 0.34-1.64); a maximum as small as observed has predictive probability 0.1%.
- The profile mean fell from +0.30 to +0.08 and its sd from 0.28 to 0.12.
- The bias on which the forecast rested is gone (HYP-9282 REFUTED), and the normalized wavelet energies at levels 2-5 declined over 2^28..2^37 (HYP-9284).

**Status.** The stationary model is wrong, in the direction of *smaller* energies. Whether a constant budget holds for all `X` is OPEN in both directions. Note that the matched halving null (random odd signs, `t(2m) = −t(m)`) has mean window energies 1.5-14 on 2^14..2^34 (above 2 in 11 of 12 replicates, maxima up to 43), so any proof of a budget must use the odd-step relation `t(n) = −t((3n+1)/2)`, not only the halving structure.

## (B) Which budgets certify the pole (CONDITIONAL, interval-certified)

The finite numerator of B3 is unchanged; only the future tail is recomputed under

`(9b)  E_X <= X² (a + b j)` for every dyadic `X = N 2^j`, `j >= 0`,

by Cauchy–Schwarz on each dyadic block, with the `j`-sum exact in Arb for `j < 400` and bounded rigorously beyond. The normalized residue `2|r| 2^(29(σ_0 − 1))` lies in:

| `(a, b)` | `(2, 0)` | `(4, 1)` | `(8, 2)` | `(16, 4)` | `(32, 8)` | `(64, 16)` |
|---|---|---|---|---|---|---|
| `N = 2^34` | [0.3215, 0.3885] | [0.2964, 0.4137] | [0.2722, 0.4379] | [0.2379, 0.4721] | [0.1895, 0.5205] | [0.1210, 0.5890] |
| `N = 2^37` | [0.3431, 0.3668] | [0.3342, 0.3757] | [0.3256, 0.3843] | [0.3134, 0.3965] | [0.2963, 0.4136] | [0.2720, 0.4379] |
| `N = 2^40` | [0.3508, 0.3591] | [0.3476, 0.3623] | [0.3445, 0.3654] | [0.3402, 0.3697] | [0.3341, 0.3758] | [0.3255, 0.3844] |

- `(2, 0)` at `N = 2^34` reproduces the published constant-budget interval. The `N = 2^37` row uses exact signed block moments through 2^37 (`collatz_boundary_moments_2p37_20261015.py`; the first 2^17 rows reproduce the stored 2^34 moments exactly; `S(2^37) = −152262364`, `H(2^37) = −67763`, `t(2^37) = −1`).
- **Nesting is not a test.** The rows are nested for every budget, but this cannot fail: by Lemma B2 the finite numerators at two census sizes differ by the exact data in between, and the census energies bound that difference (at most 1.87e-4 between 2^34 and 2^37, actually 3.75e-6) far below the difference of the radii (1.076e-3). The nesting checks the code; it is no evidence for (9a) beyond 2^37 (second referee, E1).
- **Margins** (bisection, every test in Arb). At `N = 2^34` (agreeing with the referee's independent values): largest constant budget `a = 227.385`, largest slope with `a = 2` `b = 144.393`, largest `(4c, c)` `c = 36.860`. At `N = 2^37`: `a = 1806.6`, `b = 1194.9`, `c = 292.86`. At `N = 2^40` (exact block moments through 2^40 from the streaming census; `S(2^40) = 989143994`, `H(2^40) = −55691`, `t(2^40) = +1`): `a = 14348`, `b = 9616`, `c = 2325.9`. The certified region grows by about `2^(3(2σ_0−1)) ≈ 8` per three octaves of exact data.
- In the clock-mesh formulation a coarse-mean budget `C(P_b) <= X²(a' + b' j)` with the proved fine allowance gives `E_X/X² <= 4(a' + b' j) + 4/27 + O(X^(−1/2))`; for example `(15, 4)` lies inside the certified `(64, 16)` row (minimum margin 3.85, referee).
- **Credit.** The amplitude note (sections 3 and 6.1) already proposed logarithmic budgets; the new contributions here are the interval certification, the margins, the `N = 2^37` certificate and the calibration (C). clock_duality_oct15 independently proves a sufficient coarse budget `k²/100` by comparison with B3's logarithmic pointwise allowance and propagates a linear-log budget from paid annular means (`collatz_gram_slack_20261015b.md`, `collatz_annular_overlap_20261015c.md`).

## (C) Calibration (PROVED implication)

Assume any budget in (B) at `N = 2^34` or `N = 2^37`, or (7). Then:
- `F(s)` continues meromorphically to `Re s > 1/2` and has a nonzero pole at `s_0`;
- `S(X) = Σ_(n<=X) t(n)` is not `O(X^β)` for any `β < Re s_0 = 0.99821`.

Since `|S(X)| <= B(X) = #{n <= X : n reaches 1}`, this gives `limsup log B(X)/log X >= 0.99821`. The best unconditional lower bound valid for every large `X` is `B(X) >= X^0.84` (Krasikov–Lagarias, Acta Arith. 109 (2003)). So proving any of these budgets, even the weakest certified one (about 1800 times the largest energy ever observed, beyond 2^37), would beat the known population exponent along a sequence of scales. (The weakest certified budgets are about 1000 times the largest energy ever observed for `X >= 4` beyond 2^37, and about 8000 times beyond 2^40.) The remaining coarse-mean obligation is a major open problem, not a technical gap.

## Tests

- The pre-registered 2^37..2^40 census (`collatz_census_2p40_preregistration_20261015.md`, P2): `E_X/X² <= 2` for `X = 2^37, 2^38, 2^39`. PASSED (0.069, 0.931, 0.078).
- Any proof strategy for the coarse means must fail for the halving null, whose window energies violate (9a) at most scales, so it must use the odd-step relation.
