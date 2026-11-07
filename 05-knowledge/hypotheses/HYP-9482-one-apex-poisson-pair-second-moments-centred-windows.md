---
id: HYP-9482
title: "One-apex collisions on reduced words: Poisson-pair second moments (energy 14 - 16 ln 2 up to a small positive excess, multiplicity law, smoothed slope 12 - 16 ln 2, cross-layer energy 4 - 4 ln 2 - (12 - 16 ln 2)), centred collision windows with |a - s| <= 1 + carry depth, and a halving level law"
status: "PRE-REGISTERED in 21b02648ec; outcomes at n = 22, 24, 26: P1-P5 and P8 CONFIRMED; P6 split - (a) CONFIRMED, (b) REFUTED at the boundary level a = n/2 - 3 (ratio 0.370-0.380; the n = 20 in-sample data already showed 0.372 there, MISTAKE-728), holds for 3 <= a <= n/2 - 4; P7 REFUTED (maximum multiplicity 18 at n = 26; the model itself expects 0.13 such values, the tail is about 25 percent heavier than Poisson-pair). Energy excess over 14 - 16 ln 2 decreasing: 0.032, 0.025, 0.020"
session: opus-2026-11-06-S84
date: 2026-11-06
related: [THM-4890, THM-4887, THM-4886, HYP-9479, HYP-9480, HYP-2199, THM-4885]
script: 05-knowledge/results/one_apex_collision_balance_20261106/s84_energy.py; 05-knowledge/results/one_apex_collision_balance_20261106/s84_model.py
---

# HYP-9482 - Poisson-pair second moments and centred windows for one-apex collisions (pre-registration)

## Objects

- g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2): the one-apex Hamiltonian path counts.
- Y_n = {v even : 2^(n-1) <= v < 2^n} (reduced words: top bit 1, bottom bit 0), |Y_n| = 2^(n-2).
- v* = reversed complement of the n-bit word v; it maps Y_n to Y_n and g(v*) = g(v) (converse tournament).
- E^Y_n = #{(v, v') in Y_n^2 : g(v) = g(v')}; E^Y_n(K) = #{(v, v') : |g(v) - g(v')| <= K}.
- For a collision v != v', v' != v* (non-converse pair):
  - a = common prefix length;
  - s = common suffix length;
  - w = n - a - s, the width of the window where the two words differ;
  - x, x' = the window values;
  - dX = w - bitlength|x - x'| and dXi = w - bitlength|C_w(x) - C_w(x')| (carry depths).
- Level law: L_a(n) = number of non-converse ordered collision pairs with prefix length a, per |Y_n|.
- X_n = all words of length n; X_n splits into layers Y_(n-k) followed by k trailing ones.

## Model (heuristic, extends HYP-9480 (A))

HYP-9480 (A) (kind-pasteur-2026-11-06-S13) records the Poisson law for the one-apex ("root") image: converse classes
are a Poisson process, kappa = 1/2. This gives density 1/2 of the image among odd integers; HYP-2199 observed that
density earlier.

This hypothesis applies the same model to second moments. The limit law of g(v)/2^n on Y_n is XY, with X and Y
independent and uniform on [1/2, 1). Converse classes form a Poisson process of mean mu(t) = f_Y(t)/4 per odd value,
where f_Y(t) = 4 ln(4t) on [1/4, 1/2] and 4 ln(1/t) on [1/2, 1]. The model predicts:
- E^Y_n/|Y_n| -> 14 - 16 ln 2 = 2.90965; the non-converse part is 12 - 16 ln 2 = 0.90965;
- values attained per |Y_n| -> 0.40343;
- the fraction of |Y_n| that is values with exactly 2q preimages: 0.32183, 0.06855, 0.01133, 0.00152 for q = 1..4;
- (E^Y_n(K) - E^Y_n(0))/(K |Y_n|) -> 12 - 16 ln 2 for even K;
- with independent layers, E_n/|X_n| -> 6 - 4 ln 2 = 3.22741, and the cross-layer part -> 0.31777.

In sample (n = 14..20, `insample_14_20.out`) the model fits the histogram, the attained values and the cross layers.
The energy runs above it: E^Y_20/|Y_20| = 2.93166, an excess of 0.022 that was still rising at n = 20.

## Pre-registered predictions for n = 22, 24, 26 (none of these Y_n computed before this commit)

- **P1 (energy band).** 2.925 <= E^Y_n/|Y_n| <= 2.990, and |E^Y_(n+2) - E^Y_n|/|Y| <= 0.02 for n = 22, 24.
- **P2 (multiplicity law).** Values with exactly 2q preimages, per |Y_n|, lie in:
  - q = 1: [0.312, 0.326];
  - q = 2: [0.0655, 0.0715];
  - q = 3: [0.0105, 0.0132];
  - q = 4: [0.0012, 0.0021].
- **P3 (attained values).** values/|Y_n| in [0.398, 0.406].
- **P4 (smoothed slope).** (E^Y_n(K) - E^Y_n(0))/(K |Y_n|) in [0.880, 0.930] for K = 8 and K = 16.
- **P5 (centred windows, sharp form).** Every non-converse collision pair satisfies |a - s| <= 1 + max(dX, dXi).
  - THM-4890 proves the weaker form: |a - s| <= 2 + max(dX, dXi), or a <= 4 + dXi and s <= 4 + dX.
  - The proved form is checked as a code test, not as a prediction.
- **P6 (level law).**
  - L_a(n) is within 8% of the n = 20 values 0.3506, 0.2646, 0.1540, 0.0839, 0.0436, 0.0222 for a = 1..6.
  - L_(a+1)(n)/L_a(n) is in [0.45, 0.60] for 3 <= a <= n/2 - 3.
- **P7 (multiplicity growth).** The maximum multiplicity on Y_n is at most 16.
- **P8 (cross layers on X_n).** Cross-layer energy per |X_n| is in [0.310, 0.335].

## What this would and would not show

Confirmation would mean the following:
- the Poisson-pair model describes one-apex collisions to second order, with a small positive excess;
- the collisions are centred and balanced at micro/macro scale, as proved in weaker form by THM-4890.

It would not prove linear energy. E_n = O(2^n) remains OPEN (THM-4887 proves O(n 2^(4n/3))).

## Outcomes (run after the pre-registration commit 21b02648ec; `outofsample_21_26.out`)

The program was run once, for n = 21..26. Odd n were not registered; they are reported for context.

| n | E^Y/\|Y\| | nonconv/\|Y\| | values/\|Y\| | q=1 | q=2 | q=3 | q=4 | slope K=8 | slope K=16 | max mult | cross-layer/\|X\| | E/\|X\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 22 | 2.94201 | 0.94299 | 0.40126 | 0.31759 | 0.06873 | 0.01199 | 0.00173 | 0.90297 | 0.90266 | 14 | 0.32402 | 3.26072 |
| 24 | 2.93485 | 0.93534 | 0.40143 | 0.31809 | 0.06918 | 0.01182 | 0.00164 | 0.90749 | 0.90622 | 16 | 0.32264 | 3.25743 |
| 26 | 2.92960 | 0.92984 | 0.40194 | 0.31930 | 0.06881 | 0.01173 | 0.00165 | 0.90872 | 0.90798 | 18 | 0.32149 | 3.25257 |
| model | 2.90965 | 0.90965 | 0.40343 | 0.32183 | 0.06855 | 0.01133 | 0.00152 | 0.90965 | 0.90965 | - | 0.31777 | 3.22741 |

- **P1 CONFIRMED.** All three values lie in [2.925, 2.990], and the increments are -0.0072 and -0.0052.
  - E^Y/|Y| peaked at n = 21-22 (2.9418, 2.9420) and has decreased since: 2.9328, 2.9349, 2.9303, 2.9296 for
    n = 23..26.
  - The excess over 14 - 16 ln 2 is 0.032, 0.025, 0.020 at n = 22, 24, 26.
  - The data are consistent with the excess vanishing, which would make the energy constant exactly the Poisson-pair
    value. They do not prove it.
- **P2 CONFIRMED.** All twelve fractions are inside their bands.
- **P3 CONFIRMED.** values/|Y| = 0.40126, 0.40143, 0.40194, against the model's 0.40343.
- **P4 CONFIRMED.** Slopes are 0.903-0.909 at K = 8 and K = 16, approaching 12 - 16 ln 2 = 0.90965.
- **P5 CONFIRMED.** |a - s| <= 1 + max(dX, dXi) had no violations at n = 21..26. The proved form of THM-4890 (C)
  also had no violations (code check).
- **P6 split.**
  - (a) CONFIRMED: the levels 1..6 are within 3.7 % of the n = 20 values at every registered n.
  - (b) REFUTED at the boundary level a = n/2 - 3. The ratios L_(a+1)/L_a are:
    - n = 22, a = 3..8: 0.548, 0.524, 0.517, 0.507, 0.452, 0.380;
    - n = 24, a = 3..9: 0.540, 0.531, 0.518, 0.506, 0.501, 0.465, 0.376;
    - n = 26, a = 3..10: 0.542, 0.525, 0.515, 0.514, 0.505, 0.497, 0.467, 0.370.
  - The halving law holds in [0.45, 0.60] for 3 <= a <= n/2 - 4. The last ratio is a boundary effect: windows of width
    w = n - 2a <= 6.
  - The registered range was wrong from the start. The n = 20 in-sample ratio at a = n/2 - 3 = 7 was already
    0.00401/0.01077 = 0.372, so the range should have stopped at n/2 - 4 (MISTAKE-728).
- **P7 REFUTED.** The maximum multiplicity is 18 at n = 26, from one value, so the bound of 16 fails. The bound was set
  by eye and not from the model's tail.
  - The model expects 0.13 values of multiplicity 18, 1.84 of 16, 24.2 of 14 and 282 of 12 at |Y| = 2^24. The observed
    counts are 1, 3, 30 and 356.
  - So the upper tail is about 25 % heavier than the Poisson-pair mixture. This is the over-dispersion behind the small
    energy excess.
- **P8 CONFIRMED.** The cross-layer energy is 0.32402, 0.32264, 0.32149 per |X|, against 4 - 4 ln 2 - (12 - 16 ln 2) =
  0.31777.

**Reading.** The one-apex collision statistics are Poisson-pair statistics to within 2 % at n = 26 (HYP-9480 (A) at
first order; this hypothesis at second order). The upper tail and the energy are slightly over-dispersed, and the
excess is shrinking. Linear energy remains OPEN.

## Addendum pre-registration P9-P11 (n = 28; committed before the n = 28 run)

These bands come from the committed n = 22..26 trend. The excesses of E^Y/|Y| over 14 - 16 ln 2 were 0.032, 0.025 and
0.020, with ratio about 0.75 per two levels, so the extrapolation for n = 28 is about 0.016. Following MISTAKE-728,
the tail bound uses the model's own prediction.
- **P9.** 2.918 <= E^Y_28/|Y_28| < E^Y_26/|Y_26| = 2.92960 (excess in [0.008, 0.020), and decreasing).
- **P10.** values/|Y_28| in [0.4005, 0.4040], and the multiplicity fractions for q = 1..4 are within the P2 bands.
- **P11.** The maximum multiplicity at n = 28 is at most 22. The Poisson-pair model expects 0.54 values of
  multiplicity >= 18 and 0.034 of multiplicity >= 20 at |Y| = 2^26. Allowing the observed 25 % tail inflation, the
  chance of exceeding 22 is far below 1 %.
