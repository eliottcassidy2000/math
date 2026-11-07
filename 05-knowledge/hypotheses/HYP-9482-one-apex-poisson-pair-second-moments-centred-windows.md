---
id: HYP-9482
title: "One-apex collisions on reduced words: Poisson-pair second moments (energy 14 - 16 ln 2 up to a small positive excess, multiplicity law, smoothed slope 12 - 16 ln 2, cross-layer energy 4 - 4 ln 2 - (12 - 16 ln 2)), centred collision windows with |a - s| <= 1 + carry depth, and a halving level law"
status: "PRE-REGISTERED (this commit, before any computation of Y_n for n >= 21)"
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
