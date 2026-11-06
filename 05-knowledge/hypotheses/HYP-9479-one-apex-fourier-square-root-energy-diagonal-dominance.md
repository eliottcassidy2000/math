---
id: HYP-9479
title: "Pre-registered tests for the Fourier side of the one-apex counts: (P1) square-root cancellation |S_n(theta)| ~ 2^(n/2) at five irrational theta for n = 26, 28; (P2) bounded energy E_n / |X_n| in [3.1, 3.5] for n = 26, 28; (P3) diagonal dominance of the bilinear second moment at unseen splits"
status: "P1, P2, P3 CONFIRMED out of sample (runs after the pre-registration commit 185c2b3efa); recorded in THM-4886 (D)"
session: opus-2026-11-06-S83
date: 2026-11-06
related: [THM-4886, THM-4885, HYP-9478, HYP-9474]
script: 04-computation/experiments/one_apex_fourier_20261106/
---

# HYP-9479 - square-root cancellation, linear energy and diagonal dominance for one-apex counts

## Setting

**Words and counts.**
- X_n is the set of binary words of length n with first letter 1, so |X_n| = 2^(n-1).
- H(x) = g(A(x)) is the Hamiltonian path count of the one-apex tournament of x.
- S_n(theta) = sum over x in X_n of e(theta H(x)).
- The energy is E_n = #{(x, x') in X_n^2 : H(x) = H(x')}.

**Bilinear split.** For x = u v w with |u| = a, |v| = L, |w| = l, Session 193's law gives
H(uvw) = B(u, w) + H(v) + A(u) C(v) + A(v) C(w) + boundary terms. With p = last bit of u and eps = [w begins in 0]
fixed, the bilinear second moment is

    Q = sum_{u, w} |F(u, w)|^2,   F(u, w) = sum_v e(theta (H(v) + p eps(v) + p(v) eps + A(u) C(v) + A(v) C(w))),

and its diagonal is 2^L K_u K_w, with K_u = 2^(a-2) and K_w = 2^(l-1). THM-4886 (A) proves Q <= diagonal (1 + R).

## In-sample evidence (before this commit)

- **|S_n(theta)| / 2^((n-1)/2)** at theta = phi = (sqrt5 - 1)/2, sqrt2 - 1, pi - 3, e - 2, sqrt3 - 1 for n = 14..24:
  every value lies in [0.166, 2.457] (`insample.py`).
- **E_n / |X_n|** = 2.930, 3.057, 3.188, 3.226, 3.242, 3.261, 3.257 for n = 12, 14, ..., 24. The distinct values per word
  are 0.373 at n = 20..24.
- **Q / diagonal** (`core_check.py`): 0.92-1.15 for q >= 37 at the splits (5,4,5), (6,5,6), (7,4,7), (6,6,6) and
  (8,3,8), and 1.25-4.37 at theta = 1/7 (a major arc).

## Predictions (out of sample)

- **(P1) Square-root cancellation.** For the five irrational theta above and n in {26, 28},
  |S_n(theta)| / 2^((n-1)/2) lies in [0.05, 20]. All ten values must lie in the window.
- **(P2) Linear energy.** E_n / |X_n| lies in [3.1, 3.5] for n = 26 and n = 28.
- **(P3) Diagonal dominance.**
  - Splits: (a, L, l) in {(9, 5, 9), (10, 5, 10), (9, 7, 9)}.
  - Frequencies: theta in {3/37, 5/101, 11/257, phi}.
  - Boundary bits: all four (p, eps) pairs.
  - Q / diagonal lies in [0.85, 1.2] in all 48 cases.

Confidence: P1 about 85% (the extreme ratios fluctuate); P2 about 80%; P3 about 75%.

## Results

Source: `04-computation/experiments/one_apex_fourier_20261106/outsample.out`, produced after commit 185c2b3efa.

| n | |S_n|/2^((n-1)/2) at phi, sqrt2-1, pi-3, e-2, sqrt3-1 | E_n / |X_n| | distinct / |X_n| |
|---|---|---|---|
| 26 | 1.723, 0.612, 0.413, 0.730, 1.700 | 3.2526 | 0.3734 |
| 28 | 1.979, 1.782, 1.547, 0.821, 1.443 | 3.2473 | 0.3737 |

- **(P1) CONFIRMED.** All ten ratios lie in [0.05, 20]; the actual range is 0.41-1.98.
- **(P2) CONFIRMED.** 3.2526 and 3.2473 both lie in [3.1, 3.5].
- **(P3) CONFIRMED.** Q / diagonal lies in [0.9971, 1.0686] in all 48 cases; none falls outside [0.85, 1.2].
- The vectorised block generator of g agrees with the recursion in 4,600 random checks.
