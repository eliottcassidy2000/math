---
id: THM-4837
title: "The owner's Borwein ladder: 15, 113, 80499 and 7, 56, 418 are the first non-vanishing points of the cosine transform of the odd sinc product at 2, 3, 2*pi (and 4); the first failure is always single-corner, with value pi*eps^m*(2m+1)!!/(2m)!! (2.31e-11, 2.33e-138, 4.21e-226577: the owner's 231, 233 and, to two digits, 420); the crossing index is floor(exp(2 xi - gamma)/4) outside an explicit window, and the OEIS '+1/8' formula fails at n = 17"
status: PROVED (elementary, on the cited Borwein-Borwein criterion and the GKP harmonic envelope) + FINITE-EXACT (interval arithmetic, exact rationals) + CITED; independently audited (see the note's section 8)
session: opus-2026-10-28-S67
date: 2026-10-28
note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py (+ .out)
related: [THM-4838, HYP-9436, THM-1985, THM-1990]
external:
  - "D. Borwein, J. M. Borwein, Some remarkable properties of sinc and related integrals, Ramanujan J. 5 (2001) 73-89 (criterion and first-failure value)"
  - "H. Schmid, Two curious integrals and a graphic proof, Elem. Math. 69 (2014) 11-17 (the 2cos x variant to 1/111)"
  - "R. Baillie, D. Borwein, J. M. Borwein, Surprising sinc sums and integrals, Amer. Math. Monthly 115 (2008) 888-901"
  - "Graham-Knuth-Patashnik, Concrete Mathematics, eq. (9.89) (enveloping expansion of H_n)"
  - "OEIS A092315: Song 2022 (Borwein link of a(2), a(3)), Kirchner 2020 (floor formula, check to n <= 1000), Hasler 2017 (+1/8 conjecture)"
---

# THM-4837 - the owner's Borwein ladder

**Setting.**
- P_m(x) = prod_{j=0}^{m} sinc(x/(2j+1)) and O_m = sum_{j=0}^{m} 1/(2j+1).
- U_j is uniform on [-1/(2j+1), 1/(2j+1)], independent, and S'_m = U_1 + ... + U_m.
- F_m(xi) = int_R P_m(x) cos(xi x) dx and I_m(t) = int_0^oo (sin t x / x) prod_{j=1}^{m} sinc(x/(2j+1)) dx.

**(1) Vanishing.** F_m(xi) = pi P(xi - 1 <= S'_m <= xi + 1). For xi >= 1, F_m(xi) = 0 iff O_m <= xi.

**(2) Borwein ladder.** I_m(t) = pi/2 - pi P(S'_m > t) = pi/2 - sum_{k>=0} F_m(t + 1 + 2k). So I_m(t) = pi/2 iff
O_m <= t + 1 (Borwein-Borwein).
- t = 1 is the classical integral, which breaks at sinc(x/15).
- t = 2 is the 2cos x variant (Schmid), which breaks at sinc(x/113), because 2 cos x sin x = sin 2x.

**(3) First failure.** Let m = m(xi) be least with O_m > xi and put eps = O_m - xi. Then 0 < eps <= 1/(2m+1), and
F_m(xi) = pi eps^m (2m+1)!!/(2m)!!. At t = xi - 1 this also equals pi/2 - I_m(t).
- Proof: minimality gives eps <= 1/(2m+1) <= 2 min_j a_j. So exactly one corner simplex of the box crosses the level,
  and it has volume eps^m/m!.

**(4) The owner's values.**
- At the thresholds 2, 3, 4, 2 pi:
  - (m, 2m+1) = (7, 15), (56, 113), (418, 837), (40249, 80499);
  - first values 2.3100572725e-11, 2.3324296128e-138, 4.5195214277e-1320, 4.2101759371e-226577.
- The 2 pi row: 1 + 1/3 + ... + 1/80497 < 2 pi < 1 + 1/3 + ... + 1/80499, certified by interval arithmetic. The integral
  int_R cos(2 pi x) prod_{j=0}^{m} sinc(x/(2j+1)) dx is 0 for every denominator up to 80497 and positive at 80499.
- Readings: owner's 7, 56, 418 = m(2), m(3), m(4); owner's 15, 113 and 80499 = 2m+1 at 2, 3 and 2 pi; owner's 231 and 233
  are the 3-digit mantissas at 15 and 113, and 420 matches the 2 pi mantissa to two digits (naive double-precision
  summation gives 4.208).
- The value at 15 reproduces the published rational 467807924713440738696537864469/935615849440640907310521750000 pi.

**(5) Crossing index.** With M = m + 1, O_m = (1/2) ln(4M) + gamma/2 + r(M), where
r(M) = 1/(48M^2) - 7/(1920M^4) - eps_{2M} + eps_M/2 and 0 < eps_N < 1/(252 N^6) (GKP). So 0 < r(M) < 1/(48 M^2) for all
M >= 1. Consequently, with X = exp(2 xi - gamma)/4:
- **m(xi) = floor(X) whenever frac(X) >= X/(24 floor(X)^2);**
- otherwise m(xi) is floor(X) - 1 or floor(X).

**(6) OEIS A092315 repaired.** The conjectured form floor(exp(2n - gamma)/4 + 1/8) (Hasler 2017) is false.
- At n = 17, a(17) = 81897532160124 = floor(X), while the +1/8 form gives 81897532160125.
- An independent digamma evaluation gives O_{a-1} - 17 = -5.7e-15 and O_a - 17 = +4.0e-16.
- It fails at 375 of the n <= 3000.
- The floor form (Kirchner 2020) is PROVED for 2 <= n <= 3000 by (5); see HYP-9436 for all n.

**Novelty.**
- (1)-(2) and the first-failure formula are CITED.
- NEW to the repository: the 2 pi/80499 identification, the observation that the first failure is single-corner, the
  values at 113 and 80499, the decode of the owner's six numbers and 231/233/420, the window criterion, and the
  refutation of the +1/8 conjecture.
- Not new: the Borwein link of A092315 (Song's OEIS comment, 2022).
