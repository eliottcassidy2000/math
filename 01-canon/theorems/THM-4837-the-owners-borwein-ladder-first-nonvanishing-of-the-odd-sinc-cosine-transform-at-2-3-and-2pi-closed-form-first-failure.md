---
id: THM-4837
title: "The owner's Borwein ladder (CITED core): 15, 113 and 7, 56, 418 are the Borwein breaks and crossing indices of the odd harmonic sum at 2, 3, 4, and 80499 is the Baillie-Borwein-Borwein sinc-sum/integral break at 2*pi (the k = 1 Poisson term, i.e. the cosine transform at one full turn); first-failure deficits by Borwein-Borwein's Corollary 1 (2.31e-11, 2.33e-138, one-sided 4.2101759371e-226577, whose two-sided double is Crandall's 8.42e-226577); a fully proved window criterion certifies the crossing index floor(exp(2 xi - gamma)/4) for n <= 3000 and xi = k*pi, k <= 500; the OEIS '+1/8' formula line contradicts the entry's own a(17)"
status: CITED (Borwein-Borwein 2001; Baillie-Borwein-Borwein 2008; Almkvist-Gustavsson 2014; Basel-Baillie 2016; Crandall via Baillie 2011; Schmid 2014) + PROVED (the window criterion, elementary, on the GKP/DLMF envelope) + FINITE-EXACT (interval arithmetic, re-certified by the audit); the decode of the owner's 231/233 is INTERPRETATION and of 420 SPECULATIVE; independently audited SOUND WITH CORRECTIONS (MISTAKE-692)
session: opus-2026-10-28-S67
date: 2026-10-28
note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py (+ .out)
audit: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/AUDIT.md
related: [THM-4838, HYP-9436, THM-1985, THM-1990, THM-4392]
external:
  - "D. Borwein, J. M. Borwein, Some remarkable properties of sinc and related integrals, Ramanujan J. 5 (2001) 73-89 (Theorem 2; Corollary 1 = the first-failure value)"
  - "R. Baillie, D. Borwein, J. M. Borwein, Surprising sinc sums and integrals, Amer. Math. Monthly 115 (2008) 888-901 (Theorem 2 'First Bite'; Example 1(b): the sinc sum equals the integral for N <= 40248 and fails at N = 40249, denominator 80499; Example 2: sinc products over prime reciprocals)"
  - "R. Baillie, Fun with very large numbers, arXiv:1105.3943 (Crandall's error 'about 8.42e-226577' at N = 40249; two-sided normalisation)"
  - "G. Almkvist, J. Gustavsson, More remarkable sinc integrals and sums, arXiv:1405.1265 (the 1/113 value; thresholds 5 and 7)"
  - "U. Basel, R. Baillie, Sinc integrals and tiny numbers, Elem. Math. 71 (2016) 7-20 (crossings and 40-digit deficits for thresholds 2..10)"
  - "H. Schmid, Two curious integrals and a graphic proof, Elem. Math. 69 (2014) 11-17; S. N. Majumdar, E. Trizac, PRL 123 (2019) 020201 (random-walk reading)"
  - "Graham-Knuth-Patashnik, Concrete Mathematics (9.89) / DLMF 5.11(ii) (enveloping expansions)"
  - "OEIS A092315 (Sloane 2004): Song 2022 (Borwein link), Kirchner 2020 (floor formula with an explicit-remainder test certifying n <= 1000), Hasler 2017 (+1/8 formula line)"
---

# THM-4837 - the owner's Borwein ladder

**Setting.**
- P_m(x) = prod_{j=0}^{m} sinc(x/(2j+1)) and O_m = sum_{j=0}^{m} 1/(2j+1).
- U_j is uniform on [-1/(2j+1), 1/(2j+1)], independent, and S'_m = U_1 + ... + U_m.
- F_m(xi) = int_R P_m(x) cos(xi x) dx and I_m(t) = int_0^oo (sin t x / x) prod_{j=1}^{m} sinc(x/(2j+1)) dx.

**(1) Vanishing (CITED principle).** For m >= 1, F_m(xi) = pi P(xi - 1 <= S'_m <= xi + 1). So for xi >= 1,
F_m(xi) = 0 iff O_m <= xi.
- The box-spline/frequency-shift reading is Schmid 2014 (eqs. (12), (19)) and Majumdar-Trizac 2019.
- At m = 0, xi = 1 the edge value is pi/2, so m >= 1 is needed.

**(2) Borwein ladder (CITED).** I_m(t) = pi/2 - pi P(S'_m > t) = pi/2 - sum_{k>=0} F_m(t + 1 + 2k). So I_m(t) = pi/2 iff
O_m <= t + 1.
- t = 1 is the classical integral, breaking at sinc(x/15).
- t = 2 is the 2cos x variant, breaking at sinc(x/113).
- t = 4 and t = 6 are Almkvist-Gustavsson's cosine families, with (2cos x + 2cos 3x [+ 2cos 5x]) sin x = sin 4x
  [sin 6x]. They break at m = 3091 and 168803.
- Crossings for thresholds 2..10 are tabulated by Basel-Baillie 2016, as n = m + 1.

**(3) First failure (CITED).** Let m = m(xi) and eps = O_m - xi. Then 0 < eps <= 1/(2m+1) and
F_m(xi) = pi eps^m (2m+1)!!/(2m)!!.
- This is Borwein-Borwein 2001 Corollary 1 (= BBB 2008 Theorem 2 "First Bite"). Its hypotheses say that m is the first
  failure, and for the decreasing sequence 1/(2j+1) they hold at every first failure.

**(4) The 2 pi row (CITED).** By Poisson summation, sum_{n in Z} P_m(n) - int_R P_m = 2 sum_{k>=1} F_m(2 pi k).
- BBB 2008 Example 1(b): the sinc sum equals the integral for N <= 40248 and first fails at N = 40249, denominator
  80499. Crandall estimated the error, as reported by Baillie 2011.
- The cosine transform at one full turn, F_m(2 pi), is the k = 1 term, the one-sided error.
- **New here:** F_40249(2 pi) = 4.2101759371 x 10^-226577 to 11 digits, certified. The literature's "about
  8.42 x 10^-226577" is the two-sided error sum_Z - int_R (factor 2).

**(5) Values.**
- At the thresholds 2, 3, 4, 2 pi, (m, 2m+1) = (7, 15), (56, 113), (418, 837), (40249, 80499).
- The deficits are 2.3100572725e-11 (BB 2001), 2.3324296128e-138 (Almkvist-Gustavsson), 4.5195214277e-1320 (a trivial
  instance of Corollary 1) and 4.2101759371e-226577 (one-sided).
- The 1/15 value reproduces the published rational 467807924713440738696537864469/935615849440640907310521750000 pi.

**(6) Crossing index (PROVED; a fully proved version of Kirchner's test).** With M = m + 1,
O_m = (1/2) ln(4M) + gamma/2 + r(M), where r(M) = 1/(48M^2) - 7/(1920M^4) - eps_{2M} + eps_M/2 and
0 < eps_N < 1/(252 N^6). So 0 < r(M) < 1/(48 M^2) for all M >= 1. With X = exp(2 xi - gamma)/4:
- **m(xi) = floor(X) whenever frac(X) >= X/(24 floor(X)^2)**;
- otherwise m(xi) is floor(X) - 1 or floor(X).

This certifies the floor formula for 2 <= n <= 3000 and for xi = k pi, k <= 500 (re-certified in interval arithmetic by
the audit). Kirchner's 2020 explicit-remainder test, with its remainder supported numerically, had certified n <= 1000.

**(7) OEIS observation.** The A092315 formula line floor(exp(2n - gamma)/4 + 1/8) ("conjectured", Hasler 2017)
contradicts the entry's own DATA term a(17) = 81897532160124.
- It differs from the certified value exactly when frac(X_n) >= 7/8: 375 of the n <= 3000.
- a(17) is certified by the DLMF envelope and by digamma: O_{a-1} - 17 = -5.7e-15 and O_a - 17 = +4.0e-16.
- An OEIS edit is recommended.

**(8) The owner's numbers.**
- 7, 56, 418 are m(2), m(3), m(4).
- 15 and 113 are the Borwein breaks.
- 80499 is the BBB 2 pi break.
- INTERPRETATION: 231 and 233 match the published deficit mantissas (2.31, 2.33) to three digits.
- SPECULATIVE: 420 matches the one-sided 80499 mantissa 4.2102 only to two digits. "4.20" needs a naive double sum and
  truncation, not rounding.

**Novelty.** NEW to the repository: (4) bold, (6), (7) as an explicit certified observation, and the decode (8). Every
other statement is published (see `external`) and enters the repository here for the first time.
