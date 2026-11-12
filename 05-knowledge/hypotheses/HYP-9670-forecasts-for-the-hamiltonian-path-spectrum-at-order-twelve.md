---
id: HYP-9670
title: "Forecasts for the Hamiltonian path spectrum at order 12 (registered 2026-11-11, not scorable until order 12 is enumerated): the least odd number other than 7 and 21 that is not H of a 12-vertex tournament is 488,000 (interval [476,000, 500,000]); the maximum of H is 545,000 (interval [531,205, 552,000]); the number of distinct values is 262,000 (interval [250,000, 272,000])"
status: >
  OPEN (NUMERICAL forecasts; registered in PREREGISTRATION.md, block R2, of the source packet and pushed at
  9f260a24f7). Extrapolations from the exact censuses of orders 8-11 (THM-4889 / HYP-9481, klein; A038375) and
  the proved lower bound max H(12) >= 531,205 (THM-329): L(n) 2^(n-1)/n! = 1.940, 2.034, 2.072, 2.064 for
  n = 8..11, forecast 2.085 +- 0.05; for even n, e - max H(n) 2^(n-1)/n! = 0.718, 0.620, 0.497 at n = 6, 8, 10,
  forecast next ratio in [0.72, 0.82]; |S(n)|/(max H(n)/2) = 0.968, 0.906, 0.973, 0.924 for n = 8..11.
  No model of the upper tail is claimed. The unsourced value 545,835 flagged in MISTAKE-743 lies in the interval
  for the maximum; the forecast does not use it and is not evidence for it.
  Scoring: each forecast is CONFIRMED if the exact value lies in its interval, FAILED otherwise.
  AUDIT: arithmetic of the ratios re-derived by the independent audit of the source packet (audit/AUDIT.md); the
  forecasts themselves cannot be audited before order 12 is enumerated.
session: mac-mini-2026-11-11-gauge
source: 05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md
related:
  - THM-4889
  - HYP-9481
  - THM-329
  - HYP-9357
---

# HYP-9670 - forecasts for order 12

Source: [the note](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md), section 4, and
`PREREGISTRATION.md` (R2) of the packet.

| quantity at order 12 | point forecast | interval |
|---|---|---|
| least odd number other than 7, 21 that is not `H` of a 12-vertex tournament | 488,000 | [476,000, 500,000] |
| maximum of `H` | 545,000 | [531,205, 552,000] |
| number of distinct values of `H` | 262,000 | [250,000, 272,000] |
