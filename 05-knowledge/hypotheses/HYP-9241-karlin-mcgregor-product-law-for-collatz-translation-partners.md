---
id: HYP-9241
title: "Karlin-McGregor product law for Collatz translation partners: for Haar y, the probability q_R(T) that y, y+1, ..., y+R have pairwise distinct Terras orbits up to time T satisfies q_R(T) ~ c_R q_1(T)^(R(R+1)/2), so with THM-4593's exponent 1/2 the non-coalescence exponents are R(R+1)/4 (vicious walkers: one factor per pair)"
status: >
  OPEN. NUMERICAL (20000 samples of 4096-bit random y, T <= 2048): q_2/q_1^3 = 1.11, 1.13, 1.15, 1.17, 1.19, 1.19, 1.30 at T = 16..1024;
  q_3/q_1^6 = 1.70, 1.79, 1.94, 1.93, 1.89, 1.90 at T = 16..512, while q_3 falls from 0.094 to 0.0027. Local slopes are in ratio about
  1 : 2.9 : 5.6 (vicious walkers: 1 : 3 : 6). The single-pair slope is still 0.34 at T ~ 2000 (asymptotic 1/2, THM-4593), so the
  asymptotic regime is not reached.
source: mac-mini-2026-10-07-twoanchor, 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 8.5)
related:
  - 01-canon/theorems/THM-4593-partner-coalescence-pins-the-diffusive-exponent-and-all-lags-merge-exponentially.md (q_R there is "y merges with none"; here "no two merge")
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md
  - Karlin-McGregor (1959) coincidence probabilities; Fisher (1984) vicious walkers
scripts:
  - 04-computation/experiments/twoanchor_20261007/km_exponent.py, km_analysis.py (+ km_exponent_n3000.out, km_exponent_n20000.out, km_analysis.out)
---

# HYP-9241 — a Karlin–McGregor law for Collatz partners

**Statement.**
* Let `q_R(T)` be the probability that no two of the orbits of `y, y+1, …, y+R` (Terras map, Haar y) have merged at equal time by time T.
* Conjecture: `q_R(T) ≍ q_1(T)^(R(R+1)/2)`.
* Given `q_1(T) ≍ T^(−1/2)` (THM-4581, THM-4593), the exponents are `R(R+1)/4`. These are the Karlin–McGregor / vicious-walker exponents of R+1 non-colliding diffusions.

**Difference from THM-4593.**
* THM-4593 bounds the event "y merges with none of its partners". Partners may coalesce among themselves there, and the exponent stays 1/2.
* Here every pair must stay apart. The determinantal (totally positive) structure of non-colliding walkers predicts one factor of q_1 per pair.

**Evidence.** See the status field. The ratios `q_R/q_1^(R(R+1)/2)` are stable to within about 20% over two decades of T, while the probabilities themselves fall by one to two orders of magnitude.

**Open.**
* A proof for R = 2 from THM-4581's Lyapunov weight (three coupled pair chains).
* The asymptotic exponents at T ≥ 10^5.
