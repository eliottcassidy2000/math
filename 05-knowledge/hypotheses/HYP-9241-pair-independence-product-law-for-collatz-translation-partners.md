---
id: HYP-9241
title: "Pair-independence product law for Collatz translation partners: for Haar y, the probability q_R(T) that y, y+1, ..., y+R have pairwise distinct Terras orbits up to time T is asymptotically the product over the R(R+1)/2 pairs of the pair non-merging probabilities q^(d)(T) (d = lag), so q_R(T) ~ c_R q_1(T)^(R(R+1)/2) and, with THM-4593's exponent 1/2, the non-coalescence exponent is R(R+1)/4 (the vicious-walker exponent) - but with independent pairs, not Karlin-McGregor repulsion"
status: >
  OPEN. NUMERICAL (session: 20000 samples of 4096-bit y, T <= 2048; audit B: 40000 samples of 2048-bit y, T <= 1024, and 20000 of 4096-bit y,
  T <= 2048, other seeds). q_2 = q_01 q_12 q_02 within 3% for T >= 64 (ratio 0.97-1.02); c_2 = q_2/q_1^3 = 1.11-1.22 is the lag factor
  q_02/q_01 ~ 1.20; q_3/q_1^6 = 1.66 -> 1.91 over T = 16..512 against the lag-product prediction (q_02/q_01)^2 (q_03/q_01) ~ 1.99 (q_03/q_01 ~ 1.38).
  Local-slope ratios 1 : 2.91 : 5.6-5.8. For three Brownian vicious walkers the Karlin-McGregor ratio q_2/(q_01 q_12 q_02) tends to pi/4;
  the data reject that constant. No unequal-time merges occur (impossible at these sizes: |m ln 3 - n ln 2| >= 4.4e-5 for |m| <= 2048).
  The asymptotic regime is not reached (single-pair local slope 0.30-0.35 at T ~ 2000; THM-4593 gives 1/2 at T ~ 1e6).
  Reframed after audit B (MISTAKE-586): first filed as a "Karlin-McGregor / total-positivity law".
source: mac-mini-2026-10-07-twoanchor, 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 8.5)
related:
  - 01-canon/theorems/THM-4593-partner-coalescence-pins-the-diffusive-exponent-and-all-lags-merge-exponentially.md (its q_R is "y merges with none"; here "no two merge")
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md
  - Karlin-McGregor (1959); Fisher (1984) vicious walkers (same exponent, different constant)
scripts:
  - 04-computation/experiments/twoanchor_20261007/km_exponent.py, km_analysis.py (+ km_exponent_n3000.out, km_exponent_n20000.out, km_analysis.out)
  - 04-computation/experiments/twoanchor_20261007/audit_B/km_audit.py, km_audit_summary.py (+ outputs)
---

# HYP-9241 — pair independence for Collatz partners

## Statement

* Let `q^(d)(T)` be the probability that y and y+d have not merged at equal time by Terras time T (Haar y).
* Let `q_R(T)` be the probability that no two of `y, …, y+R` have merged by time T.
* Conjecture:

      q_R(T) ~ ∏_{0 ≤ i < j ≤ R} q^(j−i)(T).

* Hence `q_R ≍ q_1^(R(R+1)/2)`. With `q_1(T) = T^(−1/2+o(1))` (THM-4581 (4): c T^(−1/2) below, C T^(−1/2)(log T)² above at sketch level), the exponent is R(R+1)/4.

## Why this is interesting

* The exponent R(R+1)/4 is the vicious-walker (Karlin–McGregor) exponent. Determinantal non-colliding walkers, however, have a different constant: π/4 for three Brownian walkers relative to the product of pair probabilities.
* The Collatz data give ratio 1. The pair failure events behave as independent, even though all partners are driven by the same parity bits of y.
* So the coalescence of translation partners looks like independent pairwise merging, not like non-crossing paths with repulsion. The walkers here are odd-step counts, which meet several times before merging.
* Contrast THM-4593. There the event "y merges with none of its partners" keeps exponent 1/2, because the partners coalesce among themselves first.

## Evidence

* See the status field.
* `q_2/q_1^3` is stable to within about 7% over T = 16..2048 in three independent runs.
* The R = 3 ratio rises toward the lag-product value of about 2.

## Open

* A proof for R = 2: asymptotic independence of the three coupled pair chains of THM-4581.
* The asymptotic regime T ≥ 10^5.
