---
id: HYP-3533
title: The LRC14 covering-floor CV criterion is an exact RESONANT-ENERGY statement CV(N_R)^2 = (sum_{N!=0}|chat(14N)|^2)/m_R^2; the sheets are SUB-BINOMIAL for coprime R but SUPER-BINOMIAL (rho up to ~2.5) for even/7-heavy R -- the explicit quantitative form of the S259 2-adic obstruction, compensated by larger m_R
status: CONFIRMED (exact identity + measured mechanism); sharpens THM-579. Not a closure.
source: mac-mini-2026-06-29-S3
theorem: THM-579
related:
  - THM-579    # the CV floor criterion this sharpens
  - HYP-3415   # kps critical-path covering floor
  - HYP-3129   # SPEC resonance lattice (this is the R-side energy)
  - THM-576    # cap lower bounds used in the clean inequality
script: 04-computation/lrc14_floor_subbinomial_sheetvar_macmini_20260629.py
result: 05-knowledge/results/lrc14_floor_subbinomial_sheetvar_macmini_20260629.out
---

# HYP-3533 — the floor CV is resonant energy; even speeds are super-binomial

## 1. Exact spectral identity (the resonant-energy form)
The Fejer/projection identity (`sum_{d=-13}^{13}(14-|d|)e(nd/14) = 196.[14|n]`) gives, exactly,
> **`CV(N_R)^2 = Var(N_R)/E[N_R]^2 = (sum_{N!=0} |chat(14N)|^2) / m_R^2`,**
i.e. the floor's coefficient of variation IS R-safe's spectral energy on the resonance lattice
`14Z\{0}`, relative to `m_R^2`.  So the THM-579 sufficient condition `CV(N_R)^2 < m_Q/(1-m_Q)`
is exactly
> **`sum_{N!=0} |chat(14N)|^2  <  m_R^2 . m_Q/(1-m_Q)`**:
the 2-adically-resonant energy of `R-safe` must be dominated by `m_R^2 . threshold`.  This puts the
whole `14 = 2.7` structure in one place -- the `14Z` lattice -- and is the cleanest crisp form of
the covering floor.

## 2. The clean 2-variable inequality (would close it under sub-binomiality)
If the 14 sheets were sub-binomial (`Var(N_R) <= 14 m_R(1-m_R)`), then `CV^2 <= (1-m_R)/(14 m_R)`,
and the criterion is implied by the ELEMENTARY
> `1 < m_R + m_Q + 13 m_R m_Q`,
which HOLDS at every worst-case cap (THM-576): r=2..6 give `>= 1.062, 1.105, 1.181, 1.219, 1.196`.
So under sub-binomiality the covering floor closes for all r with margin.

## 3. But sub-binomiality FAILS at even/7-heavy R (the 2-adic mechanism)
Measured `rho = Var(N_R)/[14 m_R(1-m_R)]`:

| R | m_R | rho | note |
|---|---|---|---|
| consec_8 {1..8} | .266 | 0.59 | sub-binomial |
| large-coprime {1,3,5,9,11,13,15} | .304 | 0.60 | sub-binomial |
| consec_7 {1..7} | .335 | 1.01 | ~binomial |
| 7-heavy (contains 7) | .382 | 1.18 | super |
| even-heavy {1,2,3,4,6,8,10,12} | .338 | 2.19 | **super** |
| even-heavy {1,2,4,6,8,10,12} | .386 | 2.48 | **super** |

Random max `rho ~ 2.0-2.2` (|R|=7..11).  So even/7 speeds CLUSTER the sheets (positive sheet
correlation = the 2-part of 14, S259's "even speeds are binding"), pushing `rho` to ~2.5.  This is
the precise quantitative form of the 2-adic obstruction.

## 4. Why the floor still holds (the compensation)
High `rho` (even-heavy R) co-occurs with high `m_R` (those R are far from the cap), so the binomial
baseline `(1-m_R)/(14 m_R)` is small and `rho x baseline = CV^2` stays below threshold.  Example:
even-heavy_7 has `rho=2.48` but `m_R=0.386`, giving `CV^2=0.283 < cap_6/(1-cap_6)=0.391`.  The
worst case for the criterion is NOT (worst m_R AND worst rho): they are ANTI-CORRELATED.  This is
why THM-579's sweep holds for all r despite super-binomiality.

## 5. Remaining obstruction (crisp)
Close the floor on the C-S branch by bounding `sum_{N!=0}|chat(14N)|^2 < m_R^2 m_Q/(1-m_Q)`
uniformly.  Two routes:
- (a) prove the joint bound `rho(R).(1-m_R)(1-m_Q) < 14 m_R m_Q` tracking the rho/m_R
  anti-correlation (the sweep says it holds, tightest at r=6, margin +0.096);
- (b) bound the resonant energy `sum_{N!=0}|chat(14N)|^2` directly via the singular series
  `chat(14N) = sum_{sum k_r r = 14N} prod ahat(k_r)`, `ahat(k)=-sin(pi k/7)/(pi k)` -- the SAME
  object as the cap (THM-576) and doublet (THM-578).
Where C-S is too lossy, splice HYP-3129's exact-low.  Caps (exact, with dips): cap_j for j=2..11 =
.7253, .6044, .4943, .3815, .2810, .1967, .1405, .0925, .0565, .0323 (cap_10 = m_P = 14249/252252).
