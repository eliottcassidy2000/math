---
id: HYP-9128
title: "AMM 12592: super-blocks with cross-annulus handoff give C* < C_* = 1 + log_5(phi^2)"
status: >
  OPEN HYPOTHESIS. Evidence is FINITE-EXACT at block level and NUMERICAL
  asymptotically. Super-blocks [N,4N) let imbalance cross the inner dyadic
  boundary through one integer handoff polynomial S_N(w). Exactly verified
  integer super-blocks reach max T(L)/L = 53/34 on [64,256), 83/53 on
  [128,512) and 157/100 on [256,1024). At those ratios no separately
  balanced block exists; this is certified exactly by Long's evaluation
  inequality. Handoff states vanishing to order about N/16 at the golden
  point w = -1 numerically converge to about 1.570 < C_* = 1.59799. The
  proved window is C* in [11/8, C_*] (THM-4467 and Long/THM-3009 constructions).
source: collatz-procgen-20260923 AMM lane (candidate H1)
depends_on:
  - 05-knowledge/results/amm12592_procgen_20260923_uniform_frontier.md
  - 01-canon/theorems/THM-4467-uniform-polya-capacity-gap-amm12592.md
---

# HYP-9128 -- super-blocks beat the golden constant

**Claim.** For `kappa = 1/16`, the states
`S_N = (1+w)^(N/16) (1-w^N)/(1-w^(N/16))` on super-blocks `[N,4N)`, `N = 4^k`,
yield integer fair extractors with `limsup T/L <= 1.573`. Hence `C* < C_*`.

**Proof plan** (from the lane):
1. explicit contours for the fold integral, checked by interval arithmetic
   uniformly in `v = i/N`;
2. Long's five-unit lattice rounding for palindromic targets;
3. finite certificates below `N_0`.

**Why it matters.** Theorem B of the lane note shows that every sub-golden
extractor must be analytic at the golden point `w = -1`. The super-blocks
achieve this through their zero there.
