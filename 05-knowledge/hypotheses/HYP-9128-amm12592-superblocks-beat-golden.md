---
id: HYP-9128
title: "AMM 12592: super-blocks with cross-annulus handoff give C* < C_* = 1 + log_5(phi^2)"
status: >
  SETTLED -- PROVED and promoted as THM-4468 (2026-09-23). Golden-zero
  super-blocks [N,4N) with S_N = (1+w)^(N/16) sum_(j<16) w^(jN/16) give an
  exactly fair extractor with T(L) <= ceil(159L/100). So
  C* <= 1.59 < C_* = 1.59799. The orchestrator audited it: fairness of the
  N = 16 and 64 blocks was re-verified from the definition, and the interval
  contour certificate for N >= 4096 was re-run.
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
