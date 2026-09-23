---
id: HYP-9127
title: "The cube-swap word: sum_k (2^10/3^9)^(k^3) is irrational in Q_2 (the smallest open instance of the Periodicity Conjecture found)"
status: >
  OPEN HYPOTHESIS. The cube-swap word Y3 puts the block 1^8 0 1 at the
  nonzero cubes and 1^9 0 elsewhere. It is explicit, supercritical
  (beta = 9/10), of zero entropy and of discrepancy width 1, and has
  Diophantine exponent 1. Its 3x+1 Bernstein number is
  -1 - 512/(3^9-2^10) - (256/3^9) sum_{k>=1} rho^(k^3) with rho = 2^10/3^9
  (exact modulo 2^30000). So "no rational has an eventually-Y3 parity
  vector" is equivalent to the 2-adic irrationality of the cubic theta
  value. Both proved mechanisms fail. Periodic approximants need
  Dio > eta (Theorem D), and here Dio = 1. The Tschakaloff-Pade argument
  needs a first-order q-difference equation (Theorem Y, squares), and
  cubes have none.
source: collatz-procgen-20260922 (mac-mini), HARD-class lane (candidate HC1)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md
  - 05-knowledge/hypotheses/HYP-9123-supercritical-strips-periodicity.md
---

# HYP-9127 -- the cube-swap word (cubic 2-adic theta value)

**Refutation form:** the 2-adic number `sum_(k>=1) (2^10/3^9)^(k^3)`
is rational.

**Why it matters.** Foundry v4 names one missing mechanism for the
Collatz divergence half, the Periodicity Conjecture, both E-SCC halves,
Mahler and Erdős: non-integrality or irrationality on words with
`Dio <= eta` and no functional equation. This is the simplest explicit
number of that kind the session found. The square-swap analogue,
`sum rho^(k^2)`, is PROVED irrational (Theorem Y).

**Context.** Even in `R`, cubic theta values at `1/(integer)` are known
only to have algebraic degree at least 4 (Ghidelli 2019; abstract only).
The partial sums have heights `2^(14.26 K^3)` against 2-adic errors
`2^(-10(K+1)^3)`, so a Liouville argument would need `log2(3^9)/10 < 1`.
