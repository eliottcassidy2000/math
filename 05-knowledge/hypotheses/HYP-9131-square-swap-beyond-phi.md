---
id: HYP-9131
title: "Square swaps beyond the Tschakaloff threshold: the 5x+1 square-swap word with 23 ones in 33 letters has an irrational Bernstein number"
status: >
  SETTLED -- PROVED (2026-09-23; theta-beyond-phi lane, orchestrator-audited).
  The 5x+1 square-swap number sum (2^33/5^23)^(k^2) is irrational in Q_2.
  Proof: Theorem H1, the 2-adic transcription of Bezivin's Hankel-determinant
  method. The tails are exponential sums, and Cauchy-Binet gives a unique
  minimal-valuation term, so the Hankel determinant has exact valuation
  L n(n+1)(2n-1)/2. This beats the archimedean clearing whenever
  mu_bar < 7/4, and here mu_bar = 1.61831. With the transcribed KRVZ
  divisibilities (cited): mu_bar < 28/11 (every 5x+1 square-swap word) and
  mu_bar < 2.87837 (every 7x+1 word). The phi threshold of Zudilin's
  one-parameter construction was not the true frontier.
source: collatz-procgen-20260923 cube-theta lane (candidates HC-CT7, HC-CT4)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260923_cube_theta.md
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md
---

# HYP-9131 -- the square tier just past phi

**Refutation form:** the Bernstein number of that word is rational, i.e. the
2-adic theta value `sum (2^33/5^23)^(k^2)` is rational.
