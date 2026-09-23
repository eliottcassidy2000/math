---
id: HYP-9131
title: "Square swaps beyond the Tschakaloff threshold: the 5x+1 square-swap word with 23 ones in 33 letters has an irrational Bernstein number"
status: >
  OPEN. This is the cheapest open instance of the PC square tier. Theorem Q
  settles every quadratic swap family when mu_bar < phi (every 3x+r map).
  For the 5x+1 square-swap word with blocks of 23 ones in 33 letters,
  mu_bar = 1.61831 = phi + 0.00028. Diagonal Pade forms beat the threshold
  at every tested size (margins +889 at mu_bar = 1.8575, n = 800), but their
  heights grow like n log n, so this is evidence, not proof (candidate
  HC-CT4 asks for a linear height bound). An improvement of Zudilin's phi
  for 2-adic theta values would settle it.
source: collatz-procgen-20260923 cube-theta lane (candidates HC-CT7, HC-CT4)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260923_cube_theta.md
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md
---

# HYP-9131 -- the square tier just past phi

**Refutation form:** the Bernstein number of that word is rational, i.e. the
2-adic theta value `sum (2^33/5^23)^(k^2)` is rational.
