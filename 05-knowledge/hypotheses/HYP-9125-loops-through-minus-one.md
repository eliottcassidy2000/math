---
id: HYP-9125
title: "Mirror of HYP-9122: loops through -1 with the least admissible number of multiplications exist for every K >= 10"
status: >
  OPEN HYPOTHESIS. FINITE-EXACT for 10 <= K <= 4000, with explicit loops
  verified by an independent checker. The exceptions K = 5..9 are PROVED by
  exhaustion. For every K >= 10 there is a loop -1 -> ... -> -1 in the
  negative E-graph with K halvings and a0(K) = ceil((K+1) log_3 2)
  multiplications, i.e. with the least conceivable ratio 2*3^(eta(K+1)).
  It reduces (PROVED) to record objects at the lower best approximations
  of log_2 3 (19/12, 84/53, 569/359, 1054/665, ...). The first three
  record cycles are the three negative Collatz cycles.
source: collatz-procgen-20260922 (mac-mini), Q1-mirror lane (candidate M1)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_q1_mirror.md
  - 05-knowledge/hypotheses/HYP-9122-loops-through-one-bounded-ratio.md
---

# HYP-9125 -- loops through -1 (Q1 mirror of HYP-9122)

**Refutation form:** a `K >= 10` for which every loop through `-1` has
more than `a0(K)` multiplications.

**Role.** The hypothesis fixes the exact escape price from the `-1` thread,
`P(m) = 3^(eta(m+1))` in `(1,3)`. That price exceeds 1 at every precision
(PROVED), so a see-saw is needed on the Q1 side too. Next records:
`25781/16266` and `50508/31867`.
