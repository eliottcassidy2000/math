---
id: HYP-9125
title: "Mirror of HYP-9122: loops through -1 with the least admissible number of multiplications exist for every K >= 10"
status: >
  OPEN HYPOTHESIS. FINITE-EXACT for 10 <= K <= 176249 (sweep lane
  2026-09-23; previously K <= 4000), with explicit, independently verified
  record objects for every record n <= 50508. For every K >= 10 there is a
  loop through -1 with a0(K) = ceil((K+1) log_3 2) multiplications. The
  exceptions K = 5..9 are PROVED by exhaustion. The analogue for E_{6 mod 8}
  is false for all K >= 35, via the rate-2/3 trap of HYP-9121's refutation.
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
