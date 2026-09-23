---
id: HYP-9126
title: "No hostile point above 3/2: every E-hostile rational has absolute value below 3/2 (both halves)"
status: >
  OPEN HYPOTHESIS. Backward: every dyadic hostile point p/2^e has value
  < 3/2. All 4 census points above 3/2 (e <= 45) and all 11 lower-clock
  generation-1 points (e <= 200) descend, with explicit certificates
  (independently replayed). Forward: every hostile -p/3^j has |x| < 3/2.
  28 of 72 candidates above 3/2 (j <= 27) are PROVED to descend; the rest
  are undetermined. With the PROVED straddle identity
  (|x|-3/2)(y-3/2) = -(rho-1)^2/(2 rho), this implies that the only
  numerator shared by two hostile generation-1 points is 1.
source: collatz-procgen-20260922 (mac-mini), Q1-mirror lane M3 and endgame lane C2
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_q1_mirror.md
  - 05-knowledge/results/collatz_procgen_20260922_q2_endgame.md
  - 05-knowledge/results/collatz_procgen_20260922_exceptional_dimension.md
---

# HYP-9126 -- the 3/2 wall

**Refutation form:** a rational hostile point (forward `-p/3^j` or backward
`p/2^e`) with absolute value `>= 3/2`.

**Why 3/2.**
* Below `3/2`, backward hostility is decidable by a finite search
  (endgame lane, Criterion H; Q1-mirror lane, section 6.1).
* Above it, certification can require near-optimal climbs to `1`, i.e.
  near-solutions of `2^A p - 2^B ~ 3^C`. The four census points above `3/2`
  descend only at depths 53 to 212.
* The wall is also where integers `m >= 2` sit, which is why the
  hypothesis bears on HYP-9124.
