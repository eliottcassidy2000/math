---
id: HYP-9124
title: "X_min: no integer m >= 2 lies in the backward hostile set of graph E (implies Q2 with no thresholds)"
status: >
  OPEN HYPOTHESIS with FINITE-EXACT support to 10^18. Every integer m >= 2
  with 3 not dividing m has a legal reverse path x -> (2^k x - 1)/3 whose
  multiplier 2^K/3^s is < 1. PROVED implication (Theorem 6.1 of the
  endgame lane): X_min implies Q2, with no thresholds and without
  HYP-9122. The canonical strategy Psi already gives a multiplicative
  descent for every m <= 10^18, via an exhaustive DFS over Psi-alive
  classes mod 3^38.
source: collatz-procgen-20260922 (mac-mini), endgame lane
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_q2_endgame.md
  - 05-knowledge/hypotheses/HYP-9120-e-scc-choice-collapse-and-hostile-rationals.md
---

# HYP-9124 -- no hostile integer (backward E)

**Refutation form:** an integer `m >= 2`, prime to 3, all of whose legal
reverse paths keep every prefix multiplier `>= 1`.

**Reductions (PROVED, endgame lane section 6).**
* Assuming HYP-9122, `X_min` is equivalent to its restriction to
  `m = 14 mod 27`.
* That restriction follows from `X_T`: for every `k >= 3` and every odd
  `w` prime to 3, the integer `2^(floor(k log_2 3)) w` has a legal path of
  multiplier `< 3^k/2^(floor(k log_2 3)+1)`.

**Obstruction.** Integers sit beyond the "3/2 wall", where the exact
hostility criterion does not terminate.
* Four census points just above `3/2` looked hostile to depth 41 but
  descend at depths 53, 65, 106 and 212.
* Integer candidates may likewise need long, near-optimal climbs.
