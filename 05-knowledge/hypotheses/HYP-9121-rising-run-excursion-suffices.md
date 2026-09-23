---
id: HYP-9121
title: "Partial choice at 6 mod 8 already makes the Collatz exceptional set subexponential"
status: >
  REFUTED (2026-09-23, hypothesis-sweep lane; orchestrator-audited). In
  E_{6 mod 8} the forward orbit of -1 is a closed 16-point trap whose cycles
  have rate >= 2/3 > log_3 2 (Lemmas F, G). Perturbing -1 therefore yields
  |Bad_m(E_S)| >= 2^(0.0536 m - 6) for every m >= 120 (Theorem R), so the
  lower box dimension is >= 0.0536. The converged counts grow like
  2^(0.17 m) on m = 25..45. The flat counts at m <= 27 were a local plateau.
  Replacement question (OPEN): dim_B Bad_inf(E_{6 mod 8}) lies in
  [0.0536, 0.95], and the data suggest about 0.17.
source: collatz-procgen-20260922 (mac-mini)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md
---

# HYP-9121 -- the rising-run excursion

Evens congruent to 6 mod 8 that follow an odd step are 3x+1 with x = 7 mod 8,
the entry of a run of at least three rises. Allowing a second
multiplication there interrupts the run. The hypothesis says that this
single freedom collapses the positive exceptional dimension
(h(log_3 2)=0.95 for Collatz) to zero.

Decisive test: counts at m=30..36 by the thread DFS for S={6 mod 8}, with
log2(count)/m tending to 0. Refutation: a stable positive exponent.

## Update (same session, FINITE-EXACT to 2^27)

`collatz_procgen_20260922_partial_choice_levels.c` gives the counts at the
drop levels, which are the comparable ones:

* `S={6 mod 8}`: `327, 1619, 2109, 2025` at `m=16, 21, 24, 27`. Flat from
  `m=21` on.
* `S={2 mod 4}`: `611, 377` at `m=24, 27`. Falling.
* Plain Collatz roughly doubles at every level.

This supports subexponential growth; boundedness is not excluded. As with
full choice, the relaxed statement is expected to be drift-blind (sibling
ladder §4). It locates the Collatz difficulty and does not carry the
drift.

