---
id: HYP-9121
title: "Partial choice at 6 mod 8 already makes the Collatz exceptional set subexponential"
status: >
  OPEN HYPOTHESIS with FINITE-EXACT support. For E_S (Collatz plus the extra
  arrow n->3n+1 at even n with n mod 8 = 6), the number of residue classes
  mod 2^m without a multiplicatively descending certificate grows
  subexponentially. At 2^22: 3,238 versus 93,222 for Collatz, 1,452 for
  S={2 mod 4}, and 782 for all evens. S={0 mod 8}, {8 mod 16} and {0 mod 16}
  give no measurable gain, and S={2 mod 8} gives little.
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

