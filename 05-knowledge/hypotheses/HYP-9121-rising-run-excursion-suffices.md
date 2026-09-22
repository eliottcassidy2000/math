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
