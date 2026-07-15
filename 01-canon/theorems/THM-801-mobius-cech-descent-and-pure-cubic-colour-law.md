---
id: THM-801
title: The seven-term tiling recursion lifts to three-face Cech descent and a pure cubic blue-colour law
status: RESERVED — set-level proof and colour derivation in hand; exact line/node audit and final scope being completed
source: codex-2026-07-15-S12
depends_on: [THM-442, THM-553, THM-796]
related: [THM-549, THM-550, HYP-2685, HYP-3234, HYP-6870, HYP-6880]
---

# THM-801 — Möbius/Čech descent for the merged metagraph

This file reserves the theorem number for the categorical lift of the old
`A+B+C-D-E-F+G` recursion.  In pin coordinates `(r,c)=(a-b-1,b)`, the three
THM-553 subtriangles define restrictions

```text
d_A: a<n,       (a,b) -> (a,b),
d_B: a-b>=3,    (a,b) -> (a-1,b),
d_C: b>=2,      (a,b) -> (a-1,b-1).
```

They cover every tile, with the three pair overlaps and one triple overlap of
THM-442.  The resulting compatible-face diagram reconstructs a tiling exactly.
The target theorem will record its complement-line descent, the new middle
`d_B` node correspondence and expanded `Omega_n` coupling tensor, and the
precise information recovered by joining this `B3` chart to the folded `B2`
mirror chart.

The same theorem will separate two incidence-algebra operations that earlier
threads sometimes blended.  Boolean Möbius inversion gives exact colour atoms;
the partition-lattice third cumulant gives the genuinely new interaction.  If
`X,Y,Z` are upper-, low-, and high-face blue and
`a=beta_n,b=beta_(n-1)`, the proof in hand is

```text
P_n(x,y,z)
 = ((1-a)+ax)((1-b)+by)((1-b)+bz)
   + ab(1-b)(x-1)(y-1)(z-1).
```

Thus every proper marginal is independent and the entire chart-change debt is
the single cubic sidecar `ab(1-b)`.  The final theorem will state the exact
Legendre-word realization and its guardrails without identifying grid
reflection with the distinct complement-line involution.

Still missing at reservation time: the checked-in audit/output, finite node
and line classification tables through `n=7`, and a final proof review of the
small-size degeneracies.
