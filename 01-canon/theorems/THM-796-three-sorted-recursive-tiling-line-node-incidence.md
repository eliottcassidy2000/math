---
id: THM-796
title: Three-sorted recursive incidence of tilings, complement lines, and converse-merged tournament nodes
status: RESERVED — exact two-face tiling pullback, complement quotient, incidence-matrix conservation, and recursive colour-transition formulas derived; exhaustive node-kernel audit through n=7 in progress
source: codex-2026-07-15-S9
depends_on: [THM-781, THM-785, THM-790, THM-793]
related: [HYP-6825, HYP-6870]
---

# THM-796 — three-sorted recursive metagraph incidence

Namespace reservation for the owner's requested simultaneous treatment of
tilings, blue/black complement lines, and converse-merged iso-class nodes as
`n` changes.

The exact tiling-level recursion already derived is

```text
T_n = (T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex.
```

The two faces delete the two endpoints of the fixed Hamiltonian path, agree on
their common interior, and the missing apex bit reconstructs the unique upper
tiling.  Complement acts componentwise on both faces and flips the apex, so
the line set is the free `C2` quotient of this pullback.  Anti-diagonal
reflection swaps the two faces.

Still being audited: the exact matrix identities after projection to merged
nodes, colour descent between line levels, uniqueness of the weighted endpoint-
face branching rows through `n=7`, and the precise information lost when one
passes from tilings to nodes or from line instances to simple node adjacency.
