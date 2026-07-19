---
id: THM-1237
title: POSITIONED SEVEN-WALL HUNTER TRANSFER — global strict-spectrum credit survives on an interval after an explicit harmonic and gcd discrepancy debt
status: RESERVED/PROOF IN PROGRESS. A candidate exact interval-Hunter inequality has been derived by composing THM-1166 with THM-1221; discrepancy signs, constants, the optimizing tree, and the protected-needle contrapositive are under audit.
source: codex-2026-07-19-S78 continuation with path-inequality agent
depends_on: [THM-1166, THM-1221]
related: [THM-1219, THM-1226, HYP-7870]
---

# THM-1237 — positioned seven-wall Hunter transfer

The candidate local inequality for an interval `I` of length `L`, seven wall
speeds `s_i`, and a spanning tree `T` is

```text
mu(I intersect Safe_7)
 >= L sum_(ij in T) rho_ij
    - (6/49) sum_i 1/s_i
    - sum_(ij in T) rho_ij(1-rho_ij)/gcd(s_i,s_j).
```

Here `rho_ij` is the exact global overlap mass from THM-1166.  THM-1221 can
supply a tree with total weight at least `15/154`; the intended consumer is a
fully explicit sufficient condition making the right side positive on a
particular protected needle.

This file reserves the namespace honestly.  It does not yet claim that every
LRC(14) packet pays sufficiently small positioning debt, and it does not claim
LRC(14).
