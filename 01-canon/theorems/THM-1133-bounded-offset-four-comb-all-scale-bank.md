---
id: THM-1133
title: Exact finite complement for all four-comb offset shapes of span at most thirty
status: RESERVED / COMPUTING — THM-1129 leaves exactly 3,539,936 legal pre-tail rows outside THM-1123; exact endpoint subtraction is running. No closure is claimed until the frozen referee reports zero failures
source: codex-2026-07-18-S73 bounded-offset continuation
depends_on: [THM-1123, THM-1129]
related: [THM-1126, THM-1127, THM-1128, MISTAKE-164]
---

# THM-1133 — bounded-offset finite complement

Namespace reservation for the exact complement left by THM-1129.  It will
test every legal triple

```text
(P,A,K),
P subset {1,...,12}, |P|=8,
A={0<a<b<c<=30},
13 max(P)<K<832,
```

not already contained in THM-1123's first-forty bank.  THM-1129's exact
ledger counts `3,539,936` such rows.  The target is

```text
7(K+c)L(P,K+A)>1.
```

This page is an honest reservation, not a theorem.  It will be promoted only
after exact endpoint subtraction, an independent replay or guard, and frozen
hashes are present.
