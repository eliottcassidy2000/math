---
id: THM-785
title: Cyclic-triangle transitivity flow on the blue/black merged metagraph
status: RESERVED — score-energy coordinate and complement-line flux identity derived; exact category-flow census and general symmetry statements still being audited
source: codex-2026-07-14-S9
depends_on:
  - HYP-6825
  - THM-781
related: [THM-550, THM-646]
---

# THM-785 — Cyclic-triangle flow on the blue/black merged metagraph

Namespace claim for the owner's requested organization of merged nodes from
the transitive class to the most balanced score classes.

Known at reservation time: if `d_i` are tournament scores and
`mu=(n-1)/2`, then

```text
sum_i (d_i-mu)^2 = n(n^2-1)/12 - 2 C_3(T),
```

so cyclic-triangle count is an exact transitivity-depth coordinate and its
maximum locus is precisely the regular/near-regular score locus.  For an
explorer tiling `t` and its all-tile complement `bar(t)`, preliminary algebra
gives

```text
C_3(bar(t))-C_3(t)=d_0(t)-d_(n-1)(t)-1.
```

Still missing at reservation time: full proof, blue/black projected-line
orientation census through `n=7`, pure-blue/mixed/pure-black flow laws,
quantification of left-right symmetry versus black drift, and separation of
general theorems from finite-atlas observations.
