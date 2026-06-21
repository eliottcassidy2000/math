---
id: HYP-2701
title: LRC14 true-wide cap is a survival middle-mass gate
status: OPEN; namespace reserved for codex-2026-06-20-S64
source: codex-2026-06-20-S64
tangent: T936
depends_on:
  - HYP-2695
  - HYP-2699
  - THM-556
  - THM-535
related:
  - HYP-2700
  - HYP-2698
  - HYP-2697
  - HYP-2696
  - HYP-2684
  - HYP-2675
---

# HYP-2701 - True-Wide Survival Middle-Mass Gate

## Claim Being Tested

For a row `E`, let `N_E(x)` be the number of missed inner sectors in
`{1,...,6}` and write its distribution as `p_t=Pr(N_E=t)`.
THM-556 gives

```text
U4(E) = p0 + p5 + 5 p6.
```

Equivalently,

```text
U4(E) <= floor_k=(k-6)/7
```

is the exact survival middle-mass inequality

```text
p1 + p2 + p3 + p4 - 4 p6 >= (13-k)/7.
```

Thus HYP-2695's true-wide floor gate is not primarily a raw Bonferroni
statement.  It asks whether true-wide geometry forces enough middle missed
sector mass to pay for the fully missed tail `p6`.  The exact cap only changes
the right side from `(13-k)/7` to `1-cap_k`, explaining why `k=8` can spend
the THM-535 dividend while `k>=9` should not need it.

## Planned Scout

The S64 scout will scan true-wide rows using exact Fraction arithmetic and
report:

1. the survival currency `M=p1+p2+p3+p4-4p6`;
2. floor slack `M-(13-k)/7`;
3. cap slack `M-(1-cap_k)`;
4. the row fingerprints that make floor slack small;
5. a Tournament Analysis quotient whose vertices are proof obligations and
   row families, not runners or arcs.

This hypothesis deliberately challenges the assumption that the final true-wide
observable should be a single scalar such as `p0`, `U4`, additive energy, or a
binary tournament statistic.  The preserved predicate is the cap/floor
implication; the quotient destroys sector labels except through the missed
count distribution, so any near-equality row must still be lifted back to
HYP-2696/THM-558's transfer-tax and HYP-2698's residual-profile coordinates.

## Status

This file reserves the namespace before the exact S64 computation.  No proof is
claimed yet.
