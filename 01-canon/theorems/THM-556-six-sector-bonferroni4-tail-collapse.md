---
id: THM-556
title: Six-sector Bonferroni level-4 collapses to the five-six missed tail
status: PROVED
source: codex-2026-06-20-S59
depends_on: []
related:
  - THM-555
  - HYP-2675
  - HYP-2691
  - HYP-2692
  - HYP-2693
---

# THM-556 - Six-Sector Bonferroni4 Tail Collapse

## Statement

Let `M` be a random subset of six labelled sectors, and write

```text
p_t = Pr(|M|=t),   0<=t<=6.
```

For `r>=0`, let

```text
S_r = E binom(|M|,r) = sum_{t=r}^6 binom(t,r) p_t.
```

Then the fourth Bonferroni upper expression for the empty-set probability is

```text
U4 = 1 - S_1 + S_2 - S_3 + S_4
   = p_0 + p_5 + 5 p_6.
```

In particular,

```text
Pr(M=empty) = p_0 <= U4,
```

with exact slack `p_5+5p_6`.

For the LRC(14) seven-sector model, taking `M(t)` to be the missed set of the
six inner sectors gives the same identity on every speed row.

## Proof

Condition on `|M|=t`.  The coefficient of `p_t` in `U4` is

```text
c_t = sum_{r=0}^4 (-1)^r binom(t,r).
```

For `t=0`, this is `1`.  For `1<=t<=4`, this is the full binomial sum
`(1-1)^t`, so it is `0`.

For `t=5`, the full sum through `r=5` is zero, hence

```text
c_5 = -(-1)^5 binom(5,5) = 1.
```

For `t=6`, the missing terms are `r=5,6`, so

```text
c_6 = -[(-1)^5 binom(6,5) + (-1)^6 binom(6,6)]
    = -[-6 + 1]
    = 5.
```

Thus only `p_0`, `p_5`, and `p_6` survive, with coefficients `1,1,5`.

## Computational Check

`04-computation/lrc14_truewide_bonferroni4_gate_codex_s59.py` asserts this
identity on every named row and every row in the exact scan boxes used for
HYP-2693/T929.
