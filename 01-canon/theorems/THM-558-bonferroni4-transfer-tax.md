---
id: THM-558
title: Bonferroni4 transfer tax for one-sector insertion
status: PROVED
source: codex-2026-06-20-S62
depends_on:
  - THM-555
  - THM-556
related:
  - HYP-2691
  - HYP-2693
  - HYP-2696
  - HYP-2675
---

# THM-558 - Bonferroni4 Transfer Tax

## Statement

Work with the six inner sectors in the LRC(14) seven-sector model.  Let `P`
be a speed prefix, let `e` be a new speed, and let

```text
M_P(t) = {inner sectors missed by P at time t}.
```

Let

```text
U4(P) = p0(P) + p5(P) + 5 p6(P),
```

where `p_j(P)=meas{|M_P(t)|=j}`.  For the one-speed insertion
`P -> P union {e}`, write `mass(a -> b)` for the measure of times whose
missed-sector count changes from `a` to `b`.

Then

```text
U4(P union {e}) - U4(P)
  = mass(1 -> 0) - mass(5 -> 4) - 4 mass(6 -> 5).
```

Equivalently,

```text
Delta U4 = Delta p0 - high_tail_tax,
high_tail_tax = mass(5 -> 4) + 4 mass(6 -> 5).
```

Thus the only positive source for the Bonferroni4 upper certificate is a
genuine one-missed-sector closure, while transitions out of five- and
six-missed states pay a signed transfer tax.

## Proof

By THM-556, the level-4 Bonferroni coefficient of a state with `r` missed
inner sectors is

```text
c(r) = 1 if r=0,
       1 if r=5,
       5 if r=6,
       0 otherwise.
```

By THM-555, inserting one speed can only keep a missed-sector state fixed or
delete one missed sector.  Hence every atom has one of the transitions

```text
r -> r     or     r -> r-1.
```

The fixed transitions contribute `c(r)-c(r)=0` to `Delta U4`.  The deleting
transitions contribute `c(r-1)-c(r)`.  Checking `r=1,...,6` gives

```text
r=1:  c(0)-c(1) =  1,
r=2:  c(1)-c(2) =  0,
r=3:  c(2)-c(3) =  0,
r=4:  c(3)-c(4) =  0,
r=5:  c(4)-c(5) = -1,
r=6:  c(5)-c(6) = -4.
```

Summing atom measures proves

```text
Delta U4 = mass(1 -> 0) - mass(5 -> 4) - 4 mass(6 -> 5).
```

Finally, THM-555 also identifies `mass(1 -> 0)` with
`p0(P union {e})-p0(P)`, so the second displayed formula follows.

## Computational Check

`04-computation/lrc14_bonferroni4_transfer_tax_codex_s62.py` verifies this
identity with exact rational arithmetic on the HYP-2691 row bank.  For every
audited insertion and schedule it asserts that the local transition formula
matches the exact difference of `U4=p0+p5+5p6` computed from the common
sector-wall refinement.
