---
source: codex-2026-05-31-S453
status: integration note
tags:
  - lonely-runner
  - tournaments
  - two-adic-grid
  - first-doubling
  - parity-seam
---

# LRC and the First-Doubling Tournament Seam

The old repo picture is a two-mode grid:

```text
odd roots:  1, 3, 5, 7, 9, ...       x+2
columns:    m, 2m, 4m, 8m, ...       x*2
```

The new point of focus is the first row step:

```text
m odd -> 2m.
```

This is the critical component.  It is not an ordinary doubling.

## The Arithmetic Seam

THM-371 records the exact arithmetic:

```text
U(2m) ~= U(m)          for m odd
phi(2m) = phi(m)
phi(4m) = 2 phi(2m)
```

The first doubling chooses a single parity sheet.  Every unit `a mod m` has
two lifts modulo `2m`, namely `a` and `a+m`, but exactly one is odd.  That is
the one that survives as a unit modulo `2m`.

So the LRC unit skeleton is not doubled at the first seam.  It is copied into
the odd sheet.

Read backward, the reduction `2m -> m` is the exact even-to-odd critical
phase.  It collapses the two residue sheets, but units remember which sheet
they came from: every surviving unit is odd.  The other sheet is not
unit-witness mass; it is quotient room.  This is why an even first-child
denominator should not be simplified to its odd parent without keeping a
ledger of what was lost in the collapse.

The tournament matching ledger does the complementary thing:

```text
floor(2m/2) = 2 floor(m/2) + 1.
```

The odd base had one unmatched vertex.  First doubling gives that vertex a
twin, creating exactly one extra pair.  Later row steps double honestly.

## SC Blowup Is the Right Tournament Analogy

There are two tournament doublings:

```text
lex blowup: old hierarchy is magnified
SC blowup:  every vertex gains strong/weak twins
```

The SC blowup is the useful model for LRC.  It keeps both `T` and `T^op` in
the lane/cross-lane data, but it erases score variation:

```text
every v_0 has score n
every v_1 has score n-1
```

This is exactly what the first seam should feel like: a hidden old object plus
a new balancing sheet.  The old tournament is still there, but ordinary score
no longer sees it.

For LRC, the old unit skeleton is still there too, but it now sits as the odd
lift inside a larger denominator.  A candidate proof should not search for a
new unit obstruction.  It should ask what happens when a repair spends the new
nonunit quotient room.

## The LRC Reading

For a first-even denominator `n=2m`, the endpoint problem has three layers:

```text
1. inherited unit layer: U(2m) ~= U(m)
2. new quotient room:   divisibility by m and even/nonunit columns
3. exported debt:       endpoints pushed below or beyond the row parent
```

This clarifies the earlier `14,16,18` split:

```text
14 = 2*7      first prime-core seam
16 = 16*1     pure dyadic row lab
18 = 2*9      first square-core seam
```

`n=14` and `n=18` should be compared as first-even children of adjacent odd
payloads.  `n=16` should not be forced into that comparison; it is what happens
after the row mode keeps doubling past the first seam.

The S453 computation gives the row-parent ladder pressure:

```text
n=14: gap/th=5/924, debt=84,  product=5/11
n=18: gap/th=1/176, debt=176, product=1
n=22: gap/th=2/561, debt=260, product=520/561
```

These are all `lpd(n)=n/2` rows.  The visible quotient ladder is the row
parent itself.

## Proof Reframe

The proof target becomes:

```text
The first seam gives no new unit witnesses.
Therefore any counterexample must use nonunit quotient room.
But nonunit quotient room exports endpoint debt to row-parent/product-depth
layers.
```

This is the LRC analogue of SC blowup:

```text
tournament: score variation disappears, old data survives in hidden lanes
LRC:        unit count does not grow, old data survives in odd lifts
```

The right invariant is not just "how many endpoints are exposed."  It is:

```text
which unit endpoints are inherited,
which repairs spend the new parity sheet,
and where the debt is exported after that spend.
```

## Practical Next Step

For `n=14`, build a local matrix with rows:

```text
odd-lift unit endpoints from U(7)
14-gate fan endpoints
row-parent quotient endpoints owned by 7
exported product-depth endpoints
```

Then search for a dual certificate where inherited unit rows force the gate,
the gate forces the fan, and the fan cannot pay the exported row-parent debt
without reopening a small denominator.  The same experiment at `n=18` should
show whether square-core payload makes the first seam easier than `n=14`.
