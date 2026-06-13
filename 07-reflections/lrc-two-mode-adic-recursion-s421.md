---
source: codex-2026-05-31-S421
status: exploratory synthesis
tags:
  - lonely-runner
  - two-adic-grid
  - column-families
  - tournament-recursion
  - bruhat-tits
---

# LRC in the Two Natural-Number Modes

The old natural-number picture is:

```text
odd roots:     1, 3, 5, 7, 9, ...        x+2
verticals:     odd, 2*odd, 4*odd, ...    x*2
```

Tournament history turned this into two modes:

```text
n -> 2n      row step, blowup;
n -> n+2     add-two / Mode-B column step on the odd top row.
```

For LRC this is more than an analogy.  It tells us not to put `14`, `16`, and
`18` on one linear shelf.

```text
14 = 2 * 7     first-even child of odd root 7
16 = 16 * 1    pure vertical dyadic column
18 = 2 * 9     first-even child of odd root 9
```

So `14` and `18` are column-mode neighbors.  `16` is row-mode depth.

## Boundary Debt Has the Row Seam

For initial segments `{1,...,n-1}`, the boundary witnesses are exactly the unit
endpoints, so their count is `phi(n)`.

Along an odd-root vertical chain:

```text
3 -> 6 -> 12:    phi = 2 -> 2 -> 4
7 -> 14:         phi = 6 -> 6
9 -> 18:         phi = 6 -> 6
```

The first doubling keeps the debt count fixed.  After that, the debt count
doubles each row step.  In density form, `phi(n)/n` halves at the first
doubling and then stabilizes.

That is the LRC version of the tournament seam: row zero to row one is not the
same as the later row steps.

## The First-Even Children

The designed first-even quotient ladders give:

```text
n=6:   gap/th=1/20,  debt=12,  product=3/5
n=10:  gap/th=1/56,  debt=48,  product=6/7
n=14:  gap/th=5/924, debt=84,  product=5/11
n=18:  gap/th=1/176, debt=176, product=1
```

The `n=14` row still has the tiniest visible gap among these, but `n=18` has a
cleaner product and larger frontier.  This matches the previous "favorite
amount of runners" instinct: `18` is not merely larger; it is the next
column-mode payload after `14`.

## What This Suggests

The proof should not ask for one universal metric too early.

Row-mode proof objects:

```text
dyadic depth;
Bruhat-Tits frontier mass;
harmonic-flow divergence;
gap-debt product.
```

Column-mode proof objects:

```text
odd payload;
product-building mass;
small-denominator sieve transfer;
odd-root x+2 transfer matrix.
```

The pure `n=16` case is the row-mode lab.  The `n=14 -> n=18` passage is the
column-mode lab.  A full LRC proof probably needs to show these two modes
commute without losing obstruction mass: row motion cannot erase debt, and
column motion can only move debt into a new odd-payload tree.
