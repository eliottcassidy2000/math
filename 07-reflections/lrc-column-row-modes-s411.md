---
source: codex-2026-05-31-S411
status: exploratory reflection
tags:
  - lonely-runner
  - natural-numbers
  - column-family
  - dyadic-row
  - product-sum
---

# LRC Column And Row Modes

The old picture was:

```text
natural numbers = odd numbers moving by +2, plus their even *2 chains
```

In the tournament language this became:

```text
Mode B / n -> n+2       move along the odd top row
blowup / n -> 2n        move down a doubling row
```

S411 asks what that means for Lonely Runner denominators.

## The Unit Seam

The first exact answer is small and clean:

```text
initial segment {1,...,n-1} leaves exactly phi(n) unprotected endpoints.
```

S411 verifies this through `n=40`, and it is theorem-shaped.  The witnesses
are the unit residues `a/n`.  So the LRC initial segment is not merely a
Dirichlet equality example; it is the unit skeleton of the denominator.

Now write:

```text
n = 2^r * b,  b odd.
```

Along a fixed odd core `b`:

```text
phi(2b) = phi(b)
phi(2^(r+1)b) = 2 phi(2^r b),  r >= 1.
```

That is the LRC analogue of the tournament pair anomaly.  The first doubling
does not double the unit skeleton.  It creates new nonunit residue room while
leaving the old unit witnesses alone.  After that, the row really doubles.

For `b=7`:

```text
n=7:   phi=6
n=14:  phi=6      first seam
n=28:  phi=12     honest doubling
```

This gives a better sentence for the fourteen-runner case:

```text
n=14 is not just composite.  It is the first even row over the prime core 7:
the same six unit witnesses, plus a new divisor channel 7.
```

## The Row Parent

For every even denominator,

```text
lpd(n) = n/2.
```

That sounds obvious, but in the LRC quotient-ladder scans it is the structural
fact that matters.  The largest-proper-divisor ladder is the first visible
multiplicative repair channel, and for even `n` it is literally the row parent
in the `2^r*b` grid.

S411 records:

```text
n=14: lpd=7,  gap/th=5/924,   debt=84
n=16: lpd=8,  gap/th=1/132,   debt=140
n=18: lpd=9,  gap/th=1/176,   debt=176
n=20: lpd=10, gap/th=1/325,   debt=180
n=24: lpd=12, gap/th=11/4080, debt=264
```

So the row-mode proof should not ask vaguely about "some divisor."  It should
ask how endpoint debt transfers from `n/2` to `n`.

This is where S398 and S399 live:

```text
S398: Archimedean gap * endpoint debt stays bounded in row repairs.
S399: dyadic endpoint law is a conservative Bruhat-Tits tree kernel.
S400/S410: Bruhat-Tits frontier/descent turns row debt into p-adic depth.
```

Those are row-mode invariants.

## The Top Row

Odd denominators do not have a row parent.  They move by `+2` along the top row.
There the relevant changes are:

```text
phi(n) changes
prime/composite status changes
product-sum factor packing changes
odd divisor channels appear or vanish
```

S411 overlays the product-sum minimum at arity `k=n-1`:

```text
n=14: k=13, min product 18, seed 1^10 + 2*3*3
n=15: k=14, min product 20, seed 1^11 + 2*2*5
n=16: k=15, min product 24, seed 1^13 + 3*8
n=18: k=17, min product 24, seed 1^14 + 2*2*6
n=22: k=21, min product 27, seed 1^18 + 3*3*3
```

This is not random decoration.  Product-sum cores are where additive slack and
multiplicative factor packing collide.  In LRC language, they indicate which
coordinates should be most fragile after scalar quotienting.

So the top-row proof should look different from the row proof.  It should
control unit skeleton changes and product-sum target coordinates as `b` moves
to `b+2`.

## The Two-Dimensional Proof Shape

The current LRC research split becomes:

```text
top-row mode:
  odd b -> b+2
  unit skeleton, product-sum cores, prime/composite jumps

row mode:
  n -> 2n
  lpd=n/2, endpoint-debt transfer, dyadic/Bruhat-Tits flow
```

The proof should probably induct on the grid, not on the line.

For example:

```text
n=14 = 2*7      first prime-core seam
n=16 = 16*1     pure dyadic row lab
n=18 = 2*9      first square-odd-core seam
```

These are adjacent in the ordinary integer line, but not adjacent in the
structural grid.  Treating them as neighboring cases may be why the problem
keeps changing costumes.  In the grid, each is a different experiment.

## Next Move

Add the following to future LRC branch-and-bound states:

```text
(row r, odd core b)
phi(n)
seam defect = 2 phi(n/2) - phi(n)
lpd(n) and lpd mode
row-parent endpoint debt
product-sum core at arity n-1
```

Then searches can branch by mode:

```text
top-row branch: product-sum/unit skeleton obstruction
row branch: endpoint-debt transfer obstruction
```

That is the closest LRC analogue so far to the tournament recursion story.
