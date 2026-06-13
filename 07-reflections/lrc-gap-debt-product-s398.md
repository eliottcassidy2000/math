---
source: codex-2026-05-31-S398
status: exploratory formalization
tags:
  - lonely-runner
  - gap-debt
  - product-law
  - dyadic
  - cayley-dickson
---

# LRC Gap Debt as an Adelic Product

The user's sentence is the right compression:

```text
gap is Archimedean size;
endpoint debt is 2-adic size;
their product is conserved.
```

That says the obstruction is not disappearing when the visible lonely interval
shrinks.  It is changing absolute values.

In the finite endpoint language, a counterexample is much sharper than "small
gap."  By THM-357, it must have full forbidden measure and every forbidden
endpoint must be strictly protected.  So it needs:

```text
Gap_A = 0
Debt_2 = 0
```

at the same time.  The product law says these two zeroes are incompatible on
the repair branch.

## The Two Zeroes

There are three ways an LRC attempt can fail to be a counterexample.

```text
Gap_A > 0      visible lonely interval;
Debt_2 > 0     exposed endpoint witness;
Gap_A*Debt_2   positive conserved obstruction.
```

The boundary-only initial segments are a useful warning.  For `{1,...,n-1}`,
the visible gap is zero, but unit endpoints survive.  In S398:

```text
n=14 initial: Gap_A=0, Debt=6
n=16 initial: Gap_A=0, Debt=8
n=18 initial: Gap_A=0, Debt=6
```

So the product lower bound is not the boundary theorem.  The boundary theorem
is simpler: debt alone wins.

The product lower bound starts after a repair has closed the obvious boundary
endpoints.  That is where the old unit witness is no longer visible, and the
question is whether it was annihilated or exported.

## The Export Rows

The export rows answer: exported.

For `n=14`, the seven-ladder to fourteen-ladder step halves the real gap and
doubles the endpoint debt:

```text
5/924  * 84  = 5/11
5/1848 * 168 = 5/11
```

For `n=18`, the mixed square-torsion row does the same thing:

```text
1/176 * 176 = 1
1/352 * 352 = 1
```

For `n=16`, the pure dyadic row is almost too informative:

```text
1/33  * 34  = 34/33
1/66  * 68  = 34/33
1/132 * 140 = 35/33
1/264 * 280 = 35/33
```

The first two rows share one product, the last two share another, and the jump
is exactly `35/34`.  That jump smells like the first correction term in the
true 2-adic debt norm.  Raw endpoint count is seeing the product, but with a
layer transition tax.

## q(x) Shadow

The hyperbola thread gives a helpful picture.  A conserved product is a
rectangular hyperbola:

```text
gap * debt = constant.
```

Moving left toward smaller real gap forces movement upward into larger
denominator debt.  This is the discrete LRC analogue of replacing a visible
triangle with a `dx/x + q` finite-part packet.  You can carve the area one way
or another, but the finite part is not erased.

Here the finite packet is endpoint debt.  The real gap is the visible triangle.
The denominator layer is the reciprocal counterterm.

## Cayley-Dickson Shadow

The Cayley-Dickson analogy becomes cleaner in this language.  Doubling does not
keep the same coordinates honest.  It keeps a norm-like obstruction honest
after the coordinates change.

At `n=16`, the `16`-gate protects the old odd unit endpoints, so it looks like
a repair.  But that repair is a zero-divisor move: it kills one coordinate and
exposes another.  The debt moves to dyadic descendant endpoints.  The product
ledger says the move was not free.

This suggests a proof strategy with a very concrete target:

```text
Define the correct weighted endpoint debt D_2.
Show every n=16 repair branch has Gap_A * D_2 >= c > 0.
Conclude that Gap_A = 0 and D_2 = 0 cannot happen.
```

The current raw-count audit may already be enough to point to the right
weights.  The `35/34` phase step is the place to stare.
