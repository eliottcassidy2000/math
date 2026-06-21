---
id: HYP-2717
title: LRC14 top-character carrier relation filter for the moderate multi-block gap
status: OPEN; analytic proof order
source: codex-2026-06-21
depends_on:
  - HYP-2716
  - HYP-2715
  - HYP-2714
  - HYP-2676
  - HYP-2645
related:
  - THM-561
  - HYP-2718
  - HYP-2708
  - HYP-2707
  - THM-546
  - THM-547
  - OPEN-Q-108
---

# HYP-2717 - Top-Character Carrier Relation Filter

## Claim

The HYP-2714 moderate-span multi-block lemma should be proved by expanding the
HYP-2716 top character in carrier Fourier modes and splitting those modes by
the carrier relation lattice.

For a split row

```text
E = {0} union_i (M_i+B_i)
```

let `F(x,theta_1,...,theta_g)` be the all-six-sector cover indicator generated
by the anchor and the coherent block shapes `B_i` at slow coordinate `x`, with
carrier phases `theta_i`.  Its carrier product is

```text
Product(E) = int_x int_{T^g} F(x,theta) dtheta dx,
```

while the actual row is

```text
p0(E) = int_0^1 F(x, M_1 x, ..., M_g x) dx.
```

Expand in carrier characters:

```text
F(x,theta) = sum_{n in Z^g} Fhat_n(x) exp(2 pi i n.theta).
```

Then the HYP-2716 top-character error has the exact formal shape

```text
p0(E)-Product(E)
  = sum_{n != 0} int_0^1 Fhat_n(x) exp(2 pi i (n.M)x) dx.
```

The proof split is:

```text
exact carrier relations:     n.M = 0
nonresonant carrier modes:   n.M != 0.
```

The exact carrier relations are unavoidable when `g>=2`: an integer vector
`M=(M_1,...,M_g)` always has a rank `g-1` integer kernel.  Thus the phrase
"the carriers become independent" is not literally equidistribution on `T^g`.
It means that the exact carrier relations have large height and therefore small
Fourier coefficients, while nonresonant modes oscillate with denominator
`|n.M|`.

## Proposed Bound

For the top character `M_6` from HYP-2716, prove a bound of the form

```text
|M_6(E)|
 <= sum_{0 != n in Lambda(M)} A_n(B_1,...,B_g)
  + sum_{n notin Lambda(M)} B_n(B_1,...,B_g)/|n.M|
 <= cap_k - Product(E),
```

where `Lambda(M)={n in Z^g : n.M=0}` and the coefficients are Fourier/BV
budgets for the cover indicator.  Low-height exact relations and small
`|n.M|` nonresonances should be routed to the HYP-2714 finite ledger; the
analytic lemma only needs the high-height tail.

This is the carrier analogue of the older relation-lattice lesson: spread is a
shadow, relation height is the object.  The HYP-2716 Krawtchouk quotient tells
which scalar to expand; HYP-2717 tells which carrier modes can contribute to
that scalar.

## Why This Is Sharper Than Raw Weyl

A naive multidimensional Weyl statement would compare the line

```text
x -> (M_1 x, ..., M_g x)
```

to the full torus.  That comparison is false at the level of exact Fourier
modes because `Lambda(M)` is nontrivial.  The right statement is a filtered
one:

1. Exact relations survive the line integral, but their height is large in the
   separated-carrier branch.
2. Nonrelations pay `1/|n.M|` after one-dimensional oscillation or BV
   integration by parts.
3. The top-character quotient from HYP-2716 is the only scalar that must fit
   under `cap-product`, so lower Krawtchouk shadows do not consume budget.

This also explains why the tested anchor-separated product can underestimate
actual cover: exact high-height carrier relations can have either sign in the
top character.  Product is a main term, not an envelope.

## Assumption Challenge and Tournament Analysis

This session considered vertices as runners, blocks, carrier phases, residual
masks, Krawtchouk weights, Fourier modes, exact relations, nonresonant modes,
and proof obligations.  For HYP-2717 the productive vertices are resonance
classes of carrier modes:

```text
exact relation modes > small-denominator nonrelations > high-denominator tail
```

The quotient preserves the top-character cover error and relation-height
budget, while destroying sector ownership and residual-mask location.  Those
lost labels must be recovered only in the finite HYP-2714 resonance ledger.

Pairwise observable: larger expected contribution to `|M_6|` before cap slack.
Switch/gauge: carrier Fourier basis after the Krawtchouk top-character quotient.
The proof-pressure tournament is transitive by construction:

```text
low-height exact relations
> small |n.M| nonrelations
> high-height exact-relation tail
> high-denominator nonrelation tail.
```

Challenged assumption: multi-block decorrelation is full-torus
equidistribution.  Corrected assumption: it is a carrier-relation tail bound.

## Next Step

Make the coefficient budgets explicit for coherent block shapes:

```text
A_n(B_1,...,B_g),   B_n(B_1,...,B_g).
```

The likely route is a bounded-variation or Selberg-Beurling smoothing estimate
for the top-character cover indicator, paired with the exact finite
low-height ledger already isolated by HYP-2714.
