---
source: oracle-2026-05-31-S21
status: proof-search synthesis
tags:
  - lonely-runner
  - n14
  - n15
  - n16
  - endpoint-debt
  - set-cover
  - hall-deficit
---

# Gate frontiers for n=14, n=15, and n=16

This pass tried to avoid another broad search and instead isolate a shared
finite obstruction across the three current Lonely Runner frontiers
`n=14,15,16`.

The setup is forced by THM-369: any counterexample must pay every small
denominator row `2 <= q <= n`, in particular with an `n`-gate. The `n`-gate then
owns `2n` endpoints. The local question is: how cheaply can those endpoints be
strictly protected, and what invoices remain unpaid?

## Local gate rows

Using only lower candidates `1..n-1`, the exact cover sizes for the `n`-gate
endpoints are:

```
n=14: 8 lower speeds
n=15: 10 lower speeds
n=16: 9 lower speeds
```

The lower covers have many private endpoints. This recovers the familiar
private-leaf flavor, especially the `n=16` nine-cover from THM-367, but it is not
the whole story because protectors can borrow room from the next residue window.

Using candidates `1..2n-1`, the local covers become much cheaper:

```
n=14: 7 speeds
n=15: 6 speeds
n=16: 5 speeds
```

That looks dangerous until the coarse rows and gcd are included.

## Minimum-window family audit

Enumerating every minimum window cover gives the useful obstruction:

- `n=14`: there are `24` minimum covers. The `6` that pay every coarse divisor
  row are all gcd-`2` even covers, so primitivity forces an extra breaker.
- `n=15`: all `32` minimum covers are primitive, but every one leaves at least
  two coarse divisor rows unpaid.
- `n=16`: the minimum cover is unique, `8,18,22,26,30` plus the `16`-gate; it is
  gcd-`2` and still misses rows `7,12,14`.

So the apparent local repair always has a second invoice: divisor debt, gcd
debt, or both.

## Quotient ladders

The quotient-ladder products line up with the existing gap-debt story:

```
n=14: scale 7/14 products = 5/11
n=15: scale 3 product = 7/11; scale 5 has smaller visible gap but product 7/10
n=16: products = 34/33,34/33,35/33,35/33
```

The `n=15` row is the new comparison point. It behaves neither like the scalar
`n=14` moat nor like the pure dyadic `n=16` row: it is a mixed odd-prime product
building, with the `3`-branch giving the best product and the `5`-branch giving
the smaller visible gap.

## Proof target

The plausible next certificate is a Hall/dual inequality. Put weights on:

1. gate-owned endpoint rows;
2. coarse denominator rows `q <= n`;
3. primitivity/gcd breaker rows.

The target inequality: every speed column pays at most one unit of dual weight,
but the total required row weight exceeds the `n-1` available columns unless a
positive gap or unprotected endpoint remains. This is the exact set-cover
version of the recent adic/product-flow picture.
