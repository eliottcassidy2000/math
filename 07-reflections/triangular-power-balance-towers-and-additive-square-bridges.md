---
source: codex-2026-06-12
status: SYNTHESIS + verified computation; HYP-2454
tags: [triangular, consecutive-sums, power-balance, squares, pell, LRC14, code72, support-gate, tournament-analysis]
---

# Triangular power-balance towers and additive/square bridges

The new prompt adds something unusually clean to the triangular bridge story:
two towers with the same interval shape, but different centers.

The ordinary tower centers at `2T_n`:

```text
1+2 = 3,
4+5+6 = 7+8,
9+10+11+12 = 13+14+15.
```

The square tower centers at `4T_n`:

```text
3^2+4^2 = 5^2,
10^2+11^2+12^2 = 13^2+14^2,
21^2+22^2+23^2+24^2 = 25^2+26^2+27^2.
```

The shared schema is:

```text
left  = C-n, ..., C
right = C+1, ..., C+n.
```

The first equality is `p=1`, the second is `p=2`.  This is a power-balance
problem, not merely a cute pair of sequences.  The defect polynomial

```text
D_p(C,n)=sum_{i=0..n}(C-i)^p - sum_{i=1..n}(C+i)^p
```

has exact integer roots `C=2T_n` for `p=1` and `C=4T_n` for `p=2`.  The script
then shows, through `p<=8,n<=40`, that for `p>=3` the positive root sits
strictly between `2pT_n` and `2pT_n+1`.  So the wild thought is: ordinary sums
and square sums may be the only exact integer-center members of this particular
left-heavy/right-light balance family.

## Why this touches the older triangular bridge

HYP-2128 says:

```text
8*C(n,2)+1=(2n-1)^2.
```

That made `2n-1` the odd-square face of the triangular pair count, exactly the
LRC worry modulus at `n=14` (`27`).  HYP-2454 adds another relation between
addition and multiplication:

```text
p=1 balance center = 2T_n,
p=2 balance center = 4T_n.
```

The factor of `2` in the center is the visible multiplication step.  But the
square tower has an ordinary additive shadow that remembers the defect:

```text
L_2(n)-R_2(n)=2T_n,
L_2(n)+R_2(n)=4S_1(n).
```

This is exactly the kind of scalar/support split the repo keeps rediscovering.
The squared equality is the scalar closure; the unsquared defect is the
retained address.

## The row that matters most: Q(3)

The prompt noticed `21+22+23+24`.  The computation makes it precise:

```text
Q_L(3) = [21,22,23,24] = F_R(4).
```

This is the only full side equality found through the search range, and it is
probably unique by a one-line length/endpoint argument.

At that same row:

```text
ordinary sum of Q_L(3) = 90 = S_1(4),
ordinary sum of Q_R(3) = 78 = C(13,2),
square sums both equal 2030.
```

That is the repo's new beacon.  The number `78` is not floating here; it is
already the `lambda_5` in the hypothetical `5-(72,16,78)` minimum-word design,
and HYP-2445 also saw it as the `D7` subgroup index.  The new tower says:
do not stare at `78` alone.  The adjacent square-balanced partner is `90`, and
the equality binds them through the side `[21,22,23,24]`.

The careful reading is not "we found the code."  It is:

```text
78 is a right-side additive shadow,
90 is the matching left-side additive shadow,
the squared equality is a scalar closure,
and any actual code route must retain support/incidence data.
```

This is perfectly aligned with the current `[72,36,16]` program, where scalar
Gleason feasibility is already solved and the hard part is the hidden
support/design lift.

## Crossover geometry

The first tower partitions square shells:

```text
F(n) = [n^2, ..., (n+1)^2-1].
```

The second tower runs through intervals centered at `2n(n+1)`, so its endpoints
cross square-shell boundaries according to Pell-like equations:

```text
2n^2+n     = m^2,
2n^2+2n   = m^2+m,
2n^2+2n+1 = m^2,
2n^2+3n   = m^2+m,
```

and their right-boundary variants.  The stored output finds endpoint
coincidences at `36,90,144,210,420,421,840,841,4900,14280,14281`, among
others.  This gives a serious next computation: solve the Pell families and
ask whether overlap types are eventually periodic in a quadratic unit.

## How to use it next

For LRC14, treat the unsquared defect as a resource ledger.  The obvious
numbers to keep together are:

```text
27 = 2*14-1, the HYP-2128 worry shell,
78 = Q_R(3) ordinary shadow and code lambda_5,
90 = Q_L(3) ordinary shadow and S_1(4),
91 = C(14,2), the HYP-2445/HYP-2444 q=91 bridge.
```

For irreducible polynomials and HYP-2452, the bridge is conceptual: a visible
coefficient vector is a boundary total; a factorization is a hidden grid lift.
Here a visible equality of square sums is a scalar total; the ordinary defect
and endpoint address are the hidden data that make the equality informative.

For the 72-code, the corresponding task is to turn:

```text
weight enumerator coefficient -> boundary total
minimum-word supports          -> hidden incidence lift
78/90 shadow                   -> candidate local address
```

into an actual constraint on the `5-(72,16,78)` design ledger.

## Assumption challenge

I initially wanted the vertices to be the numbers in the intervals.  That is
probably too literal.  The better Tournament Analysis vertices are proof
obligations:

```text
power center is integral,
endpoint aligns with a square-shell boundary,
ordinary defect is triangular,
78/90 support shadow is retained,
hidden incidence lift exists.
```

This quotient preserves the things likely to survive into a proof and destroys
the decorative interval entries.  That is the right loss function for the next
session.
