---
id: THM-358
name: lrc-initial-segment-unit-skeleton
status: PROVED
date: 2026-05-30
session: codex-2026-05-30-S362
depends_on:
  - THM-357
  - HYP-1810
---

# THM-358: Initial-Segment Lonely Runner Unit Skeleton

## Statement

Let `n >= 2` and take the initial-segment speed set

```text
V_n = {1, 2, ..., n-1}
```

with Lonely Runner threshold `1/n`.  Then the safe set is exactly

```text
{ a/n : 1 <= a <= n-1 and gcd(a,n)=1 }.
```

Equivalently,

```text
for every t in R/Z, either
  ||v t|| < 1/n for some 1 <= v <= n-1,
or
  t = a/n with gcd(a,n)=1.
```

Every point `a/n` with `gcd(a,n)=1` is an unprotected forbidden endpoint, and
there are no other unprotected endpoints.

## Proof

Consider the `n` points

```text
0, t, 2t, ..., (n-1)t
```

on the circle.  If two of these points have circular distance strictly less
than `1/n`, then their difference is `v t` for some `1 <= v <= n-1`, and

```text
||v t|| < 1/n.
```

Thus `t` is forbidden.

If no two of the `n` points have distance strictly less than `1/n`, then the
`n` circular gaps between them are all at least `1/n`.  Their sum is `1`, so
every gap is exactly `1/n`.  Hence the set of points is exactly the regular
`n`-gon

```text
{0, 1/n, 2/n, ..., (n-1)/n}.
```

Since `t` is one of these points, `t=a/n` for some `a`.  The `n` displayed
multiples are distinct, so `a` is a unit modulo `n`; otherwise the orbit of
`a/n` under multiplication by `0,1,...,n-1` would have fewer than `n` points.
Therefore `gcd(a,n)=1`.

Conversely, if `t=a/n` with `gcd(a,n)=1`, then for every `1 <= v <= n-1` the
residue `a v mod n` is nonzero, so

```text
||v t|| >= 1/n.
```

Thus these and only these points are safe.

Finally, each safe point `a/n` is a forbidden endpoint.  Since `a` is a unit
modulo `n`, choose `v` with `a v == 1 mod n` or `a v == -1 mod n`.  Then

```text
||v * a/n|| = 1/n,
```

so `a/n` lies on the boundary of a forbidden interval.  The previous paragraph
shows it lies in no open forbidden interval, so it is unprotected.  THM-357
then identifies these points as the boundary-only witnesses.

## Use

This theorem gives a proved model for the unit-boundary skeleton seen in
HYP-1810.  The standard tight example is not merely computationally tight: it
is exactly the equality case of Dirichlet's pigeonhole approximation argument.

The point `0` is part of the pigeonhole orbit, but it is not a safe forbidden
endpoint.  The endpoint skeleton itself is the nonzero unit set

```text
(Z/nZ)^* / n.
```

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- HYP-1810: unit-boundary skeleton.
- HYP-1813: Bohr-boundary descent.
- `04-computation/lonely_runner_bohr_descent_s362.py`.
