---
id: THM-397
name: lrc-n14-collective-pinch-endpoint-blocker
status: PROVED
date: 2026-06-02
depends_on:
  - HYP-2059
  - HYP-2061
  - THM-396
---

# THM-397: A collective non-shield cover of a small n=14 pinch needs an endpoint blocker

## Statement

Let `a,b` be positive speeds, set

```text
D = a+b,    g = gcd(a,b),    s = D/g.
```

Assume `2 <= s <= 14`.  Consider the pair-safe pinch times

```text
t = m/D,    1 <= m < D,    m not == 0 (mod s),
```

which are exactly the `n=14`-safe pinch times for the pair `(a,b)`.

Let `C` be any set of third speeds.  If every pair-safe pinch time is killed by
some speed in `C`,

```text
for every m not == 0 (mod s), there is c in C with ||c m / D|| < 1/14,
```

then at least one speed in `C` is an endpoint blocker modulo `D`:

```text
||c/D|| < 1/14.
```

Equivalently,

```text
c mod D in {1,...,floor((D-1)/14)}
          union {D-floor((D-1)/14),...,D-1}.
```

In particular, if `D <= 14`, no non-shield collective cover is possible.  A
failed small pair with actual sum `D <= 14` must therefore have a true
sum-multiple shield `D | c`.

## Proof

Since `s >= 2`, the residue `m=1` is pair-safe:

```text
1 not == 0 (mod s).
```

By the collective-cover assumption, some `c in C` kills the pinch `t=1/D`.
That is exactly

```text
||c/D|| < 1/14.
```

This is the endpoint-blocker condition.  Writing `r` for the least positive
residue of `c` modulo `D`, the inequality says

```text
min(r, D-r) / D < 1/14,
```

or equivalently

```text
min(r, D-r) <= floor((D-1)/14).
```

If `D <= 14`, then `floor((D-1)/14)=0`, so no nonzero residue modulo `D` can be
an endpoint blocker.  Thus any cover of all pair-safe pinches must use a speed
with `c == 0 (mod D)`, i.e. a true sum-multiple shield.

## Consequence for the HYP-2061 residual

THM-396 proves that a single universal blocker must be a sum-multiple shield.
This theorem gives the first necessary condition for the remaining collective
branch: every non-shield collective cover must expose a near-zero/near-one
residue modulo the pair sum `D`.

Therefore a HYP-2060-style n=14 counterexample cannot merely cover every small
pair-cell abstractly.  For each failed reduced-sum-`<=14` pair `(a,b)`, it must
either contain a sum-multiple shield for `a+b` or contain a speed in the endpoint
window modulo `a+b`.  This turns the collective-cover residual into an
endpoint-ledger constraint compatible with the existing endpoint-core and
pressure-owner graph machinery.
