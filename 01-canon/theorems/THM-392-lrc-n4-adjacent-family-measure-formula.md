---
id: THM-392
name: lrc-n4-adjacent-family-measure-formula
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S552
depends_on:
  - HYP-2040
  - THM-391
---

# THM-392: The adjacent n=4 LRC family has an exact safe-measure formula

## Statement

Work at LRC threshold `1/4`.  For an integer `q >= 2`, define

```text
M(q) = |{ t in [0,1) : ||t|| >= 1/4,
                         ||q t|| >= 1/4,
                         ||(q+1)t|| >= 1/4 }|.
```

Equivalently, `M(q)` is the safe measure of the adjacent triple
`(1,q,q+1)`.  Then

```text
q = 0 mod 4:  M(q) = (q+2)/(16(q+1)),
q = 1 mod 4:  M(q) = (q+3)/(16q),
q = 2 mod 4:  M(q) = (q-2)/(16(q+1)),
q = 3 mod 4:  M(q) = (q-1)/(16q).
```

In particular,

```text
M(2)=0,                        for the AP triple (1,2,3),
M(q) >= 1/28 for every q >= 3,
M(q)=1/28 iff q=6.             for the triple (1,6,7)
```

Thus the sharp positive example in the adjacent family is exactly `(1,6,7)`.

## Proof

Let

```text
I = [1/4, 3/4].
```

The speed `1` condition says `t in I`.  Write

```text
t = 1/2 + u,       -1/4 <= u <= 1/4,
x = {q t}.
```

The remaining two conditions are

```text
x in I,
x + t mod 1 in I.
```

For `u >= 0`, these are equivalent to

```text
x in [3/4 - u, 3/4],
```

an interval of length `u`.  For `u <= 0`, the corresponding interval is
`[1/4, 1/4-u]`.  The safe set is invariant under `t -> 1-t`, so the total
measure is twice the measure on `0 <= u <= 1/4`.

### Even `q`

If `q` is even, then

```text
{q(1/2+u)} = {q u}.
```

On the positive side, for an integer `j`, the condition
`{qu} in [3/4-u, 3/4]` is

```text
j + 3/4 - u <= q u <= j + 3/4,
```

or

```text
(j+3/4)/(q+1) <= u <= (j+3/4)/q.
```

This interval has length `(j+3/4)/(q(q+1))`.  It contributes positive measure
inside `[0,1/4]` exactly for

```text
j = 0, 1, ..., floor(q/4)-1.
```

So if `q=4r` or `q=4r+2`, the positive-side sum has `r` intervals and

```text
M(q) = 2/(q(q+1)) * sum_{j=0}^{r-1} (j+3/4)
     = 2/(q(q+1)) * r(2r+1)/4.
```

For `q=4r`, this becomes

```text
M(q) = (2r+1)/(8(4r+1)) = (q+2)/(16(q+1)).
```

For `q=4r+2`, this becomes

```text
M(q) = r/(4(4r+3)) = (q-2)/(16(q+1)).
```

### Odd `q`

If `q` is odd, then

```text
{q(1/2+u)} = {1/2 + q u}.
```

For an integer `j`, the positive-side condition is

```text
j + 3/4 - u <= 1/2 + q u <= j + 3/4,
```

or

```text
(j+1/4)/(q+1) <= u <= (j+1/4)/q.
```

The interval length is `(j+1/4)/(q(q+1))`.  It contributes positive measure
for `j=0,1,...,r` when `q=4r+1` or `q=4r+3`.  Hence

```text
M(q) = 2/(q(q+1)) * sum_{j=0}^{r} (j+1/4)
     = 2/(q(q+1)) * (r+1)(2r+1)/4.
```

For `q=4r+1`, this gives

```text
M(q) = (r+1)/(4(4r+1)) = (q+3)/(16q).
```

For `q=4r+3`, this gives

```text
M(q) = (2r+1)/(8(4r+3)) = (q-1)/(16q).
```

This proves the four formulas.

The lower bound follows by direct comparison.  In the `q=4r+2` branch,
`M(2)=0`, and for `r>=1`,

```text
M(4r+2) = r/(4(4r+3))
```

is increasing, so its first positive value is `M(6)=1/28`.  In the other
three residue classes, the displayed formulas are all strictly greater than
`1/28` for the allowed `q>=3`.  Therefore `M(q)>=1/28` for all `q>=3`, with
equality exactly at `q=6`.

## Verification

`04-computation/lrc_n4_adjacent_family_s552.py` verifies the formula against
exact rational breakpoint enumeration for `q=2..41`, with zero mismatches.  It
also records a Tournament Analysis fingerprint on vertices `q=2..17`: exact
safe measure is the pairwise observable, lower measure is the tighter-row
switch, and the resulting tournament is transitive with one Hamiltonian path.

## Context

This does not prove the full HYP-2040 n=4 measure gap.  It proves the main
observed near-tight adjacent family:

```text
(1,2,3) -> 0,
(1,6,7) -> 1/28,
(1,10,11) -> 1/22,
(1,14,15) -> 1/20,
...
```

The full remaining gap problem is to rule out a smaller positive measure among
non-adjacent triples.
