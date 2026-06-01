---
id: THM-394
name: lrc-n4-13q-family-measure-formula
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S553
depends_on:
  - HYP-2040
  - THM-392
  - THM-393
---

# THM-394: The n=4 family `(1,3,q)` has an exact safe-measure formula

## Statement

Work at LRC threshold `1/4`.  For an integer `q >= 1`, define

```text
M(q) = |{ t in [0,1) : ||t|| >= 1/4,
                         ||3t|| >= 1/4,
                         ||q t|| >= 1/4 }|.
```

Write `q=12m+r`, with `0 <= r < 12`.  Then

```text
M(q) = N_r(q)/(12q),
```

where

```text
r =  0: N_r(q)=q,      r =  1: N_r(q)=q+1,
r =  2: N_r(q)=q-2,    r =  3: N_r(q)=q+3,
r =  4: N_r(q)=q-2,    r =  5: N_r(q)=q+1,
r =  6: N_r(q)=q,      r =  7: N_r(q)=q-1,
r =  8: N_r(q)=q+2,    r =  9: N_r(q)=q-3,
r = 10: N_r(q)=q+2,    r = 11: N_r(q)=q-1.
```

Consequently, for distinct triples in this family:

```text
q=2 gives the AP triple (1,2,3) and M(2)=0,
q=4 gives the adjacent/additive triple (1,3,4) and M(4)=1/24,
q not in {1,2,3,4} implies M(q)>=1/18,
M(q)=1/18 first occurs at q=9, for the triple (1,3,9).
```

Thus the next small family after the adjacent and additive-return rows also
stays well above the HYP-2040 gap threshold `1/28`.

## Proof

Let

```text
I = [1/4, 3/4].
```

The first two speed conditions are especially rigid:

```text
t in I,
3t mod 1 in I.
```

Inside `[0,1)`, their intersection is, up to endpoints,

```text
J = [5/12, 7/12].
```

Therefore

```text
M(q) = |{ t in J : ||q t|| >= 1/4 }|.
```

Write

```text
t = 1/2 + u,      |u| <= 1/12.
```

Then

```text
q t = q/2 + q u.
```

If `q` is even, then `q/2` is an integer and the condition is

```text
{q u} in [1/4, 3/4].
```

If `q` is odd, then `q/2` is a half-integer and the condition is

```text
{q u} in [-1/4, 1/4] mod 1.
```

Now put `x=q u`.  Since `du=dx/q`, the measure is `1/q` times the length of
one of two periodic half-density sets inside the symmetric interval

```text
[-q/12, q/12].
```

Let `q=12m+r`.  The `2m` full unit periods contribute length `m`, because each
period has density `1/2`.  It remains only to measure the residual interval
`[-r/12,r/12]`.

For even `q`, the residual lengths for
`{x} in [1/4,3/4]` are

```text
r = 0: 0,      r = 2: 0,      r = 4: 1/6,
r = 6: 1/2,    r = 8: 5/6,    r = 10: 1.
```

For odd `q`, the residual lengths for
`{x} in [-1/4,1/4] mod 1` are

```text
r = 1: 1/6,    r = 3: 1/2,    r = 5: 1/2,
r = 7: 1/2,    r = 9: 1/2,    r = 11: 5/6.
```

Dividing `m + residual_length(r)` by `q` and simplifying gives exactly the
displayed table for `N_r(q)/(12q)`.

For the consequences, direct residue-class comparison gives:

```text
M(2)=0,
M(4)=1/24,
M(q)>=1/18 for q not in {1,2,3,4},
M(9)=1/18.
```

On each residue branch the formula is constant or monotone toward `1/12`, so
the finite check of residues `r=0,...,11` at their first admissible `q` proves
the stated minimum.

## Verification

`04-computation/lrc_n4_13q_family_s553.py` verifies the formula against exact
rational breakpoint enumeration for `q=1..80`, with zero mismatches.  It also
verifies the interval-count derivation above for the same range, again with
zero mismatches.

The script records a Tournament Analysis fingerprint on selected `q` rows.
The pairwise observable is exact safe measure, the switch is
lower-measure-as-tighter, and the selected-row tournament is transitive with
one Hamiltonian path.

## Context

After excluding the adjacent rows solved by THM-392 and the additive-return
rows solved by THM-393, a fresh scan through speeds `<=60` found `(1,3,9)` as
the next small structured row, with safe measure `1/18`.  THM-394 explains
that row as the first nontrivial minimum of the whole `(1,3,q)` family.

This still does not prove the full HYP-2040 n=4 measure gap.  It removes
another visible low-measure corridor.
