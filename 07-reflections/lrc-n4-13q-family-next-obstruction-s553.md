---
source: codex-2026-06-01-S553
status: theorem-level slice
tags:
  - LRC
  - n4
  - measure-gap
  - next-obstruction
  - exact-formula
  - tournament-analysis
---

# LRC n=4: the `(1,3,q)` corridor is solved

After THM-392 solved adjacent rows and THM-393 solved additive-return rows, a
fresh exact scan through speeds `<=60` showed the next small structured row:

```text
(1,3,9) with M=1/18.
```

THM-394 solves the full family `(1,3,q)`.  If `q=12m+r`, then

```text
M(q) = N_r(q)/(12q),
```

where

```text
r =  0: N=q      r =  1: N=q+1    r =  2: N=q-2    r =  3: N=q+3
r =  4: N=q-2    r =  5: N=q+1    r =  6: N=q      r =  7: N=q-1
r =  8: N=q+2    r =  9: N=q-3    r = 10: N=q+2    r = 11: N=q-1
```

The structure is much simpler than a generic triple.  Speeds `1` and `3`
already force

```text
t in [5/12, 7/12].
```

Writing `t=1/2+u`, `|u|<=1/12`, the third runner only asks for a periodic
half-density set inside `[-q/12,q/12]`.  Full periods contribute exactly half
their length; the entire formula is the residual length for `q mod 12`.

Consequences:

```text
q=2 -> AP (1,2,3), M=0
q=4 -> adjacent/additive row (1,3,4), M=1/24
q not in {1,2,3,4} -> M(q)>=1/18
q=9 -> first next-corridor minimum, M=1/18
```

This is not as sharp as the `1/28` global candidate; rather, it explains why
the first post-additive scan result is already far above the gap threshold.

## Next After This

A smaller exact scan through speeds `<=40`, excluding adjacent rows,
additive-return rows, and the solved `(1,3,q)` corridor, pushed the next
visible obstruction into ratio-3 territory:

```text
(2,5,15) -> 1/15
(1,7,21), (3,7,21) -> 1/14
(2,9,27) -> 2/27
(4,5,15) -> 3/40
```

The repeated signature is that one speed is three times another: `c=3b`,
`c=3a`, or `b=3a`.  This suggests the next proof target should be a
multiplicative corridor such as `(a,b,3b)` or `(a,3a,c)`.  It is a different
shape from the additive-return theorem: speeds `b` and `3b` force the same
short interval as `(1,3,q)` after the measure-preserving map `t -> bt`, but
the remaining speed `a` sees the `b` preimages.  That preimage bookkeeping is
the new difficulty.

## Tournament Analysis

The S553 `(1,3,q)` script treats selected `q` rows as vertices.  The observable
is exact safe measure, lower measure wins the switch, and selected `q` order is
the tie path.  The selected tournament is transitive, with score histogram
`{0:1,...,11:1}`, no directed 3-cycles, singleton SCCs, and one Hamiltonian
path.

## Artifacts

```text
01-canon/theorems/THM-394-lrc-n4-13q-family-measure-formula.md
04-computation/lrc_n4_13q_family_s553.py
05-knowledge/results/lrc_n4_13q_family_s553.out
```
