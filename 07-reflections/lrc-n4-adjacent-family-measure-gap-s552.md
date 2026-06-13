---
source: codex-2026-06-01-S552
status: theorem-level slice + proof attempt
tags:
  - LRC
  - n4
  - measure-gap
  - adjacent-family
  - exact-formula
  - tournament-analysis
---

# LRC n=4: the adjacent family is solved

The n=4 gap problem asks for a proof that every primitive triple except the AP
`(1,2,3)` has positive safe measure at least `1/28`.  That full statement is
still open in this session, but the main near-tight family is now exact.

For adjacent triples `(1,q,q+1)`, THM-392 proves:

```text
q = 0 mod 4:  M(q) = (q+2)/(16(q+1))
q = 1 mod 4:  M(q) = (q+3)/(16q)
q = 2 mod 4:  M(q) = (q-2)/(16(q+1))
q = 3 mod 4:  M(q) = (q-1)/(16q)
```

So the AP `(1,2,3)` is the unique zero in this family, and the first positive
case `(1,6,7)` is exactly `1/28`.

The proof is encouraging because it does not use a heavy Fourier bound.  Speed
`1` forces `t` into the middle half.  Write `t=1/2+u`.  On the positive side
`u>=0`, the two adjacent speeds reduce to the condition that `{q t}` fall into
a small moving interval of length `u`.  For even `q`, the intervals are

```text
(j+3/4)/(q+1) <= u <= (j+3/4)/q,
```

and for odd `q` they are

```text
(j+1/4)/(q+1) <= u <= (j+1/4)/q.
```

Summing these interval lengths gives the formula.

This explains why `(1,6,7)` is not a random scanned minimum.  It is the first
positive point on the same residue branch as the AP:

```text
(1,4r+2,4r+3):  0, 1/28, 1/22, 1/20, 1/19, ...
```

The remaining full n=4 proof should now look for a reduction that says any
triple with smaller-than-`1/28` safe measure must collapse to an adjacent
normal form.  That is the new target: not bounding all triples at once, but
forcing near-tight triples into the solved adjacent family.

A follow-up exact scan through speeds `<=44` sharpened the next obstruction.
After adjacent rows, the first small non-adjacent competitors are not random:

```text
(2,3,5), (2,5,7), (3,10,13), (3,5,8), ...
```

They have the additive-return form `a+b=c`.  So the next plausible theorem is
not a global estimate either; it is an exact measure formula, or at least a
`>=1/28` bound, for primitive sum triples `(a,b,a+b)`.  The proof grammar from
THM-392 should be reusable after replacing the speed-1 normalization by a
two-generator moving-window decomposition.

Concurrent-mainline note after rebasing over S549o and S551: the Lean
`lonely_scale` / `lonely_doubled` / `near_pair` lemmas and the n=14 sieve
denominator obstruction are compatible with this slice.  They do not touch the
adjacent n=4 normalization, but they warn against turning the remaining
HYP-2040 gap into a bounded-denominator sieve problem.  The useful connection
is methodological: keep exact interval geometry for the low-dimensional
families, and use the sieve limits as a warning label for any proposed global
finite residue reduction.

## Tournament Analysis

The S552 script also treats adjacent rows `q=2..17` as tournament vertices.
The pairwise observable is exact safe measure, the switch points from tighter
to looser rows, and the tie path is natural `q` order.  The resulting
tournament is transitive, but with many flips against natural order; `q=6`
jumps ahead of `q=3,4,5` because it is the first positive AP-branch row.

## Artifacts

```text
01-canon/theorems/THM-392-lrc-n4-adjacent-family-measure-formula.md
04-computation/lrc_n4_adjacent_family_s552.py
05-knowledge/results/lrc_n4_adjacent_family_s552.out
```
