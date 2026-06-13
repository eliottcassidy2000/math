---
id: THM-393
name: lrc-n4-additive-return-measure-formula
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S553
depends_on:
  - HYP-2040
  - THM-392
---

# THM-393: Primitive additive-return triples have an exact n=4 safe-measure formula

## Statement

Work at LRC threshold `1/4`.  Let `a,b` be coprime positive integers and put
`c=a+b`.  Define

```text
M(a,b) = |{ t in [0,1) : ||a t|| >= 1/4,
                         ||b t|| >= 1/4,
                         ||c t|| >= 1/4 }|.
```

Then `M(a,b)` is given by the following four cases.

If `a,b` are both odd, then

```text
a+b = 0 mod 4:  M(a,b) = (ab - 1)/(16ab),
a+b = 2 mod 4:  M(a,b) = (ab + 3)/(16ab).
```

If exactly one of `a,b` is even, let `e` be the even generator, let `o` be the
odd generator, and keep `c=a+b`.  Then

```text
e = 0 mod 4:  M(a,b) = (oc + 1)/(16oc),
e = 2 mod 4:  M(a,b) = (oc - 3)/(16oc).
```

Consequently, among primitive additive-return triples:

```text
M(a,b)=0 iff {a,b}={1,2},        giving the AP triple (1,2,3),
M(a,b)>0 implies M(a,b)>=1/28,
M(a,b)=1/28 iff {a,b}={1,6},     giving the triple (1,6,7).
```

Thus the full additive-return family `(a,b,a+b)` obeys the HYP-2040 n=4
measure gap, with the same sharp positive example as the adjacent family.

## Proof

Let

```text
I = [1/4, 3/4].
```

The speed `a` condition restricts `t` to the `a` middle-half arcs

```text
t = (2h+1)/(2a) + u,       h=0,...,a-1,       |u| <= 1/(4a).
```

On such an arc,

```text
{a t} = 1/2 + a u.
```

Write

```text
beta_h = { b(2h+1)/(2a) }.
```

The remaining two speeds are `b` and `a+b`, so they contribute

```text
y = {b t} = {beta_h + b u},
{(a+b)t} = {y + 1/2 + a u}.
```

For `u >= 0`, the two conditions `y in I` and
`y+1/2+a u mod 1 in I` are equivalent to

```text
y in [3/4 - a u, 3/4].
```

For `u <= 0`, with `w=-u >= 0`, they are equivalent to

```text
{beta_h - b w} in [1/4, 1/4 + a w].
```

So every branch reduces to summing affine intervals.  For the positive side,
an integer `j` contributes the interval

```text
(j+3/4-beta_h)/(a+b) <= u <= (j+3/4-beta_h)/b,
```

intersected with `0 <= u <= 1/(4a)`.  For the negative side, an integer `j`
contributes

```text
(beta_h+j-1/4)/(a+b) <= w <= (beta_h+j-1/4)/b,
```

again intersected with `0 <= w <= 1/(4a)`.

Because `gcd(a,b)=1`, multiplication by `b` permutes the odd residues modulo
`2a`.  Thus the branch centers `beta_h` are exactly a permutation of

```text
{ r/(2a) : r odd, 1 <= r <= 2a-1 }.
```

The preceding interval sum therefore depends only on the placement of the
quarter-points `1/4` and `3/4` in the odd-residue grid.  Evaluating the two
arithmetic sums gives:

```text
if a,b odd and a+b = 0 mod 4:  total = (ab - 1)/(16ab),
if a,b odd and a+b = 2 mod 4:  total = (ab + 3)/(16ab),
if e = 0 mod 4:                total = (oc + 1)/(16oc),
if e = 2 mod 4:                total = (oc - 3)/(16oc).
```

The only distinction between the four cases is whether the moving window hits
the quarter-point grid from the inside or outside.  The `-3` cases are the
same AP-branch deficit visible in THM-392; the `+1` and `+3` cases have a
surplus instead.

This proves the four formulas.

For the gap consequence:

1. In the `e=2 mod 4` case, the only zero is `oc=3`, hence `o=1,c=3,e=2`,
   which is `{a,b}={1,2}`.  The next smallest possible product is `oc=7`,
   attained only by `{a,b}={1,6}`, and it gives

   ```text
   (7-3)/(16*7) = 1/28.
   ```

2. In the `a,b` odd and `a+b=0 mod 4` case, the smallest possible `ab` is `3`,
   giving `(3-1)/(16*3)=1/24`.

3. The other two cases are at least `1/16`.

Therefore every positive primitive additive-return safe measure is at least
`1/28`, with equality only at `(1,6,7)`.

## Verification

`04-computation/lrc_n4_additive_return_s553.py` verifies the formula against
exact rational breakpoint enumeration for every ordered coprime pair
`a,b <= 60`, with zero mismatches.  It also verifies the moving-window branch
sum above independently, again with zero mismatches.

The same script records a Tournament Analysis fingerprint on selected
additive-return generator pairs.  The pairwise observable is exact safe
measure, the switch is lower-measure-as-tighter, the tie path is lexicographic
generator order, and the resulting tournament is transitive with one
Hamiltonian path.

## Context

THM-392 solved the adjacent family `(1,q,q+1)`.  THM-393 solves the next
observed obstruction from the S552 scan: primitive additive-return triples
`(a,b,a+b)`, including

```text
(2,3,5), (2,5,7), (3,10,13), (3,5,8).
```

This still does not prove the full HYP-2040 n=4 measure gap.  It proves that
the first non-adjacent structured obstruction family also respects the same
`1/28` lower bound.
