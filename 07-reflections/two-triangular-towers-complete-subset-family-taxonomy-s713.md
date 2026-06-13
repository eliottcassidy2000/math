---
source: monad-explorer-2026-06-12-S713
relates: T765, S712, A059270, A059255, A222716
status: exact taxonomy
---

# The complete family list: every nice subset pattern is one of 8 Pell families

## 1. Why the additive tower covers everything and the square shell does not

The additive tower is exhaustive because stage `m` is exactly

`[m^2, ..., (m+1)^2-1]`.

That is the whole square shell between consecutive squares, split into

- `A_L(m) = [m^2, ..., m^2+m]`,
- `A_R(m) = [m^2+m+1, ..., m^2+2m]`.

So the additive tower partitions all positive integers.

The raw square tower does not.  Its stage-`n` carrier is

`[2n^2+n, ..., 2n^2+3n]`,

which has length `2n+1`, and the gap before the next such carrier has length
`2n+2`.

So the raw square shell is a sparse subsector of the additive partition, not a
second partition of the naturals.

## 2. The right question

The user’s refined prompt asks for all patterns like

- `10+11+12` inside `9+10+11+12`,
- `13+14` inside `13+14+15`,
- `21+22+23+24`,
- `25+26+27` inside `25+...+30`.

The clean answer is:

every such "nice" family comes from aligning **one endpoint** of a square-shell
block with **one endpoint** of an additive half-block.

There are 4 endpoints on the square side and 4 on the additive-half side, so
there are exactly **8** families.

There is one small startup anomaly: the equation for "`B_L` prefix of `A_R`"
has the first solution `(n,m)=(1,1)`, where `B_L(1)=[3,4]` begins at `A_R(1)=[3]`
but spills one step into `A_L(2)=[4]`.  After that lone straddle, the family is
genuinely a prefix family, starting at `(3,4)`.

## 3. The 8 families

Write

- `B_L(n) = [2n^2+n, ..., 2n^2+2n]`, size `n+1`,
- `B_R(n) = [2n^2+2n+1, ..., 2n^2+3n]`, size `n`.

Then the complete family list is:

1. `B_L` is a prefix of `A_L(m)`:
   `2n^2+n = m^2`.
2. `B_L` is a suffix of `A_L(m)`:
   `2n^2+2n = m^2+m`.
3. `B_L` is a prefix of `A_R(m)`:
   `2n^2+n = m^2+m+1`.
4. `B_L` is a suffix of `A_R(m)`:
   `2n^2+2n = m^2+2m`.
5. `B_R` is a prefix of `A_L(m)`:
   `2n^2+2n+1 = m^2`.
6. `B_R` is a suffix of `A_L(m)`:
   `2n^2+3n = m^2+m`.
7. `B_R` is a prefix of `A_R(m)`:
   `2n^2+2n+1 = m^2+m+1`.
8. `B_R` is a suffix of `A_R(m)`:
   `2n^2+3n = m^2+2m`.

Each one reduces to a Pell equation, so each one is an infinite sparse family.

## 4. Where the user examples land

### `10+11+12 = 13+14` as raw pattern

This is the first hit of the Pell family

`2n^2+2n = m^2+m`,

with `(n,m)=(2,3)`.

So:

- `B_L(2)=[10,11,12]` is the **suffix** of `A_L(3)=[9,10,11,12]`;
- `B_R(2)=[13,14]` is the **prefix** of `A_R(3)=[13,14,15]`.

The host sizes are both `4` and `3`, the guest sizes are `3` and `2`, so both
deficits are exactly `1`.

### `21+22+23+24`

This is the unique exact block:

`B_L(3)=A_R(4)=[21,22,23,24]`.

It is simultaneously a prefix and a suffix of `A_R(4)`, because it fills the
whole host.  This is the only exact shared raw block between the towers.

### `25+26+27`

This is the first hit of

`2n^2+2n+1 = m^2`,

with `(n,m)=(3,5)`.

So `B_R(3)=[25,26,27]` is the **prefix** of `A_L(5)=[25,26,27,28,29,30]`.
The host has size `6`, the guest has size `3`, so the tail deficit is `3`.

### `36+37+38+39+40`

This is the first hit of

`2n^2+n = m^2`,

with `(n,m)=(4,6)`.

So `B_L(4)=[36,37,38,39,40]` is the **prefix** of `A_L(6)=[36,37,38,39,40,41,42]`.
The host size is `7`, the guest size is `5`, so the tail deficit is `2`.

## 5. The concise prediction rule

The rule is very short:

> a nice family appears exactly when a square-shell endpoint lands on an additive
> half-block endpoint.

That means one of

- `2n^2+n`,
- `2n^2+2n`,
- `2n^2+2n+1`,
- `2n^2+3n`

must equal one of

- `m^2`,
- `m^2+m`,
- `m^2+m+1`,
- `m^2+2m`.

Each of the 16 formal choices either never happens or reduces to one of the 8
families above.

## 6. Size bookkeeping

The size formulas are also clean:

- `|A_L(m)| = m+1`,
- `|A_R(m)| = m`,
- `|B_L(n)| = n+1`,
- `|B_R(n)| = n`.

So once you know which family you are in, the subset deficit is immediate:

- in `A_L`, the deficit is `(m+1)-(n+1)=m-n` for `B_L`, and `(m+1)-n=m+1-n`
  for `B_R`;
- in `A_R`, the deficit is `m-(n+1)=m-n-1` for `B_L`, and `m-n` for `B_R`.

This is the exact amount chopped off from the front or back.

## 7. What remains open

This finishes the endpoint-aligned subset problem, but not the whole overlap
problem.

What remains:

1. classify all overlaps that do **not** share an endpoint;
2. decide whether the raw-sum crossover at `90` is globally unique;
3. push the triangular analog `A222716` through the same endpoint framework.
