---
source: monad-explorer-2026-06-12-S712
relates: A059270, A059255, A222716, S586, S548
status: bounded exact computation + exact derivations
---

# Two triangular towers: the unique shared block is 21..24, and the rest is Pell

## The setup

Your first tower is the additive one:

`A_L(n) = [n^2, ..., n^2+n]`, `A_R(n) = [n^2+n+1, ..., n^2+2n]`,

so

`sum A_L(n) = sum A_R(n) = A059270(n) = n(n+1)(2n+1)/2`.

Your second tower is the square one:

`B_L(n) = [2n^2+n, ..., 2n^2+2n]`, `B_R(n) = [2n^2+2n+1, ..., 2n^2+3n]`,

so

`sum k^2 over B_L(n) = sum k^2 over B_R(n) = A059255(n)`.

The useful move was to stop staring at the sums and pass to the **endpoints**.
That is where the real structure lives.

## The first bridge: the additive tower already hides squares

One exact identity is much better than it looks:

`A059270(n) = 3 * sum_{j=1..n} j^2`.

So the "pure addition" tower is not separate from the square tower at all.  It is
already three copies of the square-pyramid numbers.  In that sense the first tower
is the additive shadow of the second.

There is then a clean multiplicative lift:

`A059255(n) = A059270(n) * ((12n^2+12n+1)/3)`.

So the square tower is the additive tower times a quadratic weight
`(2n+1)^2 - 2/3`.

This fits the repo's existing add/mult bridge language well: one tower is not
replacing the other; it is a weighted lift of it.

## The obvious examples split into two different phenomena

The user examples are not all of the same type.

### 1. The exact shared block

`21+22+23+24`

is genuinely special.  It is not merely an overlap.  It is the **unique exact raw
block shared by both towers**:

`B_L(3) = A_R(4) = [21,22,23,24]`.

This is rigid.  Equality of lengths forces `m=n+1`, and then the endpoint equation
becomes `n^2-2n-3=0`, so `n=3`.

So the `21..24` bridge is not one instance of a large exact-block family.  It is
the only one.

### 2. The truncation bridges

`10+11+12 = 13+14`

is different.  Here the square-tower raw blocks are not exact A-blocks.  They are
tails/prefixes of the additive blocks:

- `B_L(2) = [10,11,12]` is the right tail of `A_L(3) = [9,10,11,12]`;
- `B_R(2) = [13,14]` is the left prefix of `A_R(3) = [13,14,15]`.

And this is **not** isolated.  It extends to an infinite Pell family:

`m(m+1) = 2n(n+1)`,

equivalently

`(2m+1)^2 - 2(2n+1)^2 = -1`.

The first hits are

`(n,m) = (2,3), (14,20), (84,119), ...`

So `10+11+12` is the first member of an infinite "same end, shorter square-shell
tail" family.

## The 25+26+27 bridge is also infinite

The other user prompt was to think about

`25+26+27`

relative to `16+17+18+19+20` and the square tower.  The endpoint view explains the
right comparison:

`B_R(3) = [25,26,27]`

starts exactly where the next additive left block starts:

`A_L(5) = [25,26,27,28,29,30]`.

Again this is not isolated.  It is another Pell family:

`m^2 = 2n^2 + 2n + 1`,

equivalently

`2m^2 - (2n+1)^2 = 1`.

The first hits are

`(n,m) = (3,5), (20,29), (119,169), ...`

So `25..27` is the first member of an infinite "square-tower right block begins an
additive left block" family.

## What happens at 36?

The `n=4` square-left block

`B_L(4) = [36,37,38,39,40]`

starts exactly where

`A_L(6) = [36,37,38,39,40,41,42]`

starts.  This is yet another infinite family:

`m^2 = n(2n+1)`,

equivalently

`(4n+1)^2 - 8m^2 = 1`.

The first hits are

`(n,m) = (4,6), (144,204), ...`

So `36` is the first start-alignment on the left-left side.

## The triangular-number analog is the third tower

The user asked to "find others like this."  The cleanest third one is already in
OEIS and fits the repo better than higher power sums do.

The triangular-number analog is `A222716`:

`sum_{j=n^2-n-1..n^2-1} T_j = sum_{j=n^2..n^2+n-2} T_j`.

So there is indeed a **third tower**, but it is not the next power tower.  It is
the next **figurate** tower.

This matters because the OEIS crossrefs also record the stopping point:

- powers `1` and `2` admit the consecutive-run equalities (`A059270`, `A059255`);
- for powers `>2`, the analog stops (`A234319`);
- but figurate numbers continue via `A222716`.

That is a real structural split:

- the **power ladder** stops at squares;
- the **figurate ladder** keeps going.

This is exactly the kind of add-vs-mult distinction the repo already likes.

## The main lesson

The strongest new point is not a new formula.  It is the shift in viewpoint:

- sums alone hide the structure;
- interval endpoints reveal exact unique bridges plus infinite Pell truncation
  families.

So the right objects are:

1. the shell endpoints `n^2`, `n(n+1)`, `n(2n+1)`, `n(2n+3)`;
2. the exact shared block `21..24`;
3. the Pell-coded one-endpoint crossover families.

That feels like the durable frontier here.

## Open directions

1. Classify **all** nontrivial overlaps between `A`- and `B`-blocks, not just
   endpoint coincidences.
2. Determine whether the raw-sum crossover
   `21+22+23+24 = 16+17+18+19+20 = 90`
   is the only one globally, not just in the current bounded scan.
3. Push the third tower `A222716` into the same endpoint-shell language and ask
   whether its own crossover families with `A` or `B` are again Pell or a higher
   Diophantine type.
