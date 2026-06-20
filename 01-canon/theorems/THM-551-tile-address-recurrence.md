---
id: THM-551
title: The full and half tiling recurrences combine into a two-clock address for every tile
status: PROVED; verified by `tile_address_recurrence_codex_s57.py`
source: codex-2026-06-20-S57
depends_on:
  - THM-280
  - THM-442
  - THM-513
  - THM-549
  - THM-550
related:
  - HYP-2685
  - HYP-2686
  - HYP-2687
  - HYP-2688
  - HYP-2689
---

# THM-551 - Tile Address Recurrence

## Statement

In the fixed-Hamiltonian-path tiling model, a free tile is a pair

```text
(a,b),  1 <= b < a <= n,  a-b >= 2.
```

Define its two recursive clocks by

```text
beta(a,b) = a,
tau(a,b)  = a+b-1.
```

Then:

1. The map `(a,b) -> (beta,tau)` is a bijection from free tiles to integer
   pairs

   ```text
   beta >= 3,  beta <= tau <= 2beta-3.
   ```

   The inverse is

   ```text
   (beta,tau) -> (beta, tau-beta+1).
   ```

2. The full staircase count is the cumulative birth-strip count

   ```text
   F_n = binom(n-1,2) = sum_{beta<=n} (beta-2).
   ```

   Thus the full `A+B+C-D-E-F+G` recurrence of THM-442 is the cell-affine
   recursion seen through the `beta` clock.

3. With the KPS/mac-mini half-domain convention

   ```text
   a+b <= n+1
   ```

   a tile is in the half carrier at size `n` exactly when `tau<=n`.  Therefore

   ```text
   h_n = floor((n-1)^2/4) = sum_{tau<=n} floor((tau-1)/2).
   ```

   The crossing layer `tau=n` is precisely the reflection fixed line
   `a+b=n+1`, and has `floor((n-1)/2)` cells.  Thus the even and odd
   half-tiling recurrences of THM-550 are the same local address system seen
   through the parity of the `tau` clock.

4. The mirror/complement reflection at size `n` sends address

   ```text
   (beta,tau) -> (n+beta-tau, 2n-tau).
   ```

   Hence a tile is fixed exactly when `tau=n`.  If `tau<n`, the tile is already
   the canonical half-domain representative; if `tau>n`, its canonical
   representative is its reflected address above.

5. The one-flip interval-root facts of THM-513 are computable from the address
   alone:

   ```text
   b       = tau-beta+1,
   gap     = a-b-1 = 2beta-tau-2,
   c3      = gap,
   H_1flip = 1 + 2^gap,
   score defect = e_b - e_beta.
   ```

6. For `n>=4`, the local membership word in THM-442's full three-subtriangle
   inclusion-exclusion cover is also address-local.  In pin coordinates
   `r=a-b-1`, `c=b`, the three size-`n-1` subtriangles are

   ```text
   A: r+c <= n-2,
   B: r >= 2,
   C: c >= 2.
   ```

   In address coordinates these are

   ```text
   A: beta < n,
   B: 2beta-tau-2 >= 2,
   C: tau-beta+1 >= 2.
   ```

   If a tile belongs to `k` of these subtriangles, its local IE weight is

   ```text
   k - binom(k,2) + binom(k,3) = 1.
   ```

Consequently the three recurrences are not merely count recurrences.  Together
they form a complete local coordinate system for a tile: `beta` records when
the tile is born in the full tiling model, while `tau` records when its mirror
orbit crosses the half-tiling fixed line.

## Proof

For a tile `(a,b)`, `beta=a` and `tau=a+b-1`.  Since `1<=b<=a-2`,

```text
a <= a+b-1 <= 2a-3,
```

so the claimed inequalities hold.  Conversely, if
`beta <= tau <= 2beta-3`, then `b=tau-beta+1` satisfies
`1<=b<=beta-2`, so `(beta,b)` is a free tile.  This proves the bijection.

The tile first exists in the staircase for tournament size `n=a=beta`.  The
number of tiles born at `beta` is `beta-2`, one for each
`b=1,...,beta-2`.  Summing birth layers gives

```text
sum_{beta=3}^n (beta-2) = binom(n-1,2).
```

This is the layer form of the full cell-affine recursion from THM-442.

The half-domain condition `a+b<=n+1` is equivalent to

```text
tau=a+b-1 <= n.
```

For fixed `tau=t`, the possible tiles are

```text
(a,b) = (a, t-a+1),
```

with `1<=t-a+1<=a-2`.  Thus

```text
(t+3)/2 <= a <= t,
```

and the number of integer choices is `floor((t-1)/2)`.  Therefore

```text
h_n = sum_{t=2}^n floor((t-1)/2) = floor((n-1)^2/4),
```

which is THM-550 in layer form.

By THM-280/THM-549, reflection sends

```text
(a,b) -> (n+1-b, n+1-a).
```

Writing `a=beta` and `b=tau-beta+1`, the reflected tile has

```text
beta' = n+1-b = n+beta-tau,
tau'  = beta' + (n+1-beta) - 1 = 2n-tau.
```

So the address reflection formula follows.  It is fixed exactly when
`tau=n`, i.e. exactly on the half-tiling fixed line.

Finally, the inverse address formula gives `b=tau-beta+1`, so

```text
a-b-1 = beta-(tau-beta+1)-1 = 2beta-tau-2.
```

Substituting this gap into THM-513 gives the one-flip `c3`, `H`, and score
defect formulas.

For the local full inclusion-exclusion word, use `r=a-b-1` and `c=b`.  Since
`r+c=a-1=beta-1`, the first shifted-subtriangle condition
`r+c<=n-2` is equivalent to `beta<n`.  The second condition is exactly
`gap>=2`, and the third is exactly `b>=2`, giving the address tests above.
For `n>=4`, every present free tile lies in at least one of the three shifted
subtriangles.  A tile in exactly `k=1,2,3` of them contributes

```text
k - binom(k,2) + binom(k,3) = 1
```

to the local inclusion-exclusion sum, proving the local form of the full
recurrence.

## Computational Check

`04-computation/tile_address_recurrence_codex_s57.py` verifies for `n<=16`:

- address bijection and inverse;
- mirror formula on addresses;
- canonical half representatives;
- full local inclusion-exclusion weight `1` for every tile for `n>=4`, where
  the shifted-triangle cover is nondegenerate;
- birth-layer and crossing-layer count identities;
- one-flip facts recovered from `(beta,tau)`.

Stored output:

```text
05-knowledge/results/tile_address_recurrence_codex_s57.out
```

## Scope

This theorem computes local cell/root facts.  It does not by itself make
Hamiltonian-path count, OCF, or LRC sector mass cell-affine.  Its use is
address-first computation: keep the two clocks, then attach cycle-space packet
data only after the mirror quotient has been chosen.
