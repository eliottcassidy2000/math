# Triangular Towers, Moment Lifts, And Fractional Addresses

The useful reframing is that the two towers are not two unrelated sequences.
They are the first two integral members of a moment-balancing problem.

Tower A splits the square shell

```text
[n^2, (n+1)^2-1]
```

into `n+1` terms on the left and `n` terms on the right.  The first moments
match.  Tower B splits the triangular shell

```text
[T_{2n}, T_{2n+1}-1]
```

with the same unequal counts, but now the second moments match.  The
unsquared B row has defect `n(n+1)`, so squaring is not decorative; it is the
operation that repairs the additive imbalance.

That makes the examples feel less like numerology.  The block
`10+11+12 = 13+14` is the `n=2` second-moment row sitting inside the `n=3`
first-moment row after adding `9` and `15`.  The block `21..24` is stronger:
it is exactly `B_3.L` and exactly `A_4.R`.  The side-equality equations show
this whole-side hinge is unique.  It is the single place where the
second-moment triangular tower and the first-moment square tower share an
entire side.

The higher-moment extension is the important surprise.  For power `p`, the
balancing start is exact and integral for `p=1` and `p=2`:

```text
A_1(n)=n^2
A_2(n)=2n^2+n.
```

For `p>=3`, the start is almost but not quite

```text
p n^2 + (p-1)n.
```

The missing leading address is

```text
(p-1)(p-2)/(12p).
```

This is exactly the repo pattern: a scalar/product collision looks real, but
after the first two clean shadows one must attach the missing fraction/address
coordinate.  Without that coordinate the equality is only a boundary illusion.

That is why this connects to HYP-2452.  A coefficient row is a boundary total;
a factorization is a hidden convolution grid.  Here a shell interval is a
boundary support; a moment equality is a boundary total; the fractional
Faulhaber/Bernoulli address is the hidden coordinate that says whether the
next lift is honest.  The same grammar also matches the LRC14 ledgers: `q`
blocked is a scalar shadow, but the proof must remember which runner owners,
unit twists, divisor fibers, and carries realize the block.

The `n=3` B shell `[21,27]` is a good symbolic warning light.  It contains the
unique hinge `[21,24]` and ends at `27=2*14-1`, the LRC14 shell modulus.  That
does not prove a theorem, but it suggests a concrete next diagnostic: for
Q27 blockers, measure not only coverage but also moment defects across the
fibered shell.  A row that scalar-blocks every denominator may still fail to
lift if the owner/fraction address cannot be made consistent.

The follow-up overlap classifier made the pattern cleaner.  The additive tower
covers every positive integer because its rows are the square shells
`[n^2,(n+1)^2-1]`.  The square tower covers only alternating triangular shells
`[T_{2m},T_{2m+1}-1]`, skipping gaps of size `2m+2`.

There are two special exact families:

```text
B_m equation side-aligned inside A_n equation
  iff T_n = 2T_m
  iff (2n+1)^2 - 2(2m+1)^2 = -1.
```

This is the infinite Pell family `(m,n)=(2,3),(14,20),(84,119),...`; it makes
`10..14` the first nontrivial example.  In each row, the B block sits in the
middle of A with symmetric outside padding `n-m`.

The second special family is actually a singleton:

```text
B_3.L = A_4.R = [21,24].
```

So `21+22+23+24` is not the first member of a hidden exact-side ladder.  It is
the unique whole-side equality.  The looser side containments are still
infinite, but they are governed by a simple floor-square rule: for a B side
`[u,v]`, set `n=floor(sqrt(u))` and compare `[u,v]` with the A midpoint
`n^2+n` and square boundary `(n+1)^2-1`.  The resulting word is
Beatty/Sturmian rather than Pell-exact.

The new proof target I would actually trust is modest and useful:

```text
Build a moment-lift resource ledger for LRC14.
```

For each denominator/fiber, record blocked unit twists, owners, support size,
raw moment defect, and the fractional address needed to cancel that defect.
Then compare AP, V*, 2AP, the one-stranger residuals from HYP-2444, and any
multi-stranger Q27 adversaries.  If the fractional addresses are incompatible
across fibers, they may force a witness or a Bprime/owner opening.  If they
are compatible only in AP/V*/2AP, that becomes the beginning of a proof route.
