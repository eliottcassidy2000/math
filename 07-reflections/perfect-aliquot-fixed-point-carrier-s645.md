# Perfect Numbers as Aliquot Fixed Points S645

HYP-2220/S644 found the triangular/Vieta face of this story.  This is the
direct aliquot-map face: a very clean fixed-point carrier.

```text
divisor lattice D(n) -> sigma(n) -> s(n)=sigma(n)-n
```

A perfect number is a literal fixed point:

```text
s(n)=n.
```

That framing is better than saying only `sigma(n)=2n`, because the aliquot map
has neighboring dynamics.  A perfect number is a length-1 sociable cycle.
Amicable pairs are length-2 cycles.  Sociable cycles are higher-period loops.
The fixed point is not isolated in concept-space; it is the smallest possible
period of the divisor-sum flow.

The side channel is the divisor lattice.  The scalar `A(n)=sigma(n)/n` is the
visible quotient, but the proof object is the product

```text
product_{p^a || n} (1 + 1/p + ... + 1/p^a).
```

For the even perfect section, everything snaps shut:

```text
n = 2^(p-1) M
M = 2^p - 1 prime
sigma(2^(p-1)) = M
sigma(M) = 2^p
sigma(n) = 2n.
```

So Euclid-Euler is not just a construction.  It is a section of the divisor
carrier where the Mersenne prime side channel exactly balances the dyadic
geometric series.  The power of two is almost perfect by itself: its defect is
`-1`.  Adding a Mersenne prime factor pays exactly that missing unit at the
right scale.

The odd side is where the carrier language helps.  S645 finds no odd fixed
points through `100000`, but it does find odd rows extremely close to fixed.
The best in the scan is:

```text
32445 = 3^2 * 5 * 7 * 103
sigma(n)-2n = 6.
```

That is the little danger zone.  Scalar proximity to `2` does not decide the
fixed predicate.  The obstruction has to live in local prime-power weights,
congruence residues, parity, and the way the aliquot arrow enters or misses a
cycle.  This is the same lesson as the recent pi/e and LRC packets: the scalar
quotient is the public face; the proof is in the retained side channel.

The analogy stack now feels like this:

```text
pi/e:            trace/norm shadows reconstruct the hidden pair
Goldbach/Lemoine: E/O shadows reconstruct an ordered prime pair
LRC/UD:          scalar frontiers need owner/carry/spine/bulk labels
perfect numbers: sigma(n)=2n needs divisor-lattice product labels
```

One possible new technique is an aliquot side-channel jackknife.  Delete one
prime-power factor from `n`, recompute the local abundancy product, and ask how
much fixed-point mass is released or lost.  For even perfect numbers, the
dyadic factor and Mersenne factor are a perfect two-part cancellation.  For odd
near-fixed rows, the same jackknife could expose which local prime-power factor
is carrying the residual defect.

Another useful object is the defect shell:

```text
Delta(n) = sigma(n) - 2n.
```

Fixed points are `Delta=0`, powers of two sit on `Delta=-1`, quasi-perfect
would be `Delta=1`, and odd near-fixed rows are small nonzero shells.  Instead
of searching only for `Delta=0`, build the shell graph of how aliquot arrows
move defect.  The hope is to separate "near fixed because structurally close"
from "near fixed by accidental scalar cancellation."

The tournament lesson from S645 is unsurprising but useful: divisor-lattice
fixed-point data ranks first, the aliquot graph and sigma product come next,
and raw sequence listing comes last.  The good proof object should remember
the loop and the lattice at the same time.
