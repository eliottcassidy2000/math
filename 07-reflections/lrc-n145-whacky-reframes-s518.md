# Whacky n=145 Reframes (S518)

The number `145` is not just "large n."  It is `5*29`, with `144` moving
runners and threshold `1/145`.  That creates a strange exact wall:

```text
t = a/145, gcd(a,145)=1.
```

At that time every nonzero residue modulo `145` is safe.  The only blockers are
speeds divisible by `145`.  So if no speed is divisible by `145`, THM-369
already solves the system.

If a `145`-divisible speed does appear, something funny happens.  There are
only `144` moving runners.  One of them occupies residue `0`, so the speed set
cannot occupy all `144` nonzero residues.  Therefore some antipodal boundary
pair `{r,-r}` is missing at least one side.  Choose the unit wall with
`a^{-1}=r`, and one of the two LRC boundary sides is unpinned.  That is the
aperture.

So the whacky proof route is:

```text
unit wall -> one-sided aperture -> move the zero-residue embryo -> source
```

The "embryo" is the set of speeds divisible by `145`, scaled down by `145`.
It sits at the observer at every unit wall.  The local question is whether it
can move out of the forbidden cap before the pinned boundary side collapses.

This feels like the n=145 version of the observer-score story:

- source means zero blockers;
- the unit wall with zero residues is a positive blocker layer;
- the aperture tries to descend the observer score;
- failure should create endpoint-pressure labels.

The most useful thing S518 did was not solve n=145.  It found a very specific
place to attack: don't enumerate source classes, don't chase full A000568, and
don't begin with arbitrary time.  Begin at denominator `145`, where almost
everything is already safe, and study the quotient dynamics of the blockers
that remain.
