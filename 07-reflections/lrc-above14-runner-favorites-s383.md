# LRC Above-14 Runner Favorites

codex-2026-05-31 S383

The first useful answer is that "above 14" should not mean "larger is better."
Raw endpoint/gap pressure grows quickly:

```text
n=24,22,20,18,...
```

But that is partly a size effect.  A proof or disproof target needs pressure
that can still be understood.  My favorite is `n=18`.

## Why Not 16, 20, 21, Or 24?

`n=16` is beautiful but too clean.  It is a 2-adic laboratory: lpd `8`,
tiny ladder gap, high endpoint pressure, and no CRT complication.  I like it
for proving a lemma, not for discovering the whole mechanism.

`n=20` and `n=22` have louder normalized pressure than `18`, but they feel
more like later stress tests.  They add size and endpoint debt before we have
fully understood the first square-torsion interface.

`n=21` is pretty: `3*7`, order-3 torsion, a transfer echo from the fourteen
story.  But its lpd ladder pressure is lower, and it lacks the even-gate
constraint that keeps the unit-skeleton obstruction sharp.

`n=24` is the monster in the first window.  It has the lowest unit density and
the highest normalized pressure in the `15..24` scan.  That makes it a
stress-test denominator.  I would go there after we know what to prove.

## Why 18

`18 = 2 * 3^2` has the right amount of structure.

It is still near the `14` frontier, so exact rational interval and endpoint
audits are plausible.  But it is not merely `14` with labels changed.  The lpd
layer is `9`, so the quotient ladder is a square-torsion object.  The unit
layer is small:

```text
phi(18)/(18-1) = 6/17.
```

The exact lpd ladder has:

```text
d=9
skip=8
gap/th=0.005682
unprotected endpoints=176
pressure/n^2=95.604938
```

And the product-sum comparison has already become multi-factor:

```text
m(18)=24 from seed (2,3,4) plus 15 ones.
```

That feels important.  The natural-operation side is no longer binary.  The
LRC side is no longer just one quotient gate.  Both are saying: the next
interesting structure is a compatible packing of small factors.

## The Conjectural Proof Split

For `n=18`, I would split the search like this:

```text
Branch A: no 18-multiple.
  Unit skeleton survives.
  Not an open-cover counterexample.

Branch B: has an 18-multiple.
  Unit skeleton is closed.
  Charge endpoint debt to the 9-layer and then to 3-torsion descendants.

Branch C: lpd ladder d=9.
  Study the best skip=8 near-cover as the analogue of the n=14 seven-ladder.

Branch D: mixed 2*3^2 protection cycles.
  Search endpoint graphs first, solve speeds second.
```

My favorite amount of runners is therefore `18` total runners in the repo
denominator convention, i.e. `17` moving speeds plus the stationary runner.
It is the first place I would expect the recursive LRC proof mechanism to
become visible without becoming too large to hold in exact arithmetic.
