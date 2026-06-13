# Rational Shadow Carriers S635

The prompt brought in a beautiful but slightly misremembered theorem: not
"exactly one of `pi+e` and `pi*e` is rational," but "they cannot both be
algebraic."  That correction makes the connection sharper.

If

```text
S = x + y,
P = xy,
```

then `x` and `y` are roots of `T^2 - S*T + P`.  So the pair of shadows
`(S,P)` is not just two measurements.  It is the elementary symmetric carrier
of the unordered pair `{x,y}`.  If both shadows descend to a tame field, the
hidden pair descends to a quadratic over that field.

For `{e,pi}`, that cannot happen over `Qbar`, because `e` and `pi` are
transcendental.  Therefore at least one of `e+pi` and `e*pi` is
transcendental.  We do not know which, and both could be.  The uncertainty is
not a defect in the analogy; it is the whole point.  The obstruction sees only
that the missing complexity must survive in at least one of the two shadows.

S635's finite-field toy model makes this visible.  Over `F_p`, unordered pairs
have large fibers under sum alone and product alone, while the joint
`(sum, product)` map is injective.  The single shadows are cheap and lossy.
The joint carrier is the reconstruction datum.

This changes the addition/multiplication reading.  In S365 and the summand
graph thread, addition and multiplication were rival operation shadows:
addition collapses to the transitive order, multiplication retains a divisor
DAG, and product-sum witnesses are critical pairs.  Here they are not rivals.
They are complementary symmetric coordinates.  Addition gives the trace;
multiplication gives the norm.  Trace plus norm gives the quadratic carrier.

That is the bridge to odd/even.  Parity is also a field-of-definition quotient,
only at the smallest field.  Redei's theorem says a tournament's Hamiltonian
path count is always odd; the parity shadow blocks total cancellation.  The
pi/e theorem says the `Qbar` shadow blocks both elementary symmetric coordinates
from being tame.  In both cases, the proof does not tell us the full object.  It
tells us that a compressed shadow cannot lose all of the hard information.

For LRC, the immediate translation is about clocks.  A rationally
commensurable speed set has a reset period after scale is removed; the orbit is
finite-periodic, and "hit the safe box once" can be asked inside that clock.
An irrational-rank speed set has no total reset, so the correct object is the
relation lattice and the torus closure.  The reset period and the relation
lattice are the two shadows.  Treating either as the whole state risks exactly
the kind of false scalar collapse S634 warned against.

This also clarifies which clocks matter:

- common scale does not matter;
- rational ratios matter because they create finite reset periods;
- irrational rank matters because it replaces reset by torus closure;
- short additive relations matter because they create folds, shields, and
  non-random depth;
- raw elapsed time matters only after those carriers have been named.

The useful wild version is this: LRC preprocessing should be a descent audit.
Ask first what field/module the speed set lives over.  If it descends to a
finite clock, enumerate the clock.  If it does not, compute the relation
lattice and use equidistribution/discrepancy.  If it has short circuits, fold
and sieve.  The time search should be the last representation, not the first.

For unit distance, the same story is already active in HYP-2210.  The number
`57` is not an endpoint; it is

```text
20 unit-spine edges + 37 centered-hex bulk edges.
```

For rooted tournament counts, `Q_n=(P_n+F_n)/2` is not "divide by two and hope."
It is Burnside with the fixed carrier named.  For sequence work, S633 says a
hard value should be surrounded by fixed, merged, nonfixed, q-deformed, and
transporter shadows before we declare it opaque.

So the S635 principle is:

```text
If two quotient shadows together reconstruct a forbidden tame carrier,
then at least one shadow must retain the hard coordinate.
```

That principle is small enough to be algorithmic.  It can become a helper that
takes a pair of shadows, a reconstruction rule, and a forbidden descent target.
Then it returns a proof obligation: name which side channel must survive, or
prove that at least one survives without naming it.

Tournament Analysis in S635 used proof-carrier lenses as vertices.  The
tournament was transitive and ranked the joint carrier first, then the field
descent obstruction, then the LRC relation lattice.  That ranking feels
correct for the current thread: a one-shadow is useful only when the carrier
that it forgets has already been shown irrelevant.

The next hard move is to push this back into LRC algorithms.  Before asking
whether the orbit hits the safe box, classify the speed set by:

```text
scale-free ratios;
rational rank and reset period if finite;
relation lattice rank;
short circuit supports;
fold/sieve denominators;
observer-source threshold payload.
```

That is the clock analogue of `(sum, product)`.  One component tells us where
the orbit can reset; the other tells us why it may fail to behave randomly
before or without a reset.
