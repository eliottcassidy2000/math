# LRC14 Blocking-Height Dominance

The first useful answer to "does dominance grow with blocking height?" is yes,
but not in the way I expected.

If dominance means raw accumulated cover margin over every shell before the
first leak, it grows strongly.  In the one-stranger family, height versus raw
mean pair margin has Pearson correlation `0.779`; in the random primitive
sample it is `0.942`.  That is a real signal: high blocking height means a row
has built a long history of perfect unit-cover shells, and the cover carriers
do separate over that history.

But if dominance means per-shell or per-covered-unit concentration, it moves
the other way.  The normalized pair margin has correlation `-0.711` in the
one-stranger family and `-0.729` in the random sample.  The speed tournaments
also saturate to transitive orders with no directed 3-cycles in every named
hard packet.  So the high rows are not wild dominance hierarchies.  They are
balanced-cover chains with accumulated but diluted dominance.

That feels proof-useful.  A long blocker has to pay in one of two currencies:

```text
cumulative excess that marks a peelable carrier,
or balanced cover congruences that force a next-fiber leak.
```

This slots directly into the Q31 and ramified-portal work.  HYP-2471 says the
two eight-core Q27 exceptions die when the fiber ladder widens to Q31.  HYP-2480
says those exceptions are ramified 7/13 packets.  HYP-2481 adds that before the
leak, their cover-load dominance has already accumulated, but it has also
spread out enough that scalar Q27 can hide the obstruction.  That is exactly
the kind of situation where a retained side channel beats a scalar quotient.

The next experiment should be leave-one-out criticality.  For each high row,
delete one speed and measure which pre-height shells fail first.  If the top
cover carrier is also deletion-critical, we get a peelable-carrier lemma.  If
it is not, then the row is balanced in a stronger sense, and that balance should
be expressible as a congruence/valuation packet.  Either branch is better than
raw search.

The assumption challenge mattered here.  I first wanted a tournament on speeds,
because the prompt asked for dominance and speeds are the obvious carriers.
That quotient did its job, but it became transitive too quickly.  To see the
next layer of structure, the vertices probably need to be unit obligations,
denominator shells, deleted-core addresses, or proof obligations.  The speed
tournament is the shadow; the proof lives in the lift.
