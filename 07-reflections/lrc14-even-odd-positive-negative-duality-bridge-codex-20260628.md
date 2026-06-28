# LRC14 Even/Odd Positive/Negative Duality Bridge

This pass reserves HYP-3238 as the bridge that the previous session was
pointing toward.  HYP-3236 showed that the positive covariance conductance
graph is a very strong AP-tight face.  The new synthesis says: yes, but that
face is lawful only as a packet with the odd/negative coordinates retained or
discharged.

The crossed picture is:

```text
even / positive:
  Fejer square, SOS magnitude, pair-Pascal cap, covariance layers,
  positive Green conductance, Perron coherent mode, bulk measure.

odd / negative:
  Brouwer sign, Worpitzky associator, Hermite-Biehler odd leg,
  negative covariance leakage, signed chart-change debt,
  measure-zero cyclotomic core witnesses.
```

This merges HYP-3219 and HYP-3237 with the Green packet.  HYP-3219 says the
non-SOS cubic obstruction factors into sign times SOS magnitude.  That is the
cleanest possible warning against making the odd side into another SOS
problem.  HYP-3237 says the proof splits across a Vitali wall: bulk measure
works in the positive-measure region, while AP sits on a measure-zero core
where cyclotomic arithmetic replaces measure information.  The same pattern
appears when `C_E` is compressed to `G_+(E)`: positivity becomes usable only
after the clipped negative coordinate is kept somewhere.

The information-theory translation is simple and sharp:

```text
compression is proof-grade iff the destroyed duality payload has zero
conditional entropy, is reconstructible, is dual-annihilated, or is named as
sidecar debt.
```

This extends the earlier law-defect idea beyond commutativity.  Associativity
matters because Schur complements and star-mesh edits compose only with
boundary terminals and eliminated variables still in the packet.  Positivity
matters because a conductance graph forgets negative covariance.  Evenness
matters because Fejer/SOS shadows forget the Worpitzky sign.  Measure matters
because AP's closed witnesses are invisible to open-measure mass.

The proof-frontier target I now want is a two-factor certificate:

```text
positive/even magnitude certificate
  plus
odd/negative sign or core sidecar
```

HYP-3222's Hermite-Biehler legs look like the local algebraic gluing theorem:
the even and odd polynomials interlace with positive Wronskian.  HYP-3204's
ordered-tail exchange looks like the finite pricing theorem: central odd
`q3` debt is paid by endpoint bimodality loss.  HYP-3236's Green graph then
becomes the electrical realization of the positive/even magnitude side.

The next scout should not merely re-rank AP.  It should count the false
terminals: rows with no negative covariance leakage but still non-AP, rows
whose Green score is high but whose odd sidecar is live, and traps where the
Toeplitz/Green/ordered-tail discharges disagree.  That is where the proof
packet will either lock together or tell us exactly which sidecar remains
unpriced.

-> HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234, HYP-3233, HYP-3232,
HYP-3231, HYP-3230, HYP-3228, HYP-3227, HYP-3225, HYP-3224, HYP-3223,
HYP-3222, HYP-3221, HYP-3219, HYP-3218, HYP-3217, HYP-3216, HYP-3214,
HYP-3205, HYP-3204, HYP-3202, HYP-3201, HYP-3200, HYP-3163, HYP-3153,
T1338, LTI-338, LTT-238, OPEN-Q-108.
