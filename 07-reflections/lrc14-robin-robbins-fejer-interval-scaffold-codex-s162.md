# LRC14 Robin/Robbins Fejer Interval Scaffold

The divisor-function page gives two different "Robins/Robbins" lessons that
should not be collapsed.

Robin's theorem is the number-theory scalar inequality side: the sum-of-divisors
function `sigma(n)` and the extremal ratio `sigma(n)/n` are powerful precisely
because they summarize a multiplicative divisor packet.  For LRC14 that is a
warning, not the final object.  Scalar divisor data, qdiv, unit counts, lcm
payloads, and even raw `c_14` profiles can mix AP/GW, q-witness, K33, petal,
and covering routes unless the phase and endpoint labels remain attached.

Robbins' theorem is the better structural analogy.  A graph has a strong
orientation exactly when no bridge blocks bidirectional reachability.  The LRC14
certificate stack should be read the same way: a quotient is admissible only if
the certificate graph still has the bridges it needs.  The load-bearing bridges
are not just runners.  They are endpoint owners, exact q-phases, Toeplitz degree
and center, labelled packet route, K33/state-lift debt, and the AP/GW equality
atom.  Forget one of these and the quotient must carry a residual defect
certificate explaining how that forgotten coordinate is reconstructed or
annihilated.

The S162 script turns this into a first concrete Fejer certificate format.  For
each selected row it attaches the Fejer quadratic form to the labelled packet
fiber

```text
P(S) = (route, family, q_class, packet_route, state_lift, q_threshold),
```

then evaluates

```text
Q_d(x)=6/7+2*sum_{1<=k<=d}(1-k/(d+1))*c_k*cos(2*pi*k*x),
c_k=sum_{v|k, v in S} sin(pi*(k/v)/7)/(pi*(k/v)).
```

The interval arithmetic is exact-rational around a rational enclosure of `pi`
and Taylor enclosures for the trigonometric terms.  It certifies negative upper
endpoints for:

```text
near/K33 12->36                         degree 159
P10+GW                                  degree 280
covering 12->168                        degree 63
two drop(12,13)->add(14,29)             degree 41
single swap 6->63                       degree 266
```

This does not yet close the formal gap.  The hard-coded `pi` enclosure needs a
formally sourced backend, and the full HYP-2963 bank needs packet-fiber grouping
rather than one-off rows.  But the key obstruction has moved: the Fejer
certificates are no longer merely floating signs.  They can be written as
finite interval certificates attached to `P(S)` fibers, with AP/GW remaining
the only PSD-blind equality atoms in the surrounding S157 audit.

## Divisor-Function Trail

Concepts worth retaining from the one-hop reading:

- `sigma_k` is multiplicative but not completely multiplicative; its product
  formula is useful only while prime-power exponents remain visible.
- Dirichlet convolution makes `sigma = Id * 1` and `Id = sigma * mu`, reinforcing
  the rule that an averaged quotient needs an inversion side-channel.
- Ramanujan sums are exact-period primitive-root traces; `c_q(n)` is the
  Mobius-inverted primitive shell and depends on `gcd(q,n)`.
- The divisor summatory function and Dirichlet divisor problem reframe divisor
  counting as lattice-point/hyperbola geometry, matching the LRC lesson that
  scalar counts need their carrier region.
- Superabundant/Robin extremality is an ordered-maxima story for `sigma(n)/n`;
  it is analogous to the LRC tight-locus only after the divisor packet is kept
  labelled.

## Next Technical Target

The next useful script should not rescan for floating degrees.  It should take
the existing S157 first-hit table, group rows by packet key and Fejer center
type, and emit a certificate file whose rows are:

```text
packet_key, speeds or family parameters, degree, center, pi_interval_id,
sin/cos interval list, rational upper bound < 0.
```

That file is the natural input for a Lean/arb bridge.  If a packet family cannot
share a common interval template, that failure is itself a quotient guardrail:
the family is not allowed to forget the coordinate that made the interval split.
