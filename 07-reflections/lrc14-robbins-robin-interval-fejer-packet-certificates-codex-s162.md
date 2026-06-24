# LRC14 Fejer Intervals, Robbins Bridges, and Divisor Shadows

This session pushed the HYP-2974 Fejer evidence toward a proof interface.  The
floating certificates already know where they want to land: not on raw runner
sets, but on labelled packet fibers.

The useful surprise from S162 is that the interval problem appears modest in
precision.  The smallest selected margin, also the smallest full-bank S157
floating margin, is the two-swap row

```text
drop(12,13)->add(14,29),  Q_41(347/4312)=-3.36091433e-7.
```

Even with a deliberately conservative atom-L1 budget, that row asks for about
`27` bits / `9` decimal digits of relative atom enclosure.  The high-degree
packets are not the small-margin packets: `P10+GW` needs degree `280` and
`862` atoms but has margin about `9.6e-4`; `12->36` needs degree `159` but
margin about `2.4e-4`.  So the likely hard work is not numerical precision.
It is proof organization.

## Robbins as quotient discipline

Robbins' graph theorem says strong orientation is equivalent to no bridge.
That analogy feels unusually crisp here.  A Fejer certificate has a dependency
chain:

```text
exact rational center
  -> divisor atom bank
  -> trig interval enclosure
  -> signed negative margin
  -> labelled packet fiber
  -> route handoff
```

Forget a bridge and the proof stops being strongly connected.  This is the
same guardrail the repo has learned through tournament iso classes, unit
distance spines, C27/unital charts, Faulhaber moment shadows, Pollock defects,
and irreducibility/coefficient lifts: a quotient is legal only after its kernel
has been named.

There is a nice asymmetry with Robin's theorem.  Robin says a sigma-function
inequality can encode RH.  That makes scalar divisor shadows powerful, but the
project's lesson is that power is not safety.  Sigma/tau can carry global truth
while hiding the exact route label needed by the next local implication.
Dirichlet convolution and Ramanujan sums are better proof metaphors because
they keep the packet law or primitive-period trace instead of only the scalar
pushforward.

## Reset clocks versus packet clocks

The user's second-spectrum family `M_k={1,...,k-1,2k}` matters because it makes
reset-time thinking expensive.  The gap

```text
2/(2k+1) - 1/(k+1) = 1/((k+1)(2k+1))
```

is small enough that a lift-depth or total-reset-clock argument naturally pays
quadratic depth.  The Fejer route avoids measuring the whole orbit recurrence.
It asks a sharper local question:

```text
can this labelled packet fiber support a nonnegative cover function?
```

If no, one negative interval-enclosed quadratic form is enough.  That feels
like the right scale for LRC14: AP/GW are closed boundary atoms, and everything
else in the HYP-2963 bank already leaks a strict open interval.

## Creative theorem shape

The proof I would try next is not a universal bounded-degree theorem.  It is a
familywise certificate sheaf:

```text
AP/GW boundary atom
  or K33 interval certificate / state lift
  or petal interval certificate / two-block rigidity
  or covering interval certificate / boundary-moment ledger
  or few-apex lift interval certificate
  or new packet family.
```

The S162 output suggests each family may need only a small exact interval
template.  The problem becomes finite not because there are finitely many rows,
but because each surviving family has a labelled atom-bank grammar.

The most tempting implementation is:

1. Generate the Fejer atom bank as exact symbolic data:
   `(degree, center, k, v, m, rational Fejer weight)`.
2. Replace each trig value by an interval:
   `cos(2*pi*k*x)` is a root-of-unity algebraic interval,
   `sin(pi*m/7)` is one of six algebraic sine values up to sign,
   `1/pi` gets a rational interval.
3. Sum outward with rational directed rounding.
4. Store the resulting sign certificate next to the packet label.

That would turn HYP-2974 from "empirical PSD failure" into a reusable proof
module.

## Tournament-adjacent meaning

The tournament here is not a tournament of runners.  It is a tournament of
ways to remember a proof object.  The top path from S162 is:

```text
labelled interval certificate
> no-bridge assembly
> endpoint owner anchor
> Fejer divisor atom bank
> Ramanujan exact-period projector
> Dirichlet convolution packet law
> floating Fejer shadow
> Robin sigma scalar shadow
> raw divisor count
```

This is almost a philosophy of the repo in one chain.  Raw counts are useful
first derivatives.  Packet laws are second derivatives.  Interval certificates
are the point where the derivative data becomes a sign proof.

## Open speculation

There may be a clean "ear decomposition" for certificates:

```text
base cycle = AP/GW boundary equality atom
ear = one labelled family extension
orientation = interval sign / route handoff
bridge = forgotten packet coordinate
```

If this can be made precise, then the LRC14 proof would look like a strong
orientation theorem for the packet metagraph: every non-AP/GW ear can be
oriented toward a strict interval certificate or toward a forbidden state lift.
That would merge the tournament program and the Fejer program in the most
natural way I can see right now.
