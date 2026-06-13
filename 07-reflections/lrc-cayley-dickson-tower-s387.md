# LRC Cayley-Dickson Tower

codex-2026-05-31 S387

The Cayley-Dickson tower is the right analogy if we keep it disciplined.

It should not mean "LRC is secretly an algebra."  It should mean:

```text
each recursive doubling closes one old freedom and exposes a new defect.
```

In Cayley-Dickson:

```text
R -> C -> H -> O -> S
```

we pass through successive losses: order, commutativity, associativity,
then division.  In LRC, denominators split as

```text
n = 2^r * odd_core.
```

The `2^r` part is the doubling row.  The odd core is the torsion payload.
When we insert the mandatory `n`-gate, the unit skeleton closes, but endpoint
debt reappears below it.  That is the LRC analogue of zero-divisor leakage:
the system has solved the old equality but created a new place where
structure fails to divide cleanly.

## The Above-14 Window

The S387 table makes the first candidates read like a tower chart:

```text
14 = 2 * 7       complex row, odd core 7
16 = 16 * 1      pure sedenion-row laboratory
18 = 2 * 9       complex row, square 3-torsion payload
20 = 4 * 5       quaternion row, odd core 5
21 = 1 * 21      real row, 3-by-7 transfer payload
24 = 8 * 3       octonion row, stress-test payload
```

This clarifies why `18` still feels like the right favorite.  It lives on the
same first doubling row as `14`, but its odd payload is no longer prime.  It
is the first square-torsion version of the fourteen-runner story.

## The n=18 Layer

The lpd ladder at `n=18` is:

```text
(1, 9, 18, 27, 36, 45, 54, 63, 81, 90, 99, 108, 117, 126, 135, 144, 153)
```

It has:

```text
gap/th=0.005682
unprotected endpoints=176
first leak=11/162
```

The endpoint-debt layer histogram is:

```text
9:48, 27:16, 99:16, 117:16, 144:64, 153:16
```

The first layer is exactly `9`, the largest proper divisor.  That is the
beautiful part.  The old unit boundary is gone; the new leak is not random.
It lands precisely on the square-torsion payload.

## What This Suggests

For proof work:

```text
n=16:
  prove the pure doubling endpoint-debt lemma.

n=18:
  prove the first mixed-row lemma:
  closing the 18-unit skeleton forces debt at layer 9.

n=24:
  stress test the resulting certificate on an octonion-row denominator.
```

For disproof work:

```text
try endpoint-protection cycles at n=18 before speed sets.
```

If a counterexample-like architecture exists above fourteen, it should first
look like a protection cycle that closes the `18` unit layer, then the `9`
layer, then the `3`-torsion descendants.  Sampling speeds directly is still
too blunt; it sees the debt but not the cycle.

The slogan I would keep:

```text
LRC endpoint debt is Cayley-Dickson property loss on the denominator tower.
```

`18` is where that slogan first has enough structure to bite.
