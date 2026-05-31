# Hyperbola Triangle q-Gauge

codex-2026-05-31 S395

The original area is plain:

```text
int_0^1 x dx + int_1^2 dx/x = 1/2 + log(2).
```

But the user's point is not the arithmetic.  It is the two different ways of
seeing what the `1/2` is doing.

There are really two q's.

The first is a triangle-density gauge.  Pick a base length `L` and set

```text
q_L(x) = x/L^2 on [0,L].
```

Its area is always `1/2`.  Therefore

```text
FP int_0^2 (1/x + q_L(x)) dx = log(2) + 1/2.
```

This is the sense in which a triangle can move without changing the value.
The 45-degree triangle is `L=1`.  The equal-area triangle spread across
`x=0..2` is `L=2`, so its slope is `1/4`.  A literal Euclidean 30-degree
equal-area triangle has slope `1/sqrt(3)`, hence `L=3^(1/4)`.  Its log cutoff
is `log(3)/4`.

That is a nice little surprise: the 30-60-90 triangle injects a quarter-trit
of rapidity, while the base-two equal-area triangle injects one bit of
rapidity.  The visual angle is not the invariant.  The invariant is the mass
of the packet and where it sits in log coordinate.

The second q is more structural:

```text
q_ct(x) = x - 1/x.
```

Below `1`, the original line is exactly

```text
x = 1/x + q_ct(x).
```

So this q is a counterterm.  It cancels the infinite left tail of the
hyperbola.  With a cutoff,

```text
int_eps^2 dx/x + int_eps^1 (x-1/x) dx
```

converges to `1/2 + log(2)`.

This q has the feature that made the whole tangent click:

```text
q_ct(1/x) = -q_ct(x).
```

It is reciprocal-odd.  In log coordinate `x=e^t`,

```text
dx/x = dt
q_ct = 2 sinh(t)
x + 1/x = 2 cosh(t).
```

So the hyperbola is the invariant log measure, and `q_ct` is the odd
coordinate around the reciprocal fixed point `x=1`.  That is exactly the
continuous version of the anti-palindromic objects that keep appearing in the
tournament polynomial work:

```text
Q(x) = -x^d Q(1/x).
```

The triangle is the visible geometry.  The q-counterterm is the algebra.

This also explains why `1/2 + log(2)` feels like a hybrid number.  The
`log(2)` part is a length in rapidity.  The `1/2` part is a finite q-packet,
an exponential tail in rapidity:

```text
int_{-infty}^0 e^(2t) dt = 1/2.
```

So the original splice is

```text
exponential packet on (-infty,0]
+ uniform log measure on [0,log(2)].
```

That is a cleaner phrase for what is happening than "triangle plus hyperbola."
It is a packet plus a log interval.

Possible next thread: take the repo's rapidity lattice

```text
Z log(2) + Z log(3) + Z log(7)
```

and ask which constants are really log intervals, which are q-packets, and
which are cancellations between reciprocal-even and reciprocal-odd pieces.
The humble `x - 1/x` may be the continuous seed of a lot of the old
anti-palindromic structure.
