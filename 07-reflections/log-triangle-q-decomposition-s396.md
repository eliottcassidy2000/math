# Log Triangle q-Decomposition

codex-2026-05-31 S396

The starting area is clean:

```text
f(x)=x      on [0,1]
f(x)=1/x    on [1,2]

integral f = 1/2 + ln(2).
```

So the first decomposition is literally:

```text
isosceles right triangle of area 1/2
plus hyperbola strip of area ln(2).
```

The proposed second decomposition needs one correction.  If the carved piece
is a Euclidean triangle spanning `x=0` to `x=2` and the remainder is still
`ln(2)`, then the triangle must have area `1/2`.  The ramp line must be

```text
y = x/4.
```

Its angle is about `14.036` degrees, not `30` degrees.  A literal 30-60-90
triangle with area `1/2` has hypotenuse `2/3^(1/4)<2`, so it cannot span the
whole interval from `0` to `2`.  Conversely, the natural span-2 interpretations
of a 30-60-90 triangle have the wrong area.

The mathematically stable version is therefore not 30-60-90 in Euclidean
coordinates.  It is a 4:1 ramp triangle plus a signed correction `q`.

## The q Object

Subtract a general ramp `a*x` over `[0,2]`.  Then

```text
f(x)-a*x = H(x-1)/x + q_a(x),
```

where

```text
q_a(x)=(1-a)x  on [0,1]
q_a(x)=-a*x    on [1,2].
```

The area of `q_a` is

```text
integral q_a = 1/2 - 2a.
```

So the log-preserving choice is uniquely

```text
a=1/4.
```

At this value,

```text
q(x)=3x/4   on [0,1]
q(x)=-x/4   on [1,2],
integral q = 0.
```

This is the key reframing: `q` is not another positive area.  It is a
zero-mass transport term.  It moves area from the left side to cancel area on
the right side, leaving the logarithmic contribution unchanged.

An antiderivative with zero endpoints is

```text
Q(x)=3x^2/8      on [0,1]
Q(x)=(4-x^2)/8   on [1,2].
```

So `q=Q'` is a flux across the join at `x=1`.

## Mellin Tangent

The moments of the log-preserving `q` are exact:

```text
M_n = integral_0^2 x^n q(x) dx = (1-2^n)/(n+2).
```

The Mellin transform is

```text
M(s)=integral_0^2 x^(s-1) q(x) dx
    =(1-2^(s-1))/(s+1).
```

Thus

```text
M(1)=0,
M'(1)= integral q(x) log(x) dx = -ln(2)/2.
```

That is the part worth keeping.  The zero-area correction does not change the
area, but its first scale derivative remembers exactly half a bit.

The tangent I would pursue next:

```text
logarithmic remainders = zero-mass q-transports with nonzero Mellin derivative.
```

This may connect back to the repo's rapidity/log-ratio material: the geometric
piece fixes mass, while the `q` transport stores the scale information.
