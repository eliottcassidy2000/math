---
id: THM-368
name: hyperbola-triangle-q-gauge
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S395
depends_on: []
results:
  - 05-knowledge/results/hyperbola_triangle_q_gauge_s395.out
---

# THM-368: Hyperbola-Triangle q-Gauge

## Statement

The area

```text
int_0^1 x dx + int_1^2 dx/x
```

is

```text
1/2 + log(2).
```

There are two equivalent q-decompositions.

First, for any integrable `q` on `[0,2]` with

```text
int_0^2 q(x) dx = 1/2,
```

the Hadamard finite part satisfies

```text
FP int_0^2 (1/x + q(x)) dx = log(2) + 1/2.
```

In particular, the equal-area triangle family

```text
q_L(x) = x/L^2       for 0 <= x <= L,
q_L(x) = 0           otherwise,
```

has this property for every `0<L<=2`.

Second, the original ordinary splice is obtained pointwise from the
reciprocal-odd counterterm

```text
q_ct(x) = x - 1/x    for 0 < x < 1,
q_ct(x) = 0          for 1 <= x <= 2.
```

Namely,

```text
x = 1/x + q_ct(x)    for 0 < x < 1,
1/x = 1/x + q_ct(x)  for 1 <= x <= 2.
```

With a cutoff,

```text
lim_{eps->0+} [ int_eps^2 dx/x + int_eps^1 (x - 1/x) dx ]
  = 1/2 + log(2).
```

The counterterm has reciprocal parity

```text
q_ct(1/x) = -q_ct(x).
```

## Proof

The first displayed area is immediate:

```text
int_0^1 x dx = 1/2,
int_1^2 dx/x = log(2).
```

For the finite-part statement, define

```text
FP int_0^2 dx/x = lim_{eps->0+} (int_eps^2 dx/x + log(eps)).
```

Since `int_eps^2 dx/x = log(2)-log(eps)`, the finite part is `log(2)`.
Adding any integrable `q` of total mass `1/2` gives `log(2)+1/2`.

For the line family,

```text
int_0^L x/L^2 dx = L^2/(2L^2) = 1/2.
```

Thus every `q_L` belongs to the same finite-part gauge class.

For the pointwise splice, direct substitution gives

```text
1/x + (x - 1/x) = x
```

on `(0,1)` and `1/x+0=1/x` on `[1,2]`.  The cutoff computation is

```text
int_eps^2 dx/x = log(2) - log(eps),
int_eps^1 (x - 1/x) dx = (1-eps^2)/2 + log(eps).
```

Their sum is

```text
log(2) + (1-eps^2)/2,
```

which tends to `log(2)+1/2`.

Finally,

```text
q_ct(1/x) = 1/x - x = -q_ct(x).
```

## Geometry

The equal-area line family clarifies the triangle choices:

```text
base L, height 1/L, slope 1/L^2, area 1/2.
```

Thus:

```text
L=1:  slope 1,   the 45-degree isosceles right triangle;
L=2:  slope 1/4, the equal-area triangle spread across x=0..2.
```

A Euclidean 30-degree equal-area triangle has slope `1/sqrt(3)`, hence

```text
L = 3^(1/4),
height = 3^(-1/4),
log L = log(3)/4.
```

So a literal base-two 30-degree triangle does not preserve the value unless
one rescales the vertical coordinate.  The invariant is not the visual angle;
it is the area of the q-packet.

## Log Coordinate

Under `x=e^t`,

```text
dx/x = dt.
```

The hyperbola is uniform in log coordinate, and the original splice becomes

```text
int_{-infty}^0 e^(2t) dt + int_0^log(2) dt = 1/2 + log(2).
```

The reciprocal-odd counterterm becomes

```text
q_ct(x) = x - 1/x = 2 sinh(t).
```

The normalized mass `2 q_L(x) dx` is an exponential tail in `t`, supported on
`t <= log L`, with mean `log L - 1/2` and variance `1/4`.

## Mellin View

For `Re(s)>-2`,

```text
int_0^L x^s q_L(x) dx = L^s/(s+2).
```

The triangle packet moves the hyperbola's Mellin pole at `s=0` to the
packet pole at `s=-2`.  Changing `L` multiplies by `L^s`, which is translation
in log coordinate.

## Verification Record

`04-computation/hyperbola_triangle_q_gauge_s395.py` prints:

1. the original area;
2. the equal-area triangle gauge table;
3. the base-two 30-60-90 normalization warning;
4. the cutoff cancellation for `q_ct`;
5. the log-coordinate and Mellin-transform forms.

Stored output:

```text
05-knowledge/results/hyperbola_triangle_q_gauge_s395.out
```

## Related

- HYP-1863: q-gauge principle for logarithmic constants.
- HYP-1862: reciprocal q-counterterm.
- `04-computation/hyperbola_triangle_q_gauge_s395.py`.
