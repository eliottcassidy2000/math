---
id: HYP-1862
status: OPEN
source: codex-2026-05-31-S394
related:
  - HYP-1861
---

# HYP-1862: q(x)=x-1/x is the reciprocal counterterm for triangular log cutoffs

## Statement

Let

```text
p(x) = x      for 0 <= x <= 1
     = 1/x    for 1 <= x <= 2.
```

Then

```text
int_0^2 p(x) dx = 1/2 + ln 2.
```

The ordinary decomposition is triangle plus hyperbola:

```text
int_0^1 x dx + int_1^2 dx/x = 1/2 + ln 2.
```

But the more structural decomposition is

```text
p(x) = 1/x + q(x),

q(x) = x - 1/x    for 0 < x <= 1,
     = 0          for 1 <= x <= 2,
```

interpreted with a cutoff at `0`.  Here `q` is not a positive area.  It is a
reciprocal counterterm that cancels the infinite negative log-ray of `1/x` and
replaces it by the finite triangular mass `1/2`.

## Evidence

For `epsilon>0`,

```text
int_epsilon^2 dx/x = ln(2/epsilon)
int_epsilon^1 (x - 1/x) dx = 1/2 + ln(epsilon) - epsilon^2/2
```

so the combined finite part is

```text
1/2 + ln 2 - epsilon^2/2 -> 1/2 + ln 2.
```

In log-coordinate `x=e^u`, the hyperbola is flat measure:

```text
dx/x = du.
```

The patched area density is

```text
p(e^u)e^u = e^(2u)     for u <= 0,
           = 1         for 0 <= u <= ln 2.
```

Thus the triangular half-area is exactly the exponential tail
`int_{-infty}^0 e^(2u) du = 1/2`, while `ln 2` is the flat log-measure from
`0` to `ln 2`.

The bare defect

```text
q0(x)=x-1/x
```

is anti-invariant under reciprocal inversion:

```text
q0(1/x) = -q0(x).
```

It has paired even coordinate `r=x+1/x`, with

```text
r^2 - q0^2 = 4,
D q0 = r,
D r = q0,
D = x d/dx.
```

Equivalently `q0(e^u)=2sinh(u)`.

The Mellin transform of the one-sided bare counterterm is the two-pole kernel

```text
M_q(s) = int_0^1 x^(s-1)(x-1/x) dx
       = 1/(s+1) - 1/(s-1)
       = -2/(s^2-1).
```

## 30-60-90 correction

The exact `ln 2` remainder forces the carved triangle to have area `1/2`.
If a line from the origin reaches height `1` at angle `theta`, its area is
`1/(2 tan theta)`.  Therefore:

```text
theta=30 deg: area=sqrt(3)/2
theta=45 deg: area=1/2
theta=60 deg: area=1/(2sqrt(3)).
```

So a literal 30-60-90 triangle cannot keep the `ln 2` remainder while also
gluing to `(1,1)`.  It must pay a `q`-defect: an area surplus/deficit, an
endpoint jump, or a slope/gauge change.

## Predictions

1. Useful alternate decompositions of this area should be classified by their
   `q`-defect: area debt, endpoint jump, slope jump, and Mellin residue.
2. The isosceles triangle is the unique cutoff in the family `x^a` that is
   continuous at `(1,1)` and has left mass `1/2`.
3. The reciprocal counterterm `x-1/x` should recur anywhere a divergent
   logarithmic measure is replaced by a finite triangular cutoff.
4. The repository's recurring Cayley shadow `Q(Q(x))=-1/x` should be compared
   with this bare anti-reciprocal coordinate `x-1/x`.

## Sources

- `04-computation/hyperbola_triangle_q_s394.py`
- `05-knowledge/results/hyperbola_triangle_q_s394.out`
- `07-reflections/hyperbola-triangle-q-s394.md`
