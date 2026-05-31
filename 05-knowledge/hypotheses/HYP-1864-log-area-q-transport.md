---
id: HYP-1864
status: OPEN
source: codex-2026-05-31-S396
related:
  - HYP-1821
  - HYP-1822
  - HYP-1838
---

# HYP-1864: Logarithmic area remainders are zero-mass q-transports

## Statement

When a mixed polynomial/reciprocal area decomposition leaves a logarithmic
remainder, the non-logarithmic part can often be adjusted until the residual
correction is a signed zero-mass transport `q`.  The logarithm then survives
as a scale derivative of `q`, not as extra area.

In the basic two-piece model

```text
f(x)=x      on [0,1]
f(x)=1/x    on [1,2],
```

subtracting a ramp triangle `a*x` over `[0,2]` gives

```text
f(x)-a*x = H(x-1)/x + q_a(x)
q_a(x)=(1-a)x  on [0,1]
q_a(x)=-a*x    on [1,2].
```

The correction has zero mass exactly at `a=1/4`.  At that value, the removed
triangle has area `1/2`, so the remaining area is exactly `ln(2)`.

## Evidence

`log_triangle_q_decomposition_s396.py` verifies the exact ledger:

```text
area(f) = 1/2 + ln(2)
area(a*x over [0,2]) = 2a
integral q_a = 1/2 - 2a
```

Thus `q_a` is area-neutral iff `a=1/4`.

For the log-preserving correction

```text
q(x)=3x/4   on [0,1]
q(x)=-x/4   on [1,2],
```

the moments are

```text
M_n = integral_0^2 x^n q(x) dx = (1-2^n)/(n+2).
```

The Mellin transform is

```text
M(s)=integral_0^2 x^(s-1) q(x) dx = (1-2^(s-1))/(s+1),
```

so `M(1)=0` while

```text
M'(1)= integral q(x) log(x) dx = -ln(2)/2.
```

The same script also records the Euclidean obstruction to the literal
30-60-90 reading: a 30-60-90 triangle with area `1/2` has hypotenuse
`2/3^(1/4)<2`, so it cannot span `x=0` to `x=2`.  Span-2 30-60-90
normalizations have areas `2/sqrt(3)`, `2sqrt(3)`, or `sqrt(3)/2`, not
`1/2`.

## Predictions

1. Similar log constants in the repo's rapidity/natural-operation threads can
   be modeled by zero-mass `q` corrections whose Mellin derivatives recover
   `ln` ratios.
2. The correct invariant is not the carved geometric piece by itself but the
   pair `(zero-mass correction, scale derivative)`.
3. When a proposed triangle/finite geometric carve leaves a log remainder,
   checking whether the induced `q` has zero mass should identify the right
   normalization immediately.

## Sources

- `04-computation/log_triangle_q_decomposition_s396.py`
- `05-knowledge/results/log_triangle_q_decomposition_s396.out`
- `07-reflections/log-triangle-q-decomposition-s396.md`
