---
id: HYP-1863
status: OPEN
source: codex-2026-05-31-S395
related:
  - THM-368
  - HYP-1619
  - HYP-1217
  - HYP-1862
---

# HYP-1863: q-gauges expose the finite part behind logarithmic constants

## Statement

Many appearances of `log(2)` in the repo should be revisited as finite-part
statements of the form

```text
invariant measure dx/x
+ finite q-packet
= visible constant.
```

In the hyperbola-triangle model, the packet can be written either as an
equal-area triangle density `q_L`, or as the reciprocal-odd pointwise
counterterm

```text
q_ct(x) = x - 1/x.
```

The second form is likely the structural one.  It is the continuous analogue
of anti-palindromic Laurent data:

```text
q_ct(1/x) = -q_ct(x).
```

## Evidence

THM-368 proves the exact identity

```text
int_0^1 x dx + int_1^2 dx/x
= FP int_0^2 (1/x + q_L(x)) dx
= lim_{eps->0+} [int_eps^2 dx/x + int_eps^1 (x-1/x) dx]
= 1/2 + log(2).
```

The equal-area triangle gauges

```text
q_L(x)=x/L^2 on [0,L]
```

all have total mass `1/2`.  Their shape changes, but the finite part stays
fixed.  In log coordinate `x=e^t`, `dx/x=dt`, and the triangle packet becomes
an exponential tail with mean `log L - 1/2`.

The 30-60-90 tangent becomes arithmetic in the same coordinate:

```text
30-degree equal-area triangle: L=3^(1/4), log L=log(3)/4;
base-two equal-area triangle:  L=2,       log L=log(2).
```

Thus the visual triangle angle is a gauge choice, while the invariant data
are the packet area and its log-coordinate support.

## Predictions

1. Some old `log(2)` appearances should split into a log-measure part plus a
   q-packet finite part, especially in rapidity, Cayley-transform, and
   Burnside-temperature notes.
2. Anti-palindromic tournament polynomials should be read as discrete
   versions of `x-1/x`: reciprocal-odd counterterms that make a cancellation
   finite.
3. The family `q_L` should behave like a one-dimensional renormalization
   group orbit: changing `L` translates the packet in log coordinate while
   preserving its mass.
4. A useful next computation is to scan existing `log(2)`, `log(3)`, and
   `log(7)` rapidity entries for expressions that can be rewritten as
   equal-area q-packets plus `dx/x` finite parts.

## Sources

- THM-368.
- `04-computation/hyperbola_triangle_q_gauge_s395.py`
- `05-knowledge/results/hyperbola_triangle_q_gauge_s395.out`
- `07-reflections/hyperbola-triangle-q-gauge-s395.md`
