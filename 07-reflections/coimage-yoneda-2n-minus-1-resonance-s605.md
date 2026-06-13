---
source: codex-2026-06-03-S605
status: SYNTHESIS addendum; category + number theory perspective on anti-Poisson cancellation
tags: [LRC, coimage, Yoneda, category-theory, number-theory, 2n-minus-1, resonance, THM-401, anti-Poisson, cancellation]
---

# Coimage, Yoneda, and the `2n-1` Resonance Cancellation

The categorical sentence is:

```text
the master object is the coimage of the coverage observable.
```

The number-theoretic sentence is:

```text
at the LRC floor, the cancellation lives in the C = 2n - 1 resonance quotient.
```

Together they say something more precise than either sentence alone.  The
coimage tells us which quotient is allowed to be fundamental: collapse exactly
the distinctions the observable cannot see, keep exactly what still decides the
predicate.  Yoneda tells us how to recognize that quotient: probe it by every
natural observer, and if all roads represent the same answers, the object is
canonical.  The `2n-1` resonances supply the number-theory observers.

## The categorical side

Let

```text
N_delta(t) = #{v in V : ||v t|| < delta}.
```

The depth law `{p_k}` is the pushforward of Lebesgue measure through `N_delta`.
It is the coimage of `N_delta`: the minimal quotient of the clock circle that
still remembers every functional depending only on depth, especially

```text
p_0 = meas{t : N_delta(t)=0}.
```

The Yoneda reading is not ornamental.  It says that a proposed master object is
determined by its maps into all tests:

```text
depth distribution,
factorial moments,
spectral measure,
Helly certificates,
partition function,
residue/CRT witnesses,
unit-shell probes.
```

When all these probes keep reconstructing the same obstruction, the object is
not a favorite coordinate system.  It is the represented functor of the problem.
That is the categorical meaning of "fundamental."

## The number-theory side

At `delta = 1/n`, the sharp antipodal completion is

```text
C = 2n - 1.
```

At a clock `t = k/C`, the bad residues are `0,+1,-1 mod C`.  If a unit shell
`{a,-a}` is missed, multiplication by `k = a^{-1}` moves that missed shell to
`{+1,-1}` and produces a witness above the floor:

```text
M(V) >= 2/C > 1/n.
```

So a row can even hope to collapse `p_0` only after it survives this
number-theory quotient:

```text
V mod C -> antipodal shells {a,-a} -> unit action by (Z/CZ)^x.
```

This is the THM-401 / S571 resonance gate.  It is the concrete arithmetic
coimage of the endpoint witness problem.  Prime `C` makes every shell unit
visible.  Composite `C` splits the quotient into gcd strata; the nonunit shells
are where extra sporadics can hide.  For `n=14`, this is exactly `C=27=3^3`,
with unit, gcd-3, and gcd-9 lanes.

## The cancellation

THM-406 says

```text
p_0 = sum_j (-1)^j S_j,
S_j = total j-fold overlap volume.
```

The `2n-1` resonance quotient is not another low-order moment.  It is the
arithmetic compression that explains how the all-orders sum can cancel.  A
missed unit shell gives a visible witness and prevents cancellation.  A perfect
unit-shell cover removes that witness route, but still does not guarantee
collapse.  The remaining question is whether additive resonances align the
higher overlaps so the alternating sum lands exactly at zero.

This is where additive chains enter:

```text
v_k = v_i + v_j
```

means the three forbidden arcs co-resonate near the same clocks.  These
relations are not a cosmetic AP resemblance.  They are morphisms in the
resonance quotient: they force overlap packets to move together under the same
unit-character probes.  Anti-Poisson collapse is the case where enough such
relations align to cancel the whole depth polynomial at `0`.

## The combined slogan

```text
coimage  = the minimal quotient that can decide p_0;
Yoneda   = the quotient is canonical because every natural probe recovers it;
2n - 1   = the number-theory probe family at the LRC floor;
collapse = all-order inclusion-exclusion cancellation inside that quotient.
```

So the next invariant should be a resonance coimage, not a scalar moment:

```text
Res_C(V) =
  antipodal shell coverage,
  gcd strata,
  unit-character orbit data,
  additive-chain relations,
  witness-floor labels,
  depth-polynomial cancellation profile.
```

For `n=14`, compute this at `C=27`.  If `Res_27(V)` has a missed unit shell, it
exits by THM-401.  If it has only nonunit holes, the determinant/CRT and Helly
branches should decide whether the remaining higher-overlap packet cancels or
leaves positive lonely measure.  That is the promised meeting point of category
and number theory: the categorical coimage tells us what to keep, and the
number-theoretic resonance quotient tells us what the cancellation is made of.
