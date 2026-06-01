---
id: HYP-1992
name: lrc-rapidity-formal-group-bridge
status: OPEN
date: 2026-06-01
session: opus-2026-06-01-S520
---

# HYP-1992: LRC rapidity threshold via the formal group F(x,y)=(x+y)/(1+xy)

## Statement

The formal group F(x,y) = (x+y)/(1+xy) = tanh(atanh(x) + atanh(y)) is
addition in rapidity (hyperbolic arc-length) space.

Map each runner's circular distance g_i = ||v_i t|| from the observer to
a rapidity:

    r_i = atanh(1 - 2 g_i)

Then the LRC threshold g_i >= 1/n becomes:

    r_i <= atanh((n-2)/n) = 0.5 * ln(n-1)

and the observer is lonely iff all individual rapidities are at most
0.5 * ln(n-1).

**Conjecture:** The formal group controls both H(T) (via the fiber fraction
and Walsh transform) and the LRC lonely measure (via the rapidity threshold).
There should be a direct relationship between the tournament-theoretic
invariant H and the runner-theoretic lonely time.

## Evidence

- The rapidity threshold 0.5 * ln(n-1) is the Riemannian distance from
  the observer to the safe-arc boundary on the hyperbolic line.
- The total rapidity budget (n-1) * 0.5 * ln(n-1) grows as O(n log n).
- The formal group's fiber fraction f(n) = (1/2)_{n-2} / (n-2)! ~ 1/sqrt(pi n)
  involves the same Pochhammer symbol that appears in hyperbolic geometry.
- The "relativistic energy" (product of cosh(r_i)) was computed for initial
  segments and shows clean scaling.

## Why this matters

The formal group is the project's central algebraic object. If it also
controls LRC, then:
1. Techniques from tournament theory (OCF, Burnside, Krawtchouk) may apply
   to LRC.
2. The "everything is the triangle" geometric framework may extend to
   the circular runner problem via hyperbolic geometry.
3. The Cayley-Dickson tower (n=2,3,5,9,17,...) may identify special
   LRC denominators.

## Next steps

1. Compute the rapidity spectrum at known tight LRC examples.
2. Check whether the rapidity sum or product has a clean formula
   at the lonely time t = 1/(2n) for initial segments.
3. Connect to the half-turn phase tournament (THM-374) via rapidity.
4. Check whether the rapidity metric induces a natural tournament
   on runners (closer rapidity = higher competition).
