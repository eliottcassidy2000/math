---
source: kind-pasteur-2026-07-24-S168 (Opus 4.8)
status: RESULT (concrete grounding, complements HYP-9078). Verifies the radial reduction
  L(f^m) = (dm+1)! * M(g^m) for homogeneous f (g(s)=f(s,1-s), M the [0,1]-Lebesgue moment), so FC(2)-homogeneous
  IS the classical Lebesgue moment problem = Liu-Sun 2020 Thm 2.6. Places the external int_0^1 e^Q claim exactly:
  int_0^1 e^{g} = sum_m M(g^m)/m! is the generating SUM (weight 1/m!), STRICTLY WEAKER than each M(g^m)=0, so
  int e^Q != 0 does not give FC even homogeneously -- the claimed "=> FC(2)" must be NON-homogeneous, via the
  friend-cited horizontal-endomorphism rigidity, not the radial reduction.
tags: [factorial-conjecture, exponential-integral, liu-sun, moment-problem, homogeneous, grounding]
related: [HYP-9078, kps-S165, kps-S166]
---

# The radial reduction places int e^Q relative to Liu-Sun

## The reduction (verified)
For **homogeneous** `f in C[x,y]` of degree `d`, on `R_+^2` write `(x,y)=(rs, r(1-s))`, `r>=0`, `s in [0,1]`,
`dx dy = r dr ds`, `e^{-x-y}=e^{-r}`, `f = r^d g(s)` with `g(s)=f(s,1-s)`. Then
> **`L(f^m) = int int f^m e^{-x-y} = (int_0^1 g(s)^m ds)(int_0^infty r^{dm+1}e^{-r}dr) = (dm+1)! * M(g^m)`**,
> `M(g^m) = int_0^1 g^m ds` the `[0,1]`-Lebesgue moment.
Checked to 25 dps: `f=x-y` (`d=1`, `g=2s-1`) gives `L(f^m) = (m+1)! M(g^m)`: `2, 24` at `m=2,4` (odd moments 0).

## Consequence: FC(2)-homogeneous = Liu-Sun
Since `(dm+1)! != 0`, `L(f^m)=0 for all m` iff `M(g^m)=0 for all m`, which by the classical moment problem
(polynomials dense in `L^2[0,1]`; moments determine `g`) forces `g=0`, hence `f=0`. **This is exactly
Liu-Sun 2020 Thm 2.6** -- the homogeneous factorial conjecture, radially reduced to Lebesgue moments. (Cite it;
don't reprove.)

## Placing the int e^Q claim (complements HYP-9078)
The external claim's integral is the **generating sum** of these very moments:
> `int_0^1 e^{g(s)} ds = sum_{m>=0} M(g^m)/m!`  -- weight `1/m!`, versus FC's demand that **each** `M(g^m)=0`.
So `int e^Q != 0` is **strictly weaker** than FC in the homogeneous case (HYP-9078's `j!` vs `1/(j+1)` weight
gap, seen here as "each moment" vs "one weighted sum"). Therefore:
> **`int e^Q != 0` cannot imply FC(2) through the radial reduction.** The asserted implication, if real, must run
> through the **non-homogeneous** case (where no clean `r`/`s` separation exists) via the friend-cited
> **horizontal-endomorphism rigidity** -- a counterexample `f` is forced into a fibred shape whose fibre integral
> is a 1-variable `int e^Q`. That rigidity is the unverified load-bearing step (HYP-9078: "no argument supplied").

## Net
- Homogeneous FC(2) = Lebesgue moments = Liu-Sun (grounded, verified).
- `int e^Q` is the generating sum -- weaker than FC homogeneously; the `=> FC(2)` route is necessarily
  non-homogeneous + rigidity, which is the actual open theorem (not the moment reduction).
- Consistent with kps-S165's Conway-Jones reading: the integral's only zero is the `2 pi i` (transcendental)
  resonance; algebraic `g` gives the discrete/leaking version, never the full one.

Files: verify inline. Complements HYP-9078; builds on kps-S165/S166.
