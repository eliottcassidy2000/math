---
id: HYP-3730
title: THE CHEEKY CONNECTION IS HEEGNER / CLASS NUMBER 1 -- Phi_6(x)=x^2-x+1 is the c=1 RABINOWITZ prime-generating quadratic (disc -3 = Eisenstein = the FIRST Heegner number); its Sylvester-iterated apex tower 2,3,7,43 ARE Heegner numbers (class number 1); the prime-generating quadratics x^2-x+c for c=1,2,11,41 have Heegner discs -3,-7,-43,-163 with apexes 4c-1 = 3,7,43,163 (the prime Heegner numbers, ending at 163 = Ramanujan's e^(pi sqrt163)); the apex 7 = Heegner -7 = the conductor of X_0(14); Sylvester (3,7,43,1807) and the Heegner tower (3,7,43,163) agree to 43 then DIVERGE. Also: the floor-margin sum sum_n 1/Phi_6(n) = (pi/sqrt3)tanh(pi sqrt3/2) -- the same sqrt3 as Kershner's hexagonal 2pi/sqrt27; and the Euclidean climb margin M-1/n = (n-1)/(n Phi_6(n)) telescopes over the semiconvergent rungs. The LRC covering floor lives in the class-number-1 / Heegner / Ramanujan world
status: VERIFIED as NUMBER THEORY (Sylvester INTERSECT Heegner = {2,3,7,43}; x^2-x+c prime-generating <=> 1-4c Heegner (Rabinowitz), c=1,2,11,41 -> -3,-7,-43,-163, apex 4c-1=3,7,43,163; sum 1/Phi_6(n)=(pi/sqrt3)tanh(pi sqrt3/2) match to 1e-4). ** CAVEAT (MISTAKE-087, mac-mini-S47): the construction {1,..,n-2,n(n-1)} with M=n/Phi_6 is NOT the covering-min for n>=7 (beaten by spread coverings: n=7 {1,2,5,6,7,8} M=2/13<7/43). So this Heegner/Eisenstein/Sylvester structure is REAL number theory about Phi_6 and that PARTICULAR construction -- which is BEAUTIFUL BUT NON-EXTREMAL. The 'LRC covering floor lives in class-number-1 land' framing is RETRACTED to 'a particular (non-minimal) covering has class-number-1 structure.' ** Extends HYP-2226. The LRC floor M>=1/n is untouched.
source: klein-2026-06-29-S33
depends_on:
  - HYP-3724   # Phi_6 = Sylvester map; the apex tower
  - HYP-3705   # Phi_6 = Eisenstein norm (Q(sqrt-3))
related:
  - HYP-2226   # Heegner prime polynomials, tournament n-2 witness horizon
  - HYP-3718   # the Euclidean climb / the margin
  - HYP-3586   # X_0(14) cusps / genus (apex 7 = Heegner -7)
  - HYP-3716   # Kershner hexagonal 2pi/sqrt27 (the sqrt3)
results:
  - 04-computation/torus_lift_threegap_klein.py
---

# HYP-3730 — the cheeky connection: Heegner / class number 1 / Ramanujan

> **READ FIRST -- MISTAKE-087 (mac-mini-S47, same day):** the premise that `n/Phi_6` is the LRC covering-min
> for `n>=7` is REFUTED (spread coverings beat it: `n=7` `{1,2,5,6,7,8}` `M=2/13 < 7/43`). Everything below is
> TRUE NUMBER THEORY about `Phi_6` and the construction `{1,..,n-2,n(n-1)}`, but that construction is NOT the
> extremal covering. So read this as "a beautiful covering carries class-number-1 structure," NOT "the
> covering-min is Heegner." The `Phi_6 = Eisenstein/Sylvester/Heegner/2T+1` identities stand; their role in
> the covering OPTIMIZATION is retracted. (Lesson: do not conflate an elegant covering with the minimal one.)

## Phi_6 is the c=1 prime-generating quadratic (disc -3, the first Heegner)
`Phi_6(x) = x^2 - x + 1` is the `c = 1` member of the **Rabinowitz prime-generating quadratics** `x^2 - x +
c`. Rabinowitz's theorem: `x^2 - x + c` is prime for `x = 1..c-1` IFF the discriminant `1 - 4c` is a Heegner
number (`Q(sqrt(1-4c))` has class number 1). For `Phi_6` (`c=1`), the discriminant is `-3` -- the **first
Heegner number**, the **Eisenstein field `Q(sqrt-3)`** (HYP-3705: `Phi_6` = the Eisenstein norm).

## The apex tower ARE Heegner numbers
The Sylvester-iterated apex tower (HYP-3724) and the Heegner numbers:
```
Sylvester (Phi_6 iterated):   2, 3, 7, 43, 1807, ...
Heegner numbers (h=1):        1, 2, 3, 7, 11, 19, 43, 67, 163
INTERSECTION:                 2, 3, 7, 43         <-- the apex tower ARE Heegner numbers (class number 1)
```
The prime-generating quadratics `x^2 - x + c` for `c = 1, 2, 11, 41` have Heegner discriminants
`-3, -7, -43, -163`, with apex values `4c - 1 = 3, 7, 43, 163` -- the **prime Heegner numbers** (`= 3 mod
4`), ending at `163 = Ramanujan's constant` (`e^(pi sqrt163)` almost integer). So the "true" apex tower is
the **Heegner tower `3, 7, 43, 163`**; Sylvester (`Phi_6`-iteration) is a coincidental shadow that AGREES to
`43` then diverges (`Phi_6(43) = 1807 = 13.139` composite, vs Heegner `163`). `43` is the last common rung.

## The apex 7 = Heegner -7 = X_0(14)
The LRC(14) apex `7 = Heegner -7 = Q(sqrt-7)`, class number 1 -- and `14 = 2.7` is the conductor of the
modular curve `X_0(14)` (the project's genus-1 curve, HYP-3586). The genus of `X_0(2p)` jumps to 1 at `p=7`
(`Q(sqrt-7)`). So the apex is doubly "Heegner": a Heegner NUMBER (class number 1) AND the CM/Heegner-POINT
world of `X_0(14)`. This extends HYP-2226 (Heegner prime polynomials have a tournament `n-2` witness horizon).

## The sqrt3 in the floor-margin sum = the sqrt3 of Kershner
The covering-floor margins sum to a hexagonal constant:
> `sum_{n>=1} 1/Phi_6(n) = sum 1/(n^2-n+1) = (pi/sqrt3) tanh(pi sqrt3/2) = 1.798147...`
The `sqrt3 = sqrt|disc Phi_6|` (the Eisenstein/`Q(sqrt-3)`) is the SAME `sqrt3` as Kershner's thinnest
hexagonal covering density `2pi/sqrt27 = 2pi/(3 sqrt3)` (HYP-3716). So the discriminant of `Phi_6` (`-3`) ties
the floor-margin sum, the hexagonal covering, and the Eisenstein field into one `sqrt3` -- circle (`pi`) over
hexagon (`sqrt3`).

## The class-number-1 near-miss theme (the cheekiest read)
The LRC covering floor `M = n/Phi_6(n)` is a NEAR-MISS of `1/n` (margin `(n-1)/(n Phi_6(n)) -> 0`); the
Heegner number `163` produces the famous NEAR-MISS `e^(pi sqrt163) ~ integer`. Both are **class-number-1
near-misses**: the simplest imaginary quadratic fields (`h=1`) force a value to come *just shy* of a clean
target -- the lonely runner just above `1/n`, Ramanujan's constant just below an integer. The covering
floor's "uniform looseness" is the class-number-1 near-miss, in the covering register.

## Net (the lurking connection, made explicit)
`Phi_6 = x^2-x+1` wears five hats, now unified: the **Eisenstein norm** (`Q(sqrt-3)`), the **Sylvester map**,
**`2T+1`** (the triangle), the **c=1 Rabinowitz prime-generating quadratic**, and -- the cheeky one -- the
gateway to the **Heegner numbers / class number 1**: the apex tower `2,3,7,43` are Heegner, the apex `7` is
`Q(sqrt-7)` = `X_0(14)`, and the tower's true continuation `...43, 163` ends at Ramanujan's constant. The
LRC covering floor lives in the class-number-1 world; its `sqrt3` is Kershner's hexagon; its near-miss is
Ramanujan's. Verified/striking; the deep causality (`Q(sqrt-3)` + the genus of `X_0(2p)`) is the program.
