---
id: HYP-9361
title: "Open laws for the Camion function (THM-4750). (a) Split/inert equidistribution: over primes p, a non-triangle zero class (class 1, class 2 or sporadic) has 2 second-order roots (split) or 0 (inert) with density 1/2 each, independently of the class and of the other classes at the same prime. (b) The partner root of the triangle class (which is split whenever its quadratic is genuine and the root l = 0 is simple, THM-4750 (4d)) has no closed form; in particular it is equidistributed mod p. (c) Tower heights: h(p) <= 2 for all odd p, with h(p) = 2 exactly for the Camion-Wieferich primes (3, 13, 22409 below 1e5), heuristically infinitely many but only finitely many with h >= 3. (d) Degree = uniformity: for labelled structures built from 2^C(n,r) by the Session25 operations (e.g. connected r-uniform hypergraphs), the class functions interpolate at Mahler speed 1/r, are polynomials of degree r in the class index mod p^2 with leading coefficient proportional to the Fermat quotient q_2(p), and the second-order roots of a zero class follow the root statistics of a random degree-r polynomial over F_p. (e) The 2-adic zero n* of the strong groupoid count on n = 2 mod 3 is irrational."
status: >
  OPEN.
  FINITE-EXACT evidence:
  - (a): 1226 primes 7 <= p <= 10^4 give class 1 at 633 inert / 591 split, class 2 at 612 / 609, and the sporadic
    classes at 618 / 605 (second_order_*.out).
  - (b): split in all 1224 non-Wieferich cases with p <= 10^4 (as THM-4750 (4d) forces); linear at 1093 and 3511.
  - (c): Camion-Wieferich primes p <= 10^5 are 3, 13, 22409, all with h = 2 (THM-4750 (5)).
  - (d), r = 3: the speed and the cubic law are now PROVED (THM-4750 (7)). The random-cubic root statistics hold on
    60 primes <= 2200 (176 zero classes), with a forced root in class 2.
  - (e): no rational with numerator and denominator below about 2^31 matches n* mod 2^64.
  OBSERVATIONS, not claims:
  - a deficit of double roots (1 observed, about 5.5 expected; about 3% significance before multiple comparisons);
  - the sporadic-class count per prime is slightly underdispersed (variance 0.875 against the Poisson value 1;
    chi-square about 10.7 on 5 degrees of freedom).
session: mac-mini-2026-10-21-fanopsp
source: 01-canon/theorems/THM-4750-the-camion-function-strong-tournament-counts-interpolate-p-adically-at-half-speed-fermat-quotient-curvature-wieferich-collapse-second-order-zeros-in-pairs-and-the-camion-limit.md
related:
  - 05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-wieferich-type-primes-sporadic-zero-classes-prime-squares.md
scripts:
  - 04-computation/experiments/camion_wieferich_20261021/camion_second_order.c (+ second_order_*.out)
  - 04-computation/experiments/camion_wieferich_20261021/hyper3_second_order.c (+ hyper3_second_order_*.out)
  - 04-computation/experiments/camion_wieferich_20261021/two_adic_zero_fast.py (+ .out)
---

# HYP-9361 — open laws for the Camion function

## (a) Split/inert equidistribution

In a zero class, the second-order roots are the roots in `F_p` of the quadratic `g_m` (THM-4750 (4)). If its discriminant behaved like a random residue, the class would be split or inert with probability `1/2` each, and would have a double root with probability about `1/p`.

| class | inert (0 roots) | split (2) | ratio |
|---|---|---|---|
| 1 | 633 | 591 | 0.483 split |
| 2 | 612 | 609 | 0.499 |
| sporadic | 618 | 605 | 0.495 |

Classes 1 and 2 at the same prime look independent: among the 1221 primes where both have 0 or 2 roots, the (class 1, class 2) counts are (0,0) 320, (0,2) 312, (2,0) 292, (2,2) 297, with χ² = 0.14 on 1 degree of freedom. The split rate shows no dependence on `ε = (2/p)`.

Open question: is there a Chebotarev-type mechanism? The discriminant is a polynomial in `D(m)/p`, `a_1/p` and `q_2(p)X(m−2)` mod `p`; no algebraic structure is known for it.

## (b) The triangle class

`D(3) = 0` forces `l = 0` to be a root of `g_3` (base `m = 3`). So the triangle class is split whenever `g_3` is a genuine quadratic with a simple root; this part is PROVED (THM-4750 (4d)). The open part concerns the partner root `l' = 1 − 2c_1/c_2` (mod `p`), which gives the integers `n = 3 + l'(p−1)` with `p² | D(n)`. Is there a closed form for `l'`, and is it equidistributed mod `p`?

## (c) Tower heights

By THM-4750 (5), `h(p) = v_p(Λ_p)`. When `p` is Camion–Wieferich and the root is simple, `h(p) ≥ 3` iff the `p`-adic zero `ℓ` of the class-1 Camion function satisfies `ℓ ≡ p (mod p²)`.

* The heuristic is a second independent `1/p` event, so `Σ 1/p² < ∞`: finitely many `p` should have `h ≥ 3`, with expected total about `0.09` for `p ≥ 5`.
* The Camion–Wieferich primes themselves should be infinite, with density zero (`Σ 1/p = ∞`; HYP-9358 (a)).
* Up to `10^5` they are `3, 13, 22409`, with expected count about `1.87 + 1`.

## (d) Degree = uniformity

THM-4750 (1) gives Mahler speed `1/2` because the grade-`j` term of `(1+pu)^C(n,2)` has degree `2j` in `n`. With `C(n,r)` the degree is `rj`. The predictions are:

1. Mahler speed `1/r`.
2. Mod `p²`, the class sequences are polynomials of degree `r` in `l`. Their curvature coefficients are proportional to `u ≡ −q_2(p)`, so they collapse to linear at Wieferich primes.
3. The second-order roots of a zero class follow the statistics of a random polynomial of degree `r`. For `r = 3`: 0, 1, 3 roots with frequencies `1/3, 1/2, 1/6`. For tournaments (`r = 2`): 0 or 2 roots, `1/2` each.
4. Universal zero classes with an exact zero inside the polynomial range carry a forced root. Examples are the triangle at `r = 2`, and possibly `n = 2` at `r = 3`.

Test case: connected labelled 3-uniform hypergraphs `b(n)`, from `2^C(n,3) = Σ_k C(n−1,k−1) b(k) 2^C(n−k,3)`.

* At primes with `gcd(ord_p(2), 6) = 1`, the Teichmüller part `ω^C(n,3)` has period dividing `p − 1`. The tournament proof then applies with multiplier 1, so the quasi-period is `p − 1`.
* At other primes the Fourier support of `ω^C(n,3)` may need `F_(p^f)` with `f > 1`. The quasi-period can then exceed `p − 1`.
* **Now PROVED** (THM-4750 (7)) for the speed `1/r` and the degree-`r` law mod `p²` at primes with `gcd(ord_p(2), r!) = 1`.
* **Root statistics** (FINITE-EXACT, `r = 3`, 60 primes `≤ 2200`):
  * class 1: 22, 25, 1, 12 for 0–3 roots, against the random-cubic 20, 30, 0, 10;
  * class 2: 0, 31, 0, 29, with a forced root at `n = 2`;
  * sporadic classes: 21, 21, 0, 14, against 18.7, 28, 0, 9.3.
* **What remains open in (d):**
  * the statistical law (random degree-`r` root counts) itself;
  * the primes where `gcd(ord_p(2), r!) > 1` (longer periods, nodes in `F_(p^f)`).
* At the Wieferich primes the collapse to linear is proved for every `r`. But 1093 and 3511 have `ord_p(2)` divisible by 2 or 3 (364 and 1755), so for `r = 3` the clean setting does not apply to them.

## (e) The 2-adic zero

`n* ≡ 16221783116972227762 (mod 2^64)` (THM-4750 (6c)) is not a small rational. The heuristic in THM-4750 suggests `n* = log(−cζ/σ(cζ)) / log(w/σ(w))`. Here `ρ = ζw` is the root of the tournament EGF near a primitive cube root of unity, and `c` is its partial-fraction coefficient. Two problems are open:

* Is `n*` transcendental?
* Is there a combinatorial meaning for a strong "groupoid count of order `n*`"?

## What would refute

* (a): a split/inert ratio for some class type converging away from `1/2`.
* (b): a closed form for `l'`, or a bias of `l'/p` away from uniform.
* (c): a prime with `h(p) ≥ 3`. This would not refute the heuristic, only the strong form `h ≤ 2`.
* (d): a degree-`r` structure whose mod-`p²` class sequences are not polynomials of degree `r`, at primes where the period is `p − 1`.
