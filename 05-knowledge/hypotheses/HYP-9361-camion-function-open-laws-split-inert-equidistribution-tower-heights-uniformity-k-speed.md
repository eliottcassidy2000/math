---
id: HYP-9361
title: "Open laws for the Camion function (THM-4750). (a) Split/inert equidistribution: over primes p, a non-triangle zero class (class 1, class 2 or sporadic) has 2 second-order roots (split) or 0 (inert) with density 1/2 each, independently of the class and of the other classes at the same prime. (b) The partner root of the triangle class has no closed form and is equidistributed mod p. (c) Tower heights: h(p) <= 2 for all odd p, with h(p) = 2 exactly for the Camion-Wieferich primes (3, 13, 22409 below 1e5); heuristically there are infinitely many Camion-Wieferich primes but only finitely many p with h >= 3; every Hurwitz-closed structure has its own sparse list. (d) The S_r law: for structures built from 2^C(n,r) (connected r-uniform hypergraphs at primes with gcd(ord_p 2, r!) = 1), the second-order Camion zeros of a zero class are distributed like the roots of a random degree-r polynomial over F_p (root count = fixed points of a uniform permutation in S_r), with forced rational roots at exact zeros inside the analytic range. (e) The 2-adic zero n* of the strong groupoid count on n = 2 mod 3 is irrational (indeed transcendental). (f) The Camion residue rho_C(a) and the Camion lift kappa have no closed form."
status: >
  OPEN.
  AUDITED 2026-10-21 with THM-4750 (MISTAKE-654): the (b) wording ("as (4d) forces") and the (e) evidence
  (rational reconstruction) are corrected below.
  FINITE-EXACT evidence:
  - (a): 1226 primes 7 <= p <= 10^4 give class 1 at 633 inert / 591 split, class 2 at 612 / 609, and the sporadic
    classes at 618 / 605. Classes 1 and 2 are independent (chi^2 = 0.14), with no epsilon dependence. Reproduced by
    audit b4.
  - (b): the triangle class is split in all 1224 non-Wieferich cases with 7 <= p <= 10^4 (no double root occurs
    there). At p = 5 it has a double root and a ramified pair (THM-4750 (4e)).
  - (c): Camion-Wieferich primes p <= 10^5 are 3, 13, 22409, all with h = 2 (two independent searches), extended to 10^6
    in cw_flint_*.out. Base 3: 17, 31 (p <= 10^4). Connected graphs: 5, 197 (p <= 4000).
  - (d): r = 2, 3, 4 (THM-4750 §4, §6b): 1226 primes for r = 2, and 60 primes <= 2200 each for r = 3, 4. The speed and
    the degree-r polynomial law are PROVED (THM-4750 (7)); only the statistics are conjectural.
  - (e): n* mod 2^140 is known (audit b7). There is no rational of height <= 2^69.5, and the closed form
    log t / log r holds to 90 bits.
  OBSERVATIONS, not claims:
  - a deficit of double roots (1 observed, about 5.5 expected, P(<= 1) = 0.027, before multiple comparisons);
  - the sporadic-class count per prime is underdispersed against the binomial model of HYP-9358 (b) (variance 0.874,
    dispersion z = -2.47).
session: mac-mini-2026-10-21-fanopsp
source: 01-canon/theorems/THM-4750-the-camion-function-strong-tournament-counts-interpolate-p-adically-at-half-speed-fermat-quotient-curvature-wieferich-collapse-second-order-zeros-in-pairs-and-the-camion-limit.md
related:
  - 05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-wieferich-type-primes-sporadic-zero-classes-prime-squares.md
  - 05-knowledge/results/camion_function_audit_20261021.md
scripts:
  - 04-computation/experiments/camion_wieferich_20261021/camion_second_order.c (+ second_order_*.out)
  - 04-computation/experiments/camion_wieferich_20261021/hyper3_second_order.c, hyper_r_second_order.c, hyper_r_tabulate.py (+ hyper3_second_order_*.out, hyper4_second_order_5_2200.out)
  - 04-computation/experiments/camion_wieferich_20261021/cw_flint_np.py, cw_base.c, cw_connected.c (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/two_adic_zero_fast.py, two_adic_regulator.py (+ .out); audit/b7, b9
  - 04-computation/experiments/camion_wieferich_20261021/camion_residue*.{c,py}, camion_lift*.py (+ .out)
---

# HYP-9361 — open laws for the Camion function

## (a) Split/inert equidistribution

In a zero class, the second-order roots are the roots in `F_p` of the quadratic `g_m` (THM-4750 (4)). If its discriminant behaved like a random residue, the class would be split or inert with probability `1/2` each. It would have a double root (congruent zeros, possibly ramified) with probability about `1/p`.

| class | inert (0 roots) | split (2) | ratio |
|---|---|---|---|
| 1 | 633 | 591 | 0.483 split |
| 2 | 612 | 609 | 0.499 |
| sporadic | 618 | 605 | 0.495 |

* **Classes 1 and 2 look independent.** Among the 1221 primes where both have 0 or 2 roots, the (class 1, class 2) counts are (0,0) 320, (0,2) 312, (2,0) 292, (2,2) 297, with `χ² = 0.14` on 1 degree of freedom.
* **There is no dependence on `ε = (2/p)`.** The split rate is 0.483 at `ε = +1` and 0.501 at `ε = −1`.
* **Open.** Is there a Chebotarev-type mechanism? The discriminant is a polynomial in `D(m)/p`, `a_1/p` and `q_2(p)X(m−2)` mod `p`, and no algebraic structure is known for it.

## (b) The triangle class

`D(3) = 0` forces `l = 0` to be a root of `g_3` (base `m = 3`). Also `X(1) = 4`, so `g_3` is a genuine quadratic at every non-Wieferich prime. Hence the triangle class is split whenever the root `l = 0` is simple. This part is PROVED (THM-4750 (4d)).

When `l = 0` is a double root, the pair can be ramified, as at `p = 5`. No such case occurs for `7 ≤ p ≤ 10^4`.

The open part concerns the partner root `l' = 1 − 2c_1/c_2 (mod p)`. It gives the integers `n = 3 + l'(p−1)` with `p² | D(n)`.

* Does `l'` have a closed form?
* Is `l'/p` equidistributed?

## (c) Tower heights and Camion–Wieferich lists

By THM-4750 (5) and (5′), `h(p) = v_p(Λ_p) = v_p(2 − κ(2))`, where `κ(2)` is the Camion lift.

* **When `h ≥ 3`.** For a Camion–Wieferich prime with a simple root, `h(p) ≥ 3` iff the class-1 Camion zero satisfies `ℓ ≡ p (mod p²)`, a second `1/p` event.
* **Counting heuristic.** `Σ 1/p² < ∞`, so finitely many `p` should have `h ≥ 3` (expected total about `0.09` for `p ≥ 5`). The Camion–Wieferich primes themselves should be infinite with density zero, since `Σ 1/p = ∞` (HYP-9358 (a)).

Each Hurwitz-closed structure has its own lift, and hence its own sparse list:

| structure | condition | primes found | range |
|---|---|---|---|
| strong tournaments (base 2) | `p² ∣ s(p) − 2^((p−1)/2)` | 3, 13, 22409 | `≤ 10^5`; `10^6` in `cw_flint_*.out` |
| strong semicomplete digraphs (base 3) | `p² ∣ s_3(p) − 3^((p−1)/2)` | 17, 31 | `≤ 10^4` |
| connected labelled graphs | `p² ∣ cg(p) − 2^((p−1)/2) + 1` | 5, 197 | `≤ 4000` |
| Wieferich, base 2 / base 3 (classical) | `a^(p−1) ≡ 1 (mod p²)` | 1093, 3511 / 11, 1006003 | classical |

No prime appears in two lists. The base-3 Wieferich prime 11 is not base-3 Camion–Wieferich.

NUMEROLOGY, with no mechanism proposed: the base-3 list begins with the Fermat prime `17` and the Mersenne prime `31`.

## (d) The `S_r` law

THM-4750 (7) proves Mahler speed `1/r` and the degree-`r` class polynomial mod `p²` at primes with `gcd(ord_p(2), r!) = 1`. The open claim is statistical. For a zero class, the `r` second-order zeros behave like the roots of a uniformly random degree-`r` polynomial over `F_p`: Frobenius acts on them like a uniform element of `S_r`. An exact zero inside the analytic range is a forced rational root, and `S_(r−1)` governs the rest.

| `r` | class | observed (0 / 1 / 2 / 3 / 4 roots) | prediction |
|---|---|---|---|
| 2 | 1, 2, sporadic | 1:1 split/inert (table (a)) | random quadratic |
| 3 | 1 | 22 / 25 / 1 / 12 | 20 / 30 / 0 / 10 |
| 3 | 2 (forced) | 0 / 31 / 0 / 29 | 0 / 30 / 0 / 30 |
| 3 | sporadic | 21 / 21 / 0 / 14 | 18.7 / 28 / 0 / 9.3 |
| 4 | 1 | 23 / 19 / 16 / 0 / 2 | 22.5 / 20 / 15 / 0 / 2.5 |
| 4 | 2 (forced) | 0 / 18 / 32 / 1 / 9 | 0 / 20 / 30 / 0 / 10 |
| 4 | 3 (forced) | 0 / 21 / 26 / 0 / 13 | 0 / 20 / 30 / 0 / 10 |
| 4 | sporadic | 21 / 19 / 13 / 0 / 3 | 21 / 18.7 / 14 / 0 / 2.3 |

* **Heuristic.** The coefficients `a_i/p mod p` are equidistributed. No proof is in sight; it would be a Chebotarev-type statement for a non-algebraic family.
* **The forced roots** (`r = 3`, class 2; `r = 4`, classes 2 and 3) sit at the exact zeros `b_r(2) = 0` and `b_4(3) = 0`. These lie below the proved threshold, but the class polynomials already hold there (audit b8, `hyper3_class2_check.out`). That is a FINITE-EXACT observation, not a consequence of THM-4750 (7).
* **Open:**
  * the primes with `gcd(ord_p(2), r!) > 1`, which have longer periods and nodes in `F_(p^f)`;
  * the threshold at which the class polynomials actually start, which is lower than proved.

## (e) The 2-adic zero

THM-4750 (6c) gives `n* = log t / log r` in the norm-one torus of `Q_2(√−3)`. Here `ρ = ζw` is the root of the tournament EGF near a primitive cube root of unity, `C = 1/(ρT′(ρ))`, `r = σ(w)/w` and `t = −ζσ(C)/C`.

* The closed form agrees with `n*` to 90 bits (audit b9).
* `n* ≡ 1090019148193376184792522869885760632192178 (mod 2^140)`.
* There is no rational of height `≤ 2^69.5`. The `2^64` reconstruction `765490658/1598160473` is a generic fraction and fails mod `2^80` (audit b7).

Open questions:

* Is `n*` transcendental? `r` and `t` are defined by the transcendental-looking entire function `T`, so no algebraicity is expected.
* Does a strong "groupoid count of order `n*`" have any combinatorial meaning?

## (f) The Camion residue and the Camion lift

`ρ_C(a) = (σ_p(a) − 1)/p mod p` is lift-independent (THM-4750 (5′)).

* As a function on `F_p^×` (`p ≤ 19`) it has no symmetry under `a ↦ 1/a` or `a ↦ −a`, and its Lagrange degree is generic.
* For `a = 2` and `p ≤ 2000` it satisfies no small-rational relation with `q_2`, `w_p` or `B_(p−3)`.
* `κ` is not multiplicative.

Open questions:

* Is there a closed form, for instance through Gauss sums of the quadratic phase `ω^C(n,2)`?
* Do the primes with `ρ_C(2) = 0` (199, 809, 929 below 2000), where the Camion lift and the Teichmüller lift agree mod `p²`, have any meaning?

## What would refute

* **(a)** A split/inert ratio for some class type that converges away from `1/2`.
* **(b)** A closed form for `l'`, or a bias of `l'/p` away from uniform.
* **(c)** A prime with `h(p) ≥ 3`. That would refute only the strong form `h ≤ 2`, not the heuristic.
* **(d)** Root statistics that converge to a non-`S_r` law at primes where `gcd(ord_p(2), r!) = 1`.
* **(e)** A rational or algebraic `n*` of small height.
