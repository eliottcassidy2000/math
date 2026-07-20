---
id: THM-1690
title: "THE s-FORMULA, THE TUNED CANCELLATION LOCUS, AND THE FINISHING STATEMENT for the toral half of GMC(2). (1) COEFFICIENT-DEGREE FORMULA, verified on 15 cells: the order-D toral recurrence for a_m = [u^{Mm}]Rᵐ has polynomial coefficients of degree exactly s(M,N) = C(M+N,2) − gcd(M,N) + 1. My earlier min(M,N) guess is REFUTED by the decisive cell (4,6): gcd=2, min=4, and the measured s = 44 = C(10,2)−2+1, not 42. The gcd is structural — gcd(M,N) > 1 is exactly when the cyclic u ↦ ζu symmetry underlying opus THM-1625's symmetric descent is available, and each unit of gcd shaves one degree. (2) THE TUNED CANCELLATION LOCUS is the saddle-value collision locus, and it is NOT the discriminant of R. Saddles = roots of Q(u) = uR'(u) − MR(u) = Σ_j (j−M) r_j u^j (degree D, and independent of r_M — the charge-0 coefficient drops out); their values v_i = R(u_i)/u_i^M control a_m ~ Σ c_i v_iᵐ. opus's Vandermonde needs the DOMINANT v_i distinct; the residual is where they COLLIDE. Explicit: at (1,1) the collision locus is {g₀g₂ = 0} (= one-sidedness), while {disc R = 0} is the DIFFERENT locus {a saddle value = 0}; at (1,2) with R = 1+bu+cu²+u³ the collision locus is 16(c−3)³(c²+3c+9)³ = 0 — depending ONLY on c, never on b, because the saddle geometry is r_M-independent. (3) THE FINISHING STATEMENT, assembled with mac-mini THM-1645: GMC(2) = (angular CT_u) ∘ (radial Laplace L), the ANGULAR half is TNC = Duistermaat–van der Kallen = PROVED, so the ONLY remaining gap is the radial 'pointwise-nonzero ⟹ integrated-nonzero' step, blocked by ker(L) ≠ 0 (L(t−1) = 1!−0! = 0). My orthogonal-polynomial route (THM-1660/1620) BYPASSES that kernel exactly where the composite m·E_r[ψ_m] is a classical orthogonal polynomial (Hermite/Legendre, no common root) — which by THM-1670 is precisely the order-2 case (M,N)=(1,1). So the finishing statement is: GMC(2) is toral-closed and reduces to one radial Laplace-determinacy question, elementary-closed on the (1,1) low-charge stratum and open beyond, matching klein's 'sign-indefinite' and mac-mini's 'Laplace determinacy' from a third direction"
status: >
  (1) VERIFIED on 15 cells (all (M,N) with D = 2..7, plus (4,6),(3,6),(4,8)), two primes
  agreeing, >= 22 holdout rows.  The decisive discriminator is (4,6) (gcd != min).  The
  formula is an EMPIRICAL fit with a structural reading of the gcd; not proved.  It
  supersedes the min(M,N) guess of THM-1670's named-next, which was wrong.
  (2) VERIFIED/EXACT: the saddle polynomial Q = Σ(j−M)r_j u^j and its r_M-independence are
  algebraic identities; the (1,1) and (1,2) collision loci are exact symbolic computations.
  The distinction {values collide} != {disc R = 0} is proved at (1,1) by hand.
  (3) ASSEMBLY, not a new proof.  It composes mac-mini THM-1645 (TNC = DvdK, angular half
  done; radial Laplace-determinacy is the gap — all mac-mini's, cited), opus THM-1625
  (Vandermonde / collision residual), my THM-1660/1620 (the orthogonal bypass) and THM-1670
  (order = D, so the bypass is order-2 = (1,1) only).  GMC(2) is NOT proved; the honest
  remaining gap is named three ways and they agree.
source: kind-pasteur-2026-07-20-S128c123 (owner: work the s-formula, the tuned cancellation locus, and the finishing statement)
depends_on:
  - THM-1670    # order = D; whose named-next (the s-formula) this answers
  - THM-1645    # mac-mini: GMC(2) = angular o radial, angular = DvdK, radial gap
related: [THM-1620, THM-1660, THM-1625, THM-1630, THM-1640]
script: 04-computation/s_formula_and_cancellation_locus_kps_S128c123.py (+ .out)
---

# THM-1690 — the s-formula, the tuned cancellation locus, the finishing statement

Three things the owner asked for, in order.

## 1. The coefficient-degree formula

The order-`D` toral recurrence for `a_m = [u^{Mm}] R^m` (`D = M+N`) has polynomial
coefficients of degree exactly

> **`s(M,N) = C(M+N, 2) − gcd(M,N) + 1`.**

Verified on 15 cells (all `(M,N)` with `D = 2..7`, plus `(4,6), (3,6), (4,8)`), two primes,
`≥ 22` holdout. THM-1670's `named-next` guessed `min(M,N)` in place of `gcd`; that is
**refuted by `(4,6)`**: `gcd = 2`, `min = 4`, and the measured minimal degree is
`s = 44 = C(10,2) − 2 + 1`, not `42`.

| `D` | cells `(M,N): s` |
|---|---|
| 2 | (1,1): 1 |
| 3 | (1,2): 3 |
| 4 | (1,3): 6, (2,2): **5** |
| 5 | (1,4): 10, (2,3): 10 |
| 6 | (1,5): 15, (2,4): **14**, (3,3): **13** |
| 7 | (1,6): 21, (2,5): 21, (3,4): 21 |

The `M = 1` column is `C(D,2)` exactly (`gcd = 1`). The defect `C(D,2) − s` is `gcd(M,N) − 1`,
nonzero only at `(2,2), (2,4), (3,3), (4,6), (4,8), (3,6)` — exactly the `gcd > 1` cells.

**Why `gcd`.** `gcd(M,N) = g > 1` is exactly when the cyclic symmetry `u ↦ ζ u` (`ζ^g = 1`)
that drives opus THM-1625(3)'s symmetric descent is available, and each unit of `g` removes
one degree of freedom from the recurrence's coefficients. So the `gcd` in the size formula
and the `gcd` in the descent are the same `gcd`.

## 2. The tuned cancellation locus

The growth of `a_m` is `a_m ~ Σ_i c_i v_i^m`, summed over the **saddles** — the critical
points of `g(u) = R(u)/u^M`, i.e. the roots of

> `Q(u) = u R'(u) − M R(u) = Σ_j (j − M) r_j u^j`  (degree `D`).

Two facts fall out of that form: `Q` has **no `u^M` term** (the `j = M` coefficient is
`(M−M)r_M = 0`), so **the saddle geometry is independent of `r_M`** — the charge-0
coefficient of `R`, the "`b`" of the radial story. The **tuned cancellation locus** is where
the dominant saddle **values** `v_i = R(u_i)/u_i^M` collide — this is where opus's Vandermonde
argument (distinct dominant values `⟹` TNC) breaks and one is in the residual.

> **It is NOT the discriminant of `R`.** At `(1,1)`, `R = g₀+g₁u+g₂u²`: saddles at
> `u = ±√(g₀/g₂)`, values `v_± = g₁ ± 2√(g₀g₂)`. Then
> - `{v_+ = v_-} = {g₀g₂ = 0}` — the collision locus, which here *is* one-sidedness;
> - `{v_+ v_- = 0} = {g₁²−4g₀g₂ = 0} = {disc R = 0}` — a **different** locus, where a saddle
>   *value* is zero (`R` has a double root); this is the `P_0(0)` locus of THM-1670.
>
> At `(1,2)`, `R = 1+bu+cu²+u³`: the collision locus is `16(c−3)³(c²+3c+9)³ = 0` — it
> **depends only on `c`**, never on `b`, by the `r_M`-independence above.

So there are two distinct "cancellation" subvarieties, and the Vandermonde obstruction is the
value-collision one, not the discriminant. `c = 3` is the real branch (`Q = 2u³+3u²−1 =
(u+1)²(2u−1)` has a double saddle at `u = −1`).

## 3. The finishing statement

Composing with **mac-mini THM-1645**, which factors `GMC(2) = (angular CT_u) ∘ (radial
Laplace L)`:

**The toral/angular half is finished.** `CT_u[Λ_s^m]` is exactly the toral sequence of this
note; `TNC = Duistermaat–van der Kallen` is proved (THM-1630), and by THM-1645 it applies
uniformly in the radius `s`. Structurally (this note + THM-1670) that sequence is
order-`D` holonomic with coefficient-degree `C(D,2) − gcd + 1`, and its estimate-free
(Vandermonde) proof works off the tuned cancellation locus, the residual being asymmetric
value-collisions — verified non-fatal (opus THM-1625, and `u⁴−2u²−2` gives
`−2,0,16,−56,…`, not identically zero).

**The one remaining gap is radial, and it is Laplace determinacy.** By THM-1645 the whole
outstanding content of GMC(2) is the step *"`CT_u[Λ_s^m] ≠ 0` pointwise in `s` `⟹`
`∫ CT_u[Λ_s^m] e^{−s} ds ≠ 0`"*, blocked because `ker(L) ≠ 0` (`L(t−1) = 1! − 0! = 0`).

**The orthogonal-polynomial route dissolves that kernel exactly where the composite is a
classical orthogonal family.** My THM-1660/1620 never passes through "pointwise `⟹`
integrated": it shows the *integrated* quantity directly,
`m · E_r[ψ_m] = s^m He_m(b/s)` (Hermite), and Hermite polynomials have no common root. By
THM-1670 the composite is a single orthogonal family precisely in the order-2 case,
`(M,N) = (1,1)` — the `{−1,0,1}` low-charge stratum. death-star-S63/S64 push this outward
with the Legendre recurrence for a linear/higher charge-0 coefficient.

> **Finishing statement.** *GMC(2) is toral-closed: the angular half is TNC = DvdK, and
> what remains is a single radial Laplace-determinacy question — `does integrated
> non-vanishing survive `ker(L)`?` — elementary-closed (via the orthogonal-polynomial
> bypass) on the `(M,N) = (1,1)` low-charge stratum, and open beyond it.* This is the same
> remaining gap that klein-S363 reached as "sign-indefinite, positivity unavailable" and
> mac-mini reached as "Laplace determinacy, not tori" — now located a third way, from the
> toral recurrence's structure, and the three descriptions agree.

**Not claimed:** GMC(2). The radial gap is real and open; this maps it, closes its toral
half, and shows the orthogonal route is the tool that removes it where removable.

## Named next

- The radial gap in one variable: `L(g) = ∫₀^∞ g(s) e^{−s} ds` with `ker(L)` spanned by the
  Laguerre-orthogonality relations — the composite's determinacy is a **moment-problem**
  question on `[0,∞)`, and the right object is likely the `Laguerre`/`Charlier` side of the
  same Askey scheme whose `Hermite`/`Legendre` side is the toral layer. If the radial
  functional is Laguerre-diagonal, its kernel is explicit and the descent is a second
  no-common-root.
- Prove the `s`-formula. The `gcd` defect strongly suggests a `Riemann–Hurwitz` count for the
  branched cover `z^M = t R(z)`: the coefficient-degree is the number of finite singularities
  of the order-`D` ODE, and the `gcd` symmetry merges `gcd − 1` of them.
