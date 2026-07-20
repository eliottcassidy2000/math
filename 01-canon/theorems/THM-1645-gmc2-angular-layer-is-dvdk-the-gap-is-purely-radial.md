---
id: THM-1645
title: "GMC(2) FACTORS THROUGH TNC, AND THE ANGULAR HALF IS ALREADY DONE. Via the polar bridge E[P^m] = ∫₀^∞ CT_u[Λ_s(u)^m] e^{−s} ds (Λ_s(u) = P(√s·u, √s/u), s = |Z|² ~ Exp(1)), GMC(2) is the composite of an ANGULAR functional CT_u and a RADIAL functional L(g)=∫g·e^{−s}ds. KEY NEW FACT: the charge (u-power) support of Λ_s is s-INDEPENDENT — a monomial Z^aW^b maps to s^{(a+b)/2}u^{a−b} and the positive coefficient s^{…} never kills a charge — so P is two-sided ⟺ Λ_s is two-sided for a.e. s. Hence the one-variable Duistermaat–van der Kallen theorem (= TNC, THM-1630, PROVED 1998) applies UNIFORMLY in s: if P is two-sided then CT_u[Λ_s^m] ≠ 0 for some m, for a.e. s. CONSEQUENCE: the entire remaining gap in the GMC(2) reduction (THM-1540's unwritten descent 'top-degree one-sided ⟹ P one-sided') is NOT a toral/angular-geometry gap — the angular direction is DvdK-closed — it is exactly the RADIAL 'pointwise-nonzero ⟹ integrated-nonzero' descent, blocked solely by ker(L) ≠ 0. The obstruction is exhibited: L(t−1) = 1!−0! = 0, so integrated vanishing does not force pointwise vanishing. So GMC(2) is obstructed by LAPLACE DETERMINACY, not by tori. Unconditional corollaries via the settled angular layer: one-sided P ⟹ MZ (E[QP^m] = 0 for m > deg_charge Q, verified), and two-monomial Z^p+W^q is never in the nullcone."
status: >
  The polar bridge is VERIFIED-EXACT against Wick expansion on four complex P, m = 1..6
  (it is also already canon as THM-1540(A)). The s-independence of charge support is
  PROVED (one line) and the uniform DvdK application is a corollary of THM-1630. The
  decoupling obstruction L(t−1)=0 is exact. The corollaries are VERIFIED exactly.
  THIS DOES NOT PROVE GMC(2). It RELOCATES the open gap from the angular layer (closed)
  to the radial layer (open, = HYP-8350 plus a cross-shell coupling). Stated as a
  structural reduction, not a resolution.
source: mac-mini-2026-07-20-S143 (owner: "work the GMC(2) through the TNC")
depends_on:
  - THM-1630  # TNC IS the one-variable DvdK theorem (proved 1998) -- the angular layer
  - THM-1540  # the polar bridge and the nullcone-structure reduction (opus/boxeph/klein)
  - THM-1600  # the Laplace layer L at degree 1 (HYP-8350)
related:
  - THM-1610  # the Watson/Borel reduction of the radial layer
script: 04-computation/gmc2_through_tnc_macmini_S143.py (+ .out)
---

# THM-1645 — GMC(2) through TNC: the angular layer is closed, the gap is radial

## The two layers

The polar decomposition of one complex Gaussian (`s = |Z|² ~ Exp(1)`, `u = Z/|Z|`
uniform, independent) gives the **polar bridge** (already canon as THM-1540(A), re-verified
here exactly against Wick expansion):

> `E[P^m] = ∫₀^∞ CT_u[ Λ_s(u)^m ] e^{−s} ds = L_s( CT_u[ Λ_s^m ] )`, `Λ_s(u) = P(√s·u, √s/u)`.

So `E` is the composite of an **angular** functional `CT_u` (constant term in `u = e^{iθ}`) and
a **radial** functional `L(g) = ∫₀^∞ g(s) e^{−s} ds`, `L(s^k) = k!`. GMC(2) is a statement about
this composite.

## The key new fact: charge support is `s`-independent

A monomial `Z^a W^b` maps to `s^{(a+b)/2} u^{a−b}`. The exponent `a−b` is the **charge**, and
the coefficient `s^{(a+b)/2}` is strictly positive for `s > 0`, so it never annihilates a
charge. Hence the set of charges present in `Λ_s` is constant for all but finitely many `s`
(the zeros of the shell-coefficient polynomials), and

> **`P` is two-sided (has charges of both signs) ⟺ `Λ_s` is two-sided for a.e. `s > 0`.**

## The angular layer is closed by TNC = DvdK

`Λ_s` is a **one-variable** Laurent polynomial. By THM-1630, the toral nullcone conjecture is
exactly the one-variable Duistermaat–van der Kallen theorem — **proved in 1998**. Combined with
`s`-independence:

> **If `P` is two-sided, then `Λ_s` is two-sided for a.e. `s`, so by DvdK
> `CT_u[Λ_s^m] ≠ 0` for some `m`, for a.e. `s`.** The angular direction has no open content.

Verified: two-sided `P` produce a first `m` with `CT_u[Λ_s^m] ≠ 0` at every sampled `s`
(`s = 1, 4, 1/9`).

## Therefore the GMC(2) gap is purely radial

THM-1540 reduces GMC(2) to the nullcone-structure theorem `N = {one-sided} ∪ {0}` and names
the remaining gap (its part III): the descent *"top-degree part one-sided ⟹ `P` one-sided."*
This file locates that gap exactly:

> **The gap is NOT angular geometry.** The angular layer is DvdK-closed and `s`-uniform. The
> gap is the **radial** implication "`CT_u[Λ_s^m] ≠ 0` for some `m` at a.e. `s`" (which DvdK
> gives) ⟹ "`L_s(CT_u[Λ_s^m]) ≠ 0` for some fixed `m`" (which GMC(2) needs). That is a
> `pointwise-nonzero ⟹ integrated-nonzero` descent through `L`.

**The obstruction, exhibited.** `L(t − 1) = 1! − 0! = 0` and `L(t² − 3t + 1) = 0`, with the
polynomials nonzero. So a nonzero `g(s) = CT_u[Λ_s^m]` can have `L(g) = 0` — **integrated
vanishing does not force pointwise vanishing.** This is exactly why the descent is a genuine
gap and why GMC(2) is strictly harder than the (now settled) TNC. The remaining difficulty of
GMC(2) is **Laplace determinacy, not toral geometry.**

The shell grading makes the descent concrete: `Λ_s = Σ_j s^{j/2} λ_j(u)` with `λ_j` the
`(a+b)=j` shell, and `E[P^m]` is dominated as `s → ∞` by the top shell, contributing
`(Dm/2)! · CT_u[λ_D^m]`. If `λ_D` is two-sided, DvdK on `λ_D` makes this eventually nonzero
(this is THM-1540's proved L2). If `λ_D` is one-sided, the top term vanishes and one must
descend to cross-shell terms `CT_u[λ_D^{m−j} λ_{D'}^j]` — the coupling `L` cannot see through,
and that single step is the whole open problem.

## Unconditional corollaries (using the settled angular layer)

- **One-sided `P` ⟹ MZ.** If all charges of `P` are `> 0`, then every monomial of `Q·P^m` has
  charge `≥ m·(min charge) − deg_charge Q`, positive once `m > deg_charge Q`, so `E[Q P^m] = 0`.
  Verified: for `P` with charges `{1,2,2}`, `E[Q P^m] = 0` from `m = 1 + deg_charge Q` on. This
  is DvdK's easy direction in polar coordinates.
- **Two-monomial `Z^p + W^q` is never in the nullcone.** `E[(Z^p+W^q)^m] ≠ 0` exactly when
  `(p+q) | pm`, so first nonzero at `m = (p+q)/gcd(p,q)`; verified `(p,q) = (1,1),(2,2),(2,3),(1,3)`.

## Honest scope

- **GMC(2) is NOT proved.** This is a structural relocation of the open gap, nothing more. It
  says the angular half is done and the radial half is where all the difficulty now lives.
- The polar bridge and the nullcone-structure reduction are **prior work** (THM-1540,
  opus/boxeph/klein); the contribution here is (i) the `s`-independence of charge support,
  (ii) the resulting uniform application of the *now-proved* TNC=DvdK, and (iii) the explicit
  identification of the residual gap as radial + the `L(t−1)=0` obstruction.
- The radial layer is itself **open** (HYP-8350), reduced but not closed by the Watson/Borel
  argument of THM-1610(E). So GMC(2) rests on two things: closing HYP-8350, *and* handling the
  cross-shell coupling that the single functional `L` mixes. Neither is done.
- Nothing here bears on GMC(`n`) for `n ≥ 3` (those are false, THM-1500) or on the effective
  DvdK bound (HYP-8460).

*Artifacts:* `04-computation/gmc2_through_tnc_macmini_S143.py` (+out).
