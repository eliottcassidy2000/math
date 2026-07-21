---
id: THM-1810
title: "THE BOSONIC/FERMIONIC (PERMANENT/DETERMINANT) DICHOTOMY EXPLAINS WHY GMC(2) IS HARD AND TRANSITIVITY IS EASY. The determinant/discriminant/Vandermonde side is FERMIONIC and ALTERNATING: the sign lets the 3-cycle-reversing involution cancel all intransitivity in one step, so only the n! transitive tournaments survive (THM-1805, boxeph THM-1800) — transitivity is a clean classical fact. Its symmetric partner is the PERMANENT/HAFNIAN, which is the GAUSSIAN moment: E[|Z|^{2a}] = a! = per(J_a) (bosonic Wick = permanent of the covariance), whereas the fermionic moment VANISHES on a repeated state, det(J_a)=0 for a≥2 (Pauli). The Laplace moment engine E[P^m] = L_s(CT_u[Λ_s^m]) is this permanent/hafnian functional: NO sign, so NOTHING cancels, so there is no one-step or finite-moment-cutoff detection — which is exactly the EMP floor (detection depth ≥ d+1, growing with radial degree, THM-1790). GMC(2) is hard for the same structural reason the permanent is #P-hard and the determinant is P: it lives on the bosonic side where the alternating collapse is unavailable. The polar bridge crosses the det/per divide — angular DvdK = the fermionic GIT nullcone, radial EMP = the bosonic Laplace moment — which is why coupling them cannot be finite-degree-uniform (THM-1770)."
status: >
  VERIFIED-EXACT and structural. (1) E[|Z|^{2a}] = a! = per(J_a) checked a = 1..6; the general
  complex-Gaussian Wick = permanent of the covariance matrix is classical (bosonic Wick). (2)
  det(J_a) = 0 for a ≥ 2 (fermionic vanishing / Pauli), checked. (3) det(Vandermonde) = signed
  tournament sum (n! transitive survivors) vs per(Vandermonde) = unsigned (all 2^{C(n,2)}
  tournaments), verified n = 3,4. The DICHOTOMY and its consequence for GMC(2) hardness is a
  structural reading, not a new theorem — it unifies THM-1805 (fermionic/transitivity),
  THM-1790 (the bosonic EMP floor), and THM-1770 (the bridge is not degree-uniform) under one
  frame, and it is complementary to boxeph's fermionic THM-1800.
  Explains a difficulty; proves no open problem. GMC(2) remains OPEN.
source: klein-2026-07-20-S387 (owner: think more about discriminants and determinants and other similar concepts)
depends_on:
  - THM-1805  # the Vandermonde signed tournament sum (the fermionic side, transitivity)
  - THM-1790  # the EMP floor (the bosonic side: no cancellation ⟹ growing detection depth)
related:
  - THM-1800  # boxeph: binary forms ↔ tournaments, the determinant/discriminant (fermionic) dictionary
  - THM-1475  # Pf(S(T)) is the odd function — the tournament's fermionic (Pfaffian) moment
  - THM-1770  # the bridge is not finite-degree-uniform (it crosses the det/per divide)
attribution: >
  The fermionic side (Vandermonde/discriminant = transitivity) is THM-1805 (mine) and boxeph's
  THM-1800. This file adds the BOSONIC complement (permanent/hafnian = the Gaussian moment engine)
  and the resulting hardness reading. The permanent-is-#P-hard / determinant-is-P contrast is
  Valiant's; its use here to explain GMC(2)'s difficulty is the new framing.
script: 04-computation/bosonic_fermionic_klein_S387.py (+ .out)
---

# THM-1810 — the bosonic side is why GMC(2) is hard

## Two moments of the same Gaussian

Complex Wick has two incarnations, the same combinatorics with and without signs:

| | combinatorics | moment | on a repeated state |
|---|---|---|---|
| **bosonic** (Gaussian) | matchings, **all `+`** | `E[Z^aZ̄^b] = δ_{ab}·a!` = **`per(J_a)`** (hafnian/permanent) | `a! ≠ 0` — rich |
| **fermionic** (Grassmann) | matchings, **signed** | Pfaffian / **determinant** | `det(J_a) = 0` (`a≥2`) — **Pauli** |

Verified: `E[|Z|^{2a}] = a! = per(J_a)` for `a = 1..6`; `det(J_a) = 0` for `a ≥ 2`. The Gaussian
moment is literally the **permanent** of the covariance; the fermionic partner **vanishes** on
any repeated state.

## The same split on tournaments

The Vandermonde matrix `V_{ij} = x_i^{j}` gives both invariants, and both are tournament sums:

```text
det V = Σ_T sgn(T)·x^{score(T)}   (FERMIONIC)   — only the n! TRANSITIVE tournaments survive
per V = Σ_T       x^{score(T)}    (BOSONIC)      — ALL 2^{C(n,2)} tournaments contribute, +
```

Verified `n = 3,4`: the signed sum has exactly `n!` nonzero-coefficient scores (the transitive
survivors, THM-1805); the unsigned sum totals `2^{C(n,2)}` — every tournament, no cancellation.

## The load-bearing reading — why GMC(2) is hard

- **Fermionic / determinant side** (boxeph THM-1800, THM-1805, THM-1475's Pfaffian). A sign
  twists the sum, and the 3-cycle-reversing involution **cancels all intransitivity in one
  step**. Transitivity is therefore a clean, finite, classical fact, and its GIT nullcone
  (angular DvdK) is a theorem. The determinant is `P`.
- **Bosonic / permanent side** (this file). The Gaussian moment engine
  `E[P^m] = L_s(CT_u[Λ_s^m])` is the **permanent/hafnian** functional — **no sign, so nothing
  cancels**; every Wick matching contributes with `+`. With no alternating collapse available,
  the nullcone cannot be detected in one step or at a finite moment cutoff. This is exactly the
  **EMP floor**: detection depth `≥ d+1`, growing with the radial degree (THM-1790). The
  permanent is `#P`-hard.

> **GMC(2) is hard for the same structural reason the permanent is `#P`-hard and the determinant
> is `P`: it lives on the bosonic side, where the alternating cancellation that trivialises
> transitivity is unavailable.**

## The bridge crosses the divide

The polar bridge is exactly a crossing of the `det/per` divide:

```text
   angular layer   CT_u   = the FERMIONIC GIT nullcone  = DvdK        (cancellation available — proved)
   radial layer    L_s    = the BOSONIC Laplace moment  = EMP         (no cancellation — depth grows)
```

Coupling a cancelling (fermionic, degree-blind) layer to a non-cancelling (bosonic,
degree-growing) one is why the coupling **cannot be finite-degree-uniform** (THM-1770), and why
the EMP/Laplace asymptotic — not elimination — is the tool that survives the degree (THM-1790).
The det/per language names the obstruction: you cannot alternate-cancel a permanent.

## Scope

A structural synthesis, exact where checked (`n ≤ 4`, `a ≤ 6`), unifying the fermionic
(transitivity, THM-1805/1800/1475) and bosonic (moment engine, THM-1790) threads and reading
GMC(2)'s difficulty off Valiant's permanent/determinant divide. Proves no open problem; it
explains why the open problem is open, which constrains the search: any GMC(2) proof must be a
genuinely bosonic (permanent-side) argument — the fermionic tricks that work for transitivity do
not transfer.

*Files: `04-computation/bosonic_fermionic_klein_S387.py` (+ `.out`).*
