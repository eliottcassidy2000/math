---
id: THM-1345
title: The two-variable Jacobian conjecture in the ℂ*-equivariant category — JC₂ is a THEOREM there (equivariant Keller ⟹ invertible for EVERY ℂ*-action: hyperbolic ⟹ linear [boxeph-S144], elliptic ⟹ triangular [completed here]); the Poisson reframing det JF = {P,Q} with the ℂ*-action = Hamiltonian flow of xy, so equivariant-JC₂ = the classical shadow of eigen-graded DC₁; and the weight-additive bracket makes the leading-form obstruction {P_top,Q_top}=0 the recursive engine locating full JC₂'s hardness in the descent through the weight filtration.
status: >
  PROVED/VERIFIED (category-restricted) — (A) equivariant Keller ⟹ invertible:
  hyperbolic weights (a,−b) ⟹ LINEAR (boxeph-S144's leading-coefficient
  argument, re-derived: det = (s·fg)' telescopes at a=b=1, the single top-degree
  term (1+a·deg g+b·deg f)·f_top·g_top can't vanish for all weights — verified
  (a,b) ≤ 3); elliptic weights (w₁,w₂)>0 ⟹ TRIANGULAR/linear (properness from
  F⁻¹(0)={0} + finite étale over simply-connected ℂ² ⟹ degree 1 ⟹ automorphism;
  formal inverses terminate on samples). (B) det JF = {P,Q} and the bracket is
  weight-additive (verified); the ℂ*-action is the Hamiltonian flow of xy
  (X_{xy} = x∂_x − y∂_y). (P5) the leading-form instrument validated (rediscovers
  {P_top,Q_top}=0 on a genuine tame automorphism). NOT PROVED — full JC₂
  (OPEN since 1939; the notorious false-proof graveyard — this file makes NO
  claim on it). The propagation of {P_top,Q_top}=0 down the filtration is exactly
  the classical Abhyankar–Moh–Suzuki-hard open content.
source: death-star-2026-07-20-S59q (HYP-8155; owner: work creatively on JC₂). Builds on boxeph-S144 (hyperbolic equivariant case — CREDITED).
depends_on:
  - THM-1300, THM-1305  # the 3D counterexample and its equivariant anatomy
related:
  - boxeph-S144 (hyperbolic equivariant Keller ⟹ linear; the dim-2 no-go)
  - the surviving ladder DC₂ ⟹ JC₂ ⟹ DC₁ (S59m)
scripts:
  - 04-computation/jc2_equivariant_and_poisson_deathstar_S59q.py -> 05-knowledge/results/jc2_equivariant_and_poisson_deathstar_S59q.out
---

# THM-1345 — JC₂ in the equivariant category, and where the real difficulty lives

**Honesty banner.** Full JC₂ is OPEN and this file does not touch it. Everything
below is either (i) a theorem in the restricted ℂ*-equivariant category, or
(ii) a reframing that *locates* the difficulty. JC₂ is the graveyard of false
proofs; the value here is structural, not a resolution.

## 1. det JF is the Poisson bracket (the right language)

For F = (P, Q): ℂ² → ℂ², det JF = P_x Q_y − P_y Q_x = **{P, Q}**, the Poisson
bracket of the standard symplectic form dx∧dy (verified). So JC₂ reads:

  **a polynomial canonical pair ({P,Q} = const ≠ 0) is a coordinate system.**

This is the classical (ℏ→0) shadow of Dixmier's DC₁: [P,Q] = 1 in the Weyl
algebra A₁ ⟹ (P,Q) generate. The repo's surviving ladder is
**DC₂ ⟹ JC₂ ⟹ DC₁** — proving JC₂ would settle Dixmier's original 1968
question DC₁.

## 2. The ℂ*-action is a Hamiltonian symmetry; equivariant = ad_{xy}-graded

The action λ·(x,y) = (λx, λ⁻¹y) is the Hamiltonian flow of **H = xy**:
X_{xy} = {xy, ·} = x∂_x − y∂_y. A Keller pair is *equivariant* iff
{xy, P} = P and {xy, Q} = −Q, i.e. P, Q are **eigenvectors of ad_{xy}** of
eigenvalue ±1 (weights ±1). So "equivariant Keller pair" = "canonical pair
diagonalizing the Hamiltonian xy." THM (this file): such pairs are LINEAR.
Equivalently: **the classical shadow of eigen-graded DC₁ holds** — which is
exactly the sub-case Dixmier-type homogeneity arguments have always reached.

## 3. Equivariant JC₂ is a theorem — for every ℂ*-action

- **Hyperbolic** weights (a, −b), a,b>0 (boxeph-S144, re-derived). Equivariance
  forces P = x·f(s), Q = y·g(s), s = x^b y^a, and
  **det = fg + s(a·fg′ + b·f′g)**. At a=b=1 this telescopes to (s·fg)′, so
  det = c ⟹ fg = c ⟹ f,g constant. For general (a,b) the top s-degree
  coefficient is the *single* term (1 + a·deg g + b·deg f)·f_topg_top ≠ 0, so
  det constant forces deg f = deg g = 0: **F is linear.** (Verified (a,b) ≤ 3.)
- **Elliptic** weights (w₁, w₂) > 0 (completed here). P, Q weighted-homogeneous
  with w(P)+w(Q) = w₁+w₂. Keller ⟹ F⁻¹(0) is discrete (étale) and ℂ*-invariant
  ⟹ = {0} ⟹ F finite; finite étale over simply-connected ℂ² ⟹ trivial cover
  ⟹ degree |F⁻¹(0)| = 1 ⟹ **F is an automorphism** (a weighted-triangular one,
  e.g. (x, y + x^k); formal inverses terminate on all samples).

**Net: in the ℂ*-equivariant category, JC₂ is TRUE.** This is the exact 2D
mirror of the 3D non-descent (boxeph-S144): the machinery that BUILT the
3D counterexample cannot even pose a problem in 2D — it produces only
automorphisms.

## 4. Why 2D has no room (W1 made exact)

The 3D counterexample is a cascade of coupled units (THM-1305: A=v³, B, C, D,
E₀ with the free unit v = 1+t surviving the c₁∧c₀ system). Its 2D equivariant
analog collapses to the *single* telescoped constraint (s·fg)′ = const ⟹
fg = const — one product equation, no coupling, so no nonconstant unit can
appear and the unit-crossing that powers the collision cannot form. The
counterexample's engine needs ≥ 2 coupled units; **dim 2 supplies exactly
one.** This is the precise content of the earlier "no room at n=2" intuition.

## 5. Where full JC₂ actually lives: descent through the weight filtration

The bracket is **weight-additive**: {P_a, Q_b} has weight a+b (verified). Give
ℂ[x,y] the (1,−1) grading; write P = Σ_w P_w, Q = Σ_w Q_w. For a Keller pair,
the top-weight part of {P,Q} is {P_A, Q_B} (A, B the top weights) and it must
VANISH: **the leading weight-forms Poisson-commute**, hence are algebraically
dependent. (Instrument validated in P5: for the genuine automorphism
(x + (y+x²)³, y + x²) the top ordinary-degree forms satisfy {P_top,Q_top}=0;
a random non-Keller pair fails it.) This is the classical Abhyankar–Moh
leading-form obstruction, recovered from the bracket's grading.

So JC₂ decomposes cleanly:
- **base case** (single weight = equivariant): SETTLED (§3);
- **inductive step** (propagate {P_A,Q_B}=0 down the filtration to force the
  whole pair triangularizable): the OPEN, AMS-hard content.

The equivariant theorem removes the base case from the mystery and names the
remaining difficulty precisely: **can a canonical pair evade being organized by
any Hamiltonian ℂ*-symmetry?** A JC₂ counterexample would be exactly a
canonical pair whose leading-form dependence never propagates to a global
coordinate — a pair "anti-equivariant at every weight." Whether that is
possible is the whole question; this file does not decide it.

## 6. The DC₁ / observer bridge (graded honestly)

Under the repo's observer dictionary A_k ↔ tournaments on n = 2k+1 vertices
(S59m §4), DC₁ ↔ A₁ ↔ **n = 3, the 3-cycle** — the atomic odd cycle where
Rédei parity is born, and the surviving floor of the whole doubling tower.
So the last open descendant of the JC/Dixmier collapse sits, under this
dictionary, on the repo's most-studied atom. This is an ANALOGY (graded
speculative), but it says the natural repo tool for the surviving frontier is
its 3-cycle / atomic-parity technology — and it continues the S59p
"conservation of +1" thread: the +1 that the classical reduction exiles is the
same +1 (u = 1+xy, the formal-group denominator, the observer) whose absence
in 2D leaves "no room."
