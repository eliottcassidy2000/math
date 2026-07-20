---
id: THM-1305
title: The equivariant anatomy of the Jacobian counterexample — the six-function normal form (zA+y^kB, yC+x^{k−1}zD, x(E₀−E₁s)) is fully general for weights (1,−1,−k); the s-graded determinant law (c₂, c₁, c₀ derived and verified); c₂ = 0 forces A = v^{k+1}, D = δv² at EVERY k (the cube of u is theorem, not choice); every rational in the owner's expression is DERIVED (Φ = tC+E₀D = 4t+6 linear by det-cancellation, root t* = −3/2; Ψ = E₀A+t²B = (1+t)(2+t), Ψ(t*) = −1/4; E₀(t*) = 13/2; unit-crossing −1/2); and the counterexample is INFINITESIMALLY RIGID modulo equivariant automorphisms (moduli tangent 0; each row unique up to scalar given the other two).
status: >
  PROVED/VERIFIED-EXACT — the normal form (weight-monomial classification, three
  lines); the decode identities (I)–(VI) (exact script); the s-grading formulas
  for det J (hand-derived, engine-verified at k=2: c₂ = 0, c₁ = 0, c₀ = −2);
  the universal relation c₂ = −2A′D + (k+1)AD′ (E₁ = const), hence A = v^{k+1},
  D = δv²; the rigidity computations (nullity 9 = orbit-tangent 9 within the
  22-parameter chart; per-row solution spaces are exactly the scalar lines
  through the witness). EVIDENCE (numeric, labeled) — the k=3 rung: the reduced
  2-parameter slice is EXACTLY inconsistent ((e₁+1)² = 0 forced, then
  conflicting irrationals over ℚ(√13) vs ℚ(√305)); the full 7-parameter
  homogeneous quadratic system + δ hunt is running at close (verdict to be
  appended). OPEN — higher degree caps, deg(v) ≥ 2, s-dependent E beyond
  linear, non-equivariant deformations at higher ambient degree.
source: death-star-2026-07-19-S59n (HYP-8080; owner prompt "sporadic or families + rationals archaeology")
depends_on:
  - THM-1300  # the verified counterexample and its C* structure
related:
  - HYP-8080; kind-pasteur S128c97 (Campbell stratum); mac-mini S129 (equivariant-family question, filed); the three S59n archaeology sweeps (session record)
scripts:
  - 04-computation/jacobian_decode_identities_deathstar_S59n.py -> 05-knowledge/results/jacobian_decode_identities_deathstar_S59n.out
  - 04-computation/jacobian_equivariant_family_deathstar_S59n.py -> 05-knowledge/results/jacobian_equivariant_family_deathstar_S59n.out
  - 04-computation/jacobian_orbit_tangent_deathstar_S59n.py -> 05-knowledge/results/jacobian_orbit_tangent_deathstar_S59n.out
  - 04-computation/jacobian_weight_ladder_solver_deathstar_S59n.py -> 05-knowledge/results/jacobian_weight_ladder_solver_deathstar_S59n.out
  - 04-computation/jacobian_k3_proper_solve_deathstar_S59n.py -> 05-knowledge/results/jacobian_k3_proper_solve_deathstar_S59n.out
  - 04-computation/jacobian_k_ladder_hunt_deathstar_S59n.py -> 05-knowledge/results/jacobian_k_ladder_hunt_deathstar_S59n.out
---

# THM-1305 — anatomy, decode, and rigidity of the equivariant Keller counterexample

## 1. The normal form (PROVED, fully general)

For the ℂ*-action with weights w(x,y,z) = (1, −1, −k), k ≥ 2, the invariant
ring is ℚ[t, s], t = xy, s = x^k z, and EVERY equivariant polynomial map with
component weights (−k, −1, +1) is exactly

  F₁ = z·A(t) + y^k·B(t),  F₂ = y·C(t) + x^{k−1}z·D(t),  F₃ = x·E(t, s)

(the weight-monomial classification is three lines). The owner's map is the
case k = 2, E = E₀ − s: **six univariate functions**
A = u³, B = u(4+3t), C = 1+12t+9t², D = 3u², E₀ = 2−3t, E₁ = 1 (u = 1+t).

## 2. The determinant law (derived + verified)

det JF is weight-0, hence ∈ ℚ[t,s]; with E = E₀ − E₁s it is quadratic in s:
det = c₀(t) + c₁(t)s + c₂(t)s². Hand-derivation (verified exactly by engine at
k = 2, reproducing 0, 0, −2):

- c₂ = −(k+1)E₁(A′D − AD′) + (k−1)E₁A′D − tE₁′(…)  |_{E₁ = 1} = **−2A′D + (k+1)AD′**
- c₁ = E₀(A′D−AD′) + (k−1)E₀′AD − 2t^kB′D − 2kt^{k−1}BD + kt^kBD′ + (k+1)AC + (k+1)tAC′ − tA′C
- c₀ = (E₀+tE₀′)(t^kB′D + kt^{k−1}BD − AC − tAC′) − E₀′(t^{k+1}B′D − t²AC′) − t^{k+1}(B′C − kBC′)

**Keller ⟺ c₂ = 0, c₁ = 0, c₀ = nonzero const.**

**The cube is forced:** c₂ = 0 ⟺ (k+1)AD′ = 2A′D ⟺ D³⁺ᵏ⁻¹…, concretely
**A = v^{k+1}, D = δv²** for a polynomial v and constant δ — at EVERY k the
D-row is quadratic in v (they coincide with "D = A′-shaped" only at k = 2,
which is what makes the witness look deceptively differentiable: D = A′ there).
For the witness: v = u, δ = 3.

## 3. The rationals decode (VERIFIED — nothing in the expression is arbitrary)

On the invariant plane, the collision curve is {E = 0} = {s = E₀(t)}, and:

- **Φ(t) := tC + E₀D** (the F₂ = 0 condition on the curve) computes to **4t + 6**
  — cubic and quadratic terms cancel IDENTICALLY (the same cancellation that
  makes det constant), leaving a LINEAR polynomial whose unique root is
  **t\* = −3/2**. A root must exist unless Φ is a rootless constant: within
  this class, *collision is the generic fate*, which is why "Keller ⟹
  injective" had to fail here if it failed anywhere (kind-pasteur's
  Campbell-minimal stratum, seen from inside).
- **Ψ(t) := E₀A + t²B** (the collision value) collapses to **(1+t)(2+t)** —
  roots −1, −2 — and **Ψ(t\*) = −1/4** is the collision target.
- The two weight-0 units CROSS at t\*: u(t\*) = (4+3t)(t\*) = **−1/2**.
- **E₀(t\*) = 2 − 3t\* = 13/2** = s\* = the z-coordinate at λ = 1. So
  **13 = 2·Φ-root arithmetic (2 + 9/2), a DERIVED value** — it varies across
  any family and is not a structural constant. (The archaeology sweep's
  verdict: the repo's load-bearing 13s — Φ₆(14)/14 = 13 + 1/14, 13 = 14² mod
  183, D/(13D+k) — do NOT bridge to this 13; the one shared verb is
  "squaring": 13 = 14² mod 183 vs λ ↦ λ². Graded speculative.)
- The invariant coordinates of the collision orbit are **(t\*, s\*) =
  (−3/2, 13/2)** — the owner's two "mystery rationals" are exactly the
  collision's coordinates on the quotient plane.

Small naturals (1,2,3,4) = weights + unit coefficients: STRUCTURAL.
Big rationals (13/2, −1/4) = values of derived polynomials at the forced root:
DERIVED. That is the complete answer to "which naturals appear."

**Amendment (S59o, THM-1320):** the −2 row of the decode table is now DERIVED: det = c₀(0) = −E₀(0)·A(0)·C(0) = −(2·1·1) — a polynomial identity at every k. Nothing in the expression remains arbitrary.

## 4. Rigidity: sporadic, not a family (within the equivariant chart)

Three independent computations at k = 2:

- **Moduli tangent = 0**: in the 22-parameter coefficient chart the det-
  constancy linearization has nullity 9, and the tangent space of the
  (equivariant Aut)∘F∘(equivariant Aut) orbit accounts for ALL 9 directions.
- **Per-row rigidity**: fixing any two rows, the third row's solution space is
  exactly the 1-dimensional scalar line through the witness (verified for all
  three rows).
- The 9-dim kernel decomposes as reparametrization only — **no genuine
  deformation exists at these degrees**: the counterexample is infinitesimally
  SPORADIC modulo automorphisms, not a member of a positive-dimensional family.
  (Honest scope: degree caps as in the scripts; E linear in s; deformations
  outside the weight class untested at close.)

## 5. The k-ladder question (k = 3: obstruction found)

The natural "next rung" — weights (1,−1,−3), predicted fiber 1 + 3 (orbit
λ³ = 1) — meets real resistance: with A = v⁴, D = δv², the c₁-system is
consistent (7-parameter affine solution space), but imposing c₀-constancy:

- the natural 2-parameter slice is EXACTLY inconsistent: t⁷ forces
  (e₁+1)² = 0, then t¹ demands e₀ ∈ ℚ(√13) while t² demands e₀ ∈ ℚ(√305);
- **the full system is EMPTY — two validated instruments agree** (appended at
  close): (a) exact mod-p projective scans over ALL δ and the whole c₁-kernel
  space: k=2 validation finds Keller points at every prime (exactly one
  δ-orbit, matching rigidity), while **k=3 has ZERO points at p = 5, p = 7,
  and δ-sampled p = 11**; (b) the Keller-forced complex Newton hunt (the
  reciprocal-variable equation c₀(0)·ν = 1 deletes the degenerate basin):
  k=2 validation rediscovers the witness orbit from random starts (coefficient
  ratios (1, 1.75, 0.75) and (1, −1.5) — the witness modulo the complex
  equivariant torus, rigidity seen numerically), while **k=3 scores 0/80**.
  VERDICT (evidence-grade, natural degrees): the k = 2 witness is not the
  first rung of a ladder but an **ISOLATED SPECIES** — sporadic at the
  deformation level AND at the weight-type level. W1 sharpens: only the
  ℤ/2-protected rung exists; the even-fiber species the parity lemma cannot
  protect appears not to exist at all (at these degrees).

## 5b. W3 tested (Φ-linearity is NOT universal)

On the k=3 c₁-solution space, Φ = tC + E₀D retains nonzero coefficients at
t²..t⁵ (exact λ-linear functionals, script appended) — c₁ = 0 alone does NOT
force Φ linear beyond k = 2. The collision-forcing linearity of Φ at the
witness is a property of its rung, not of the grading calculus.

## 6. What this hands the fleet

- The normal form + grading is a complete calculus for weight-(1,−1,−k)
  Keller maps: any future equivariant counterexample search is 3 linear
  solves + 1 quadratic system, not a blind det expansion.
- The TRAP reading (polysemous-constants template): the collision is two
  units crossing at one point — a codimension-1 coincidence, and Φ linear
  means the coincidence is UNAVOIDABLE in this class. "Keller counterexamples
  = maps whose unit-crossing cannot be dodged" is the bold one-line summary.
- Deg(v) ≥ 2, s-quadratic E, and mixed-weight classes remain unexplored — the
  three cheapest continuations, each one session with this file's machinery.
