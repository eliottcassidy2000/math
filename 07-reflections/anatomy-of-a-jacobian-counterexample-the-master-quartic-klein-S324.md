# Anatomy of a Jacobian counterexample: the master quartic, the forced S₃, and why the "3" is not a choice

**Instance:** klein-2026-07-19-S324 (owner: long session, deepest possible understanding
of JC counterexamples — sporadic vs families, explicit examples, rational-thread sweep,
bold pattern-hunting licensed). Complements opus-S414's claimed angles (torus normal
form, quotient contraction, general coefficient moduli — NOT duplicated here).
**Everything computational below is exact** (sympy; the elimination, discriminant,
fiber counts, and the k=2 ansatz are reproducible from the S323/S324 scripts).

## 1. The structure theorem: one quartic governs everything

Eliminating (x, y, z) from F = (a, b, c) together with W = 1 + xy gives the
**universal fiber cubic**

> **D(a,b,c)·W³ + (b² − 12a)·W − 4a = 0,  D := 27a²c² − 18abc + 16a + b³c − b²**

- **No W² term, identically** — the three fiber w-values sum to 0 at EVERY target
  (equivalently Σ xᵢyᵢ = −3 over every fiber): the first universal trace law.
- At the collision target (−1/4, 0, 0): D = −4 and the cubic is **4W³ − 3W = 1 =
  T₃(W)** — Chebyshev on the nose; the fiber's w-values are cos(2πj/3) = {1, −½, −½}.
- **disc_W = −4·S²·D with S = ∂D/∂c** (54a²c − 18ab + b³). The Galois closure's
  quadratic subfield is ℚ(√(−D)): the cover is a T₃-shape twisted by the square
  class of its own leading coefficient. One quartic D and its c-derivative generate
  the entire geometry.
- On S = 0 (D ≠ 0) the double root is **r³ = 2a/D exactly** (sampled instances
  rational: r = 6/11, 3/4, 6/13 — the 11³, 13³ in the D-values are forced by this
  identity, not numerology).
- **The exhibited collision target lies ON S = 0** (a = −1/4, b = c = 0 ⟹ S = 0,
  D = −4): the cubic factors −(W−1)(2W+1)² — the collision is a fully-rational
  point of the double-root stratum, its twisted pair sharing w = −1/2 = r with
  r³ = 2(−1/4)/(−4) = 1/8 ✓. The famous three points are the S = 0 geometry made
  rational.

## 2. The fiber-count spectrum is {3, 1} — never 2 — and F is plausibly surjective

Exact fiber counts: 3 whenever D ≠ 0 (even on S = 0, where two fiber points share a
w-value but stay distinct); **1 on D = 0** (the cubic degenerates to linear — two
preimages escape to infinity; the target stays IN the image). No 0-fiber point found;
no fiber-2 point exists off the escape locus. Two independent instruments agree:
the mod-p Chebotarev census measured **N₂ ≡ 0 at every prime** (S323) — the same
structural zero seen now from the ℂ side. CONJECTURE (S324): **F is surjective** —
a surjective, everywhere-étale, 3-to-1-generic polynomial self-map of ℂ³.
(Étale + non-proper is exactly how this coexists with simple connectivity: F is
not a covering map; the D = 0 hypersurface is where properness fails.)

## 3. Why degree 3, and why S₃: the minimality mechanism

- A degree-2 étale self-cover of ℂ³ would be Galois: its deck transformation is a
  **fixed-point-free algebraic involution of ℂ³**. If finite-order automorphisms
  of ℂ³ are linearizable (the dim-3 linearization theorems; citation to pin),
  every such involution has a fixed point — **degree 2 is obstructed**.
- A ℤ/3-Galois triple cover needs a fixed-point-free order-3 automorphism — same
  obstruction. **The measured S₃ (not ℤ/3) is forced.**
- A **non-Galois degree-3 cover has trivial deck group** — it needs NO global
  automorphism at all and dodges the obstruction entirely. The counterexample's
  configuration (degree 3, monodromy S₃, √(−D) twist) is thus the MINIMAL shape a
  JC counterexample cover can have.
- Consistent negative, computed: the mechanical **k = 2 transplant of the formula
  fails** — in the 7-parameter ansatz (w²z + p·wy² + qy², y + 2xwz + rxy² + vxy³,
  sx + tx²y + ux²z) the ONLY constant-Jacobian point is the degenerate s = t =
  u = 0, det = 0. The "3" in the formula is structural, not a dial.

## 4. Sporadic or families? The refined answer

Counterexamples form (at least) three lawful multiplication structures, and the
interesting moduli live BELOW them:

1. **The conjugation orbit** ψ∘F∘φ, ψ, φ ∈ Aut(ℂ³) — an infinite-dimensional
   "trivial family" (the tame group alone is infinite-dimensional). Essential
   moduli = orbits under this action.
2. **The composition monoid**: F∘F is Keller (det 4) and 9-to-1-generic; F∘τFσ,
   padded products F × id, the ℂ⁶ cotangent lift (S323, symplectic), the Weyl
   endomorphism φ_F (S323) — counterexamples breed by composition. Degrees
   realize 3^m; mixed compositions with automorphisms give the full monoid.
3. **Within-shape moduli**: opus-S414's coefficient-moduli scan (their claim) will
   measure the local dimension around F. The S324 data predict it is SMALL along
   the degree axis (k = 2 dies) and nontrivial along conjugation directions.
4. Higher-dimensional families are FORCED by the classical reductions run in
   reverse: cubic-homogeneous and Drużkowski counterexamples must exist in some
   larger dimension (Bass–Connell–Wright/Yagzhev machinery applied to F ⊕ id) —
   a named construction lead, untouched today.

Prediction (wild, testable): **every minimal-degree counterexample in dim 3 has
fiber spectrum {3, 1}, monodromy S₃, and a master quartic with disc = −4S²D
shape** — the mechanism of §3 leaves no alternative configuration at degree 3.

## 5. The rational-thread sweep (owner's directive), with confidence tags

- **−1/4 (HIGH, mechanism):** −1/4 = −∏ⱼ cos(2πj/3) = −e₃ of the collision fiber's
  w-values = the constant coefficient story of T₃(W) = 1. The SAME 1/4 that is the
  LRC 3-runner floor (equal spacing of three points on ℝ/ℤ). Both are elementary
  shadows of ℤ/3 acting on the circle — the one genuine cross-project mechanism
  found in the sweep.
- **3/2 (MED):** |xy| = 3/2 at the two twisted preimages = 1 − (−½) = the cosine
  spread cos(0) − cos(2π/3); also w = −½ there: the twisted points sit AT the
  nontrivial third-root cosine.
- **13/2 (LOW, numerology):** z = 13/2 with 13 = 9 + 4 = 3² + 2² — the session's
  two primes squared and summed; also LRC(14)'s omnipresent 13. No mechanism.
- **{2,3}-smoothness (HIGH, pattern):** EVERY constant in the map, the cubic, the
  discriminant (27, 18, 16, 12, 4, 3, 2, 54, 108, 216 = 6³) is 2^a·3^b. The
  example is defined over ℤ[1/6], |S₃| = 6. WILD CONJECTURE (H2): minimal JC
  counterexamples are {2,3}-smooth — the smallest-primes wall, the same shape as
  LRC(14)'s apex-7 (there the wall is the smallest prime NOT under control; here
  the entire object is built from the two smallest primes).
- **27, 12, 16 ↔ repo constants** (2/27 rung, the 12m ladder, 2⁴): numeral
  coincidences, cataloged per the directive, NO mechanism claimed.

## 6. Explicit examples produced or identified this arc

1. F itself (verified S323; anatomy above).
2. F∘F and all monoid words (9:1, 27:1, …; det powers of −2).
3. τ∘F∘σ and the conjugation orbit (trivial moduli).
4. λ-scalings: F_λ = λ⁻¹·τ_λ∘F∘σ_λ (torus conjugates — opus-S414's normal form).
5. The padded F × id on ℂⁿ, n ≥ 3.
6. The symplectic ℂ⁶ lift Φ (S323; det 1, Φ*ω = ω, non-injective).
7. The Weyl endomorphism φ_F of A₃ (S323; the DC(3) witness).
8. NEGATIVE: no k = 2 sibling in the natural 7-parameter shape (this session).

## 7. Leads filed (backlog-grade)

(i) Pin the dim-3 linearization citations → upgrade §3 to THEOREM ("no étale
degree-2 self-cover of ℂ³"); (ii) prove/refute surjectivity of F (the D = 0
escape analysis is nearly a proof — one lifting lemma missing); (iii) the D-as-
discriminant near-miss: D vs disc(at³ + bt² + t + c) differ in exactly two
coefficients by factors of 4 — find the correct hidden cubic (a weighted/
regularized resultant?); (iv) construct the S₄/degree-4 analog (non-Galois
quartic cover — the next rung of the minimality ladder: does a 4:1 étale
self-map exist, or does 3 remain the only possibility?); (v) the
cubic-homogeneous/Drużkowski descendants of F in higher dimension;
(vi) mod-p: deficiency δ(p) asymptotics (Lang–Weil error terms measured in the
S323 census — frac3 approaches 1/6 from below as ~c/p; extract c).

## 8. Cross-links

T1547 + the S323 reflection (verification, φ_F, ℂ⁶ lift, S₃ census) ·
opus-S414 (concurrent: torus normal form, moduli — complementary) ·
everything-is-the-triangle (the {2,3}/smallest-prime motif; Cayley–Dickson
doubling = the JC↔DC 2n functor) · the repo's structural-zero epistemology
(N₂ ≡ 0 here; blue/black and SC-NS zeros there) · CONSTANTS-INDEX (2/27, 12m —
numeral-only flags).
