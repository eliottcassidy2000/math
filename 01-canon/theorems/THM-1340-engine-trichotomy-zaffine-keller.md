---
id: THM-1340
title: THE ENGINE TRICHOTOMY for z-affine Keller maps — the projectivized top coefficient lands in a RATIONAL plane curve; degree-1 engines give units, degree-2 engines are impossible, and the minimal viable engine is a rational cubic ∈ {nodal, cuspidal}; F and ALL known variant-radical counterexamples share THE cuspidal cubic γ(s) = (1, 3s, −s³) (cusp at s = ∞, hyperflex at s = 0) — projectively unique — so Mount-Everest-in-the-z-affine-class reduces to [exclude the nodal engine] + [classify B over the unique engine]
status: PROVED for the structural core (det split; engine-curve degeneracy; rationality via Lüroth; the cuspidal identification with the corrected cusp location; the unification across radicals 1+λxy — all exact, frozen out); DERIVED for the degree transfer (engine degree = generic fiber degree — verified on F (3=3), stated as the bridging lemma with the line/conic exclusions conditional on it plus Keller-1939/Campbell); the n=4 (z,w)-affine framework + first box-bounded negative included
source: mac-mini-2026-07-20-S136 (owner: work on Mount Everest + find counterexamples in dimension 4)
depends_on: [THM-1310 (z-linear anatomy lead), THM-1315, THM-1330, Keller 1939, Campbell]
related: [THM-1305 (rigidity), opus-S415/klein-S327-S328 (variant examples — unified here), HYP-8140]
output: 05-knowledge/results/jc_engine_curve_theorem_macmini_S136.out
---

# THM-1340 — the engine trichotomy

## 1. The det split (exact)

For z-affine F = A(x,y)·z + B(x,y) (A, B: ℂ² → ℂ³), det JF is quadratic in z with

    [z²]: det[A_x, A_y, A],   [z¹]: det[A_x, B_y, A] + det[B_x, A_y, A],   [z⁰]: det[B_x, B_y, A].

Keller ⟺ the first two vanish identically and the third is a nonzero constant.
(Verified on F: 0, 0, −2.)

## 2. The engine is a rational plane curve (proved)

det[A_x, A_y, A] ≡ 0 says A, A_x, A_y are everywhere dependent: the projectivization
P(A): ℂ² → ℙ² has rank-≤1 differential, so its image lies in a CURVE (or point), and A
factors as A = λ(x,y)·γ(t(x,y)) — the ENGINE. Since γ is dominated by a polynomial map from
ℂ², the engine curve is unirational, hence (Lüroth) RATIONAL: **elliptic and higher-genus
engines are impossible.** For F: A = u³·γ(x/u) with γ(s) = (1, 3s, −s³), exactly.

## 3. The trichotomy (degree transfer marked DERIVED)

Bridging lemma (DERIVED; equality verified on F where 3 = 3): the generic fiber degree of a
z-affine Keller map equals the degree of its engine curve (the fiber resolvent is the
engine-curve equation traced along the syzygy parametrization). Granting it:

- **engine degree 1 (line):** fiber degree 1 ⟹ birational ⟹ invertible (Keller 1939) — units only;
- **engine degree 2 (conic):** fiber degree 2 ⟹ Galois ⟹ invertible (Campbell) — impossible;
- **engine degree 3:** the minimal viable case; rational plane cubics are NODAL or CUSPIDAL
  (smooth excluded by §2). F's engine is the CUSPIDAL cubic: in the chart at [0:0:1],
  (p,q) = (−1/s³, −3/s²) satisfies 27p² + q³ = 0 — a (2,3)-cusp AT s = ∞ (s = 0 is a
  hyperflex; an earlier in-session line mislocating the cusp at s = 0 is corrected here).

## 4. Unification (exact): all known variants share THE engine

For every radical u_λ = 1 + λxy (λ = 1 the original; λ = −2, 3 = opus-S415's variants):
A_λ = (u_λ³, 3x u_λ², −x³) satisfies det[A_x, A_y, A] ≡ 0 and P(A_λ) = [1 : 3s : −s³] with
s = x/u_λ — THE SAME cuspidal cubic, reparametrized. Cuspidal plane cubics form a single
PGL₃-orbit (classical), so **the A-part of every known degree-3 counterexample is one
projectively rigid object.** Mount Everest in the z-affine class therefore reduces to:
(i) EXCLUDE the nodal-cubic engine (named open); (ii) classify the B-solutions of the two
remaining conditions over the fixed cuspidal engine (kps-THM-1310's design equations are
this classification's first chart; their in-box uniqueness becomes engine-level).

## 5. The n = 4 program: the (z,w)-affine framework and its first negative

For F = A z + B w + C: ℂ⁴ → ℂ⁴ (A, B, C: ℂ² → ℂ⁴), det JF is quadratic in (z,w) with SIX
conditions: det[A_x,A_y,A,B] = det[B_x,B_y,A,B] = mixed = 0 (the PENCIL degeneracy — the
engine surface analog), two linear-in-C conditions, and det[C_x,C_y,A,B] ≡ c ≠ 0. Fibers =
intersections of two plane curves: Bezout allows degree 4 — a route to a NEW species (which
THM-1330 makes automatically irreducible). FIRST RESULTS (frozen out): the shifted double
cusp-engine A = (u³, 3xu², −x³, 0), B = (0, u³, 3xu², −x³) PASSES the pencil degeneracy
(z², zw, w² ≡ 0), as does the monomial-partner variant; but on the deg-≤4 C-ansatz the two
linear conditions leave a 4-dimensional solution space on which det[C_x,C_y,A,B] is
IDENTICALLY ZERO — **no Keller map from this ansatz (box-bounded negative, deg ≤ 4)**. The
n=4 hunt is now a structured program: find (A,B) pencils that are degenerate WITHOUT killing
the C-determinant — the first named obstruction.

## 6. Referee

Frozen out: the det split on F; A = u³γ(x/u) exact; the cusp-at-infinity chart identity;
the λ = −2, 3 unification; the n=4 six-condition framework; candidates 1–2 passing the
pencil block; the 214×60 exact linear system, nullspace dim 4, and the identically-zero
determinant form.
