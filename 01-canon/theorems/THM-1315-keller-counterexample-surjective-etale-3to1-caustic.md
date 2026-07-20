---
id: THM-1315
title: THE KELLER COUNTEREXAMPLE IS A SURJECTIVE, EVERYWHERE-ÉTALE, GENERICALLY 3-TO-1 SELF-MAP OF ℂ³ — with fiber cubic φ = −2y³ + 3by² − 18ay + (18ab − b³ − 27a²c), caustic Δ = −2²3⁶·a²·K, and ramification entirely at infinity (Res ∝ the same K, forced by étaleness)
status: PROVED (syzygy reduction + explicit fiber cubic, symbolic; surjectivity by complete case analysis with referee probes; S₃ monodromy by certified specialization + van der Waerden — the independent syzygy-route confirmation of the THM-1310 target; all computations exact sympy, frozen out)
source: mac-mini-2026-07-19-S131 (owner: work the S3 pin and the curve)
depends_on: [the verified counterexample (HYP-8070 kps / HYP-8095 mac-mini / HYP-8075 death-star), THM-1300 (context)]
related: [THM-1305, THM-1310 (kps conic route, in progress), opus-S414 moduli lane, HYP-8085 (collision family), HYP-8100]
script: 04-computation/jc_syzygy_s3_asymptotic referee blocks (inline in the frozen out)
output: 05-knowledge/results/jc_syzygy_s3_asymptotic_macmini_S131.out
---

# THM-1315 — the counterexample's global geometry

Let F = (F₁,F₂,F₃): ℂ³ → ℂ³ be the verified Keller counterexample (det JF ≡ −2).

## 1. The syzygies and the fiber cubic

With u = 1+xy and G = u²z + y²(1+3u): F₁ = uG, F₂ = y + 3xG, and

- **Syzygy 1:** 3x·F₁ ≡ (1+xy)(F₂ − y)  (identically in ℂ[x,y,z]);
- **Syzygy 2:** x³·G + u²·F₃ ≡ x(xy+2)  (z-free).

Hence over a target (a,b,c) with a ≠ 0 the fiber solves x = (b−y)/(3a − y(b−y)) and one
explicit cubic

    φ_{a,b,c}(y) = −2y³ + 3by² − 18ay + (18ab − b³ − 27a²c) = 0,

whose leading coefficient −2 NEVER vanishes: every fiber cubic has exact degree 3. The
parameter c enters only as the vertical shift −27a²c — the family is a versal cubic family,
which is the structural reason full S₃ monodromy is unavoidable.

## 2. Monodromy = S₃ (independent syzygy-route confirmation)

At (a,b,c) = (1,0,1): φ = −2y³ − 18y − 9·… (see frozen out) is irreducible over ℚ with
discriminant −125388 < 0 (a fortiori non-square) ⟹ Galois group S₃ at the specialization;
by van der Waerden specialization (specialized ≤ generic ≤ S₃), the generic monodromy is S₃.
Certified at four further specializations. (kps-c98's THM-1310 conic-pair route is the
first-claimed lane; this route is its independent confirmation — double-construction norm.)

## 3. The caustic, and ramification-at-infinity coherence

    Δ(a,b,c) = −2916·a²·K,   K = 27a²c² − 18abc + b³c + 16a − b²   (−2916 = −2²·3⁶),
    Res_y(φ, 3a − y(b−y)) = 27·a²·K  — THE SAME K.

Since F is étale (det ≡ −2), no finite ramification is possible; coherently, wherever φ has
a double root (Δ = 0) that root simultaneously lies on the x-denominator conic (Res = 0):
**the double root always escapes to x = ∞**. The would-be branch surface {K = 0} is the
non-properness (Jelonek-type) locus: over it the finite fiber count drops 3 → 2 by a sheet
leaving through infinity, never by two finite points colliding.

## 4. Surjectivity (complete case analysis)

- a ≠ 0, c ≠ 0: φ and the denominator conic share at most ONE root (the matching
  φ = −2(y²−by+3a)(y−r) has NO solution (r,c) — verified symbolically), and y = b (the x=0
  root) satisfies φ(b) = −27a²c ≠ 0; hence ≥ 2 honest finite preimages.
- a ≠ 0, c = 0: the plane x = 0 maps onto {c = 0} as (y,z) ↦ (z+4y², y, 0) — covered.
- a = 0, b ≠ 0: φ(0,b,c) = −(y−b)²(2y+b); the root y = −b/2 gives x = 2/b, z determined.
- a = 0, b = 0: (c/2, 0, 0) ↦ (0,0,c) directly.

**F(ℂ³) = ℂ³.** So the first Keller counterexample is a surjective, étale, generically
3-to-1 endomorphism of affine 3-space — non-injective with no finite branching and no
omitted values: all its failure of invertibility is carried through infinity (and, per the
2-adic ladder of THM-1300's formal inverse, through unbounded dyadic denominators).

## 5. Referee

Frozen out `jc_syzygy_s3_asymptotic_macmini_S131.out`: syzygies (identically zero), the
cubic and its degree/leading coefficient, five S₃ specializations, Δ and Res factored, the
unsolvable degenerate matching, and exact Gröbner probes at (0,1,1), (0,0,1), (0,1,0),
(0,0,0) — 1, 1, 3, 1 preimages respectively, matching the case analysis.
