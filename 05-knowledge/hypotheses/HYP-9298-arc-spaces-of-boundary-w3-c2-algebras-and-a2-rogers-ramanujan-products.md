---
id: HYP-9298
title: "Arc-space A_2 Rogers-Ramanujan identities: with wt x = 2, wt y = 3, the arc algebra of the C_2-algebra R_q of the boundary W_3 minimal model W_3(3,q) (q >= 7, 3 not dividing q) has Hilbert series prod_{n != 0,+-1 mod q} (1-t^n)^{-1} * prod_{n != 0,+-1,+-2 mod q} (1-t^n)^{-1}; explicitly HS J_infty(C[x,y]/(xy, y^2 - x^3)) = (mod 7 product) and HS J_infty(C[x,y]/(y^2 - x^3, x^2 y)) = (mod 8 product)"
status: >
  OPEN (conjecture). Evidence is FINITE-EXACT and pure commutative algebra (arc_hs.py, independent of the vertex-algebra
  code; ranks mod 2^31-1). The arc Hilbert series equals the product:
  - through weight 26 for q = 7 (ideal (xy, y^2 - a x^3), any a != 0 up to rescaling);
  - through weight 26 for q = 8 ((y^2 - a x^3, x^2 y));
  - through weight 24 for q = 10 ((x y^2 - b x^4, y^3 - c x^3 y));
  - through weight 24 for q = 11 ((y^3 - d x^3 y, x^2 y^2 - e x^5)).
  The ideals are the Li-degree-0 parts of grI computed by the singular-support engine (w3_c2_relations.py); the
  coefficients for q = 10, 11 are the engine's exact rationals.
  Vertex-algebra side: the W_3(3,q) vacuum characters equal the same products through weights 20-24
  (plethystic-exponent check of the engine's dimensions, which agree with the assumption-free maximal-submodule
  computation). So the conjecture is equivalent to classical freeness of W_3(3,q) (HYP-9297, N = 3) together with
  the product character formula. Product formulas for boundary W_3 characters are standard in spirit
  (A_2 Rogers-Ramanujan identities, Andrews-Schilling-Warnaar 1999; Feigin-Foda-Welsh 2008) but were not
  re-verified here.
  Literature (2026-10-18): Bruschek-Mourtada-Schepers and Mourtada's survey (Bull. AMS 62, 2025) treat C[x]/(x^k)
  (Andrews-Gordon) and the node C[x,y]/(xy). No two-variable weighted fat point with an A_2 product was found.
session: kind-pasteur-2026-10-18-S1
source: 05-knowledge/results/singular_support_arc_spaces_classical_freeness_20261018.md
scripts:
  - 04-computation/voa_singular_support/arc_hs.py (commutative algebra), w3_c2_relations.py (C_2 ideals), run_w3.py
related: [HYP-9297, THM-4671, Bruschek-Mourtada-Schepers 2011, Mourtada (Bull. AMS 2025), Andrews-Schilling-Warnaar 1999]
---

# HYP-9298 — arc spaces of the boundary `W_3` `C_2`-algebras and `A_2` products

## Statement

Grade `C[x, y]` by `wt x = 2` and `wt y = 3`. Give the jets `x_j, y_j` weights `2 + j` and `3 + j`. Let
```
P_q(t) = prod_(n ≢ 0, ±1 (mod q)) (1 - t^n)^(-1)  ·  prod_(n ≢ 0, ±1, ±2 (mod q)) (1 - t^n)^(-1).
```
Let `R_q` be Zhu's `C_2`-algebra of the boundary `W_3` minimal model `W_3(3, q)`, the chiral algebra of the 4d
Argyres-Douglas theory `(A_2, A_(q-4))`. Then `HS(J_infty Spec R_q)(t) = P_q(t)`.

Explicit cases, up to weighted rescaling of `x` and `y`:
- `q = 7`: `R_7 = C[x,y]/(xy, y^2 - x^3)`.
- `q = 8`: `R_8 = C[x,y]/(y^2 - x^3, x^2 y)`.
- `q = 10`: `R_10 = C[x,y]/(x y^2 - (2/123) x^4, y^3 - (8/41) x^3 y)`.
- `q = 11`: `R_11 = C[x,y]/(y^3 - (11/138) x^3 y, x^2 y^2 - (11/1035) x^5)`.

The `q = 10` and `q = 11` coefficients are as found. A weighted rescaling can normalize two of them.

## Context

- **`A_1` case.** It is the theorem of Bruschek-Mourtada-Schepers: `HS(J_infty C[x]/(x^k))` with `wt x = 2` is the
  Andrews-Gordon product, which is the `Vir_(2,2k+1)` character.
- **`A_2` case.** This is the `A_2` analogue. The modulus `q` and the two-factor product mirror the `A_2`
  Rogers-Ramanujan identities.
- **Vertex-algebra reading.** A proof would give classical freeness of every boundary `W_3` model, via the known
  character, and hence the Schur indices of `(A_2, A_(q-4))` as arc-space Hilbert series.
- **Curiosity (OBSERVATION).** The cusp alone has
  `HS J_infty(C[x,y]/(y^2 - x^3)) = prod_(n>=2)(1-t^n)^(-1) · 1/((1-t^3)(1-t^4)(1-t^5))` through weight 24.

## Modulus rigidity (FINITE-EXACT, 2026-10-18)

For `q >= 10` the `C_2`-ideal carries a modulus that weighted rescaling cannot remove.
- `q = 10`: `R = (x y^2 - beta x^4, y^3 - gamma x^3 y)`, where only `beta/gamma` matters. The vertex algebra gives
  `beta/gamma = 1/12`.
- `q = 11`: `R = (y^3 - d x^3 y, x^2 y^2 - e x^5)`, where only `d/e` matters. The vertex algebra gives `d/e = 15/2`.

The arc Hilbert series equals the product **only at the vertex-algebra value**, through weight 22.

`q = 10`:
- ratio `1/12`: equal;
- ratio `1`: first differs at weight 11 (45 vs 44);
- ratios `1/6` and `-1/12`: first differ at 12 (65 vs 66).

`q = 11`:
- ratio `15/2`: equal;
- ratio `1`: first differs at 13 (93 vs 92);
- ratios `7` and `8`: first differ at 14 (131 vs 132).

**Reading.** The rank of the arc ideal is lower semicontinuous, so the generic modulus has the smallest arc Hilbert
series. That generic value falls below the product. The vertex-algebra modulus is a special jump point where the arc
space is just large enough to carry the character.

So classical freeness of `W_3(3,q)` is not a property of the shape of the `C_2`-ideal. It pins the modulus. This is a
rigidity phenomenon with no analogue in the one-variable case, where `C[x]/(x^k)` has no modulus.

## Tests

- **Push and extend.** `q = 7, 8` to weight 34; `q = 13, 14` (needs the engine's `C_2` ideals).
- **Hostile control, DONE.** The monomial degenerations fail:
  - `(xy, y^2)` for `q = 7` first differs at weight 8 (12 vs 11);
  - `(y^2, x^2 y)` for `q = 8` first differs at weight 10 (27 vs 26).

  `(xy, y^2 + x^3)`, in the rescaling orbit, equals the product. So the identity needs the genuine cusp relation
  `y^2 ~ x^3`, i.e. the vertex-algebra ideal.
