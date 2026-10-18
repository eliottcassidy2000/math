---
id: HYP-9298
title: "Arc-space A_2 Rogers-Ramanujan identities: with wt x = 2, wt y = 3, the arc algebra of the C_2-algebra R_q of the boundary W_3 minimal model W_3(3,q) (q >= 7, 3 not dividing q; conjecturally the (A_2,A_{q-4}) Argyres-Douglas chiral algebra, Cordova-Shao) has Hilbert series prod_{n != 0,+-1 mod q} (1-t^n)^{-1} * prod_{n != 0,+-1,+-2 mod q} (1-t^n)^{-1}; explicitly HS J_infty(C[x,y]/(xy, y^2 - x^3)) = (mod 7 product) and HS J_infty(C[x,y]/(y^2 - x^3, x^2 y)) = (mod 8 product)"
status: >
  OPEN (conjecture). Evidence is FINITE-EXACT and pure commutative algebra (arc_hs.py, independent of the vertex-algebra
  code; ranks mod 2^31-1, second prime 2^31-19 where stated). The arc Hilbert series equals the product:
  - through weight 26 for q = 7 (ideal (xy, y^2 - a x^3), any a != 0 up to rescaling); the independent referee
    pushes it to weight 34 (mod 2^31-1) and 30 (mod 2^31-19);
  - through weight 26 for q = 8 ((y^2 - a x^3, x^2 y)); the referee gives 34 and 30;
  - through weight 24 for q = 10 ((x y^2 - b x^4, y^3 - c x^3 y));
  - through weight 24 for q = 11 ((y^3 - d x^3 y, x^2 y^2 - e x^5));
  - through weight 22 for q = 13, at both primes, for
    (x y^3 - (13/448) x^4 y, y^4 - (13/56) x^3 y^2 + (169/150528) x^6). There are two moduli: a/b = 1/8 and
    c/b^2 = 1/48 in the notation (x y^3 - a x^4 y, y^4 - b x^3 y^2 + c x^6).
    CORRECTION (audit 2, MISTAKE-633): an earlier version stated the x^6 coefficient as -28529/3278. That value was a
    single-prime rational-reconstruction artefact, and the ideal as written fails at weight 16 mod 2^31-19.
    Perturbing the x^6 coefficient (e.g. to -1) breaks the identity at weight 16.
  The ideals are the Li-degree-0 parts of grI computed by the singular-support engine (w3_c2_relations.py). All
  coefficients have been reconstructed by CRT from two primes.
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
Let `R_q` be Zhu's `C_2`-algebra of the boundary `W_3` minimal model `W_3(3, q)`. This is conjecturally the chiral
algebra of the 4d Argyres-Douglas theory `(A_2, A_(q-4))` (Cordova-Shao conjecture). Then
`HS(J_infty Spec R_q)(t) = P_q(t)`.

Explicit cases, up to weighted rescaling of `x` and `y`:
- `q = 7`: `R_7 = C[x,y]/(xy, y^2 - x^3)`.
- `q = 8`: `R_8 = C[x,y]/(y^2 - x^3, x^2 y)`.
- `q = 10`: `R_10 = C[x,y]/(x y^2 - (2/123) x^4, y^3 - (8/41) x^3 y)`.
- `q = 11`: `R_11 = C[x,y]/(y^3 - (11/138) x^3 y, x^2 y^2 - (11/1035) x^5)`.
- `q = 13`: `R_13 = C[x,y]/(x y^3 - (13/448) x^4 y, y^4 - (13/56) x^3 y^2 + (169/150528) x^6)`.

The coefficients are as found and confirmed at two primes. A weighted rescaling acts on them through a single
character, so it can normalize **one** coefficient. One modulus remains for `q = 10, 11` and two for `q = 13`.

## Context

- **`A_1` case.** It is the theorem of Bruschek-Mourtada-Schepers: `HS(J_infty C[x]/(x^k))` with `wt x = 2` is the
  Andrews-Gordon product, which is the `Vir_(2,2k+1)` character.
- **`A_2` case.** This is the `A_2` analogue. The modulus `q` and the two-factor product mirror the `A_2`
  Rogers-Ramanujan identities.
- **Vertex-algebra reading.** A proof would give classical freeness of every boundary `W_3` model, via the known
  character, and hence the Schur indices of `(A_2, A_(q-4))` as arc-space Hilbert series.
- **Curiosity (OBSERVATION, corrected).** For the cusp alone,
  `HS J_infty(C[x,y]/(y^2 - x^3))` equals `prod_(n>=2)(1-t^n)^(-1) · 1/((1-t^3)(1-t^4)(1-t^5))` through weight 33 and
  **fails at weight 34** (13436 vs 13435; block weight 34, Li 11: 56 vs 55). Audit 2 established this exactly over Z.
  - The candidate formula is the complete-intersection Euler characteristic of the arc ideal.
  - Failure is expected: the cusp is not log canonical (lct = 5/6), so its jet schemes have excess dimension.
  - An earlier version claimed the formula through weight 24 only, as an unqualified identity.

## Modulus rigidity (FINITE-EXACT, 2026-10-18)

For `q >= 10` the `C_2`-ideal carries a modulus that weighted rescaling cannot remove.
- `q = 10`: `R = (x y^2 - beta x^4, y^3 - gamma x^3 y)`, where only `beta/gamma` matters. The vertex algebra gives
  `beta/gamma = 1/12`.
- `q = 11`: `R = (y^3 - d x^3 y, x^2 y^2 - e x^5)`, where only `d/e` matters. The vertex algebra gives `d/e = 15/2`.

The arc Hilbert series equals the product **only at the vertex-algebra value**, through weight 22. The independent
referee computed the full rank-drop locus at the first separating weight (`modulus_locus.py`, FINITE-EXACT mod p):
- `q = 10`, weight 12: the locus is exactly `{1/12, 1}`;
- `q = 11`, weight 14: the locus is exactly `{15/2, 1}`.

Every other ratio is generic, including `∞`. Generic ratios give the same Hilbert series and fall below the product:
65 vs 66 at weight 12 for `q = 10`, and 131 vs 132 at weight 14 for `q = 11`.

Ratio 1 is a different special point. It is not a fat point:
`(x y^2 - beta x^4, y^3 - beta x^3 y) = (x, y)·(y^2 - beta x^3)` contains a cusp curve. That is why its arc Hilbert
series lies above the product (45 vs 44 at weight 11 for `q = 10`; 93 vs 92 at weight 13 for `q = 11`).

**Reading.** The rank of the arc ideal is lower semicontinuous, so generic moduli give the smallest arc Hilbert
series. Among the fat points of the family, the vertex-algebra modulus is the unique jump point where the arc space
is just large enough to carry the character.

So classical freeness of `W_3(3,q)` is not a property of the shape of the `C_2`-ideal. It pins the modulus. This is a
rigidity phenomenon with no analogue in the one-variable case, where `C[x]/(x^k)` has no modulus.

## Derived reading (OBSERVATION, elementary)

In every case `R_q` is a complete intersection, with generators in weights `q-2` and `q-1`. So the Koszul complex of
the arc ideal (the derived arc space) has graded Euler characteristic
```
1 / ( prod_(m=2)^(q-3) (1-t^m) · prod_(m=3)^(q-2) (1-t^m) ).
```
The factor multiset is `{2..q-2} + {3..q-3}`, which is exactly the `n < q` truncation of `P_q`. The `A_1` analogue
holds too: the Euler characteristic of the derived arc space of `x^N` is the `n < 2N+1` truncation of `AG_N`. So the
conjecture says that `H_0` alone carries the full product, while the higher Koszul homology accounts for the factors
`n >= q`.

## Proof route (lead, not a result)

Bruschek-Mourtada-Schepers count the standard monomials of an initial ideal of the arc ideal, as partitions with
difference conditions. `arc_initial_ideal.py` computes, for `q = 7`, the minimal generators of the initial ideal under
lex with the largest variable of highest weight (`out/arc_initial_ideal_q7_rwtx_W15.txt`). There are 49 generators up
to weight 15. They begin with:
- the infinite families `x_0 y_j` and `y_0 y_j`;
- `x_0^4` and `x_0^3 x_j`;
- `x_1 y_1^2`, `y_1^3` and further degree-3 families.

Turning these into a two-colour partition rule whose generating function is `P_7` would prove the `q = 7` case. That
would be an `A_2` analogue of BMS.

## Tests

- **Push and extend.** `q = 7, 8` reached weight 34 (audit 2). `q = 14` needs the engine's `C_2` ideals, which should
  be reconstructed from two primes.
- **Hostile control, DONE.** The monomial degenerations fail:
  - `(xy, y^2)` for `q = 7` first differs at weight 8 (12 vs 11);
  - `(y^2, x^2 y)` for `q = 8` first differs at weight 10 (27 vs 26).

  `(xy, y^2 + x^3)`, in the rescaling orbit, equals the product. So the identity needs the genuine cusp relation
  `y^2 ~ x^3`, i.e. the vertex-algebra ideal.
