---
id: HYP-9368
title: "The genus-zero inner automorphisms of the Leech lattice VOA are classified by the Niemeier lattices, twice: for v in Lambda (x) Q, Tr_{V_Lambda}(e^{2 pi i v(0)} q^{L_0 - 1}) - 24 is a normalised Hauptmodul iff, up to Co_0 (x) Lambda, v = 0, v is a deep hole (23 classes; non-Fricke umbral Hauptmoduln, THM-4759 (3)), or v is a Fricke point c''_N (23 classes, one per Niemeier lattice N, of order h(N) + 1; Fricke Hauptmoduln, THM-4759 (10)). A Fricke point is a point of order n with a unique nearest Leech point at squared distance 2 - 2/n whose proper multiples are deep holes or lie in Lambda; the Fricke points are exactly the c''_N, and they are classified by their orbifold neighbour N"
status: >
  OPEN (conjecture). The 'if' directions are PROVED in THM-4759 (3), (9) and (10). Evidence for 'only if':
  - Order 2: no Fricke point exists, since v = -v mod Lambda pairs the minimal vectors (PROVED). This matches the
    absence of a Niemeier lattice with h = 1.
  - Order 3: a Fricke point is v = y/3 with y a norm-12 vector having no lambda in Lambda_4 at y.lambda = 6.
    order3_fricke_census.py found two orbit invariants among 2977 sampled norm-12 vectors. The k = 0 class is that
    of the M24-fixed vector (2^24)/sqrt 8, and the two orbit sizes |Co_0|/|M24| + |Co_0|/|U6(2):2| account exactly
    for all 34,417,656,000 norm-12 vectors. So, given the classical two-orbit decomposition, there is exactly one
    order-3 Fricke class: the (h+1)-point of A1^24 (= the D4^6 shallow multiple 2c). This matches the prediction
    (one Niemeier lattice with h = 2) (VERIFIED).
  - Prime order p: the Z/p-orbifold at a Fricke point has dim V_1 = 24p, because each twisted sector contributes the
    constant term 24 of U_p F. So it is a Niemeier lattice with Coxeter number p - 1. Fricke points of prime order p
    therefore exist only for p - 1 in {2, 4, 6, 10, 12, 16, 18, 22, 30, 46}; all of these are realised (PROVED
    argument).
  OPEN: the injectivity statement, that two Fricke points with the same orbifold neighbour N are equivalent; and the
  composite orders, where (h+1)-points already realise every Fricke level h + 1. A natural route is the analogue for
  Fricke points of Moller-Scheithauer's uniqueness of the inverse orbifold automorphism.
  History: a first draft (S53, before fricke_h_plus_one.py) restricted the Fricke points to the eleven bad-prime
  multiples of deep holes, with orders 3, 4, 5, 6, 7, 9, 10, 11, 15, 23. The (h+1)-points REFUTED that version, by
  adding orders 8, 13, 14, 17, 19, 26, 31 and 47, before it was registered.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_h_plus_one.py (+ _q2.out, .out)
  - 04-computation/experiments/deep_holes_theta_20261022/leech_h_plus_one.py (+ .out), weyl_h_plus_one.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/order3_fricke_census.py (+ .out)
related:
  - THM-4759 (deep-hole twining, Fricke points, Fricke descent)
  - THM-4758 (swap theorem)
---

# HYP-9368 — genus-zero inner automorphisms of V_Λ

## Statement

Let v ∈ Λ ⊗ Q have order n modulo Λ. Write Z_v = Tr_{V_Λ}(e^{2πi v(0)} q^{L_0 − 1}) = Θ_Λ^{[v]}/η²⁴.

**Conjecture.** Z_v − 24 is the normalized Hauptmodul of a genus-zero group if and only if, up to Co_0 ⋉ Λ, one of
the following holds:
1. v = 0, giving J;
2. v is a deep hole of type X (23 classes), giving T^X + k with the umbral lambency group Γ^X (THM-4759 (3));
3. v is a Fricke point c″_N (23 classes), of order h(N) + 1, giving the Fricke Hauptmodul of Γ_0(h+1) + (h+1)
   (THM-4759 (10)).

Equivalently, the Fricke points of Λ are exactly the 23 points c″_N. Here a Fricke point is a point with a unique
nearest Leech point at squared distance 2 − 2/n, all of whose proper multiples are deep holes or lie in Λ.

## Predictions

- No Fricke point of order 2, 12, 16, 18, 20, 21, 22, 24, 25, 27–30, 32–46, or any order ≥ 48.
- Exactly two classes of order 7, 11, 13, 19 and 31, and one class for each order 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 23,
  26 and 47.
- Each Fricke class is determined by its orbifold neighbour N. For example, the D16E8 and E8³ shallow multiples of
  order 15 both have neighbour D8³, so they should be equivalent; they do have equal theta series of Λ + v
  (fricke_multiples.py).

## Why the "only if" is plausible

Suppose v is neither in Λ nor a deep hole. Then Λ + v has vectors of norm < 2, which gives Z_v a pole at the cusp 0.

A Hauptmodul has a single pole on its curve, so the invariance group Γ ⊇ Γ_0(n) must identify the cusps ∞ and 0.

If Γ is contained in the normaliser of Γ_0(n), the element doing this is of Fricke type. Then the pole at 0 must be
simple in the local parameter of width n, so the minimal norm is exactly 2 − 2/n with a single minimal vector. There
must be no further poles, so every proper multiple of v is a deep hole or lies in Λ.

If instead Γ contains SL_2(Z)-elements outside the normaliser, then Z_v = J + const. That forces θ_{Λ+v} = Δ(J + c),
so v ∈ Λ.

A complete argument also has to handle the remaining intermediate groups. Granting that, the "only if" reduces to the
geometric classification of Fricke points.
