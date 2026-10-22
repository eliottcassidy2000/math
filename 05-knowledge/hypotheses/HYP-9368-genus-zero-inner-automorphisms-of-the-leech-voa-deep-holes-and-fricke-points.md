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
  - Order 4, the first composite order: a Fricke point is y/4 with y of norm 24 in a type-4 class mod 2 Lambda and
    max_{Lambda_4} y.mu <= 7. order4_fricke_census.py found 2046 such samples, all with a single orbit invariant,
    that of the E6^4 shallow multiple 3c (= the A2^12 (h+1)-point). This matches the predicted single class
    (VERIFIED by sampling).
  - Prime order p, classification (PROVED as a sketch; the argument is below). Up to Co_0 x Lambda, the Fricke
    points of order p correspond bijectively to the Niemeier lattices with h = p - 1. The ingredients are:
    - (i) the orbifold at a Fricke point is V_N' with dim V_1 = 24p, so h(N') = p - 1;
    - (ii) the dual automorphism of V_N' is inner, e^{2 pi i u(0)} with u in (1/p)N', and regular, because its
      fixed lattice is rootless;
    - (iii) Kac: a regular element of order h + 1 has Kac coordinates (2 at one special node; 1 elsewhere), and the
      special nodes form one Omega-orbit;
    - (iv) gcd(h + 1, |P/N'|) = 1, so u is congruent to w(rho'/p) mod N' with w in W;
    - (v) orbifold duality.
  OPEN:
  - Composite orders. Steps (ii)-(v) work for any order. The missing step is dim V'_1 = 24e at a Fricke point of
    composite order e. The twisted sectors whose twist is a deep hole need the hole vertices and the g-action, and
    the computations agree with 24e in all cases tested.
  - A fully referenced write-up of (i)-(v).
  History: a first draft (S53, before fricke_h_plus_one.py) restricted the Fricke points to the eleven bad-prime
  multiples of deep holes, with orders 3, 4, 5, 6, 7, 9, 10, 11, 15, 23. The (h+1)-points REFUTED that version, by
  adding orders 8, 13, 14, 17, 19, 26, 31 and 47, before it was registered.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_h_plus_one.py (+ _q2.out, .out)
  - 04-computation/experiments/deep_holes_theta_20261022/leech_h_plus_one.py (+ .out), weyl_h_plus_one.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/order3_fricke_census.py (+ .out), order4_fricke_census.py (+ .out)
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

**Geometric form.** Let v be a Fricke point of order n, with nearest Leech point λ. Put c = v + (v − λ)/(n − 1), the
extrapolation away from λ. Then (n − 1)c ∈ Λ and |c − λ|² = 2 + 2/(n − 1). The conjecture says that c is a deep hole
of a type with Coxeter number h = n − 1, and so v = c + (λ − c)/(h + 1) is its (h+1)-point (THM-4759 (10),
fricke_nudge.py).

Equivalently, the 24 norm-2 vectors μ_i of Λ + v satisfy (v − μ_i)·(v − λ) = −1/n, and they extend to the 24 + k
vertices of c. In short: every Fricke point is a nudged deep hole.

## Predictions

- No Fricke point of order 2, 12, 16, 18, 20, 21, 22, 24, 25, 27–30, 32–46, or any order ≥ 48.
- Exactly two classes of order 7, 11, 13, 19 and 31, and one class for each order 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 23,
  26 and 47.
- Each Fricke class is determined by its orbifold neighbour N. For example, the D16E8 and E8³ shallow multiples of
  order 15 both have neighbour D8³, so they should be equivalent; they do have equal theta series of Λ + v
  (fricke_multiples.py).

## The classification at prime order (sketch)

Let v be a Fricke point of prime order p, and g = e^{2πi v(0)}.

1. **Dimension.** Since gcd(i, p) = 1 for i ≠ 0, every twisted sector V_Λ(g^i) is a Galois and modular transform of
   the g-twisted one. Its g-invariant part has character U_p F = Σ c_F(pn) qⁿ, whose constant term is c_F(0) = 24.
   So V′ = V_Λ^{orb(g)} has dim V′_1 = 24 + 24(p − 1) = 24p.
2. **The orbifold is a lattice VOA.** The orbifold of a lattice VOA by an inner automorphism is the lattice VOA of
   the Kneser neighbour, so V′ = V_{N′} with N′ Niemeier and |R(N′)| = 24(p − 1). Hence h(N′) = p − 1, which forces
   p ∈ {3, 5, 7, 11, 13, 17, 19, 23, 31, 47}.
3. **The dual automorphism.** The dual automorphism τ of V_{N′} fixes V_{Λ_0} ⊃ Heisenberg, so τ = e^{2πiu(0)} with
   u ∈ (1/p)N′. Its fixed lattice is Λ_0, which is rootless, so exp(2πiu) is regular.
4. **Kac coordinates.** In each component, a regular element of order p = h + 1 lies in the open alcove with Kac
   coordinates s_i ≥ 1 and Σ a_i s_i = h + 1. So exactly one node has s_i = 2, and it has label a_i = 1, i.e. it is
   special. Ω permutes the special nodes transitively. Hence u_i ∈ W·(ρ_i/p) + P_i.
5. **Glue.** u − w(ρ′/p) lies in P(X′) and p times it lies in N′. Since gcd(p, |P/N′|) = 1, it lies in N′. So
   u ≡ w(ρ′/p) mod N′.
6. **Duality.** Orbifold duality recovers (V_Λ, g) from (V_{N′}, τ). Conjugacy of inner automorphisms of V_Λ is
   Co_0 ⋉ Λ-equivalence of v.

So the order-p Fricke classes are indexed by the Niemeier lattices with h = p − 1:
- two classes each at p = 7, 11, 13, 19 and 31;
- one class each at p = 3, 5, 17, 23 and 47.

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
