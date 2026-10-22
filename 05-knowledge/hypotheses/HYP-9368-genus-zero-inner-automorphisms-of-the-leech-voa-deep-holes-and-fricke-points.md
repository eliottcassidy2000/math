---
id: HYP-9368
title: "Fricke points of the Leech lattice and genus-zero inner automorphisms of V_Lambda. (A) The strict Fricke points (order n, a unique nearest Leech point at squared distance 2 - 2/n, all proper multiples deep holes or in Lambda) are exactly the 23 points c''_N of THM-4759 (10), one per Niemeier lattice N, of order h(N) + 1, classified by their orbifold neighbour N. (B) The genus-zero inner automorphisms (Tr_{V_Lambda}(e^{2 pi i v(0)} q^{L_0 - 1}) - 24 a normalised Hauptmodul) are 0, the 23 deep holes, the strict Fricke points, and a further family of Fricke-type points whose shallow proper multiples cancel in the twining. Its known members are the shell-1 points of the D-type Niemeier lattices (A1^24 = D2^12, D4^6, D6^4, D8^3, D10E7^2, D12^2, D16E8; orders 4, 8, 12, 16, 20, 24, 32), which twine by 4A, 8A, 12H, 16C, ...; the full classification of (B) is open"
status: >
  OPEN (conjecture), with the 'if' directions PROVED in THM-4759 (3), (9), (10).
  Part (A), evidence:
  - Order 2: no strict Fricke point exists, since v = -v mod Lambda pairs the minimal vectors (PROVED).
  - Order 3: a Fricke point is y/3 with y a norm-12 vector having no lambda in Lambda_4 at y.lambda = 6.
    order3_fricke_census.py finds two orbit invariants among 2977 sampled norm-12 vectors; the k = 0 one is that of the
    M24-fixed vector (2^24)/sqrt 8. The orbit sizes |Co_0|/|M24| + |Co_0|/|U6(2):2| account exactly for all
    34,417,656,000 norm-12 vectors. So there is exactly one order-3 class, given the classical two-orbit decomposition
    (VERIFIED).
  - Order 4 (2v deep): order4_fricke_census.py finds 2046 samples, all with the invariant of the A2^12 (h+1)-point
    (VERIFIED by sampling).
  - Order 6 (v = a - b with a an A1^24 hole and b an A2^12 hole): order6_fricke_census.py finds 1399 samples, all with
    6-neighbour A4^6, i.e. dim V'_1 = 24e (VERIFIED by sampling).
  - Prime order p: classification by the Niemeier lattices with h = p - 1 (PROVED as a sketch, below). It uses:
    dim V'_1 = 24p via U_p F; the dual automorphism is inner and regular; Kac coordinates (2;1,...,1);
    gcd(h + 1, |P/N'|) = 1; orbifold duality.
  - Composite orders: the missing step is dim V'_1 = 24e. It holds in all samples.
  Part (B), REFUTATION OF THE FIRST REGISTERED VERSION, and the new family. The first registered version said
  'Z_v - 24 Hauptmodul iff 0 / deep hole / strict Fricke point', arguing that the absence of further poles forces every
  proper multiple to be deep. That is wrong: the twined contributions of shallow vectors at an intermediate cusp can
  cancel. swap_shell_test.py and the shell-1 Fricke test (fricke_shell_test.py) find the following.
  - The Leech points c_w of the shell-1 cusps (m = h + 2) all have a unique minimal vector at the Fricke norm
    2 - 2/(h+2).
  - For the seven D-type lattices the twining is Fricke-invariant through q^2: order 4 (A1^24), 276q + 2048q^2 = 4A;
    order 8 (D4^6), 36q + 128q^2 = 8A; order 12 (D6^4), 14q + 36q^2 = 12H (A058486); order 16 (D8^3), 8q + 16q^2 = 16C
    (A058516); order 20 (D10E7^2), 5q + 10q^2; order 24 (D12^2), 4q + 6q^2; order 32 (D16E8), 2q + 4q^2. In each case
    some proper multiple is shallow.
  - The other sixteen types fail Fricke invariance, D24 (order 48, X_0^+(48) not genus zero) among them.
  - Their orbifold neighbour is N itself, of Coxeter number e - 2. So at orders 4, 8, 12, 16 (20, 24, 32) there are
    inner automorphisms with Fricke Hauptmodul twining and dim V'_1 = 24(e - 1), not 24e.
  - Confirmation through q^4 for A1^24 and D4^6 is running (fricke_shell_test_q4.out).
  OPEN:
  - The complete list in (B).
  - A conceptual reason for the D-type family (the shell-1 cusps correspond to Kac type (3;1,...,1)-like points; why D?).
  - A referenced write-up of the prime-order sketch.
  History:
  - The first S53 draft allowed only the eleven bad-prime multiples; the (h+1)-points refuted it.
  - The second version (registered as c31e9a6a3f..23dc8edd9e) claimed the 'only if' for all orders; the shell-1
    D-type points refuted it. Both refutations happened within the session.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_h_plus_one.py (+ _q2.out, .out)
  - 04-computation/experiments/deep_holes_theta_20261022/leech_h_plus_one.py (+ .out), weyl_h_plus_one.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/order3_fricke_census.py (+ .out), order4_fricke_census.py (+ .out), order6_fricke_census.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_shell_test.py (+ .out) (the shell-1 points; the (B) family)
related:
  - THM-4759 (deep-hole twining, Fricke points, Fricke descent)
  - THM-4758 (swap theorem)
---

# HYP-9368 — Fricke points and genus-zero inner automorphisms of V_Λ

## Statement

Let v ∈ Λ ⊗ Q have order n modulo Λ. Write Z_v = Tr_{V_Λ}(e^{2πi v(0)} q^{L_0 − 1}) = Θ_Λ^{[v]}/η²⁴.

Call v a **strict Fricke point** if it satisfies both of these:
- Λ + v has a unique vector of norm 2 − 2/n;
- every proper multiple jv (j | n, 1 < j < n) is a deep hole or lies in Λ.

THM-4759 (10) proves that then Z_v − 24 is the Fricke Hauptmodul of Γ_0(n) + n, and that the 23 points c″_N (orders
h(N) + 1) are strict Fricke points.

**Conjecture (A).** The strict Fricke points are, up to Co_0 ⋉ Λ, exactly the 23 points c″_N. Each is determined by its
orbifold neighbour N, which has h(N) = n − 1.

**Geometric form of (A).** Let v be a strict Fricke point with nearest Leech point λ. Then c = v + (v − λ)/(n − 1)
should be a deep hole of Coxeter number n − 1, and v = c + (λ − c)/n its (h+1)-point: every strict Fricke point is a
nudged deep hole (fricke_nudge.py).

**Conjecture (B), revised.** Z_v − 24 is a normalized Hauptmodul exactly when v is one of the following:
- v = 0;
- a deep hole (23 classes);
- a strict Fricke point;
- a Fricke-type point with shallow but cancelling multiples.

Known members of the last family are the shell-1 points of the seven D-type Niemeier lattices, of orders 4, 8, 12,
16, 20, 24 and 32. The complete list is open.

## The family that refuted the first registered (B)

The cusps w of type N adjacent to z satisfy −w·ρ = (h/2)|y − ρ_N/h|². Shell 0 gives the Fricke cusp z″ (m = h + 1),
and shell 1 gives m = h + 2.

For every type the shell-1 Leech point c_w has a unique minimal vector of norm 2 − 2/(h + 2), the Fricke norm. For the
seven D-type lattices
- A1²⁴ = D2¹² (order 4),
- D4⁶ (8),
- D6⁴ (12),
- D8³ (16),
- D10E7² (20),
- D12² (24),
- D16E8 (32),

the twining Z_{c_w} is Fricke-invariant (exact q-series identity through q²) and matches a Monster Fricke series:
4A, 8A, 12H (OEIS A058486), 16C (A058516), and three more at levels 20, 24 and 32 still to be identified. Yet some proper
multiple of c_w is shallow. So the absence of extra poles comes from cancellation of the twined shallow states, not
from deep multiples.

The orbifold of V_Λ at c_w is V_N, of Coxeter number e − 2, so dim V′_1 = 24(e − 1). These are therefore not strict
Fricke points, and they do not contradict (A).

The other sixteen types fail Fricke invariance at shell 1. This includes D24, of order 48, consistent with
X_0^+(48) not having genus zero.

## Evidence for (A)

- **Order 2.** Impossible: v ≡ −v pairs the minimal vectors.
- **Order 3.** One class, the M24-point (2²⁴)/(3√8). This rests on the two-orbit decomposition of the norm-12 vectors
  (order3_fricke_census.py).
- **Order 4, with 2v deep.** One invariant in 2046 samples, that of the A2¹² point (order4_fricke_census.py).
- **Order 6.** All 1399 samples have neighbour A4⁶ (order6_fricke_census.py).
- **Prime orders.** The sketch below.

## The classification at prime order (sketch)

Let v be a strict Fricke point of prime order p, and g = e^{2πi v(0)}. Prime order has no proper multiples, so (A)
and (B) agree here.

1. **Dimension.** Since gcd(i, p) = 1 for i ≠ 0, every twisted sector V_Λ(g^i) is a Galois and modular transform of
   the g-twisted one. Its g-invariant part has character U_p F, whose constant term is c_F(0) = 24. So V′ = V_Λ^{orb(g)}
   has dim V′_1 = 24p.
2. **The orbifold is a lattice VOA.** V′ = V_{N′} with N′ Niemeier and |R(N′)| = 24(p − 1). Hence h(N′) = p − 1.
3. **The dual automorphism.** τ = e^{2πiu(0)} with u ∈ (1/p)N′. It is regular, since its fixed lattice Λ_0 is
   rootless.
4. **Kac coordinates.** A regular element of order h + 1 has Kac coordinates 2 at one special node and 1 elsewhere.
   Ω permutes the special nodes transitively, so u_i ∈ W·(ρ_i/p) + P_i.
5. **Glue.** gcd(p, |P/N′|) = 1, so u ≡ w(ρ′/p) mod N′.
6. **Duality.** Orbifold duality recovers (V_Λ, g) from (V_{N′}, τ).

So the order-p Fricke classes are indexed by the Niemeier lattices with h = p − 1:
- two classes each at p = 7, 11, 13, 19 and 31;
- one class each at p = 3, 5, 17, 23 and 47.

## Predictions for (A)

- No strict Fricke point of order 2, 12, 16, 18, 20, 21, 22, 24, 25, 27–30, 32–46, or any order ≥ 48. At 12, 16,
  20, 24 and 32 the (B) family supplies genus-zero twinings instead.
- Exactly two classes at orders 7, 11, 13, 19 and 31, and one at each of 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 23, 26
  and 47.
- The D16E8 and E8³ shallow multiples of order 15 are equivalent: both have neighbour D8³ and equal theta series.
