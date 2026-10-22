---
id: HYP-9368
title: "Fricke points of the Leech lattice and genus-zero inner automorphisms of V_Lambda. (A) The strict Fricke points (order n, Z_v - 24 the normalized Gamma_0(n)+n Hauptmodul - by THM-4760 equivalently rational, a unique nearest Leech point at squared distance 2 - 2/n and exactly 24 at squared distance 2 - and all proper multiples deep holes or in Lambda) are exactly the 23 points c''_N of THM-4759 (10), one per Niemeier lattice N, of order h(N) + 1, classified by their orbifold neighbour N. (B) The genus-zero inner automorphisms (Tr_{V_Lambda}(e^{2 pi i v(0)} q^{L_0 - 1}) - 24 a normalised Hauptmodul) are 0, the 23 deep holes, the strict Fricke points, and a further family of Fricke-type points whose shallow proper multiples cancel in the twining. Its known members are seven shell-1 points of Borcherds' 2-targets (A1^24, D4^6, D6^4, D8^3, D10E7^2 with E7 carrier, D12^2, D16E8 with E8 carrier; orders 4, 8, 12, 16, 20, 24, 32), PROVED in THM-4760 to twine by 4A, 8A, 12H, 16C, 20F, 24I, 32A. Conjecture (B'): there are no others, so the Fricke orders realised by translations are exactly 3-17, 19, 20, 23, 24, 26, 31, 32, 47"
status: >
  OPEN (conjecture). The 'if' directions are PROVED in THM-4759 (3), (9), (10) and THM-4760 (3), (4); THM-4760 also
  supplies the rationality input of THM-4759 (10).
  Part (A), REFINED after the S53 audit (order-5 census, audit_S53/order5_strict_census.out). The first wording
  defined a strict Fricke point by geometry alone: a unique vector of norm 2 - 2/n plus deep proper multiples. At
  order 5 the audit samples 45,890 such points in two classes.
  - One class is the A3^8 (h+1)-point (n_2 = 24).
  - The other has n_2 = 23, and its conjugate 2v has the shallower minimum 6/5. So its twining is irrational and not
    a Hauptmodul.
  So geometry alone does not suffice. The definition now requires the Fricke Hauptmodul property, which by THM-4760
  adds rationality and n_2 = 24. The prime-order sketch below silently used exactly this, through the claim that all
  twisted sectors are Galois transforms of one another.
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
  cancel.
  - The members found are PROVED in THM-4760 by its cusp criterion: rationality (vector-valued Sturm), a unique
    Fricke-norm vector, 24 norm-2 vectors, and Phi_c | P at the intermediate cusps. They are seven shell-1
    points (carriers A1, D4, D6, D8, E7 in D10E7^2, D12, E8 in D16E8), of orders 4, 8, 12, 16, 20, 24 and 32,
    twining by 4A, 8A, 12H, 16C, 20F, 24I and 32A (Ligozat certificates for the eta quotients; q^4 agreement for 4A
    and 8A).
  - Their only shallow multiple is 2v. It carries a pair +-mu with opposite eigenvalues, and its orbifold is
    Borcherds' 2-descent source N_1 (h = 2h_1 + 2).
  - The orbifold of v is N itself, of Coxeter number e - 2, so dim V'_1 = 24(e - 1), not 24e.
  - All other shell-1 classes fail, each provably. The D10 (in D10E7^2), D16 (in D16E8), E8^3 and D24 points
    pass (i)-(iii) and fail only rationality (fields contain Q(sqrt 5), Q(sqrt 2), Q(sqrt 2), Q(sqrt 3)). The selection
    is Galois acting on component Kac classes, THM-4760 (4'): E7 at order 20 has its two classes swapped by d = 3, which
    the swap of the two E7's undoes. The even-h lattices outside the nine fail at the single
    shallow vector of 2v (a pole at 1/2), and odd h fails at n_2 = 23.
  - Shell 2 of every type, shells 2-6 of the 14 types with h <= 14, and shells 2-8 of the eight types with h >= 16
    give nothing (sampled where large).
  - (B') below is the sharpened conjecture.
  OPEN:
  - (B') for cusps not adjacent to a deep-hole cusp. The cusps at distance 2 and 3 have been sampled and give only the
    30 known classes.
  - A referenced write-up of the prime-order sketch.
  - The glue compatibility in THM-4760 (4'). The exact Kac-class Galois action is computed there: E8's carrier class
    is fixed and E7's classes are swapped by d = 3. Rationality itself is proved by the certificate.
  History:
  - The first S53 draft allowed only the eleven bad-prime multiples; the (h+1)-points refuted it.
  - The second version (registered as c31e9a6a3f..23dc8edd9e) claimed the 'only if' for all orders; the shell-1
    D-type-lattice points refuted it. Both refutations happened within the session.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_h_plus_one.py (+ _q2.out, .out)
  - 04-computation/experiments/deep_holes_theta_20261022/leech_h_plus_one.py (+ .out), weyl_h_plus_one.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/order3_fricke_census.py (+ .out), order4_fricke_census.py (+ .out), order6_fricke_census.py (+ .out)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_shell_test.py (+ .out) (the shell-1 points; the (B) family)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_cusp_certificate.py (+ .out, _deep.out), fricke_monster_classes.py (+ .out), fricke_cusp_orbifolds.py (+ .out) (THM-4760)
related:
  - THM-4759 (deep-hole twining, Fricke points, Fricke descent)
  - THM-4758 (swap theorem)
  - THM-4760 (cusp criterion; proves the (B) members and the 2-descent table)
---

# HYP-9368 — Fricke points and genus-zero inner automorphisms of V_Λ

## Statement

Let v ∈ Λ ⊗ Q have order n modulo Λ. Write Z_v = Tr_{V_Λ}(e^{2πi v(0)} q^{L_0 − 1}) = Θ_Λ^{[v]}/η²⁴.

Call v a **strict Fricke point** if it satisfies both of these:
- Z_v − 24 is the normalized Hauptmodul of Γ_0(n) + n. By THM-4760 (1) this means: Z_v is rational, Λ + v has a
  unique vector of norm 2 − 2/n, and Λ + v has exactly 24 vectors of norm 2. Condition (iii) is vacuous when the
  multiples are deep.
- every proper multiple jv (j | n, 1 < j < n) is a deep hole or lies in Λ.

The geometric conditions alone (a unique Fricke-norm vector and deep multiples) are not enough. The audit's order-5
census finds an irrational class with n_2 = 23 among them.

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

The members of the last family found so far are seven shell-1 points (carriers A1, D4, D6, D8, E7 in D10E7², D12,
E8 in D16E8), of orders 4, 8,
12, 16, 20, 24 and 32 (PROVED in THM-4760).

**Conjecture (B′), sharpened.** There are no further Fricke-type points. So the Fricke Hauptmodul twinings of V_Λ by
translations are exactly the 30 classes of THM-4760 (3) and (4), with the 23 orders 3–17, 19, 20, 23, 24, 26, 31, 32
and 47. In particular none exists at the remaining genus-zero Fricke levels 18, 21, 25, 27, 29, 35, 36, 39, 41, 49,
50, 59 and 71.

Evidence:
- Primes are covered by (A).
- On the cusps adjacent to a deep-hole cusp, shells 0–2 of all types and shells 2–8 of the eight types with h ≥ 16
  give nothing new. These shells reach every missing level from 18 to 50.

## The family that refuted the first registered (B) (now THM-4760 (4))

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

the twining Z_{c_w} − 24 is the Fricke Hauptmodul of Γ_0(e) + e. This is PROVED in THM-4760 by the cusp criterion,
with the Monster names 4A, 8A, 12H (A058486), 16C (A058516), 20F (A058555), 24I (A058579) and 32A (A058629) certified
by Ligozat. Yet 2c_w is shallow, carrying a pair ±μ whose twined contributions cancel. So the absence of extra poles
comes from cancellation of the twined shallow states, not from deep multiples. Rationality is the remaining filter:
the D10, D16, E8³ and D24 shell-1 points cancel too, but are irrational. In D10E7² and D16E8 the rational carrier
is the E-component, not the D-component (THM-4760 (4), (4′)).

The orbifold of V_Λ at c_w is V_N, of Coxeter number e − 2, so dim V′_1 = 24(e − 1). These are therefore not strict
Fricke points, and they do not contradict (A).

The other shell-1 classes fail. This includes D24, of order 48, consistent with X_0^+(48) not having genus zero.
D24 cancels at every cusp, but it is irrational.

## Evidence for (A)

- **Order 2.** Impossible: v ≡ −v pairs the minimal vectors.
- **Order 3.** One class, the M24-point (2²⁴)/(3√8). This rests on the two-orbit decomposition of the norm-12 vectors,
  which is classical: orbits O6a and O6b, ATLAS p. 181; see Sun–Wang, arXiv 2111.10999, §2.3, whose Lemma 2.11 gives
  the unique minimal vector of the class (order3_fricke_census.py).
- **Order 4, with 2v deep.** One invariant in 2046 samples, that of the A2¹² point (order4_fricke_census.py).
- **Order 6.** All 1399 samples have neighbour A4⁶ (order6_fricke_census.py).
- **Prime orders.** The sketch below.

## The classification at prime order (sketch)

Let v be a strict Fricke point of prime order p, and g = e^{2πi v(0)}. Prime order has no proper multiples, so (A)
and (B) agree here.

1. **Dimension.** Since gcd(i, p) = 1 for i ≠ 0 and Z_v is rational, every twisted sector V_Λ(g^i) is a Galois and
   modular transform of the g-twisted one. Its g-invariant part has character U_p F, whose constant term is
   c_F(0) = 24. So V′ = V_Λ^{orb(g)} has dim V′_1 = 24p. Without rationality this step fails: in the audit's
   irrational order-5 class, Λ + 2v has minimum 6/5, unlike Λ + v.
2. **The orbifold is a lattice VOA.** V′ = V_{N′} with N′ Niemeier and |R(N′)| = 24(p − 1). Hence h(N′) = p − 1.
3. **The dual automorphism.** τ = e^{2πiu(0)} with u ∈ (1/p)N′. It is regular, since its fixed lattice Λ_0 is
   rootless.
4. **Kac coordinates.** A regular element of order h + 1 has Kac coordinates 2 at one special node and 1 elsewhere.
   Ω permutes the special nodes transitively, so u_i ∈ W·(ρ_i/p) + P_i. This is Kac, LNM 848 (1981), Thm 1: the
   regular elements of order h + 1 form one conjugacy class, credited to Macdonald (1972).
5. **Glue.** gcd(p, |P/N′|) = 1, so u ≡ w(ρ′/p) mod N′.
6. **Duality.** Orbifold duality recovers (V_Λ, g) from (V_{N′}, τ).

So the order-p Fricke classes are indexed by the Niemeier lattices with h = p − 1:
- two classes each at p = 7, 11, 13, 19 and 31;
- one class each at p = 3, 5, 17, 23 and 47.

## Predictions for (A)

- No strict Fricke point of order 2, 12, 16, 18, 20, 21, 22, 24, 25, 27–30, 32–46, or any order ≥ 48. At 12, 16,
  20, 24 and 32 the (B) family supplies genus-zero twinings instead (THM-4760).
- Exactly two classes at orders 7, 11, 13, 19 and 31, and one at each of 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 23, 26
  and 47.
- The D16E8 and E8³ shallow multiples of order 15 are equivalent: both have neighbour D8³ and equal theta series.
