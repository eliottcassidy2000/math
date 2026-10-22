---
id: THM-4759
title: "Deep-hole twining of the Leech VOA is umbral. For every Niemeier root system X (Coxeter number h, k components, glue N/Q, Coxeter Frame shape pi^X): theta(N + rho_N/h) = |N/Q| eta(tau)^24 eta_{pi^X}(tau/h) and theta(Lambda + c) = Delta(tau) (24 + k + |N/Q| eta_{pi^X}(tau/h)); the principal automorphism e^{2 pi i rho_N(0)/h} of V_N has graded trace T^X = 1/eta_{pi^X}, the Cheng-Duncan-Harvey principal modulus of lambency X, and the deep-hole automorphism e^{2 pi i c(0)} of V_Lambda has graded trace T^X + 24 + k; the two traces agree in every weight except weight one. The poles of T^X away from infinity are the shallow multiples d c of the deep hole, which occur exactly at the bad primes d of every component with d || h, and the lambency is h + {h/p : p such a prime}. Those multiples v = pc (Fricke points) twine V_Lambda by 24 + the Fricke Hauptmodul of Gamma_0(h/p) + (h/p): 3A, 4A, 5A, 6B, 7A, 9A, 10D, 11A, 15C and 23AB from the nine A-free types. More generally, every Niemeier lattice gives a Fricke point c'' = c + (lambda - c)/(h+1) of order h+1, one Leech wall away. Its Fricke swap theorem theta(Lambda + c'') - theta(N + rho/(h+1)) = 24 Delta is PROVED, and it twines by 24 + T_{(h+1)+} (3A ... 47AB, 26B). The orbifold of V_Lambda at pc is the Niemeier lattice of the Borel-de Siebenthal descent X_p (delete the label-p nodes): D24 -> D12^2 -> D6^4 -> A3^8, D16E8 -> D8^3 -> D4^6 -> A1^24, E8^3 -> D8^3, A8^3, A4^6. So every power of a deep-hole automorphism twines V_Lambda by 24 + a Hauptmodul, and all the groups are monstrous except Gamma_0(25)"
status: >
  PROVED modulo one classical input, plus VERIFIED structure.
  - (1) The per-component identity theta_{Q + rho/h} = eta^l eta_pi(tau/h) is CLASSICAL: Kostant's conjugacy of the
    principal element with a Coxeter element, applied to the lattice VOA V_Q (equivalently Lepowsky and
    Kac-Kazhdan-Lepowsky-Wilson, the principally specialised basic representation; Garvan-Kim-Stanton for type A).
    Independently VERIFIED through q^29 for A_1..A_8, A_12, D_4..D_10, E_6, E_7, E_8 (principal_cores.py).
  - (2) and (3) are PROVED from (1), the glue-coset lemma proved below, and THM-4758.
  - (4) T^X is a Hauptmodul for the CDH lambency group. This is CITED (Cheng-Duncan-Harvey 2014, Prop. 2.1, proved
    by inspection) and independently re-certified here (hauptmodul_cert.py).
  - (5) The shallow-multiple law has three parts.
    - The bounds (min >= 2 - 2d/h, and min >= 2 outside the common bad primes) and the attainment are PROVED, from
      Kostant's conjugacy, the vacuum anomaly of lattice twisted modules, and the exponent inequality
      sum_i (m_i mod e)(e - (m_i mod e)) >= l e (e - 1)/6. The exponent inequality is proved here for A_n and D_n,
      and checked exactly for E_6, E_7, E_8 and for 2815 pairs (R, e) (exponent_law.py).
    - The uniqueness of the minimal vector and the norm-2 counts are VERIFIED for all 23 types by exact enumeration
      on both sides (multiple_cosets.py on the Niemeier side, deep_hole_multiples.py on the Leech side).
    - The component rule for Atkin-Lehner invariance is PROVED for A_n and D_n and checked for E (360 pairs).
  - The Monster corner (6) is CITED: CDH Section 2.3 and Carnahan, "51 constructions of the Moonshine module" (2017),
    Cor. 4.6.
  - (7) Fricke twining is PROVED, given the VERIFIED uniqueness of the shallow vector, by a principal-part argument
    on X_0(e). The q-series are checked through q^3 for all eleven Fricke points (fricke_multiples.py).
  - (8) Fricke descent is VERIFIED for all eleven cases by exact root enumeration (fricke_neighbors.py), and agrees
    with the divisor-lattice table of THM-4758 (4).
    - It is identified as the subsystem of roots of height divisible by p (height_descent.py).
    - It is PROVED for prime e = h/p by the inclusion R(L_e) >= X_p together with the root count 24(e - 1).
  - (9) is PROVED from (3), (5) and (7).
  - (10) The (h+1)-points are covered as follows.
    - VERIFIED for all 23 types: the unique minimal vector of norm 2 - 2/(h+1), the 24 norm-2 vectors and the deep
      proper multiples (leech_h_plus_one.py, fricke_h_plus_one.py), and the nudge geometry (fricke_nudge.py).
    - PROVED from those: the Fricke twining.
    - PROVED: the height equidistribution mod h+1.
  - (11) The Fricke swap theorem is PROVED by the same vector-valued Sturm certificate as THM-4758, run on the 23 pairs
    (rho, z'') (swap_fricke_proof.py). It also supplies the Leech-side facts of (10) from the N side.
  - The M24-point remark is VERIFIED (sampling plus exact orbit-size arithmetic). It is not a proof that the k = 0
    vectors form a single orbit.
  Prior-art search (agent, 2026-10-22): no published statement of (2), (3) or (5) was found. The search predates
  (7)-(11), which have not yet been searched. CDH Remark 2.2 asks for
  a conceptual proof that T^X is genus zero; (2) identifies T^X as a VOA trace but does not by itself prove genus zero.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/principal_cores.py (+ .out) (per-component identity, glue lemma, anti-self-duality)
  - 04-computation/experiments/deep_holes_theta_20261022/coxeter_frame_twined.py (+ .out) (Theta_N^[rho/h] = eta^24/eta_pi, all 23 types, by enumeration)
  - 04-computation/experiments/deep_holes_theta_20261022/twined_hauptmodul.py (+ .out, _b.out) (Leech twined traces, 8 pure A types)
  - 04-computation/experiments/deep_holes_theta_20261022/hauptmodul_cert.py (+ .out) (Hauptmodul certificate, lambencies, component rule)
  - 04-computation/experiments/deep_holes_theta_20261022/multiple_cosets.py (+ .out), deep_hole_multiples.py (+ .out) (shallow multiples)
  - 04-computation/experiments/deep_holes_theta_20261022/exponent_law.py (+ .out) (exponent inequality, twisted-module bound)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_multiples.py (+ .out) (Fricke twining, eleven Fricke points)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_neighbors.py (+ .out) (Fricke descent: neighbour root systems)
  - 04-computation/experiments/deep_holes_theta_20261022/order3_fricke_census.py (+ .out) (the order-3 Fricke point is the M24-point)
  - 04-computation/experiments/deep_holes_theta_20261022/height_descent.py (+ .out) (Fricke descent = roots of height divisible by p)
  - 04-computation/experiments/deep_holes_theta_20261022/swap_fricke_proof.py (+ .out) (the Fricke swap theorem: certificates for the 23 pairs (rho, z''))
  - 04-computation/experiments/deep_holes_theta_20261022/swap_shell_test.py (+ .out) (the next shell m = h + 2 fails)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_nudge.py (+ .out) (Fricke point = deep hole nudged 1/(h+1) toward a second-shell point; z'' = z + Leech wall)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_shell_test.py (+ .out) (shell-1 points: Fricke-invariant twining for the seven D-types despite shallow multiples)
related:
  - THM-4758 (swap theorem and deep-hole theta identity)
  - HYP-9363 (proved in THM-4758)
  - HYP-9368 (genus-zero inner automorphisms of V_Lambda: deep holes and Fricke points)
---

# THM-4759 — deep-hole twining of the Leech VOA is umbral

## Notation

Notation is as in THM-4758, plus the following.
- X is the root system of the Niemeier lattice N, with Coxeter number h and k components.
- Q = Q(X) is its root lattice and |N/Q| the order of the glue code.
- c is the corresponding deep hole of Λ, and ρ_N is the Weyl vector of N.
- For an irreducible simply-laced root system R, the **Coxeter Frame shape** π_R = Π_a a^{m_a} is defined by
  det(x − c_R) = Π_a (x^a − 1)^{m_a} for a Coxeter element c_R:

      A_n: (n+1)·1^{−1};   D_n: (2n−2)·2·(n−1)^{−1}·1^{−1};
      E_6: 12·3·2·6^{−1}·4^{−1}·1^{−1};   E_7: 18·3·2·9^{−1}·6^{−1}·1^{−1};
      E_8: 30·5·3·2·15^{−1}·10^{−1}·6^{−1}·1^{−1}.

- π^X is the product over the components, and η_π(τ) = Π_a η(aτ)^{m_a}. Every π^X has degree Σ a m_a = 24 and
  Σ m_a = 0.
- **T^X := 1/η_{π^X}(τ) = q^{−1} − k + O(q).** This is the principal modulus of Cheng–Duncan–Harvey (CDH 2014,
  Prop. 2.1).
- σ_N = e^{2πiρ_N(0)/h} ∈ Aut(V_N) is the principal (Kostant) automorphism. g_c = e^{2πic(0)} ∈ Aut(V_Λ) is the
  deep-hole automorphism. Both have order h.
- Z_V(g) := Tr_V(g q^{L_0 − 1}) is the graded trace.

## Statement

**(1) Components (classical).** For every irreducible simply-laced root lattice Q of rank ℓ:

    Σ_{x ∈ Q} e^{2πi x·ρ/h} q^{|x|²/2} = η(τ)^ℓ / η_π(τ),     θ_{Q + ρ/h}(τ) = η(τ)^ℓ η_π(τ/h).

Every coset Q + w (w ∈ Q^*) has the same ρ/h-shifted theta series. π is anti-self-dual: a ↦ h/a sends π to π^{−1}.
Moreover Π_a a^{m_a} = det(1 − c_R) = |Q^*/Q|.

**(2) Niemeier lattices.**

    θ_{N + ρ_N/h}(τ) = |N/Q| · η(τ)^{24} · η_{π^X}(τ/h),       Z_{V_N}(σ_N) = Θ_N^{[ρ_N/h]}/η²⁴ = T^X.

So CDH's principal modulus is the graded trace of the principal automorphism of the Niemeier lattice VOA.

**(3) The Leech lattice.**

    θ_{Λ + c}(τ) = Δ(τ) · (24 + k + |N/Q| · η_{π^X}(τ/h)),
    Z_{V_Λ}(g_c) = T^X(τ) + 24 + k,
    Tr_{V_Λ(g_c)}(q^{L_0 − 1}) = T^X(−1/τ) + 24 + k   (the g_c-twisted module V_{Λ+c}).

Hence **Tr(g_c | V_Λ(n)) = Tr(σ_N | V_N(n)) for every n ≠ 1**. In weight 1 the traces are 24 (the Heisenberg algebra,
fixed pointwise) and −k (each simple factor of the Lie algebra V_N(1) contributes −1).

For X = ∅ (c = 0) this reads Z_{V_Λ}(1) = J + 24. Twining by a deep hole replaces J by the umbral principal modulus
T^X + k = q^{−1} + O(q).

**(4) Genus zero (CDH Prop. 2.1, re-certified).** T^X is a Hauptmodul for Γ_0(h) + S_X. Its poles on X_0(h) are
exactly the cusps 1/(h/e) for e ∈ S_X, all simple. Γ_0(h) + S_X is the CDH lambency:
- h for the 14 types with an A-component;
- 6+3, 10+5, 14+7, 18+9, 22+11, 30+15 and 46+23 for D4⁶, D6⁴, D8³, D10E7², D12², D16E8 and D24;
- 12+4 for E6⁴, and 30+6,10,15 for E8³.

**(5) Bad primes, shallow multiples and Atkin–Lehner involutions.** Let d | h with 1 < d < h, and put e = h/d.

- (a) Suppose d is a bad prime of every component of X (bad primes: D_n has 2; E_6 and E_7 have 2, 3; E_8 has 2, 3,
  5; A_n has none), and gcd(d, e) = 1. Then:
  - W_e fixes T^X;
  - the coset N + dρ_N/h has minimal norm 2 − 2d/h, attained by exactly one vector;
  - it has φ(d)·T(d) vectors of norm 2, where T(d) is the number of components with bad prime d (THM-4758);
  - on the Leech side, dc lies at distance √(2 − 2d/h) from a unique Leech point. So dc is an interior point of a
    Voronoi cell, not a hole.
- (b) In every other case, W_e does not fix T^X. Both N + dρ_N/h and Λ + dc have minimal norm ≥ 2, and dc is again a
  deep hole, of order e and hence of Coxeter number e. In all 39 instances its type contains A_{e−1}: A1²⁴, A2¹²,
  A3⁸, A4⁶, A5⁴D4 (5 instances, never D4⁶), A7²D5², A8³ or A9²D6 (never D6⁴).

  Partial reason. The orbifold of V_Λ by g_c^d is the divisor lattice L_e, and L_e ⊇ {x ∈ N : d | x·ρ_N}. So
  R(L_e) contains X_d = {α : d | ht α}. On an A_{h−1} component, X_d = {e_i − e_j : i ≡ j mod d} = A_{e−1}^d.

**(5′) The exponent law behind (5).** Let R be irreducible with rank ℓ, Coxeter number h and exponents m_1, …, m_ℓ.
For a divisor e of h put J_e(m) = r(e − r), where r = m mod e. Then:

    min_{x ∈ P(R)} |x + ρ/e|²  ≥  (1/2e²) Σ_i J_e(m_i)  ≥  (ℓ/12)(1 − 1/e).

Equality holds on the right exactly when e = 1, or e = h/p for a bad prime p of R with p ∥ h. In those cases the
minimum on the left attains the bound. The case e = h of the middle term is the strange formula:
Σ_i m_i(h − m_i) = ℓh(h+1)/6.

Consequently **S_X = {1} ∪ {h/p : p a bad prime of every component of X with p ∥ h}**. The component rule holds for
every irreducible R: π_R is W_e-invariant (e ∥ h, e > 1) if and only if e = h/p for a bad prime p ∥ h(R). Concretely:
- D_n: e = n − 1 when n is even;
- E_6: e = 4;
- E_7: e = 9;
- E_8: e = 6, 10, 15;
- A_n: never.

**(7) Fricke twining by the shallow multiples.** In case (5a) put v = dc. This is a point of order e = h/d; call
it a **Fricke point**. Then

    Tr_{V_Λ}(e^{2πi v(0)} q^{L_0 − 1}) = 24 + T_{e+e}(τ),

where T_{e+e} = q^{−1} + O(q) is the normalized Hauptmodul of the Fricke group Γ_0(e) + e. Moreover Λ + v contains
exactly 24 vectors of norm 2 and exactly one of norm 2 − 2/e.

The eleven Fricke points come from the nine A-free types. These are exactly the types with a Fricke point, because
A_n has no bad primes.

| X | d | e | Z − 24, through q³ | Γ |
|---|---|---|---|---|
| D4⁶ | 2 | 3 | q⁻¹ + 783q + 8672q² + 65367q³ | 3+ (T_3A) |
| E6⁴ | 3 | 4 | q⁻¹ + 276q + 2048q² + 11202q³ | 4+ (T_4A, A107080) |
| D6⁴ | 2 | 5 | q⁻¹ + 134q + 760q² + 3345q³ | 5+ (T_5A) |
| E8³ | 5 | 6 | q⁻¹ + 78q + 364q² + 1365q³ | 6+6 (T_6B, OEIS A007255) |
| D8³ | 2 | 7 | q⁻¹ + 51q + 204q² + 681q³ | 7+ (T_7A) |
| D10E7² | 2 | 9 | q⁻¹ + 27q + 86q² + 243q³ | 9+ (T_9A, A007266) |
| E8³ | 3 | 10 | q⁻¹ + 21q + 62q² + 162q³ | 10+10 (T_10D, A058100) |
| D12² | 2 | 11 | q⁻¹ + 17q + 46q² + 116q³ | 11+ (T_11A) |
| D16E8, E8³ | 2 | 15 | q⁻¹ + 9q + 19q² + 42q³ | 15+15 (T_15C, A058510) |
| D24 | 2 | 23 | q⁻¹ + 4q + 7q² + 13q³ | 23+ (T_23AB) |

The squarefree composite levels have the eta-quotient Hauptmoduln (η(pτ)η(qτ)/(η(τ)η(pqτ)))^{24/((p−1)(q−1))} minus
their constant terms:
- 6+6: exponent 12, constant 12;
- 10+10: exponent 6, constant 6;
- 15+15: exponent 3, constant 3.

For 9+ the Hauptmodul is x + 3 + 27/x with x = (η(τ)/η(9τ))³.

There is no Fricke point of order 2: an order-2 point v satisfies v ≡ −v mod Λ, so its minimal vectors come in pairs.

**(8) Fricke descent (Borel–de Siebenthal).** The cyclic orbifold V_Λ^{orb(e^{2πiv(0)})} at a Fricke point v = pc is
the lattice VOA of a Niemeier lattice N′ with Coxeter number e − 1. Its root system is X_p: delete a node of label p
from the extended Dynkin diagram of every component, keeping the maximal-rank subsystem of Coxeter number h/p − 1.
- p = 2: D_n ↦ D_{n/2}², E7 ↦ A7, E8 ↦ D8.
- p = 3: E6 ↦ A2³, E8 ↦ A8.
- p = 5: E8 ↦ A4².

All eleven cases:
- D4⁶ → A1²⁴, D6⁴ → A3⁸, D8³ → D4⁶, D12² → D6⁴, D24 → D12²;
- D10E7² → A7²D5²;
- D16E8 → D8³;
- E6⁴ → A2¹²;
- E8³ → D8³ (p = 2), A8³ (p = 3), A4⁶ (p = 5).

So p = 2 halves the D-index: D24 → D12² → D6⁴ → D3⁸ = A3⁸, and D16E8 → D8³ → D4⁶ → D2¹² = A1²⁴. The corresponding
Fricke levels are 23, 11, 5 = 3·2^k − 1 and 15, 7, 3 = 2^k − 1.

*Mechanism.* X_p is exactly the set of roots of N whose height α·ρ_N is divisible by p, i.e. the fixed-point subsystem
of the principal element exp(2πiρ_N/p) of order p (height_descent.py, all eleven cases).

Its size is 24(h/p − 1). Per component this is the exponent identity Σ_i (m_i mod p) = ℓp/2, since
#{α : p | ht α} = 2Σ_i ⌊m_i/p⌋ and Σ_i m_i = ℓh/2:
- D_n with n even, p = 2: all exponents are odd;
- E6, p = 3: residues 1, 1, 2, 1, 2, 2;
- E7, p = 2: all exponents are odd;
- E8, p = 2, 3, 5: likewise.

**Conceptual proof for prime e = h/p.** The descent lattice is L_e. Its classes S_e meet the Niemeier line in the
classes of height ≡ 0 mod p, so R(L_e) ⊇ X_p. The Fricke twining (7) gives dim V′_1 = 24e, i.e. |R(L_e)| = 24(e − 1)
= |X_p|. Hence R(L_e) = X_p. For composite e this is verified only.

**The Coxeter numbers move by the affine maps H ↦ pH + 1.** Write H = h + 1, the order of the Fricke point in (10).
Then h = p(h′ + 1) becomes H = pH′ + 1. Descent at the bad prime p is H ↦ (H − 1)/p:
- p = 2, the Cunningham chain 5 → 11 → 23 → 47: A3⁸ → D6⁴ → D12² → D24;
- p = 2, the Mersenne chain 3 → 7 → 15 → 31: A1²⁴ → D4⁶ → D8³ → D16E8 or E8³; also 9 → 19: A7²D5² → D10E7²;
- p = 3, the Collatz map 3n + 1: A2¹² (H = 4 = 3¹ + 1) → E6⁴ (H = 13), and A8³ (H = 10 = 3² + 1) → E8³ (H = 31);
- p = 5: A4⁶ (H = 6) → E8³ (H = 31).

E8³ is reached in three ways: 31 = 2·15 + 1 = 3·10 + 1 = 5·6 + 1.

This also explains the divisor lattices of THM-4758 (4). L_d is the orbifold of V_Λ by g_c^{h/d}:
- Coxeter number d − 1 exactly when h/d is a common bad prime ∥ h (Fricke descent);
- Coxeter number d otherwise; then (h/d)c is a deep hole of an A-containing type.

All 23 rows of swap_theorem_proof.out agree.

**(9) Genus zero for every power of a deep-hole automorphism.** For every deep hole c and every m ∈ Z,
Tr_{V_Λ}(g_c^m q^{L_0 − 1}) − 24 is the normalized Hauptmodul of a genus-zero group:
- J, if mc ∈ Λ;
- T^{X′} + k′, if mc is a deep hole of type X′;
- T_{e+e}, if mc is a Fricke point.

All of these groups are monstrous except Γ_0(25), which arises from the A24 hole and its powers prime to 5. The
composite Fricke levels are the Monster classes 6B, 10D and 15C, identified through OEIS A007255, A058100 and
A058510.

**(10) Every Niemeier lattice has a Fricke point: the (h+1)-construction.** Let z″ be the cusp of type N with
z·z″ = −1 and z″·ρ = −(h + 1). It exists because ρ = (ρ_N, h, h + 1) in II_{25,1} = N ⊕ U. Let
c″ = proj_Λ(z″)/(h + 1), a point of order e = h + 1. Then for all 23 types:
- Λ + c″ has a unique vector of norm 2 − 2/(h+1) and exactly 24 vectors of norm 2. Every proper multiple jc″
  (j | e, 1 < j < e) is a deep hole. The same holds on the N side: N + ρ_N/(h+1) has a unique minimal vector,
  ρ_N/(h+1) itself, of norm 2h/(h+1).
- The (h+1)-neighbour of Λ at c″ is N. Equivalently, Λ is the (h+1)-neighbour of N at ρ_N/(h+1): a second holy
  construction, next to the classical one at ρ_N/h.
- **Tr_{V_Λ}(e^{2πic″(0)} q^{L_0 − 1}) = 24 + T_{(h+1)+(h+1)}**, the normalized Fricke Hauptmodul.

| e = h + 1 | N | Z − 24 = q⁻¹ + … |
|---|---|---|
| 3 | A1²⁴ | 783q + 8672q² (3A) |
| 4 | A2¹² | 276q + 2048q² (4A) |
| 5 | A3⁸ | 134q + 760q² (5A) |
| 6 | A4⁶ | 78q + 364q² (6B) |
| 7 | A5⁴D4, D4⁶ | 51q + 204q² (7A) |
| 8 | A6⁴ | 36q + 128q² (8A, A007265) |
| 9 | A7²D5² | 27q + 86q² (9A) |
| 10 | A8³ | 21q + 62q² (10D) |
| 11 | A9²D6, D6⁴ | 17q + 46q² (11A, A058205) |
| 13 | A11D7E6, E6⁴ | 12q + 28q² (13A, A034319) |
| 14 | A12² | 10q + 24q² + 51q³ (14C, A058504) |
| 15 | D8³ | 9q + 19q² (15C) |
| 17 | A15D9 | 7q + 14q² + 29q³ (17A, A058530) |
| 19 | A17E7, D10E7² | 6q + 10q² + 21q³ (19A, A058549) |
| 23 | D12² | 4q + 7q² + 13q³ (23AB, A058570) |
| 26 | A24 | 3q + 6q² + 9q³ (26+26 = T_26B, A058597; (η(2τ)η(13τ)/(η(τ)η(26τ)))² − 2) |
| 31 | D16E8, E8³ | 3q + 3q² + 6q³ (31AB, A058628) |
| 47 | D24 | q + 2q² + 3q³ (47AB, A058690) |

*Geometry: a Fricke point is a deep hole nudged towards the second shell.* fricke_nudge.py checks all 23 types.
- **The cusp.** z″ = z + r, where r is a Leech root (a wall of the Leech chamber, r·ρ = −1) with r·z = −1. So z″ is
  the reflection of the deep-hole cusp in a wall adjacent to it. Such walls are the Leech roots r_λ with λ in the
  second shell of the hole, |λ − c|² = 2 + 2/h (THM-4758's r_λ·z = h(1 − |λ − c|²/2)).
- **The point.** c″ = c + (λ − c)/(h + 1): the deep hole moved 1/(h+1) of the way towards a second-shell Leech point
  λ. Then |c″ − λ|² = h²(2 + 2/h)/(h + 1)² = 2 − 2/(h + 1), and λ is the unique nearest Leech point.
- **The vertices.** Exactly 24 of the 24 + k hole vertices stay at distance √2 from c″. The other k, one per
  component, move away. In the hole diagram, which is the extended Dynkin diagram, the finite simple roots stay and
  the k affine nodes leave.

This is the geometric meaning of the defects 24 + k and 24 in (11).

*Worked example: A1²⁴ in standard Leech coordinates.*
- c = (1²⁴)/√8 is an A1²⁴ deep hole, fixed by M24. Its 48 vertices are (1²⁴ ∓ 4e_i)/√8. Both shapes, (−3, 1²³) and
  (5, 1²³), are Leech vectors.
- The origin lies in its second shell: |c|² = 3 = 2 + 2/h.
- So c″ = c + (0 − c)/3 = (2/3)c = (2²⁴)/(3√8), which is the M24-point of the order-3 census.
- From c″, the vertices (1²⁴ − 4e_i)/√8 stay at squared distance (23/9 + 121/9)/8 = 2, while the partners
  (1²⁴ + 4e_i)/√8 move to (23/9 + 169/9)/8 = 8/3.
- The deep hole and the Fricke point both lie on the M24-fixed line through (1²⁴).

*Heights mod h + 1 (PROVED).* The dual automorphism on V_N is e^{2πiρ_N(0)/(h+1)}. A root α lies in its ζ^j-eigenspace
iff ht(α) ≡ j mod (h + 1). Per component, the count is n_j + n_{h+1−j}, where n_j = #{exponents ≥ j}. Under
m ↦ h − m this equals ℓ. So the roots of N are equidistributed: exactly 24 in each nonzero class of ht mod h + 1, and
none in the zero class.

This gives two consequences:
- dim V′_1 = 24 + 24h = 24(h+1), and V_Λ^{orb} = V_N;
- the 24 norm-2 vectors of Λ + c″ correspond to the 24 simple roots of N.

So the "vertices" of a Fricke point are the simple roots.

*Relation to generalised deep holes.* Möller–Scheithauer realise each V_N from V_Λ by a generalised deep hole, unique
up to conjugacy. Their twisted sectors all have conformal weight ≥ 1 (extremality, R(g) = 0), and for V_N that
automorphism is g_c.

The Fricke point g_{c″} also has orbifold V_N, so it must fall outside their class. It does: its twisted sector has
the single state of weight 1 − 1/(h+1) < 1. So for each N the two inner automorphisms are the extremal one (the deep
hole) and the Fricke one, which is maximally non-extremal: one shallow state, Fricke-symmetric.

Consequences:
- **Supersingular primes.** The prime values of h + 1 are 3, 5, 7, 11, 13, 17, 19, 23, 31 and 47. Each of them is
  therefore a prime p for which X_0(p)/W_p has genus zero: Ogg's supersingular primes, here forced by a lattice. The
  missing supersingular primes are 2, 29, 41, 59 and 71, for which h = p − 1 is not a Niemeier Coxeter number.
- **The shallow multiples of (7).** Each is the (h+1)-point of the corresponding Fricke descent X_p of (8). For
  example, the D4⁶ point 2c is the (h+1)-point of A1²⁴. This is why its order is e = h(X_p) + 1.
- **The A1²⁴ point is fixed by M24.** Its stabiliser contains the umbral group G^X = M24, and by the order-3 census it
  is the point (2²⁴)/(3√8).

**(11) The Fricke swap theorem (PROVED).** The certificate of THM-4758 applies verbatim to the cusp pair (ρ, z″), with
M″ = {ρ, z″}^⊥ and D = (Z/(h+1))². It succeeds for all 23 types (swap_fricke_proof.py, 8 s): Θ_{M″} − σ*Θ_{M″} = Δ·v″,
with line sum 24. Hence

    θ_{Λ+c″} − θ_{N+ρ_N/(h+1)} = 24 Δ,     Θ_Λ^{[c″]} − Θ_N^{[ρ_N/(h+1)]} = 24 Δ,

and so

    Tr_{V_N}(e^{2πiρ_N(0)/(h+1)} q^{L_0 − 1}) = T_{(h+1)+(h+1)}   exactly.

The weight-one trace is 24 + 24(ζ + … + ζ^h) = 0, by the height equidistribution. The Leech trace is this plus 24.
N + ρ_N/(h+1) has no vectors of norm 2, and Λ + c″ has 24.

The defect vector is always the pure one, v″ = 24(1_{S_{h+1}} − 1_{S_1}), with no intermediate terms for any
composite h + 1 (swap_fricke_proof.out). Unlike THM-4758, no bad primes enter.

The two theorems sit side by side:

| | point of Λ | order | defect | what the defect counts |
|---|---|---|---|---|
| THM-4758 | deep hole c | h = Σ_{i≥0} a_i | 24 + k | nodes of the extended Dynkin diagrams (the hole vertices) |
| (11) | Fricke point c″ | h + 1 | 24 | nodes of the finite Dynkin diagrams (the simple roots) |

The twining sides are parallel too. V_N twined by the principal element of order h gives the non-Fricke umbral T^X.
Twined by the element of Kac type (2;1,…,1) and order h + 1, it gives the Fricke T_{(h+1)+}.

*Levels 0 and 1 (ANALOGY with a precise core).* At level k the Weyl point of the level-k alcove is
(ρ̂ + kΛ_0)/(k + h^∨), with Kac coordinates (k + 1; 1, …, 1).
- k = 0 gives ρ/h, Kostant's principal element: the deep hole.
- k = 1 gives ρ/(h + 1), Kac type (2; 1, …, 1): the Fricke point. This is the point at which the level-1 (WZW)
  modular data S_{0λ} ∝ Π sin(π(λ + ρ, α)/(1 + h)) live.

So the two swap theorems are the levels 0 and 1, the level of V_N itself. The next shell would be level 2, and there
the swap fails.

The (h+2)-neighbour of N at ρ/(h+2) is still rootless, i.e. it is Λ, for A1²⁴, A2¹², D4⁶, A24, E8³ and D24. But its
coset minimum is |ρ/(h+2)|² = 2h(h+1)/(h+2)², which is below the Fricke value 2 − 2/(h+2). The twining therefore has a
pole of higher order at 0, and is not of Fricke type. Equivalently, the deep hole sees the 24 + k nodes of the affine diagrams and the Fricke point sees
the 24 finite simple roots.

*Only the two nearest cusp types.* The cusps w of type N adjacent to z satisfy −w·ρ = (h/2)|y − ρ_N/h|², where y runs
over the shells of N + ρ_N/h. Shell 0 gives m = h + 1. The next shell, m = h + 2, FAILS the certificate for all 23
types (swap_shell_test.py): the coset theta series already differ below norm 2. So the swap holds
exactly for m = h and m = h + 1.

*But the Fricke twining extends one shell further for the D-types.* The shell-1 Leech points (m = h + 2) still have
a unique minimal vector at the Fricke norm 2 − 2/(h + 2). For the seven D-type lattices A1²⁴ = D2¹², D4⁶, D6⁴, D8³,
D10E7², D12² and D16E8, their twining is Fricke-invariant through q², and the series are 4A, 8A, 12H, 16C, ... at
orders 4, 8, 12, 16, 20, 24 and 32.

In each case some proper multiple is shallow. For A1²⁴, 2c_w is λ/2 with |λ|² = 4, with two minimal vectors ±λ/2
whose twined contributions cancel. These points are not strict Fricke points and not swap pairs. They are the (B)
family of HYP-9368, and they refute the first registered form of that hypothesis (fricke_shell_test.py).

**(6) The Monster corner (CITED).** For X ≠ A24:
- π^X is the Frame shape of a fixed-point-free element g_X ∈ Co_0, so Z_{V_Λ}(ĝ_X) = 1/η_{π^X} = T^X as well (CDH
  §2.3);
- Γ^X is a monstrous non-Fricke group;
- Carnahan's orbifold duality (2017, Cor. 4.6) applies. His hypotheses are fixed-point-freeness and positive twisted
  vacuum energy E_0 = −(1/24) Σ_a m_a/a, which holds componentwise here; for example D_n contributes (1/2 + 1/h)/24.
  It gives V_Λ^{orb(ĝ_X)} ≅ V^♮, with dual element g'_X ∈ M of McKay–Thompson series T^X + k. His appendix lists the
  22 reciprocal eta products with the classes 2B, 3B, 4C, 5B, 6E, 6C, 7B, 8E, 9B, 10E, 10B, 12I, 12B, 13B, 14B,
  16B, 18D, 18C, 22B, 30G, 30A and 46AB. These class names were read through a summarising tool and should be checked
  against the source before being relied on.

So each X ≠ A24 has four graded traces that agree in every weight n ≥ 2:

| VOA, automorphism | Z(τ) | weight-1 trace |
|---|---|---|
| V_N, σ_N (inner, principal) | T^X | −k |
| V_Λ, g_c (inner, deep hole) | T^X + 24 + k | 24 |
| V_Λ, ĝ_X (lift of g_X ∈ Co_0) | T^X | −k |
| V^♮, g'_X ∈ M | T^X + k | 0 |

**A24 is the exception, but only on the deep-hole side.** Its Frame shape 25·1^{−1} is the only Coxeter Frame shape
of a Niemeier lattice that is not a Conway-group Frame shape, and Γ_0(25) is the only lambency group that is not
monstrous (CDH). So the A24 hole twines the Leech VOA by the Γ_0(25) Hauptmodul with no Monster partner. Its Fricke
point of order 26 = 25 + 1, by contrast, twines by T_26B + 24 (10).

## Proof

**(1).** Kostant (1959) showed that the principal element exp(2πiρ/h) of the adjoint group is conjugate to a lift ĉ
of a Coxeter element in the normaliser of the torus. All lifts of c are T-conjugate, because 1 − c is invertible on
Lie(T). These elements act on the level-1 lattice VOA V_Q ≅ L_ĝ(Λ_0) by inner automorphisms, so

    Θ_Q^{[ρ/h]}/η^ℓ = Tr_{V_Q}(e^{2πiρ(0)/h} q^{L_0 − ℓ/24}) = Tr_{V_Q}(ĉ q^{L_0 − ℓ/24}).

Since c fixes no nonzero vector, only the Fock space of e^0 contributes to the right side. That contribution is
q^{−ℓ/24} Π_n det(1 − c qⁿ)^{−1} = 1/η_π, using det(1 − cx) = Π_a (1 − x^a)^{m_a}, Σ a m_a = ℓ and Σ m_a = 0.

Poisson summation and η(−1/τ) = √(τ/i) η(τ) then give θ_{Q^* + ρ/h} = |Q^*/Q|^{1/2} Π_a a^{m_a/2} η^ℓ Π_a η(τ/a)^{−m_a}.
Anti-self-duality, Π a^{m_a} = det(1 − c) = |Q^*/Q| and the glue lemma turn this into θ_{Q+ρ/h} = η^ℓ η_π(τ/h).

*Glue lemma.* A special automorphism of the extended Dynkin diagram acts on the fundamental alcove by
x ↦ w x + ω_i. It permutes the affine simple roots, so it fixes ρ/h, the unique point of the alcove at which every
affine simple root takes the value 1/h. Hence w(ρ/h) = ρ/h − ω_i, so the isometry w carries Q + ρ/h onto
Q − ω_i + ρ/h. As i runs over the special nodes, −ω_i runs over Q^*/Q.

Anti-self-duality is a finite check on the five families. principal_cores.py checks it, together with the char
polynomials (sympy), the coset equalities and the identity, through q^29.

**(2).** N + ρ_N/h is the disjoint union of the |N/Q| glue cosets Π_i (Q_i + g_i + ρ_i/h). By the glue lemma each
has theta series Π_i θ_{Q_i + ρ_i/h} = η²⁴ η_{π^X}(τ/h). The twined form follows by Poisson summation, since N is
unimodular and the anti-self-duality holds componentwise. coxeter_frame_twined.py confirms the twined form by direct
enumeration of N through q³ for all 23 types.

**(3).** Combine (2) with THM-4758 (2), θ_{Λ+c} = θ_{N+ρ_N/h} + (24 + k)Δ, and with its Poisson dual
Θ_Λ^{[c]} = Θ_N^{[ρ_N/h]} + (24 + k)Δ. Divide by Δ = η²⁴. The twisted-module form uses
T^X(−1/τ) = |N/Q| η_{π^X}(τ/h), which is the computation in (1) for N. twined_hauptmodul.py enumerates the Leech
vectors weighted by e^{2πiλ·c} for the 8 pure A types and finds 24 + T^X + k through q⁴.

**(4).** By Ligozat's formula, the order of the eta quotient T^X at every cusp 1/d of Γ_0(h) is

    (h/24) Σ_a gcd(d,a)² (−m_a) / (gcd(d, h/d) d a).

This shows that the poles are exactly the S_X-orbit of ∞, all simple, and that S_X acts freely on that orbit.
Invariance under W_e holds exactly when the exponents satisfy m_{a∗e} = m_a, where a ∗ e = ae/gcd(a,e)². The sign is
fixed numerically at 30 digits. Then T^X descends to X_0(h)/S_X with a single simple pole, so the quotient has genus 0
and T^X is a Hauptmodul. hauptmodul_cert.py runs this certificate for all 23 types and recovers the CDH lambencies.

**(5′).**

*Left inequality.* By Kostant, σ^d = e^{2πi dρ(0)/h} is conjugate in Aut(V_Q) to ĉ^d. The element c^d has
eigenvalues e^{2πi m_i/e} on the Cartan. Every irreducible σ^d-twisted V_Q-module V_{Q + w + ρ/e}, w ∈ P/Q, is
therefore a ĉ^d-twisted module. Its conformal weight is at least the vacuum anomaly
(1/4e²) Σ_i J_e(m_i) (Bakalov–Kac, Dong–Lepowsky).

*Right inequality (the exponent inequality).*
- **A_n.** The exponents 1, …, n run d = h/e times through every residue mod e, except that 0 is hit once fewer, and
  J(0) = 0. So Σ J = d·e(e² − 1)/6, which exceeds (de − 1)e(e − 1)/6 by (e − 1)(d + 1)e/6 > 0 for e > 1.
- **D_n, h = 2(n − 1), e odd.** Then e | n − 1. The odd exponents 1, 3, …, 2n − 3 hit every residue (n − 1)/e times,
  and the extra exponent n − 1 is ≡ 0. So Σ J = (n − 1)(e² − 1)/6. Its excess over ne(e − 1)/6 is
  (e − 1)(n − 1 − e)/6 ≥ 0, with equality iff e ∈ {1, n − 1}.
- **D_n, e = 2e′ even.** The odd exponents hit every odd residue (n − 1)/e′ times, giving
  (n − 1)(2e′² + 1)/3, plus e′² if (n − 1)/e′ is odd. This exceeds ne′(2e′ − 1)/3 by
  (e′/3)(t + e′t + 1 − 2e′ + 3[t odd]e′) > 0, where t = (n − 1)/e′. (For t ≥ 2 the bracket is ≥ 3, and for t = 1 it
  is 2 + 2e′.)
- **E_6, E_7, E_8.** Finite check.
- exponent_law.py checks (A_1..A_200, D_4..D_200, E) and checks the left inequality by exact minimisation over every
  coset (A_1..A_8, D_4..D_10, E).

**(5).**

*Lower bounds.* Summing (5′) over the components (Σ ℓ_i = 24) gives min|N + dρ_N/h|² ≥ 2 − 2d/h. The inequality is
strict unless d is a bad prime ∥ h of every component.

Norms in N + dρ_N/h lie in (2d/h)Z, since x·ρ_N ∈ Z, x² ∈ 2Z and |ρ_N|² = 2h(h + 1) (Venkov). So strictness forces
min ≥ 2. This proves (b) on the Niemeier side, and the bound in (a).

*Attainment in (a).* By the component rule, W_e fixes T^X, so T^X has a simple pole at W_e(∞) = 1/d. By (2), the
expansion of T^X at 1/d is Σ_{y ∈ N + dρ_N/h} e(y·ρ_N/h) q^{|y|²/2}/η²⁴, up to a constant. When gcd(d, e) = 1 the
cusp 1/d has width e, so a simple pole means exponent −1/e: some y has |y|² = 2 − 2d/h.

*Leech side.* The swap σ maps the vertical line of D through the class of dρ_N/h to the horizontal line through the
class of d′c, with gcd(d′, h) = d. By THM-4758 (1), the coset theta series of N + dρ_N/h and Λ + d′c agree below q¹.
The covering radius √2 of Λ forces min|Λ + dc|² ≤ 2, so in case (b) dc is a deep hole.

*Component rule.*
- A_n: h·1^{−1} ↦ (h/e)·e^{−1} under W_e, which is never π for e > 1.
- D_n: h ↦ h/e forces e = n − 1. That e is a Hall divisor iff n − 1 is odd, and then 2 ↔ h and n − 1 ↔ 1.
- E: finite check.
- hauptmodul_cert.py checks the rule for 360 pairs (R, e).

*Verified parts.* The uniqueness of the minimal vector, the counts φ(d)T(d) and the types in (b) are verified by exact
enumeration up to norm 2 + 2/h, on both sides:
- N + dρ_N/h (multiple_cosets.py);
- Λ + dc (deep_hole_multiples.py).

The type of dc in (b) is read off from the vertex count 24 + k′ and the Coxeter number e.

**(7).**

*Setting.* Put F = Θ_Λ^{[v]}/η²⁴. It is Γ_0(e)-invariant: e·q(v) ∈ Z, and the coefficients are rational, so the
diamond operators act trivially. Its expansion at a cusp a/c of Γ_0(e), c | e, is governed by the twisted sector
Λ + cv:
- at c = e this is the vacuum, giving a simple pole at ∞;
- at c = 1 it is the Fricke point itself, with a single vector of norm 2 − 2/e. Since 0 has width e, this is a
  simple pole;
- for 1 < c < e, gcd(cd, h) = cd is not prime, so cv is a deep hole by (5b) and there is no pole.

*Fricke invariance.* G := F∘W_e = θ_{Λ+v}(eτ)/Δ(eτ) has the same two simple poles. Its leading coefficient is 1
because the shallow vector is unique. So F − G is holomorphic at ∞. Since (F − G)∘W_e = −(F − G), it is also
holomorphic at 0, hence constant and odd, hence 0.

*Conclusion.* F descends to X_0(e)/W_e with a single simple pole. So that curve has genus 0 and F − 24 is its
normalized Hauptmodul. The constant term of G then gives the 24 norm-2 vectors.

*Check.* fricke_multiples.py compares Θ^{[v]}/Δ with θ_{Λ+v}(eτ)/Δ(eτ) through q³ for all eleven cases. The series in
the table are its output.

**(10).**

*On the Leech side,* leech_h_plus_one.py builds z″ from the closest-vector normalisation of ρ_N, and checks that z″ is
null with −z″·ρ = h + 1. It enumerates Λ + c″ exactly. fricke_h_plus_one.py checks the Fricke identity and the proper
multiples, through q³ for all 23 types (fricke_h_plus_one.out; through q² in fricke_h_plus_one_q2.out).

*On the N side,* weyl_h_plus_one.py enumerates N + ρ_N/(h+1) and the (h+1)-neighbour, which is rootless and hence Λ.

*Twining.* The Fricke identity follows from the principal-part argument of (7). Its inputs are the unique shallow
vector and the deep proper multiples, both verified.

*Root-lattice part of uniqueness.* ρ/(h+1) has Kac coordinates (2; 1, …, 1), so it lies strictly inside the
fundamental alcove: θ(ρ/(h+1)) = (h − 1)/(h + 1) < 1. That alcove lies in the Voronoi cell of 0 in Q, so 0 is the
unique closest point of Q. The glue cosets are covered by the enumeration only.

**(8).** fricke_neighbors.py constructs Λ′ = Λ_0 + Zw, with Λ_0 = {x ∈ Λ : x·v ∈ Z} and w ∈ v + Λ of even norm. It
enumerates its roots exactly as the norm-2 vectors x ∈ Λ + jv with (x − jw)·v ∈ Z, and splits them into irreducible
components. The result is X_p in all eleven cases. For prime e the identification is PROVED by the inclusion
R(L_e) ⊇ X_p plus the root count 24(e − 1) (see (8)); for composite e it is verified.

**(9).** This combines (3), (5) and (7). The types of the deep-hole multiples come from deep_hole_multiples.py.

## Remarks

- **The order-3 Fricke point is the M24-point.** Let y_0 = 3x_0, where x_0 is the unique norm-4/3 vector of Λ + 2c
  for D4⁶.
  - y_0 is a norm-12 vector with no λ ∈ Λ_4 at inner product 6.
  - Its inner products with Λ_4 have exactly the distribution of the M24-fixed vector (2²⁴)/√8. Computed by hand:
    5: 24, 4: 759, 3: 6072, 2: 21528, 1: 43056, 0: 53682.
  - In a sample of 2977 norm-12 vectors (order3_fricke_census.py) only two such distributions occur.
  - The two orbit sizes |Co_0|/|M24| = 33,965,568,000 and |Co_0|/|U6(2):2| = 452,088,000 add up exactly to the
    34,417,656,000 norm-12 vectors.

  So, up to Co_0 ⋉ Λ, the D4⁶ Fricke point is (2²⁴)/(3√8) and is fixed by M24, and it is the only order-3 point that
  twines V_Λ by T_3A + 24 (HYP-9368).

- **What is new.** (2) gives a VOA meaning to CDH's principal moduli, in the direction asked for by CDH Remark 2.2: T^X
  is the trace of Kostant's principal element on V_N, and of the deep-hole automorphism on V_Λ up to 24 + k.

  (5) and (5′) then go most of the way towards a structural reason for genus zero:
  - the poles of T^X away from ∞ sit at the twisted sectors of weight < 1;
  - by the exponent law those sectors occur exactly at the common bad primes p ∥ h;
  - the Atkin–Lehner involutions W_{h/p} at the same primes fix π^X.

  So T^X descends to X_0(h)/S_X with a single simple pole.

  What remains is case-dependent: the E-type checks, and the fact that S_X contains no element beyond
  ∩_R S_R (checked for the 23 types). A uniform, case-free proof of the exponent law and of the component rule is
  open.
- **27, 3n + 1, {25, 1}.**
  - For A8³ (h = 9 = 3², k = 3), |N/Q| = 27 = 24 + k, so θ_{Λ+c} = 27Δ(1 + η_{9³1^{−3}}(τ/9)) is the balanced case.
    Its Hauptmodul is T_{9B}, and 3c is an A2¹² hole.
  - For A2¹², T^X + 12 = T_{3B}, and θ_{Λ+c} counts 12-tuples of 3-cores. The size n of a 3-core and its abacus
    vector x ∈ A_2 satisfy 3n + 1 = |3x + ρ|²/2, the Eisenstein norm of a primary integer. So Collatz's 3n + 1 is the
    norm form on the deep-hole coset A_2 + ρ/3.
  - {25, 1} is literally the Frame shape 25¹·1^{−1} of A24, the exceptional corner of (6).
- **Least resistance.** ρ/h is the point of the alcove where every affine simple root takes the same value 1/h; it is
  the barycentre only in type A. Kostant's theorem says that the "straight" flow exp(2πitρ), stopped at t = 1/h, is
  conjugate to the product of the reflections in all simple walls. The swap theorem pairs this point of N with the deep
  hole of Λ, the point of Λ ⊗ R farthest from Λ.

  The two cusp types for which the swap holds, m = h and m = h + 1, are reached from z by crossing zero walls and one
  wall of the Leech chamber.
- **Collatz (NUMEROLOGY, recorded honestly).** The Collatz orbit of 3 is 3 → 10 → 5 → 16 → 8 → 4 → 2 → 1. Together
  with the doublings 6, 12 and 24 of 3, every term is the index of a D-component of some Niemeier root system (D3 = A3,
  D2 = A1²).
  - Every halving step in it is a Fricke descent (8): D24 → D12, D12 → D6, D6 → D3, D10 → D5, D16 → D8, D8 → D4,
    D4 → D2.
  - The two odd steps 3 → 10 and 5 → 16 have no lattice meaning that we know of.
  - The D-indices 7 and 9 occur only next to A-components, which have no Fricke points.
