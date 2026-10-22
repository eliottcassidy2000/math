---
id: THM-4758
title: "The swap theorem for the Leech cusp, and the deep-hole theta identity: for each of the 23 Niemeier lattices N with roots (Coxeter number h, Weyl vector rho_N, k components), the vector-valued theta function of M = {rho, z}^perp in II_{25,1} is symmetric under the swap of its two isotropic lines up to Delta times the invariant vector sum_d (T(d) - T(h/d)) 1_{S_d} (T = number of hole vertices at d = h, number of components with the bad prime p at d = p); consequently theta(Leech + deep hole) - theta(N + rho_N/h) = (24 + k) Delta (HYP-9363), and for the pure A-types theta(Leech + deep hole) = h^(k/2) (eta(tau)^h / eta(tau/h))^k + (24 + k) Delta, counting Leech points around the hole by h-core partitions"
status: >
  PROVED for each of the 23 types by a finite exact certificate (swap_theorem_proof.py).
  - The inputs are Borcherds' theorem that the vector-valued theta function of a lattice is a modular form for the
    Weil representation (Duke 1998, Thm 4.1), and a vector-valued Sturm lemma (proved below from the valence formula).
  - The cusp vectors come from Borcherds' table of norm-0 vectors of II_{25,1}
    (https://math.berkeley.edu/~reb/lattices/norm0); its N_0, N_1, N_2 columns are reproduced.
  - Independently corroborated by exact enumeration of both sides through q^3 for all 23 types (all_holes.py) and, for
    A_24, A_1^24, A_2^12, deep into the expansion (S52 hyp9363_holes.py).
  - The pure-A eta-quotient form uses Garvan-Kim-Stanton (Invent. Math. 101 (1990)): the theta series of
    A_{h-1} + rho/h is the h-core generating function.
  - The defect formula c_d = T(d) - T(h/d) is PROVED by exhaustion over the 23 types (FINITE-EXACT, complete). Its
    conceptual explanation is OPEN.
  Prior-art search: see the note. Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/swap_theorem_proof.py (+ .out) (the proof certificates)
  - 04-computation/experiments/deep_holes_theta_20261022/all_holes.py (+ all_holes_J2.out) (independent enumeration)
  - 04-computation/experiments/deep_holes_theta_20261022/swap_reciprocity.py, swap_dimension.py (structure, dimensions)
related:
  - HYP-9363 (stated in S52; proved here)
  - THM-4752 (the A_24 cusp nu = (1/2^25 | 5/2) in tournament coordinates)
  - THM-4759 (Coxeter Frame shapes, umbral twining, Fricke points, Fricke swap theorem)
  - HYP-9368 (genus-zero inner automorphisms of V_Lambda)
---

# THM-4758 — the swap theorem for the Leech cusp

## Setting

Let II = II_{25,1} and let ρ be the Weyl vector (the Leech cusp). Fix a Niemeier lattice N with roots, and let z be
its cusp from Borcherds' table: a primitive null vector in the closed Leech chamber with z·ρ = −h and z^⊥/z ≅ N.
Write:
- R for the root system of N, k for its number of components, and ρ_N for its Weyl vector;
- c for the corresponding deep hole of the Leech lattice Λ; it has 24 + k vertices;
- M = {ρ, z}^⊥, a positive definite even lattice of rank 24.

Then D = M^*/M ≅ (Z/h)² is hyperbolic. In coordinates (a, b),
- q(a, b) = ε·ab/h for a unit ε;
- the two isotropic lines are H_Λ = Λ_M/M = {(a, 0)} and H_N = N_M/M = {(0, b)}, where Λ_M and N_M are the
  projections of ρ^⊥ and z^⊥;
- the maximal isotropic subgroups are S_d = dZ × (h/d)Z for d | h, with S_1 = H_Λ and S_h = H_N. Each S_d defines
  an even unimodular overlattice L_d of M, a Niemeier lattice, with L_1 = Λ and L_h = N.

The swap is σ(a, b) = (λb, a/λ). It is an isometry of D exchanging H_Λ and H_N. The unit λ is chosen so that σ maps
the Leech line A (projections of {x : x·ρ = −1}) onto the Niemeier line B (projections of {x : x·z = −1}).

Θ_M = Σ_γ θ_{M+γ} e_γ is the vector-valued theta function.

## Statement

**(1) Swap theorem.**

    Θ_M − σ*Θ_M = Δ · Σ_{d | h} c_d 1_{S_d},     c_d = T(d) − T(h/d),

where T is defined on the divisors of h by:
- T(1) = 0;
- T(h) = 24 + k, the number of vertices of the deep hole;
- for a prime p | h with p < h, T(p) is the number of components of R for which p is a bad prime, i.e. divides a
  coefficient of the highest root. Equivalently it is a torsion prime of the simply connected group: D_n has 2,
  E_6 and E_7 have 2 and 3, E_8 has 2, 3 and 5, and A_n has none;
- T(d) = 0 for every other divisor.

In particular θ_{M+γ} = θ_{M+σγ} for every non-isotropic class γ.

**(2) Deep-hole theta identity (HYP-9363).**

    θ_{Λ+c}(τ) − θ_{N+ρ_N/h}(τ) = (24 + k) Δ(τ).

Equivalently, for every j ≥ 0,

    #{λ ∈ Λ : |λ − c|² = 2 + 2j/h} − #{x ∈ N : |x − ρ_N/h|² = 2 + 2j/h} = (24 + k) τ(1 + j/h) [h | j].

In terms of the roots of II_{25,1}:

    #{r : r·ρ = −1, r·z = −j} − #{r : r·z = −1, r·ρ = −j} = (24 + k) τ(1 + j/h) [h | j].

**(3) Pure A-types.** For R = A_{h−1}^k with k(h − 1) = 24, so h ∈ {2, 3, 4, 5, 7, 9, 13, 25}:

    θ_{Λ+c}(τ) = h^{k/2} (η(τ)^h / η(τ/h))^k + (24 + k) Δ(τ).

In q^{1/h}, the function η(τ)^h/η(τ/h) is the generating function of h-core partitions. Hence:
- **A24.** The number of Leech points at squared distance 2 + 2j/25 from the A24 hole is 5p(j−1) for 1 ≤ j ≤ 24,
  because every partition of fewer than 25 is a 25-core. At j = 25 it is 5p(24) − 600 = 7275.
- **A1²⁴.** The count is 4096 × (the number of ways to write j − 1 as a sum of 24 triangular numbers), plus 48τ.
- **A2¹².** The count is 729 × (the 12-fold convolution of a₃(n) = d_{1,3}(3n+1) − d_{2,3}(3n+1)), plus 36τ.
  So the Leech lattice seen from its A2¹² deep hole is governed by the divisors of 3n + 1 (Granville–Ono). Here
  729 = 3⁶ is the size of the ternary Golay code, the glue code of N(A2¹²).

**(3′) Twined form: a Hauptmodul.** Write Θ_L^{[x]}(τ) = Σ_{y∈L} e^{2πi y·x} q^{|y|²/2}. Applying τ → −1/τ to (2)
gives Poisson duals, since Λ and N are unimodular:

    Θ_Λ^{[c]} − Θ_N^{[ρ_N/h]} = (24 + k) Δ    (all 23 types).

For the pure A-types, the S-transform of Garvan–Kim–Stanton gives Θ_N^{[ρ/h]} = η(τ)^{24+k}/η(hτ)^k. Hence

    Θ_Λ^{[c]}(τ)/η(τ)²⁴ = 24 + T_h(τ),   T_h(τ) = (η(τ)/η(hτ))^{24/(h−1)} + 24/(h−1).

Here T_h = q^{−1} + 0 + O(q) is the normalized Hauptmodul of Γ_0(h). For h = 2, 3, 4, 5, 7, 9 and 13 it is the
Monster McKay–Thompson series of class 2B, 3B, 4C, 5B, 7B, 9B and 13B (Conway–Norton). For h = 25 it is the Γ_0(25)
Hauptmodul.

In vertex-algebra terms, the deep-hole automorphism g_c = e^{2πi c(0)} of the Leech lattice VOA has graded trace

    Tr_{V_Λ}(g_c q^{L_0 − 1}) = 24 + T_h(τ),

the Leech character 24 + J with J replaced by T_h. Equivalently, its twined trace coincides with that of a
fixed-point-free Conway element of Frame shape h^k/1^k, plus the constant 24 + k, which lives on the weight-one
space.

This is checked independently by direct enumeration of Leech vectors weighted by e^{2πiλ·c}, through q⁴, for A1²⁴,
A2¹², A4⁶, A6⁴ and A24 (twined_hauptmodul.py).

**(4) Divisor lattice.** The root counts R(L_d) = 24h_d of the intermediate Niemeier lattices are tabulated in
swap_theorem_proof.out. Often h_d = d. For the D- and E-types it drops by one; for example D24 gives L_23 with h = 22,
and E8³ gives L_15 with h = 14.

This is explained in THM-4759 (8). L_d is the orbifold of V_Λ by g_c^{h/d}:
- h_d = d − 1 exactly when h/d is a bad prime of every component, with h/d ∥ h (Fricke descent: the roots of height
  divisible by h/d);
- otherwise h_d = d, and (h/d)c is a deep hole of an A-containing type.

## Proof

**Lemma (vector-valued Sturm bound).** Let ρ be a unitary representation of SL_2(Z) with finite image, and let
F ∈ M_k(ρ). If every component of F vanishes to order > k/12 at ∞, then F = 0.

*Proof of the lemma.* Each component of F|_k g is a linear combination of components of F. For each j put

    N_j = Π_{g ∈ Γ_ρ \ SL_2(Z)} F_j|_k g,

a scalar modular form for SL_2(Z) of weight k·[SL_2(Z) : Γ_ρ]. Its order at ∞ exceeds its weight divided by 12, so
the valence formula gives N_j = 0. Some factor then vanishes, and so F_j = 0.

**Proof of (1).** Θ_M lies in M_12(ρ_M) (Borcherds 1998, Thm 4.1; the signature is 24, so there is no phase).
Because σ ∈ O(D) commutes with ρ_M, σ*Θ_M is in M_12(ρ_M) too. Let v_γ be the q¹ coefficient of
E = Θ_M − σ*Θ_M, i.e. the difference of the numbers of norm-2 vectors in M + γ and M + σγ.

For each type, swap_theorem_proof.py checks exactly:
- (i) E has no coefficient at exponents < 1, i.e. no difference at norms < 2. This uses all vectors of M^* of norm
  ≤ 2, between 49 and 1149 per type.
- (ii) v is ρ_M-invariant. For T, v is supported on isotropic classes. For S, the identity
  (1/h) Σ_δ e(−b(γ, δ)) v_δ = v_γ holds in Z[ζ_h], verified by reduction modulo Φ_h.

Then F = E − Δ·v lies in M_12(ρ_M) and vanishes to order > 1, so F = 0 by the lemma. The decomposition of v into the
1_{S_d} is computed exactly and matches c_d = T(d) − T(h/d) for all 23 types.

**Proof of (2).** The affine hyperplane {x ∈ II : x·ρ = −1} = r_0 + ρ^⊥, with r_0 = e_0 − e_1, projects onto the coset
A = π_M(r_0) + Λ_M. Through the Leech roots r_λ = λ + ρ′ + (λ²/2 − 1)ρ, which satisfy
r_λ·z = h(1 − |λ − c|²/2), this coset is isometric to Λ + c. So θ_{Λ+c} = Σ_{γ ∈ A} θ_{M+γ}.

Likewise {x : x·z = −1} projects onto B = σ(A), and θ_{N+ρ_N/h} = Σ_{γ ∈ B} θ_{M+γ}; the line positions are
checked by the script. Hence

    θ_{Λ+c} − θ_{N+ρ_N/h} = Σ_{γ ∈ A} E_γ = (Σ_{γ ∈ A} v_γ) Δ.

The script finds Σ_{γ ∈ A} v_γ = 24 + k for every type. Conceptually this is −c_1, the number of hole vertices.

**Proof of (3).** It follows from (2) together with three facts:
- Garvan–Kim–Stanton: θ_{A_{h−1}+ρ/h} = η(τ)^h/η(τ/h);
- the rotation of coordinates (an element of W(A_{h−1})) maps A_{h−1} + g + ρ/h to A_{h−1} + (g − [1]) + ρ/h, so
  every shifted glue coset has the same theta series;
- the glue code has order h^{k/2}.

All eight cases are also checked against exact enumeration (all_holes_J2.out).

## Remarks

- Representation theory alone forces (1) only for prime h. For prime h, the torus-invariant, σ-antisymmetric weight-12
  cusp space has dimension 1. For h = 4, 9 and 25 it is larger, of dimension 2, 3 and 5 (swap_dimension.py). So the
  vanishing of the extra components is arithmetic, and is certified case by case.
- For prime h there is also a two-line proof. W(R) realizes the torus on D, and the sum rule
  Σ_v R_v = Σ_u C_u = θ_{M^*} gives R(1) − C(1) = (θ_N − θ_Λ)/(h − 1) = 24hΔ/(h − 1) = (24 + k)Δ.
- Dividing (2) by η²⁴ gives a statement about the deep-hole orbifold of the Leech VOA (Möller–Scheithauer): the
  characters of the c-twisted V_Λ-module and of the (ρ_N/h)-twisted V_N-module differ by the constant 24 + k.
- **Conceptual reading of the proof.** The coset characters ch_{(i,j)} = θ_{M+γ(i,j)}/η²⁴ are the sector characters
  of the Z_h × Z_h orbifold that connects V_Λ and V_N. A weight-0 vector-valued modular function is determined by its
  principal part up to an invariant constant vector. The principal parts are exactly the vectors of M^* with norm < 2.
  So the swap theorem says: the principal parts are swap-symmetric, and therefore the characters are swap-symmetric up
  to constants. The constants are the weight-1 dimension differences, and they are what the bad-prime formula
  computes.
- **The Leech cusp is special.** The same certificate applied to arbitrary pairs of cusps of II_{25,1} passes for 275
  of the 276 pairs of Borcherds' representatives. It fails for A4⁶ × E8³ (m = 7), whose M^* has asymmetric vectors of
  norm 12/7 (cusp_pairs.py).

  Around the Leech cusp it holds for exactly the two nearest cusp types of each N:
  - m = h (this theorem);
  - m = h + 1, the cusp z″ = z + r one Leech wall away. This is the Fricke swap theorem, THM-4759 (11), with
    defect 24.

  It fails at the next shell, m = h + 2 (swap_shell_test.py).
