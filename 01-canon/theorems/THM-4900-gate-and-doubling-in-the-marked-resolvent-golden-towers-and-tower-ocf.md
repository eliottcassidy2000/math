---
id: THM-4900
title: "The strong gate and the skew-Sylvester doubling in the marked resolvent of a tournament, and deep gate towers: the gate acts linearly on (N, P) = (u^T adj(zI-M) u, det(zI-M)) with eigenvalues (sqrt(z^2+4) +- 1)^2 (the padding E is the same with sqrt(z^2+4) -> z); nested gates have closed-form characteristic polynomials, an explicit spectrum over one vertex, and a universal resolvent limit (2 + sqrt(z^2+4))/z, the metallic mean of index 4/z; z = 1 is the self-similar point of the doubling, where THM-1950's inverse response s(T) = F_T(1) = u^T (I + M^T M)^(-1) u is an exact doubling invariant and the gate is the golden map s -> (5s+1)/(s+1), whose orbit from one vertex is F_(2k+2)/F_(2k-1) -> phi^3; a tower form of the odd-cycle formula, H(Gamma^r T) = sum_j sw_j(T) H(Gamma^r TT_(n-2j)); and deep towers forget the core except through the cyclic triangles: rho_(d,r) = H(Gamma^r TT_(d-2))/H(Gamma^r TT_d) <= 1/(r+1), so H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) + O(rho_(n,r)/r)"
status: >
  PROVED: statements 1-7 (exact linear algebra, Mobius dynamics, the cover-polynomial operator and an induction on its
  tower weights).
  FINITE-EXACT companions: every identity checked exactly in 04-computation/experiments/microcosm_gate_resolvent_20261106/
  (sympy on random cores of orders 1-4, numerical spectra to k = 10, tower values to total order 1201; statement 7 on
  all r <= 120 and on all pairs of cores of orders 4-6).
  KNOWN ingredients, credited:
  - the marked pair (p, F), the padding E as a Mobius boundary transfer, the doubling D as the substitution
    w^2 = (z^2+1)/2, the Cayley linearization of E and its all-depth law C_(E^k) = ((z+1)/(z-1))^(2k) C
    (session35_dimensions_dualities_20261023/recursions/RESULT.md);
  - the gate Gamma, its recognizability, HC(Gamma T) = H(T), the cover-polynomial operator
    P_Gamma(q) = q[P(q+1) + qP(q) - P(q-1)] and H(Gamma^r T) = (L^r P_T)(1) (Session210, graded-spectrum-gates/PROOF.md);
    the operator was first written in Session150 transfer/TRANSFER.md and Session159 profile_operator/NOTE.md (1)-(2)
    (2026-11-02/03), with the one-vertex family values 3, 15, 123, ..., 496479;
  - the size-weighted odd-cycle form of P_T at every order, sw_j = sum 2^|K| over collections K of vertex-disjoint
    directed odd cycles with sum (|C|-1) = 2j, so sw_1 = 2 c3 and sw_2 = 2 c5 + 4 (disjoint triangle pairs)
    (Session159 cycle_profiles/NOTE.md (1), 2026-11-03; THM-078 is the odd-n case);
  - the transitive-core towers H(Gamma^r TT_n) = sum_q c_(r,q) q^n with
    c_(r+1,q) = (q-1)c_(r,q-1) + q^2 c_(r,q) - (q+1)c_(r,q+1) (Session159 lagrange/NOTE.md (5), (7));
  - s = F_T(1) is THM-1950's total inverse response 1^T(I+K)^(-1)1, with 0 < s <= n, s = n iff regular, the
    velocity-addition law under series composition, and the strong minima for n <= 8 (THM-1950, klein-2026-07-21-S400;
    THM-4051 at order 7; THM-3729: det(I+K) s = odd rooted-Pfaffian square sum);
  - the series (order-join) law at every z, i.e. P +- zN multiply (THM-3322, in the coordinates det(I - zC), C = 2A - J);
  - the skew-Sylvester double D (THM-447);
  - D(Gamma(T)) has no gate, i.e. the tournament-level non-commutation in statement 5 (Session211
    recursion-bridge/PROOF.md section 5; CONCURRENT, pushed earlier the same day; re-proved here in one line).
  New here:
  - the gate's linear transfer on the marked pair and its eigen-coordinates;
  - the closed-form gate towers, their universal (metallic) limit and the explicit one-vertex spectrum;
  - z = 1 as the self-similar point of the doubling: s(D T) = s(T), and the golden gate map with its Fibonacci orbit;
  - the scalar shadow of the gate/doubling interaction: commutation at z = 1, failure of the BS relation there, and
    P_(D Gamma T) = P_(D E T);
  - the tower odd-cycle formula (6) as stated (a one-line consequence of the KNOWN ingredients);
  - statement 7: the operator form of the tower weights (KNOWN in content, see 7(a)), their monotonicity, the
    forgetting bound rho <= 1/(r+1), and the first-order law with an explicit error.
  AUDIT: statements 1-6 audited by audit H1 (05-knowledge/results/microcosm_audit_H1_20261106.md): PASS WITH
  CORRECTIONS, all applied (credits to THM-1950 and Session159; MISTAKE-729). Statement 7 was added after audit H1
  and checked by its follow-up (report section 7): PASS, no mathematical errors, three cosmetic fixes applied; no
  earlier occurrence of rho <= 1/(r+1), the first-order law or the sandwiches was found in the repository.
source: mac-mini-2026-11-06-microcosm (owner directive 2026-11-06: micro/macrocosm; the duality of x+2 and x*2 recursion on tournaments; an n-vertex tournament as a strong (n+2)-tournament with a recognizable gate)
depends_on: [THM-447, THM-291, THM-078, THM-1950, THM-3322, THM-1370, session35 recursions RESULT.md, session210 graded-spectrum-gates PROOF.md, session211 recursion-bridge PROOF.md, session159 cycle_profiles + profile_operator + lagrange NOTE.md]
verification:
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/gate_resolvent_exact.py (+ .out): statements 1, 3, 4 exactly (sympy), Fibonacci orbit
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_spectrum.py (+ .out): statement 2 (closed form, explicit spectrum)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_limits.py (+ .out): statements 2, 3, 5 numerically
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_ocf.py (+ .out): statement 6 exactly
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_weights_diffusion.py (+ .out): statement 7(a), 7(b)1, 7(d) sandwich
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_forgetting_bound.py (+ .out): statement 7(b)2 and 7(c)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/run_all.sh: reproduces every output
---

# THM-4900 — the gate and the doubling in the marked resolvent, and deep gate towers

## Setting

- **Tournament data.** For a tournament T on n vertices, let M = A − Aᵀ be its dominance matrix (skew-symmetric,
  ±1 off the diagonal) and u the all-ones vector. The *marked pair* is

      P_T(z) = det(zI − M),    N_T(z) = uᵀ adj(zI − M) u,    F_T = N_T/P_T = uᵀ(zI − M)⁻¹u.

  This is Session35's (p, F); this file uses the letter P for p. In statements 6 and 7, P_T(q) instead denotes the
  ordered path-cover polynomial; the variable (z or q) tells them apart.
- **Operations.**
  - **Series composition** A ⇒ B: every vertex of A beats every vertex of B. This is THM-3322's "order join".
  - **Padding** E(T): s → T → t with s → t. This is THM-291's canonical Mode-B section.
  - **Gate** Γ(T): the same with t → s. This is Session210's a → b → T → a: the arc t → s is the recognizable gate,
    d⁺(t) = d⁻(s) = 1. By Session210's convention Γ(∅) = TT_2.
  - **Skew-Sylvester double** D(T), with dominance [[M, M + I], [M − I, −M]] (THM-447).
- **The x+2 / x·2 grid.** On orders, E and Γ act as n ↦ n + 2 and D as n ↦ 2n. These satisfy the Baumslag–Solitar
  relation d∘a = a²∘d (Session35 §1).
- **Notation.** φ = (1 + √5)/2 is the golden ratio. The tournament invariant of statement 3 is called s(T), as in
  THM-1950.

Write s_z = √(z² + 4) in statements 1–2 (not to be confused with the invariant s(T)).

## Statement

**1. The gate is linear on the marked pair.**

    (N, P)_Γ(T) = G(z)·(N, P)_T,    G(z) = [[z²+9, 2z], [2z, z²+1]].

The padding E has G_E(z) = [[z²+1, 2z], [2z, z²+1]] (Session35). Consequences:
- P_Γ(T) = P_T·(z² + 1 + 2zF_T) = P_E(T): flipping the apex arc does not change the characteristic polynomial.
- F_Γ(T) = ((z²+9)F_T + 2z)/(z² + 1 + 2zF_T).
- **Eigen-coordinates.** G(z) is symmetric, with eigenvalues (s_z ± 1)² and eigen-coordinates

      Ψ± = zN + (−2 ± s_z)P,    Ψ±(ΓT) = (s_z ± 1)²·Ψ±(T).

  Setting s_z ↦ z and dropping the −2 gives E's coordinates N ± P, multiplied by (z ± 1)². That is Session35's Cayley
  linearization C = (1+F)/(1−F) ↦ ((z+1)/(z−1))²·C.

**2. Nested gates.** For every core T and k ≥ 0:

    P_(Γ^k T)(z) = [(s_z+1)^(2k)·Ψ+(T) − (s_z−1)^(2k)·Ψ−(T)] / (2s_z).

- **Universal limit.** For every real z > 0 and every core T,

      F_(Γ^k T)(z) → (2 + s_z)/z,

  with error O(((s_z−1)/(s_z+1))^(2k)). At z = 4/m the limit is the metallic mean (m + √(m²+4))/2:
  - the golden ratio φ at z = 4;
  - silver 1 + √2 at z = 2;
  - 2 + √5 = φ³ at z = 1.

  The macro limit forgets the core.
- **Over one vertex.** For k ≥ 1,

      P_(Γ^k K1)(z) = z(z²+3)·[(s_z+1)^(2k−1) + (s_z−1)^(2k−1)] / (2s_z).

  The spectrum of Γ^k(K1), which has 2k+1 vertices, is {0, ±i√3} together with

      ±i·√(4 + cot²((2j+1)π / (2(2k−1)))),    j = 0, …, k−2.

  Every non-core eigenvalue has modulus above 2. For example, the 5-vertex tower has spectrum {0, ±i√3, ±i√7}.

**3. The self-similar point of the doubling.** The substitution w² = (z²+1)/2 fixes z = ±1. So

    F_(D T)(±1) = F_T(±1),    P_(D T)(1) = 2ⁿ P_T(1)².

Put s(T) := F_T(1) = uᵀ(I − M)⁻¹u = uᵀ(I + MᵀM)⁻¹u. This is THM-1950's total inverse response 𝟙ᵀ(I + K)⁻¹𝟙 with
K = M (a scalar equals its transpose).
- 0 < s(T) ≤ n, with equality iff T is regular (KNOWN: THM-1950).
- **s(D T) = s(T) exactly** (new).
- **Towers, for real z > 0.**
  - D^k: F_(D^k T)(z) = (z/w_k)·F_T(w_k), where w_k² − 1 = (z² − 1)/2^k. Hence F_(D^k T)(z) → z·s(T): the doubling
    tower remembers exactly one number of its core.
  - E^k: F_(E^k T)(z) → 1, for z ≠ 1.
  - Γ^k: the limit (2 + s_z)/z of 2. The gate tower forgets its core entirely.

  The D- and E-tower limits are one-line corollaries of Session35's F_D formula and its all-depth law for E.

**4. The algebra at z = 1.**
- **Series composition is relativistic addition.**

      F_(A⇒B) = (F_A + F_B)/(1 + F_A F_B)    at every z,

  equivalently, the Cayley transforms multiply (THM-3322's join law). In particular
  s(A ⇒ B) = (s_A + s_B)/(1 + s_A s_B), which is THM-1950's velocity-addition law for s.
  - Since s(K1) = 1 is absorbing, s(T) = 1 whenever T has a one-vertex strong component: transitive tournaments, every
    E(T), every T with a source.
  - s is the relativistic sum of the values of s on the strong components. For example s(C3 ⇒ C3) = 3/5.
- **The gate is the golden map** (new).

      s(Γ T) = (5s(T) + 1)/(s(T) + 1).

  This is hyperbolic, with fixed points φ³ = 2 + √5 (attracting) and −φ⁻³, and multiplier φ⁻⁴.
  - In the coordinate y = (s − 1)/2 (the normalization of q = (H − 1)/2 for path counts), it is y ↦ (2y+1)/(y+1).
    That is the square of the continued-fraction map y ↦ 1 + 1/y.
  - From the single vertex (y = 0), the orbit runs, for k ≥ 1, through the convergents of φ that lie below φ:

        y(Γ^k K1) = F_(2k)/F_(2k−1),    s(Γ^k K1) = F_(2k+2)/F_(2k−1)  (1, 3, 4, 21/5, 55/13, 72/17, …) → φ³.

**5. x+2 against x·2.**
- At z = 1, D acts trivially. So s(DΓT) = s(ΓDT) for every T, but s(DΓT) ≠ s(Γ²DT) for every T. The integer relation
  d∘a = a²∘d fails here and commutation holds instead.
- At the tournament level, D(Γ(T)) never has a gate. Hence D(Γ(T)) ≇ Γ(X) for every X, and in particular
  DΓ ≠ Γ²D (compare Session35's DE ≠ E²D). KNOWN: Session211 `recursion-bridge` PROOF §5, pushed earlier the same
  day.
- P_(DΓT) = P_(DET).
- On the marked pair, the +2 operations are linear (G_E, G_Γ), while ×2 is quadratic with a substitution:
  P_D = 2ⁿP(w)² and N_D = 2ⁿ(z/w)N(w)P(w).

**6. The tower odd-cycle formula.** Write P_T(q) = Σ_j sw_j(T)·q^(n−2j) for the ordered path-cover polynomial
(Session210 Lemma 3). Session159 cycle_profiles/NOTE.md (1) identifies sw_j at every order as Σ 2^|K| over
collections K of vertex-disjoint directed odd cycles with Σ_(C∈K)(|C|−1) = 2j; THM-078 is the odd-n case.
- sw_0 = 1, and sw_1(T) = 2c3(T) for every n. This is KNOWN (Session159 (1)); an independent elementary proof is
  below.
- For every r ≥ 0,

      H(Γ^r T) = Σ_j sw_j(T)·H(Γ^r TT_(n−2j)),

  where TT_m is the transitive tournament on m vertices and P_(TT_m) = q^m. The core enters only through its
  odd-cycle statistics. The tower enters only through universal weights U_(d,r) := H(Γ^r TT_d) (for d = 0, with the
  convention Γ(∅) = TT_2).
- At r = 0 every weight is 1 and this is the odd-cycle formula H(T) = Σ_j sw_j(T).

**7. Deep gate towers forget the core, except through the cyclic triangles.** Let L f(q) = q[f(q+1) + q f(q) − f(q−1)]
be Session210's gate operator, so that U_(d,r) = (L^r q^d)(1). Put

    ρ_(d,r) = U_(d−2,r)/U_(d,r)    (d ≥ 2).

- **(a) Operator form of the weights** (KNOWN in content; new as an operator). Let

      M = ∂² + 2 sinh(y)·∂ = e^(−2 cosh y)·∂·e^(2 cosh y)·∂    (∂ = d/dy),

  and pair polynomials with entire functions by ⟨P, Φ⟩ = [P(∂)Φ](0). Then ⟨LP, Φ⟩ = ⟨P, MΦ⟩, so for every tournament
  T and r ≥ 0

      H(Γ^r T) = [P_T(∂)·M^r e^y](0),    U_(d,r) = (d/dy)^d [M^r e^y](0).

  - The macro part of a tower is one Sturm–Liouville (diffusion-generator) operator applied r times to e^y. The core
    enters only as the constant-coefficient operator P_T(∂), read at the origin.
  - On exponentials, M e^(ky) = k e^((k+1)y) + k² e^(ky) − k e^((k−1)y). So M^r e^y = Σ_q c_(r,q) e^(qy) with
    Session159's coefficients (lagrange/NOTE.md (7)), and U_(d,r) = Σ_q c_(r,q) q^d is its formula (7) for d ≥ 1 (for d = 0 keep the q = 0
    coefficient, which Session159 discards, with 0^0 = 1).
  - Session159 profile_operator/NOTE.md (2) wrote the gate on the profile side as a second-order operator in z. The
    operator M is the transposed form, on the evaluation side.
- **(b) Monotonicity and the forgetting bound.** For all e ≥ 2 and r ≥ 0:
  1. U_(e,r) ≥ U_(e+2,r−1) (r ≥ 1): at fixed total order N = e + 2r, a deeper tower carries more paths than a
     bigger transitive core.
  2. **U_(e,r) ≥ (r+1)·U_(e−2,r)**, i.e. ρ_(e,r) ≤ 1/(r+1), uniformly in e.
- **(c) First-order law.** Let T, T′ have the same order n, and σ(X) = Σ_(j≥1) |sw_j(X)|. For r + 1 > σ(T′),

      H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r) + E,
      |E| ≤ ρ_(n,r)·[σ(T) + σ(T′) + 2|c3(T) − c3(T′)|·σ(T′)] / (r + 1 − σ(T′)).

  - So deep towers forget their core, except through the cyclic-triangle count, to first order in ρ.
  - The error is O(ρ/r), and it is not 2Δc3·ρ·o(1). When c3(T) = c3(T′), the towers still differ at second order.
    For example, the order-5 classes with P = q⁵ + 8q³ + 2q and q⁵ + 8q³ + 6q both have c3 = 4, but their towers
    differ by −2.2·10⁻⁴ at r = 20 (audit H1).
  - Γ^((n−1)/2)(K1) has order n and c3 = ((n−1)/2)². So the comparison with the one-vertex tower of the same total
    order is the case T′ = Γ^((n−1)/2)(K1).
- **(d) Universality.**

      1 ≤ U_(d,r)/U_(d+2,r−1) ≤ 1 + 2^d/r,    1/(1 + 2^d/r) ≤ ρ_(d,r)/ρ_(d+2,r−1) ≤ 1 + 2^(d−2)/r.

  So the weights, and the forgetting rate, depend asymptotically only on the total order N = d + 2r.

The bound in (b) is ρ ≤ 1/(r+1) ≈ 2/N. The data suggest the true rate is about (log N)²/N² (HYP-9490, OPEN).

## Proofs

**1.** Order the vertices s, core, t. Then

    zI − M_Γ = [[z, −uᵀ, 1], [u, B, −u], [−1, uᵀ, z]],    B = zI − M.

- **Determinant.** The Schur complement of B is [[z + F, 1 − F], [−1 − F, z + F]], with determinant
  z² + 1 + 2zF. This gives P_Γ, and the same computation with the apex sign reversed gives P_E.
- **Response.** Solve (zI − M_Γ)(x, v, y) = (1, u, 1). The core rows give v = (1 − x + y)B⁻¹u, so uᵀv = (1 − b)F with
  a = x + y and b = x − y. The two boundary rows give za − b = 2 and zb + a = 2(1 − b)F. Hence

      b = 2(zF − 1)/(z² + 1 + 2zF),    F_Γ = a + (1 − b)F = ((z²+9)F + 2z)/(z² + 1 + 2zF).

  Multiplying by P_Γ gives N_Γ = (z²+9)N + 2zP and P_Γ = 2zN + (z²+1)P.
- **Eigenvectors.** G(z)·(z, s_z−2)ᵀ = (z³ + 5z + 2zs_z, s_z(z²+1) − 2)ᵀ = (z² + 5 + 2s_z)·(z, s_z−2)ᵀ, and
  z² + 5 + 2s_z = (s_z+1)². The other eigenvector is (z, −s_z−2), with eigenvalue (s_z−1)². Since G is symmetric,
  these are also left eigenvectors.

**2.**
- **Closed form.** By 1, Ψ±(Γ^k T) = (s_z±1)^(2k)·Ψ±(T), and P = (Ψ+ − Ψ−)/(2s_z).
- **Limit.** Let R = Ψ+/Ψ− = (zF − 2 + s_z)/(zF − 2 − s_z). Then R(Γ^k T) = ((s_z+1)/(s_z−1))^(2k)·R(T).
  - For real z > 0 we have s_z > 2, so |(s_z+1)/(s_z−1)| > 1.
  - R(T) ≠ 0, because F_T(z) = z·uᵀ(z²I + MᵀM)⁻¹u > 0 while the root of Ψ+ is (2 − s_z)/z < 0.
  - So R → ∞, i.e. zF − 2 − s_z → 0, at the stated rate.
  - At z = 4/m: s_z = (2/m)√(m²+4) and (2 + s_z)/z = (m + √(m²+4))/2.
- **One vertex.** K1 has N = 1 and P = z, so Ψ+ = z(s_z−1) and Ψ− = −z(s_z+1). Substituting and using
  s_z² − 1 = z² + 3 gives the formula.
  - Its bracket vanishes iff ((s_z+1)/(s_z−1))^(2k−1) = −1, i.e. s_z = −i·cot((2j+1)π / (2(2k−1))), so
    z² = −4 − cot²(…).
  - The value s_z = 0 (j = k−1) is not a root: the bracket divided by 2s_z equals 2k−1 there.
  - Pairs j ↔ 2k−2−j give the same square. This leaves k−1 values and 2(k−1) eigenvalues; with 0 and ±i√3 that is
    2k+1, all of them.

**3.**
- **Fixed point.** w(1) = 1, so the two Session35 formulas give the claims.
- **The form of s.** (I − M)⁻¹ = (I + M)(I − M²)⁻¹, and uᵀM(I − M²)⁻¹u = 0 since M(I − M²)⁻¹ is skew. Hence
  s = uᵀ(I + MᵀM)⁻¹u > 0.
- **Bounds** (as in THM-1950). (I + MᵀM)⁻¹ ⪯ I gives s ≤ uᵀu = n, with equality iff MᵀMu = 0, i.e. Mu = 0 (regular).
- **D-tower.** Iterating Session35's F_D(z) = (z/w)F(w) telescopes to (z/w_k)F(w_k), and w_k → 1. F_T is continuous at
  1, since its poles are imaginary.
- **E-tower.** By Session35, C_E = ((z+1)/(z−1))²·C, and |(z+1)/(z−1)| > 1 for real z > 0 with z ≠ 1. Since
  F_T(z) > 0 (proof of 2), C_T(z) ≠ 0, so C → ∞ and F → 1.

**4.**
- **Series law.** Let A ⇒ B have dominance [[M_A, J], [−Jᵀ, M_B]]. Solving (zI − M)x = u blockwise gives
  x_A = (1 + σ_B)R_A u, where σ_B = uᵀx_B and R_A = (zI − M_A)⁻¹. Likewise x_B = (1 − σ_A)R_B u. So
  σ_A = (1 + σ_B)F_A and σ_B = (1 − σ_A)F_B. Hence

      F = σ_A + σ_B = (F_A + F_B)/(1 + F_A F_B).

  THM-3322 states this law as "P ± zN multiply" in its coordinates. Its f = zN/P equals F_M at (1+z)/z. At z = 1 it is
  THM-1950's velocity-addition law.
- **Absorbing vertex.** F_(K1) = 1/z, so s(K1) = 1, and 1 ⊕ x = 1.
- **Gate map.** Statement 1 at z = 1 gives (10s + 2)/(2s + 2).
- **Conjugation.** With s = 2y + 1: (5(2y+1) + 1)/(2y + 2) = 2·(2y+1)/(y+1) + 1. So the map is y ↦ (2y+1)/(y+1), which
  is the continued-fraction map applied twice.
- **Fibonacci orbit.** By induction,

      y_(k+1) = (2F_(2k) + F_(2k−1))/(F_(2k) + F_(2k−1)) = F_(2k+2)/F_(2k+1).

  These are the convergents F_(m+1)/F_m of φ with m odd, which lie below φ, and y → φ, i.e. s → 2φ + 1 = φ³.
- **Dynamics.** The Möbius matrix [[5, 1], [1, 1]] has eigenvalues 3 ± √5 = 2φ^(±2), so the multiplier is φ⁻⁴.

**5.**
- **Commutation at z = 1.** D is the identity there, so s(DΓT) = g(s(T)) = s(ΓDT), where g is the golden map.
- **No BS relation at z = 1.** g is injective, so g(x) = g(g(x)) forces x to be a fixed point of g, i.e. φ³ or −φ⁻³.
  These are irrational, while s(T) is rational.
- **No gate after doubling.** THM-447(4): the scores of D(X) are 2s_i + 1 and |X| − 1. Since Γ(T) is strong, every
  s_i ≥ 1, so every score is at least 2. A gate needs a vertex of out-degree 1.
- **Characteristic polynomials.** P_D depends only on P, and P_Γ = P_E.

**6.**
- **Coefficient sw_1.** P_T(q) = Σ_π C(q + f(π), n), where f(π) counts forward adjacencies (Session210 (7)).
  - The coefficient of q^(n−2) is Σ_π e_2(f, f−1, …, f−n+1)/n!
    = [n(n−1)·Σ_π f² − n(n−1)²·Σ_π f + K_n] / (2·n!).
  - Σ_π f = (n−1)·n!/2 does not depend on T.
  - In Σ_π f², adjacent position pairs count directed 2-paths, and there are C(n,3) + 2c3 of those. Each contributes
    2(n−2)(n−3)!. All other position pairs give amounts that do not depend on T.
  - So the coefficient is a(n) + 2c3. Since P_(TT_n) = qⁿ (each vertex chooses its bin, and the order inside a bin is
    forced), a(n) = 0.
- **Tower formula.** L is linear and H(Γ^r T) = (L^r P_T)(1) (Session210 Theorem 4).

**7.**
- **(a)** Take P = q^e. Since (q+1)^e − (q−1)^e = 2Σ_(i odd) C(e,i) q^(e−i),

      L q^e = q^(e+2) + Σ_(k≥0) 2C(e, 2k+1)·q^(e−2k),

  so ⟨L q^e, Φ⟩ = Φ^(e+2)(0) + Σ_k 2C(e,2k+1)·Φ^(e−2k)(0). By Leibniz,
  ∂^e(Φ″ + 2 sinh(y)Φ′)(0) = Φ^(e+2)(0) + 2Σ_i C(e,i)·sinh^(i)(0)·Φ^(e−i+1)(0), and sinh^(i)(0) = [i odd]. The two
  agree. Also ⟨P, e^y⟩ = P(1). As a check, r = 1 gives M e^y = e^y + e^(2y) − 1, hence
  H(Γ T) = P_T(1) + P_T(2) − P_T(0) = H(T) + P_T(2), which is Session210 (5).
- **(b)** The displayed L q^e has nonnegative coefficients. Splitting off the first step gives

      U_(e,r) = U_(e+2,r−1) + Σ_k 2C(e,2k+1)·U_(e−2k,r−1),    U_(e,0) = 1,        (R)

  with nonnegative terms, which proves 1. Item 2 is by induction on r.
  - At r = 0 both sides are 1.
  - Suppose U_(h,r−1) ≥ λ·U_(h−2,r−1) for every h ≥ 2, where λ = r. Pascal's rule twice gives
    C(e,j) − C(e−2,j) = 2C(e−2,j−1) + C(e−2,j−2). So (R) splits as

        U_(e,r) = U_(e+2,r−1) + Σ_k 2C(e−2,2k+1)·U_(e−2k,r−1) + S,
        S = Σ_(k≥0) [4C(e−2,2k) + 2C(e−2,2k−1)]·U_(e−2k,r−1).

  - Apply the hypothesis to U_(e+2,r−1), and to every U_(e−2k,r−1) with C(e−2,2k+1) ≠ 0. These are all at heights
    ≥ 2. By (R) at e − 2, this gives U_(e,r) ≥ λ·U_(e−2,r) + S.
  - (R) at e − 2 reads U_(e−2,r) = U_(e,r−1) + Σ_(k≥1) 2C(e−2,2k−1)·U_(e−2k,r−1). Compare it with S term by term:
    - the k = 0 term of S is 4U_(e,r−1) ≥ U_(e,r−1);
    - for k ≥ 1, S contains 2C(e−2,2k−1)·U_(e−2k,r−1).

    So S ≥ U_(e−2,r), and U_(e,r) ≥ (λ + 1)·U_(e−2,r) = (r+1)·U_(e−2,r).
- **(c)** By 6, H(Γ^r X) = U_(n,r) + 2c3(X)·U_(n−2,r) + R_X, with |R_X| ≤ Σ_(j≥2) |sw_j(X)|·U_(n−2j,r). By (b),
  U_(n−2j,r) ≤ U_(n−4,r) ≤ U_(n−2,r)/(r+1) for j ≥ 2, and R_X = 0 for n < 4. So:
  - H(Γ^r T) − H(Γ^r T′) = 2Δc3·U_(n−2,r) + δ, with |δ| ≤ (σ(T) + σ(T′))·U_(n−2,r)/(r+1);
  - H(Γ^r T′) = U_(n,r)(1 + ε′), with |ε′| ≤ σ(T′)·ρ_(n,r) ≤ σ(T′)/(r+1).

  Divide, writing U_(n−2,r) = ρ_(n,r)·U_(n,r). The error is E = (δ/U_(n,r) − 2Δc3·ρ_(n,r)·ε′)/(1 + ε′), which gives the
  bound.
- **(d)** Σ_k 2C(d,2k+1) = 2^d for d ≥ 1 (it is 0 for d = 0, where the ratio is exactly 1), and by (b) every U_(d−2k,r−1) ≤ U_(d,r−1) ≤ U_(d+2,r−1)/r. Substituting in (R) gives
  the first sandwich. For the second, write

      ρ_(d,r)/ρ_(d+2,r−1) = [U_(d−2,r)/U_(d,r−1)]·[U_(d+2,r−1)/U_(d,r)],

  and apply the first sandwich at d − 2 and at d. (For d = 2 the first factor is exactly 1, since L·1 = q².) ∎

## Remarks

- **Micro/macrocosm.**
  - **The condensation multiplies.** The marked response of any tournament is the relativistic sum (product of Cayley
    transforms) of its strong components' responses, just as H is the product of their path counts.
  - **The gate nests.** A whole arbitrary tournament sits inside a strong shell. Iterating forgets it:
    - the response tends to the metallic limit, at rate φ⁻⁴ at z = 1;
    - the path count is a universal tower combination of the core's odd-cycle statistics (6), and its first-order
      trace is the cyclic-triangle count alone (7).
  - **Doubling is invisible at the self-similar point z = 1.** It remembers exactly s.
  - **Two limits of one family.** Session159 lagrange (7) gives H(Γ^r TT_n) ~ r!(r+1)^n for fixed depth r as n → ∞:
    a big core dominates a shallow tower. Statement 7 is the opposite limit, fixed core and r → ∞, where the tower
    dominates and the core is forgotten.
- **Numerology, typed as such.**
  - φ³ = 1 + 2φ is also THM-224's exceptional eigenvalue: COINCIDENCE, no mechanism.
  - The 7-vertex tower has s = 21/5, and the 5-vertex tower has H = 15 and spectrum {0, ±i√3, ±i√7}: NUMEROLOGY with
    respect to {7, 21} and {3, 15}.
- **What does not follow.** The marked pair does not determine H, nor isomorphism. Statements 6 and 7 use the cover
  polynomial, not the resolvent.
