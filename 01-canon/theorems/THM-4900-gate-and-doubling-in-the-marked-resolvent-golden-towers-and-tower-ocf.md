---
id: THM-4900
title: "The strong gate and the skew-Sylvester doubling in the marked resolvent of a tournament: the gate acts linearly on (N, P) = (u^T adj(zI-M) u, det(zI-M)) with eigenvalues (sqrt(z^2+4) +- 1)^2 (the padding E is the same with sqrt(z^2+4) -> z); nested gates have closed-form characteristic polynomials, an explicit spectrum over one vertex, and a universal resolvent limit (2 + sqrt(z^2+4))/z, the metallic mean of index 4/z; z = 1 is the self-similar point of the doubling, where phi(T) = F_T(1) = u^T (I + M^T M)^(-1) u is an exact doubling invariant, series composition is relativistic addition with the one-vertex tournament absorbing, and the gate is the golden map whose orbit from one vertex runs through Fibonacci convergents to phi^3; and a tower form of the odd-cycle formula, H(Gamma^r T) = sum_j sw_j(T) H(Gamma^r TT_(n-2j)) with sw_1 = 2 c3"
status: >
  PROVED: statements 1-6 (exact linear algebra, Mobius dynamics, and the cover-polynomial operator).
  FINITE-EXACT companions: every identity checked exactly in 04-computation/experiments/microcosm_gate_resolvent_20261106/
  (sympy on random cores of orders 1-4, numerical spectra to k = 10, tower values to total order 1201).
  KNOWN ingredients, credited:
  - the marked pair (p, F), the padding E as a Mobius boundary transfer, the doubling D as the substitution
    w^2 = (z^2+1)/2, and the Cayley linearization of E (session35_dimensions_dualities_20261023/recursions/RESULT.md);
  - the gate Gamma, its recognizability, HC(Gamma T) = H(T), the cover-polynomial operator
    P_Gamma(q) = q[P(q+1) + qP(q) - P(q-1)] and H(Gamma^r T) = (L^r P_T)(1) (Session210, graded-spectrum-gates/PROOF.md);
  - the series (order-join) law, i.e. P +- zN multiply (THM-3322, in the coordinates det(I - zC), C = 2A - J);
  - the size-weighted odd-cycle form of P_T for odd n (THM-078);
  - the skew-Sylvester double D (THM-447);
  - D(Gamma(T)) has no gate, i.e. the tournament-level non-commutation in statement 5 (Session211
    recursion-bridge/PROOF.md section 5; CONCURRENT, pushed earlier the same day; re-proved here in one line).
  New here:
  - the gate's linear transfer and its eigen-coordinates;
  - the closed-form gate towers, their universal (metallic) limit and the explicit spectrum;
  - the self-similar point z = 1, the invariant phi and its golden gate map;
  - the doubling and padding towers;
  - the scalar shadow of the gate/doubling interaction: commutation at z = 1, failure of the BS relation there, and
    P_(D Gamma T) = P_(D E T);
  - the tower odd-cycle formula and sw_1 = 2 c3 for every order.
  Self-audited only (PRE-AUDIT).
source: mac-mini-2026-11-06-microcosm (owner directive 2026-11-06: micro/macrocosm; the duality of x+2 and x*2 recursion on tournaments; an n-vertex tournament as a strong (n+2)-tournament with a recognizable gate)
depends_on: [THM-447, THM-291, THM-078, THM-3322, THM-1370, session35 recursions RESULT.md, session210 graded-spectrum-gates PROOF.md, session211 recursion-bridge PROOF.md]
verification:
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/gate_resolvent_exact.py (+ .out): statements 1, 3, 4 exactly (sympy), Fibonacci orbit
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_spectrum.py (+ .out): statement 2 (closed form, explicit spectrum)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_limits.py (+ .out): statements 2, 3, 5 numerically
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_ocf.py (+ .out): statement 6 exactly
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/run_all.sh: reproduces every output
---

# THM-4900 — the gate and the doubling in the marked resolvent

## Setting

- **Tournament data.** For a tournament T on n vertices, let M = A − Aᵀ be its dominance matrix (skew-symmetric,
  ±1 off the diagonal) and u the all-ones vector. The *marked pair* is

      P_T(z) = det(zI − M),    N_T(z) = uᵀ adj(zI − M) u,    F_T = N_T/P_T = uᵀ(zI − M)⁻¹u.

  This is Session35's (p, F); this file uses the letter P for p.
- **Operations.**
  - **Series composition** A ⇒ B: every vertex of A beats every vertex of B. This is THM-3322's "order join".
  - **Padding** E(T): s → T → t with s → t. This is THM-291's canonical Mode-B section.
  - **Gate** Γ(T): the same with t → s. This is Session210's a → b → T → a: the arc t → s is the recognizable gate,
    d⁺(t) = d⁻(s) = 1.
  - **Skew-Sylvester double** D(T), with dominance [[M, M + I], [M − I, −M]] (THM-447).
- **The x+2 / x·2 grid.** On orders, E and Γ act as n ↦ n + 2 and D as n ↦ 2n. These satisfy the Baumslag–Solitar
  relation d∘a = a²∘d (Session35 §1).

Write s = √(z² + 4).

## Statement

**1. The gate is linear on the marked pair.**

    (N, P)_Γ(T) = G(z)·(N, P)_T,    G(z) = [[z²+9, 2z], [2z, z²+1]].

The padding E has G_E(z) = [[z²+1, 2z], [2z, z²+1]] (Session35). Consequences:
- P_Γ(T) = P_T·(z² + 1 + 2zF_T) = P_E(T): flipping the apex arc does not change the characteristic polynomial.
- F_Γ(T) = ((z²+9)F_T + 2z)/(z² + 1 + 2zF_T).
- **Eigen-coordinates.** G(z) is symmetric, with eigenvalues (s ± 1)² and eigen-coordinates

      Ψ± = zN + (−2 ± s)P,    Ψ±(ΓT) = (s ± 1)²·Ψ±(T).

  Setting s ↦ z and dropping the −2 gives E's coordinates N ± P, multiplied by (z ± 1)². That is Session35's Cayley
  linearization C = (1+F)/(1−F) ↦ ((z+1)/(z−1))²·C.

**2. Nested gates.** For every core T and k ≥ 0:

    P_(Γ^k T)(z) = [(s+1)^(2k)·Ψ+(T) − (s−1)^(2k)·Ψ−(T)] / (2s).

- **Universal limit.** For every real z > 0 and every core T,

      F_(Γ^k T)(z) → (2 + s)/z,

  with error O(((s−1)/(s+1))^(2k)). At z = 4/m the limit is the metallic mean (m + √(m²+4))/2:
  - golden φ at z = 4;
  - silver 1 + √2 at z = 2;
  - 2 + √5 = φ³ at z = 1.

  The macro limit forgets the core.
- **Over one vertex.** For k ≥ 1,

      P_(Γ^k K1)(z) = z(z²+3)·[(s+1)^(2k−1) + (s−1)^(2k−1)] / (2s).

  The spectrum of Γ^k(K1), which has 2k+1 vertices, is {0, ±i√3} together with

      ±i·√(4 + cot²((2j+1)π / (2(2k−1)))),    j = 0, …, k−2.

  Every non-core eigenvalue has modulus above 2. For example, the 5-vertex tower has spectrum {0, ±i√3, ±i√7}.

**3. The self-similar point of the doubling.** The substitution w² = (z²+1)/2 fixes z = ±1. So

    F_(D T)(±1) = F_T(±1),    P_(D T)(1) = 2ⁿ P_T(1)².

Put φ(T) := F_T(1) = uᵀ(I − M)⁻¹u = uᵀ(I + MᵀM)⁻¹u.
- 0 < φ(T) ≤ n, with equality iff T is regular.
- φ(D T) = φ(T) exactly.
- **Towers, for real z > 0.**
  - D^k: F_(D^k T)(z) = (z/w_k)·F_T(w_k), where w_k² − 1 = (z² − 1)/2^k. Hence F_(D^k T)(z) → z·φ(T): the doubling
    tower remembers exactly one number of its core.
  - E^k: F_(E^k T)(z) → 1, for z ≠ 1.
  - Γ^k: the limit (2 + s)/z of 2. The gate tower forgets its core entirely.

**4. The algebra at z = 1.**
- **Series composition is relativistic addition.**

      F_(A⇒B) = (F_A + F_B)/(1 + F_A F_B)    at every z,

  equivalently, the Cayley transforms multiply. In particular φ(A ⇒ B) = (φ_A + φ_B)/(1 + φ_A φ_B).
  - Since φ(K1) = 1 is absorbing, φ(T) = 1 whenever T has a one-vertex strong component: transitive tournaments, every
    E(T), every T with a source.
  - φ is the relativistic sum of the φ's of the strong components. For example φ(C3 ⇒ C3) = 3/5.
- **The gate is the golden map.**

      φ(Γ T) = (5φ(T) + 1)/(φ(T) + 1).

  This is hyperbolic, with fixed points φ³ = 2 + √5 (attracting) and −φ⁻³, and multiplier φ⁻⁴.
  - In the coordinate y = (φ − 1)/2 (the normalization of q = (H − 1)/2 for path counts), it is y ↦ (2y+1)/(y+1).
    That is the square of the continued-fraction map y ↦ 1 + 1/y.
  - From the single vertex (y = 0), the orbit consists of the even Fibonacci convergents:

        y(Γ^k K1) = F_(2k)/F_(2k−1),    φ(Γ^k K1) = F_(2k+2)/F_(2k−1)  (1, 3, 4, 21/5, 55/13, 144/34, …).

**5. x+2 against x·2.**
- At z = 1, D acts trivially. So φ(DΓT) = φ(ΓDT) for every T, but φ(DΓT) ≠ φ(Γ²DT) for every T. The integer relation
  d∘a = a²∘d fails here and commutation holds instead.
- At the tournament level, D(Γ(T)) never has a gate. Hence D(Γ(T)) ≇ Γ(X) for every X, and in particular
  DΓ ≠ Γ²D (compare Session35's DE ≠ E²D). KNOWN: Session211 `recursion-bridge` PROOF §5, pushed earlier the same
  day.
- P_(DΓT) = P_(DET).
- On the marked pair, the +2 operations are linear (G_E, G_Γ), while ×2 is quadratic with a substitution:
  P_D = 2ⁿP(w)² and N_D = 2ⁿ(z/w)N(w)P(w).

**6. The tower odd-cycle formula.** Write P_T(q) = Σ_j sw_j(T)·q^(n−2j) for the ordered path-cover polynomial
(Session210 Lemma 3; THM-078 identifies sw_j for odd n as the size-weighted counts of independent odd-cycle sets).
- sw_0 = 1, and sw_1(T) = 2c3(T) for every n.
- For every r ≥ 0,

      H(Γ^r T) = Σ_j sw_j(T)·H(Γ^r TT_(n−2j)),

  where TT_m is the transitive tournament on m vertices, TT_0 is empty, and P_(TT_m) = q^m. The core enters only
  through its odd-cycle statistics; the tower enters only through universal weights U_(d,r) := H(Γ^r TT_d).
- At r = 0 every weight is 1 and this is the odd-cycle formula H(T) = Σ_j sw_j(T).

## Proofs

**1.** Order the vertices s, core, t. Then

    zI − M_Γ = [[z, −uᵀ, 1], [u, B, −u], [−1, uᵀ, z]],    B = zI − M.

- **Determinant.** The Schur complement of B is [[z + F, 1 − F], [−1 − F, z + F]], with determinant
  z² + 1 + 2zF. This gives P_Γ, and the same computation with the apex sign reversed gives P_E.
- **Response.** Solve (zI − M_Γ)(x, v, y) = (1, u, 1). The core rows give v = (1 − x + y)B⁻¹u, so uᵀv = (1 − b)F with
  a = x + y and b = x − y. The two boundary rows give za − b = 2 and zb + a = 2(1 − b)F. Hence

      b = 2(zF − 1)/(z² + 1 + 2zF),    F_Γ = a + (1 − b)F = ((z²+9)F + 2z)/(z² + 1 + 2zF).

  Multiplying by P_Γ gives N_Γ = (z²+9)N + 2zP and P_Γ = 2zN + (z²+1)P.
- **Eigenvectors.** G(z)·(z, s−2)ᵀ = (z³ + 5z + 2zs, s(z²+1) − 2)ᵀ = (z² + 5 + 2s)·(z, s−2)ᵀ, and z² + 5 + 2s = (s+1)².
  The other eigenvector is (z, −s−2), with eigenvalue (s−1)². Since G is symmetric, these are also left eigenvectors.

**2.**
- **Closed form.** By 1, Ψ±(Γ^k T) = (s±1)^(2k)·Ψ±(T), and P = (Ψ+ − Ψ−)/(2s).
- **Limit.** Let R = Ψ+/Ψ− = (zF − 2 + s)/(zF − 2 − s). Then R(Γ^k T) = ((s+1)/(s−1))^(2k)·R(T).
  - For real z > 0 we have s > 2, so |(s+1)/(s−1)| > 1.
  - R(T) ≠ 0, because F_T(z) = z·uᵀ(z²I + MᵀM)⁻¹u > 0 while the root of Ψ+ is (2 − s)/z < 0.
  - So R → ∞, i.e. zF − 2 − s → 0, at the stated rate.
  - At z = 4/m: s = (2/m)√(m²+4) and (2 + s)/z = (m + √(m²+4))/2.
- **One vertex.** K1 has N = 1 and P = z, so Ψ+ = z(s−1) and Ψ− = −z(s+1). Substituting and using s² − 1 = z² + 3 gives
  the formula.
  - Its bracket vanishes iff ((s+1)/(s−1))^(2k−1) = −1, i.e. s = −i·cot((2j+1)π / (2(2k−1))), so
    z² = −4 − cot²(…).
  - The value s = 0 (j = k−1) is not a root: the bracket divided by 2s equals 2k−1 there.
  - Pairs j ↔ 2k−2−j give the same square. This leaves k−1 values and 2(k−1) eigenvalues; with 0 and ±i√3 that is
    2k+1, all of them.

**3.**
- **Fixed point.** w(1) = 1, so the two Session35 formulas give the claims.
- **The form of φ.** (I − M)⁻¹ = (I + M)(I − M²)⁻¹, and uᵀM(I − M²)⁻¹u = 0 since M(I − M²)⁻¹ is skew. Hence
  φ = uᵀ(I + MᵀM)⁻¹u > 0.
- **Bounds.** (I + MᵀM)⁻¹ ⪯ I gives φ ≤ uᵀu = n, with equality iff MᵀMu = 0, i.e. Mu = 0 (regular).
- **D-tower.** Iterating Session35's F_D(z) = (z/w)F(w) telescopes to (z/w_k)F(w_k), and w_k → 1. F_T is continuous at
  1, since its poles are imaginary.
- **E-tower.** By Session35, C_E = ((z+1)/(z−1))²·C, and |(z+1)/(z−1)| > 1 for real z > 0 with z ≠ 1. So C → ∞ and
  F → 1.

**4.**
- **Series law.** Let A ⇒ B have dominance [[M_A, J], [−Jᵀ, M_B]]. Solving (zI − M)x = u blockwise gives
  x_A = (1 + σ_B)R_A u, where σ_B = uᵀx_B and R_A = (zI − M_A)⁻¹. Likewise x_B = (1 − σ_A)R_B u. So
  σ_A = (1 + σ_B)F_A and σ_B = (1 − σ_A)F_B. Hence

      F = σ_A + σ_B = (F_A + F_B)/(1 + F_A F_B).

  THM-3322 states this law as "P ± zN multiply" in its coordinates. Its f = zN/P equals F_M at (1+z)/z.
- **Absorbing vertex.** F_(K1) = 1/z, so φ(K1) = 1, and 1 ⊕ x = 1.
- **Gate map.** Statement 1 at z = 1 gives (10φ + 2)/(2φ + 2).
- **Conjugation.** With φ = 2y + 1: (5(2y+1) + 1)/(2y + 2) = 2·(2y+1)/(y+1) + 1. So the map is y ↦ (2y+1)/(y+1), which
  is the continued-fraction map applied twice.
- **Fibonacci orbit.** By induction,

      y_(k+1) = (2F_(2k) + F_(2k−1))/(F_(2k) + F_(2k−1)) = F_(2k+2)/F_(2k+1).

- **Dynamics.** The Möbius matrix [[5, 1], [1, 1]] has eigenvalues 3 ± √5 = 2φ^(±2), so the multiplier is φ⁻⁴.

**5.**
- **Commutation at z = 1.** D is the identity there, so φ(DΓT) = g(φ(T)) = φ(ΓDT), where g is the golden map.
- **No BS relation at z = 1.** g(x) = g(g(x)) forces x to be a fixed point of g, i.e. φ³ or −φ⁻³. These are irrational,
  while φ(T) is rational.
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

## Remarks

- **Micro/macrocosm.**
  - **The condensation multiplies.** The marked response of any tournament is the relativistic sum (product of Cayley
    transforms) of its strong components' responses, just as H is the product of their path counts.
  - **The gate nests.** A whole arbitrary tournament sits inside a strong shell. Iterating forgets it:
    - the response tends to the metallic limit, at rate φ⁻⁴ at z = 1;
    - the path count is a universal tower combination of the core's odd-cycle statistics (6). HYP-9490 conjectures that
      its first-order trace is the cyclic-triangle count alone.
  - **Doubling is invisible at the self-similar point z = 1.** It remembers exactly φ.
- **Numerology, typed as such.**
  - φ³ = 1 + 2φ is also THM-224's exceptional eigenvalue: COINCIDENCE, no mechanism.
  - The 7-vertex tower has φ = 21/5, and the 5-vertex tower has H = 15 and spectrum {0, ±i√3, ±i√7}: NUMEROLOGY with
    respect to {7, 21} and {3, 15}.
- **What does not follow.** The marked pair does not determine H, nor isomorphism. Statement 6 uses the cover
  polynomial, not the resolvent.
