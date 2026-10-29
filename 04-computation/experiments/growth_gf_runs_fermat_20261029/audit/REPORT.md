# Audit record: THM-4845, THM-4846, THM-4847, HYP-9445 (session mac-mini-2026-10-29-growth)

*Independent adversarial audit by a subagent of mac-mini-2026-10-29-growth, 2026-10-29. The producer saved this report
verbatim from the auditor's final message. The auditor's scripts and outputs are in this directory.*

**Auditor:** independent adversarial audit, 2026-10-29.

**Scope.** Worktree `math-wt-pyth-cuboids`. The files were audited as they stand at HEAD `03c05b4762`: checkpoint `a37e12900d` plus the producer's later commits `447812773b`, `e02f4f5b13`, `5047208384`, `92514b16e2` and `03c05b4762`. The additions made after the audit started were also checked: THM-4846 (F) and (D′), the palette family and A002785 in THM-4845 (B), and the order-9 data in HYP-9445.

**Method.** Every number was recomputed by an algorithm different from the producer's:
- direct Burnside over all of S_n;
- the explicit divisor realisation of the Gaussian field, plus Gauss–Hermite quadrature;
- core kernels by brute force over labelled permutations;
- run observers by a subset dynamic programme that builds paths by adding start vertices, never by cutting permutations; it was cross-checked against brute-force path-cover enumeration for all k and n ≤ 5;
- Smith forms, ranks, determinants and nullspaces with FLINT (python-flint 0.9) instead of PARI;
- star kernels as nullspaces of column-restricted matrices;
- strong counts mod p by power-series inversion of S = 1 − 1/T instead of Wright's convolution.

Class lists come from `gentourng`; the class counts 2, 4, 12, 56, 456, 6880 equal A000568. Everything is in `04-computation/experiments/growth_gf_runs_fermat_20261029/audit/`.

## Summary verdict

The mathematics is sound and almost every number reproduces exactly. Per-class agreement with the producer covers the Smith forms of M and the unglued N, spur₂/spur₃ and the boundary-killing flag:
- (n, k) = (5,2), (5,3), (6,2), (6,3);
- all 456 classes at (7,3);
- HC, Smith form of M₂ and spur₂ for all 456 classes at order 7 and all 6880 at order 8.

No result needs withdrawal. What must be fixed:

1. **FAIL — THM-4847 (C3) and its title.**
   - The claim is "at a Fermat prime p ≥ 5, s₂(p) ≡ 1", and the title says "a Fermat prime has s₃(p) = −1 and s(p) = +1".
   - At p = 5: s(5) = 544 ≡ −1 (mod 5), because (2/5) = −1.
   - At p = 3: s(3) = 2 ≡ −1, and s₃(3) ≡ 0.
   - Corrected: s₂(p) ≡ (2/p), which is +1 exactly for the Fermat primes p ≥ 17; s₃(p) ≡ −1 holds for every Fermat prime p ≥ 5.
2. **FAIL (definition) — THM-4846 (E) and results note §2.**
   - The "glue loss" counts are counts of classes with spur₂ + spur₃ &gt; 0, which only sees the first p-adic level.
   - THM-4846 §0 defines spurious factors as "invariant factors that N has and M lacks". By that definition more classes lose: the table below shows the gap.
   - "The first S₂ losses" at order 6 is false. At order 5, class `1101110111` (scores 12223, H = 13, two Hamiltonian cycles) has M₂ with Smith form 1¹⁶ 2² 48 0, while N has 1¹⁶ 2² 96 0.
   - A nonsingular example at order 6, k = 2: `111011110111111` has M = 1²⁴ 2² 4² 48 12288 and N = 1²⁴ 2² 4² 96 6144, with the same determinant.

   | k | order 5 | order 6 | order 7 | order 8 |
   |---|---|---|---|---|
   | 2: Smith form of N ≠ M | 1 | 24 | 273 | 5255 |
   | 2: spur &gt; 0 (producer's count) | 0 | 18 | 192 | 4224 |
   | 3: Smith form of N ≠ M | 4 | 31 | 373 | — |
   | 3: spur &gt; 0 (producer's count) | 4 | 30 | 360 | — |

3. **Citation and novelty corrections**, itemised below:
   - Hubbard–Stratonovich for the Gaussian trick;
   - Oberschelp, Harary–Palmer and Wright for the classical correction expansions;
   - Garey for the Fibonacci count;
   - Lidl–Niederreiter for (B), which is classical, not new;
   - OEIS A002785 already carries the twisted-Burnside formula.
4. **Minor gaps and wording**, listed per claim.

## THM-4845

- **(A) Gaussian representation: PASS.**
  - The covariance identity E[X_pX_q] = ln2·Σ_{d|gcd} φ(d) = ln2·gcd holds.
  - With the explicit realisation, a_n = Σ_m 2^{−|m|/2+½Σ_d φ(d)M_d²}/z_m equals Davis for n ≤ 30.
  - Gauss–Hermite quadrature gives E[Z_n(W)] = 1, 1, 1, 2, 4, 12, 56 for n = 0..6, to 9 decimals.
  - e(m) = (mᵀGm − |m|)/2 was re-derived from pair orbits, and the THM-516 form holds. The core-kernel generating function (3) has coefficient B[m,t] of THM-514.
  - **CITE-NEEDED:** this is the Hubbard–Stratonovich (Gaussian linearisation) identity applied to THM-516's quadratic form; it should be named.
  - **OVERCLAIM (mild):** "answering OPEN-Q-105(a)/HYP-2538".
    - The representation does diagonalise the coupling, but only inside an expectation that exists coefficientwise: the averaged series has radius 0, and no closed evaluation is given.
    - Recommended wording: "reformulates / partially addresses". OPEN-Q-105(a) should stay open for a closed form.
    - Likewise, "the Gaussian-conditional Newton recursion is the recursion that gives the generating function" is framing: the recursion runs on random polynomials, and evaluating it is the Burnside sum.
  - **Title wording:** "the core kernel B[m,t] is E[…]" should say "Σ B[m,t]x^m v^t = E[…]".
- **(B) Field family: PASS.**
  - Direct Burnside over all σ ∈ S_n equals A000568, A000088, A000273 and A001174 for n ≤ 8.
  - The divisor-realised field reproduces all four families for n ≤ 8, and tournaments to n = 12.
  - The palette exponent q^u t^v holds algebraically, with u = (mᵀGm − Σ_odd m_j)/2 − Σ_even m_j.
  - A002785: a direct twisted Burnside count, using the orbits of τ(a,b) = (σb, σa), gives 1, 1, 2, 2, 8, 12, 88, 176 for n = 1..8. The producer's formula gives those values, and its values for n ≤ 16 equal OEIS A002785.
  - **CITE:** OEIS A002785 notes "also self-converse tournaments" (McKay) and already has an equivalent explicit formula (Pab Ter, 2005). The Gaussian version is a rewrite of it.
- **(C) Cameron–Martin is the inherited PDE: PASS.**
  - The exponent identity e(m) = e(m − e_p) + (p−1)/2 + Σ_j gcd(p,j)(m−e_p)_j holds for every base q: 2364 checks to degree 25.
  - It is exactly session98's PROJECTIVE_RESPONSE (1) at q = 2.
- **(D) Meromorphic normalized series: PASS.** Priority belongs to session98's `ALL_POLES_CORE_EXPANSION`, which is now credited.
  - Brute-force kernels for m ≤ 9: K[3,1] = 4, K[5,1] = 96, K[6,2] = 1280, K[7,1] = 5760, K[8,2] = 43008, K[9,1] = 645120, K[9,3] = 9175040.
  - P₁ to P₄ are identical to the producer's, and the leading-coefficient formula holds for k ≤ 7.
  - The exact expansion holds for n ≤ 60.
  - At y = 0.3 the core sum equals the power series: 1.447300039618993.
  - (1 − y/4)⁴F(y) tends to 32 (31.9904 at y = 4(1 − 10⁻⁴)).
  - Truncations m ≤ 40, 60, 90 of (5) agree to 20 digits at y = 20, −50 + 7i, 1000 + i and 3 + 2i.
  - The bound max e(μ)/C(m,2) = 1/3 over cores with m ≤ 40 is sharp at the 3-cycle.
  - The convergence and pole-order proofs are correct.
  - **FAIL (notation):** the title and (A)/(D) call μ "permutations" while weighting by |μ|!/z_μ. The sum must run over cycle *types*; otherwise each type is overcounted by its class size.
  - **GAP (fixable):** the remainder claim "4^{−(K+1)n}P_{K+1}(n)(1 + o(1))" has no proof.
    - One-line fix: F − Σ_{k≤K+1} H_k is holomorphic on |y| &lt; 4^{K+2}, so its coefficients are O(R^{−n}) for every R &lt; 4^{K+2}.
    - Measured remainder over next term at n = 12: 1.0061, 1.0554, 1.1326, 1.0071 for K = 0..3. At n = 20: 1.00000067, 1.0000139, 1.00031, 1.0057.
  - **Graph statement: PASS.**
    - The expansion is exact for n ≤ 30, and the brute-force fixed-point-free kernel agrees for m ≤ 8.
    - Q₁ = 4C(n,2) and Q₂ = 32C(n,3) + 192C(n,4) = (8/3)(n)₃(3n − 7).
    - The revised convergence bound e(μ) ≤ (C(m,2) + m)/2 is valid.
  - **CITE-NEEDED:** these correction terms are classical:
    - graphs: Pólya; Oberschelp, Math. Ann. 174 (1967); Harary–Palmer, *Graphical Enumeration*, §9.1;
    - tournaments: Moon 1968; Wright, Glasgow Math. J. 11 (1970) 97–101, as session98 cites.

    Only the global meromorphy and the exact pole orders are new, and they are session98's.
- **(E) Not P-recursive: PASS.** Both arguments are correct and standard.
- **(F) Borel–Cantelli remark: GAP.**
  - W_p^{1/p} → 1 almost surely is right, since P(|X_p| &gt; 2√(p·ln2·ln p)) ≤ 2/p².
  - "Radius exactly 1" does not follow: the exponential of a radius-1 series can have a larger radius.
  - Fix: the coefficients are ≥ 0, so Pringsheim's theorem applies once Σ W_p/p = ∞ almost surely. Along primes, X_p = √ln2 (ξ₁ + √(p−1) ξ_p) ≥ ln p infinitely often (conditional second Borel–Cantelli), giving infinitely many terms ≥ 2^{−1/2}.

## THM-4846

- **(A) One-arc observer: PASS.** The vacancy-walk form has priority in `four_vertex_glue` §1, now credited.
  - For every class at n = 4, 5, 6:
    - M + Mᵀ is the adjacency matrix of ST_n (24/36, 120/240 and 720/1800 vertices/edges);
    - each vertex pair labels exactly (n−1)! edges;
    - M is S_{n−1}-equivariant;
    - the trivial block is Aᵀ and the sign block is (−1)^{|x−y|−1}A[x][y].
  - ST₄ is isomorphic to G(12,5), the Nauru graph. The Nauru graph is classically known to be this Cayley graph (Wikipedia); cite it.
  - Determinants match the stated lists exactly:

    | n | singular classes | nonzero determinants |
    |---|---|---|
    | 4 | 3 | 1 (strong) |
    | 5 | 7 | 1, 1, 1, 2⁸, 3⁶ |
    | 6 | 33 | 1 (×8: 7 strong, 1 not), 2⁸, 2⁸3¹²5¹⁰7¹⁰, 2²⁴3²⁸5¹⁰, 2³⁶11¹⁰, 2⁴⁰13¹⁰ (×4), 2⁴² (×2), 3³⁰ (×2), 7¹⁰ (×2), 23¹⁰ |

  - The multiplicity n·f^λ in the Schur factorization is correct.
  - **Vacuous:** "every exponent is a combination of the S₅-dimensions 1, 4, 5, 6" holds for every integer, since 1 is in the list.
- **(B) Two arcs: PASS.** The proof is correct, and per A[β,α] matches the direct plus sector for every class at n ≤ 7.
- **(C) T5 stars: PASS.**
  - rank Pc = 6; the five stars lie in the kernel and span it; the relation is Σ(−1)^w v_w = 0.
  - T5 − w has scores 1122, a unique Hamiltonian cycle and diagonals w+1→w+3 and w+2→w+4; c_w is ± on the diagonal ends.
  - M₃ kills every lifted tetrahedron boundary; rank 56; Smith form 1⁵² 3² 6² 0⁴, which is (R7).
  - Uniqueness: no other boundary killer exists for (n,k) = (5–7, 3), (4–8, 2), (6,4), (7,4) or (7,5).
  - **OVERCLAIM (interpretation; §3 and results note):** "each spurious Z/6 is the shadow of one kernel direction / one per independent star" depends on a basis choice. The rotation of Z/5 acts on the span of the stars through Φ₅, which is irreducible over F₂ and F₃. Only the 4-dimensional correspondence is canonical.
- **(D) Ghost lemma: PASS.**
  - The proof was rechecked: Z₍₂₎[S₃] splits, the image is Λ ⊕ M₂, and the block form and rank formula follow.
  - spur₂ equals the connecting rank for every class at k = 2 through order 8 (6880 at order 8), every class at k = 3 through order 7, and the 239 order-8 classes with at most 3 Hamiltonian cycles.
  - At T5 the connecting map has rank 4, onto coker P̄.
  - The lemma is first-level only; see item 2 of the summary.
- **(D′) Mod-3 ghost lemma, post-launch: PASS.**
  - The derivation was redone by hand: the 3-local image is W′ ≡ [[u,∗],[0,v]], with pieces I and II as stated, and the reduction when S̄ is invertible holds.
  - The formula equals the direct mod-3 rank difference on every class for n = 4..7 (0 mismatches). Classes with spur₃ &gt; 0: 0, 3, 21, 271.
- **(E) Census: PASS, except the glue-loss item above.**
  - Orders 4 and 5: as stated.
  - Star law fails for 12/168/2808 classes at orders 6/7/8, with matching breakdowns.
  - **Harmless bug:** `compound_kernel.py` imposes C(u) = 0 for every u at a sink, contradicting its own docstring. The true star kernel at a sink is all of Q^{n−1}. This mis-reports per-vertex dimensions but changes no failure count.
  - Torsion against Hamiltonian cycles (HC), all as claimed:
    - k = 2: classes with HC ≥ 2 are 3, 27, 332 and 5953 at orders 5–8, all with torsion;
    - unique-HC torsion-free counts are 1, 1, 3, 8, 16, 36 for orders 3–8;
    - k = 3: the same through order 7;
    - at order 8, k = 3: all 60 classes with HC = 2 and all 124 with HC = 3 have torsion, and the 55 unique-HC classes split 36/19 with the same 19 classes as at k = 2;
    - k = 4, order 6: 27/27;
    - the one-arc exceptions are as stated.
  - **CITE fix:** F_{2n−6} is Garey's count (J. Combin. Theory Ser. B 13 (1972) 266–269). Douglas (1970) characterised these tournaments. "Douglas's Fibonacci numbers" is a misattribution.
- **(F) Sectors as compounds, post-launch: PASS.**
  - Plus sector = [x_V] per P[α,β], sign sector = [x_V] det P[α,β].
  - The Jacobi convention (−1)^{Σα+Σβ}, with rows V−β and columns V−α, is correct for every k.
  - Verified for k = 2, 3, 4 on every class with n ≤ 6: 13944 + 23664 + 12600 entries, with P built from path tables.
  - **Wording:** non-principal minors expand over path-and-cycle systems, not over permutations.
- **Data-file note:** in `fast_m2_n*.out`, the field "unglued=" is the concatenation of the plus and minus Smith forms, not the Smith form of N. spur₂ is unaffected.

## THM-4847

- **(A1) Coprimality: PASS** (classical; restates S58 §1.2).
- **(A2) Periodicity: PASS.** The certificate type should read "a common prime factor".
- **(A3) Algebra: PASS on the facts.**
  - **OVERCLAIM:** "F_n is not the value of a nontrivially factored polynomial at X = 2" is false as worded, since N = ab is the value of (X + a − 2)(X + b − 2).
  - It should say that the family X^{2^n} + 1 has no algebraic or Aurifeuillean factorization.
- **(A4) Π⁰₁ classification: PASS.**
- **(B) Function-field analogue: PASS as mathematics, but classical.**
  - Source: Lidl–Niederreiter Thm 2.47; that X⁴ + 1 is reducible modulo every prime is a textbook fact. "NEW HERE" should be dropped, and the comparison typed ANALOGY.
  - Checked: X^{2^n} + 1 over F_q (q prime ≤ 23, n ≤ 5) is irreducible only for n = 1 with q ∈ {3, 7, 11, 19, 23}, and every factor has degree ord_{2^{n+1}}(q).
- **(C1) Tournament counts are blind: PASS** for n = 2..14.
- **(C2) Three-state counts are Pépin: PASS.** The (⇒) proof is correct; checked for n = 1..14.
- **(C3) Strong counts at a prime: FAIL**, as in item 1.
  - Direct values (Wright with exact binomials mod p², no vanishing assumed): s₂(p) ≡ −1, −1, +1, +1 and s₃(p) ≡ 0, −1, −1, −1 for p = 3, 5, 17, 257.
  - The producer's "FINITE-EXACT semicomplete Pépin form at 5, 17, 257" only evaluates 3^{C(p,2)}, so it is circular.
- **F₅ numbers: PASS** (by series inversion):
  - p = 641: r = 257, q even, D(257) ≡ 254, s₃(N) ≡ 150;
  - p = 6700417: r = 641, q odd, D(641) ≡ 4217483, s₃(N) ≡ 5459519;
  - the quasi-periodicity laws were re-verified for b = 2, 3 and primes 5..37.
- **(D) Symmetry versus counting: PASS.**
- **Citations: PASS.** F₅–F₃₂ composite and F₃₃ unknown (Wikipedia); the open problems; Boklan–Conway arXiv:1605.01371, whose 4/2³² = 2⁻³⁰ bound is for n ≥ 33.
- **"Odd-order units" row: PASS.**

## HYP-9445

- **(1): PASS** for all counts.
  - The sign sector S₂ has torsion for all 3, 27, 332 and 5953 classes with HC ≥ 2.
  - For unique-HC classes, S₂ is torsion-free in 32 of 32 through order 7 and 54 of 55 at order 8. The exception `1010111111111111111101111101` has S₂ = 1²¹ 2⁴ 0³.
  - Odd-only torsion occurs in 8 classes at order 7 and 225 at order 8.
  - Not re-checked: the order-9/10 random samples, and the full order-9 claim for HC ≥ 2 (177 989 classes).
- **(2): PASS.**
  - Nonsingular unique-HC counts are 1, 1, 3, 7, 16, 36.
  - Every nonsingular M₂ of a unique-HC class is unimodular, and every torsion class is singular.
  - Order 9, recomputed from the producer's 144 recorded classes: all are distinct and each has one HC, so by Garey's 144 the list is complete. They split 81 unimodular and 63 torsion, all 63 singular. The pre-registered prediction of 80 failed, as the producer reports.
  - "3·12³ and 2·6²" are the k = 3 Smith forms; at k = 2 the same classes have 4², 3³·6 and 2².
  - **FAIL (minor):** the growth constant 2.24698 is the largest root of y³ − 2y² − y + 1, equal to 1 + 2cos(2π/7). It is *not* a root of y³ − y² − 2y + 1, whose largest root is 1.80194.
  - The (n+1)·2^{n−6} fit and the A033303 fit are NUMEROLOGY (exact fits).
- **(3): PASS.** At order 6 two strong classes have no star vector at all, with kernel dimensions 2 and 3.
- **(8a): partly wrong.**
  - A single failing prime factor already decides the test.
  - For F₇ and F₈ the smallest known factors give r ≈ 4.1·10¹⁶ and 2.3·10¹⁴, far beyond "about 10⁷" and beyond series inversion.

**Cross-references.** Every path cited in the five files exists. The bare "THM-1430" and "THM-589" references are ambiguous because the repository has duplicate numbers.

## Upgrades from this audit (FINITE-EXACT)

1. **F₆ fails both strong rotation tests at both factors.**
   - At 274177 (r = 65537): D ≡ 215736 and s₃ ≡ 123651.
   - At 67280421310721 (r = 274177): D ≡ 32409983229131 and s₃ ≡ 35318451187950.
2. **F₉ fails** at 2424833 (r = 786433): D ≡ 522107, s₃ ≡ 623467.
3. **F₁₁ fails** at 319489 (r = 90113: D ≡ 61029, s₃ ≡ 273957) and at 974849 (r = 761857: D ≡ 350016, s₃ ≡ 751119).

So F₅, F₆, F₉ and F₁₁ are neither Camion pseudoprimes nor pass the semicomplete Pépin test. The method was cross-checked against Wright's recurrence modulo a 47-bit prime. Also contributed: the proofs closing the (D) remainder gap and the (F) radius gap, (D′) through order 7, and the order-9 unique-HC split.

**Reproduction.** In the `audit/` directory:
- `audit_gf.py`, `mero_check.py`, `audit_selfconverse.py`;
- `audit_census.py n k [--slice a:b]`, `compare_census.py`, `merge_n7k3.py`, `audit_fast.py`;
- `audit_star.py`, `audit_onearc.py`, `audit_t5.py`, `audit_ghost3.py --n7`, `audit_jacobi.py`, `audit_s2.py`, `n9_unique_check.py`;
- `audit_fermat.py --f6`, `fermat_extend.py`, `fermat_primes_direct.py`, `check_series_bigp.py`.

Outputs are in `audit/out/`. No existing file was edited and nothing was committed.

Sources:
- [Boklan–Conway, arXiv:1605.01371](https://arxiv.org/abs/1605.01371v2)
- [Wikipedia: Fermat number](https://en.wikipedia.org/wiki/Fermat_number)
- [Wikipedia: Nauru graph](https://en.wikipedia.org/wiki/Nauru_graph)
- [Garey 1972 cited in Springer s00373-006-0641-8](https://link.springer.com/article/10.1007/s00373-006-0641-8)
- [Demaria–Kiihl, enumeration of the Douglas tournament](https://repositorio.unicamp.br/acervo/detalhe/1385918)
- [UNF SOARS poster, 2020](https://digitalcommons.unf.edu/soars/2020/spring_2020/39)
- [OEIS A002785](https://oeis.org/A002785)
