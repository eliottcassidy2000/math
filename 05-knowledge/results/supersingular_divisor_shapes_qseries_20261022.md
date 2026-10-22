# Supersingular counts, the divisor shapes p^2, pq, p^2 q, and six Rogers–Ramanujan-type families

**Session:** opus-2026-10-22-ssj (temporary identity: this worktree has no `.machine-id`), 2026-10-22.
**Brief (relayed from the session working on HYP-9371):** merge an arXiv source and a set of classical identities; compute F, U, S at n = p^2, pq, p^2 q; relate them to the {p^2 q r, p, p^3} work and to the supersingular count / F_p-F_{p^2} split; look for genuine connections with Bailey mod 9, Dyson mod 27, Göllnitz–Gordon, Jackson–Slater, Rogers mod 14 and Rogers–Selberg.
**Status:** pre-audit at the checkpoint; the audit record is in section 10. Labels are per item:
* PROVED (elementary; classical in substance, no priority claimed);
* CITED (Deuring, Eichler, Deligne–Rapoport, Katz–Mazur, Delfs–Galbraith, Ogg, Tate, Kubert);
* FINITE-EXACT (scripts in [`04-computation/experiments/supersingular_fus_qseries_20261022/`](../../04-computation/experiments/supersingular_fus_qseries_20261022/README.md));
* NUMEROLOGY / REFUTED for the integer coincidences.
* SOURCE MISMATCH for the arXiv link.

No HYP was opened. Every conjectural formula that came up was either proved here or refuted (section 9).

## 0. The answers in one page

1. **The arXiv link is not a paper on supersingular curves.** arXiv:2610.05474v1 is Bizeul–Klartag–Lehec, *Presenting a proof of the Kannan–Lovász–Simonovits conjecture* (math.FA, 4 Oct 2026). Section 1 records what it proves and what it uses.
   * The supersingular facts in the brief are standard. They were checked computationally against standard formulas, not against that paper, for all primes 5 ≤ p < 3000 by three independent routes.
   * Which paper was meant must be confirmed by the user.
2. **F, U, S are counts of divisors, and at the three shapes they are constants.** Here D(N) = {d : d | N, 1 < d < N}, F = |D|, S = #squarefree elements of D, U = #prime elements of D (DB1).

   | N | D(N) | F | U | S | F − U − S |
   |---|---|---|---|---|---|
   | p^2 | {p} | 1 | 1 | 1 | −1 |
   | pq | {p, q} | 2 | 2 | 2 | −2 |
   | p^2 q | {p, p^2, q, pq} | 4 | 2 | 3 | −1 |

   These values do not depend on the primes (PROVED; trial division over all primes < 400 as a control).
   * The three shapes are exactly the solution sets of F = U + S − 1 (p^2, p^2 q) and F = U + S − 2 (pq). These are the first two levels below the balanced set {p, p^3, p^2 q r}.
   * A new prime-adjunction recursion, δ(Nr) = 2δ(N) + U(N) − [N squarefree] − [N prime], builds the chain p^2 → p^2 q → p^2 q r (defects −1, −1, 0). Raising the exponent instead gives p^2 → p^3 (−1 → 0).
3. **F, U, S cannot see supersingular counts.** No identity between them is possible.
   * F, U, S are prime-blind, while every supersingular count is affine in p with slope ψ/12.
   * The only genuine meeting point is the bad reduction of X_0(N). The divisors of the ℓ-part of N index the components of the fibre at ℓ, and the supersingular points glue them.
   * At pq: g(X_0(pq)) = 2g(X_0(q)) + H(p,q) − 1 = 2g(X_0(p)) + H(q,p) − 1.
   * At p^2 q: g(X_0(p^2 q)) = 2g(X_0(p^2)) + H(q,p^2) − 1.
   * At p^2: three components (divisors 1, p, p^2) meet at the ss(p) supersingular points. The middle component, the unique element of D(p^2), has multiplicity p − 1.
4. **The six identity families are three theta families, and their products are coordinates on genus-zero modular curves.**
   * The three theta families sit at levels 7, 8, 9: Rogers mod 14 is Rogers–Selberg at q^2, Dyson mod 27 is Bailey at q^3, and Jackson–Slater is Göllnitz–Gordon at −q.
   * Ratios of their product sides are Hauptmoduln of X_1(7), X_1(8), X_1(9) (explicit j-maps of degrees 24, 24, 36, verified to O(q^520)).
   * This is the genuine link to supersingular curves:
     * **mod 12:** X_1(N) has exactly (p−1)d_N/12 supersingular points, that is 2(p−1), 2(p−1), 3(p−1), with no mod-12 correction, because Γ_1(N) is torsion-free. Among X_0(7), X_0(8), X_0(9) only X_0(7) keeps a correction, mod 3, from its two elliptic points of order 3.
     * **F_p / F_{p^2}:** on the supersingular locus the p^2-Frobenius is the diamond operator ⟨p⟩. It permutes the identities of each family through the class of p in (Z/N)^×/±1, which is the identities' own index set.
     * **class numbers:** #X^ss(F_p) = (h(−4p) s(floor) + [p ≡ 3 mod 4] h(−p) s(surface))/2. The local numbers s are explicit; for example (1 + (−p/7)) n_p on X_0(7), and on X_1(8) a rule that depends on p mod 16.
   * The moduli 14 and 27, the class-number-one reading of 7, 9, 14, 27, and the shapes of 9, 14, 12 (= p^2, pq, p^2 q) are typed NUMEROLOGY, each with a hostile test.

## 1. The source: arXiv 2610.05474v1 (SOURCE MISMATCH; CITED, unaudited)

**Identification.**
* The arXiv API and the HTML page (sha256 prefix `ebf2bf23e635b236`, fetched 2026-10-22) both give *Presenting a proof of the Kannan–Lovász–Simonovits conjecture*, by Pierre Bizeul, Boaz Klartag and Joseph Lehec, arXiv:2610.05474v1 [math.FA], 04 Oct 2026.
* The full page text contains no "supersingular", "elliptic curve" or "j-invariant"; the only "elliptic" is elliptic regularity.
* Sixteen neighbouring identifiers were checked: ±1, ±2 and digit transpositions of 05474, and 05474 in other months. None is about elliptic curves.
* The most recent arXiv papers on supersingular curves are 2609.11037 (supersingularity testing on twisted Hessian curves), 2609.03249 (memetic search for supersingular curves over F_p) and 2608.13371 (Gross vectors mod 2). They were listed, not read.
* The brief's description fits any introduction to the subject, so the intended source is unknown.

**What the KLS paper proves.**
* Theorem 1.1: there is a universal C such that for every n and every isotropic log-concave probability measure μ on R^n, the Poincaré constant satisfies C_P(μ) ≤ C. Equivalently, ‖Cov μ‖_op ≤ C_P(μ) ≤ C ‖Cov μ‖_op.
* This is the KLS conjecture. By Cheeger, Buser, Ledoux, Gromov–Milman and E. Milman it is equivalent to dimension-free concentration and isoperimetry in convex bodies.
* The paper's history section records the earlier bounds √n (Lee–Vempala), n^{o(1)} (Chen), polylog (Klartag–Lehec), log n (Klartag), √(log n), and C·16^{log* n} (Song–Zhang, arXiv:2610.01447).

**What it uses.**
* The third-moment bound |E⟨X,u⟩X⊗X|^2 ≤ 8|u|^2. This comes from the sharp thin-shell inequality Var|X|^2 ≤ 8n (Chen–Klartag) in the anisotropic form of Bizeul and Letwin.
* An induction through Eldan's stochastic localization, with inverse covariance as the quadratic variation, giving cumulant bounds |κ_m(u,·,…,·)|^2 ≤ K^{m−1}((m−1)!)^2 |u|^2 (Theorem 4.1).
* A suspension construction (Proposition 7.1) that realises derivatives of tilt-averages F_f(z) as mixed cumulants in a larger dimension, giving |∇^d F_f(0)/d!| ≤ R^d ‖f‖_2.
* The Song–Zhang criterion ([34, Theorem 5.1], reproved as Theorem 3.1), which turns this into C_P ≤ C R^2.
* The acknowledgements state: "Most proofs and mathematical ideas in this paper were found by ChatGPT; a notable exception is the idea to use suspension which was suggested by the authors."

**Repository status.** CITED and unaudited. No map to an open repository problem was found or sought.

## 2. The supersingular facts in the brief, checked (CITED + FINITE-EXACT)

Script: `supersingular_census.gp`, run on all 428 primes 5 ≤ p < 3000 with 0 failures.

| Fact (standard) | Check |
|---|---|
| ss_p(j) = ∏ (j − j_E) over supersingular E has degree ⌊p/12⌋ + {0,1,1,2} for p ≡ 1, 5, 7, 11 mod 12, and equals g(X_0(p)) + 1 | degree law on 428 primes; g from PARI `mfdim` |
| every root lies in F_{p^2} | every irreducible factor of ss_p over F_p has degree ≤ 2, and ss_p is squarefree |
| #{ss j ∈ F_p} = (h(−4p) + [p ≡ 3 (4)] h(−p))/2 (Delfs–Galbraith), i.e. h(−4p)/2, h(−p), 2h(−p) for p ≡ 1 (4), 7 (8), 3 (8) | linear factors against PARI `qfbclassno`, which is cross-checked by an independent reduced-form count |
| #{conjugate pairs in F_{p^2} \ F_p} = g(X_0(p)^+) = (2g + 2 − h(−4p) − [p ≡ 3 (4)]h(−p))/4 | quadratic factors on 428 primes |
| all ss j lie in F_p iff p ∈ {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71} (Ogg 1975) | reproduced for p < 3000 |

Three routes compute ss_p.
* **Route A (Hasse invariant):** E_{p−1} ≡ 1 (mod p) as a q-expansion. Hence ss_p(j) = j^δ (j − 1728)^ε f(j), where f is the polynomial part in j of 1/(Δ^m E_4^δ E_6^ε) and p − 1 = 12m + 4δ + 6ε.
* **Route B (Legendre line):** the squarefree part of Res_λ(H_p(λ), 256(λ^2−λ+1)^3 − jλ^2(λ−1)^2). It equals route A on 123 primes p < 700.
* **Route C (point counting):** #{j ∈ F_p : a_p(E_j) = 0} equals the number of linear factors on 301 primes p < 2000.

This inherits and extends [ogg_triangular_chowla_20260927 §2](ogg_triangular_chowla_20260927.md), which checked the list for p < 200 (p < 1000 in its audit).

## 3. F, U, S at p^2, pq, p^2 q and the {p, p^3, p^2 q r} work (PROVED + FINITE-EXACT)

**Definitions recovered first.** This is DB1 of [arithmetic braids](arithmetic_braids_20260917_divisors.md), restated in [predecessor_prime_shapes_20261014b §1](predecessor_prime_shapes_20261014b.md) and [trijection_reset_20261020/arithmetic.md §1](trijection_reset_20261020/arithmetic.md).
* For N ≥ 2: F(N) = #{d | N : 1 < d < N}, S(N) the squarefree ones, U(N) the prime ones.
* With r = ω(N): F = τ(N) − 2, S = 2^r − 1 − [N squarefree], U = r − [N prime].
* U ⊂ S, so F = U + S is an equality of counts, not a disjoint union.
* The "Chowla reading" σ(n) − n − 1 in [ogg_triangular_chowla_20260927 §5](ogg_triangular_chowla_20260927.md) is a different object (a sum of divisors), and it is not used here.

**Values** (table in section 0). They follow from the formulas: τ = 3, 4, 6 and r = 1, 2, 2.
* The check `fus_shapes.py` part 1 covers all 78 squares p^2 and all 6006 ordered pairs of distinct primes below 400 for pq and p^2 q, plus a slow trial-division path on primes < 60.
* Hostile controls: the collision q = p sends pq to p^2 and p^2 q to the balanced p^3 (defect −1 → 0). The literal exponent reading p^q r, at 40 = 2^3·5, gives (6, 2, 3) with defect +1, a different family.

**Relation to the balanced shapes.**
1. **Level sets** (DB2, re-run exhaustively for N ≤ 2·10^6):
   * δ = 0 exactly at p, p^3, p^2 q r;
   * δ = −1 exactly at p^2, p^2 q;
   * δ = −2 exactly at pq;
   * every other negative defect is a squarefree composite with δ = −ω (−3 to −7 occur below 2·10^6).

   So the asked shapes are the complete solution sets of F = U + S − 1 and F = U + S − 2.
2. **Prime-adjunction recursion** (PROVED; checked on 59,997 pairs N ≤ 20,000).
   * For a prime r ∤ N:
     ```text
     δ(N r) = 2 δ(N) + U(N) − [N squarefree] − [N prime].
     ```
   * Proof: F(Nr) = 2F(N) + 2, S(Nr) = 2S(N) + 1 + [N squarefree], U(Nr) = U(N) + 1 + [N prime]. Subtract.
   * Chains:
     * p^2 →(·q) p^2 q →(·r) p^2 q r has defects −1 → −1 → 0. At the last step the balance comes from U(p^2 q) = 2 cancelling 2δ = −2.
     * p →(·q) pq has defects 0 → −2.
     * p^2 →(·p) p^3 has defects −1 → 0 (τ: 3 → 4).
   * So **both balanced composites are reached from the δ = −1 shape p^2**: p^3 by raising the exponent, and p^2 q r by two prime adjunctions through the other δ = −1 shape p^2 q.
   * This is the precise content found for the owner's "p^2 q r is p and p^3 combined". [procgen_tcpc_20261001 §4.3](procgen_tcpc_20261001_tournament_clock_prime_collatz.md) typed that phrase as ANALOGY via the exponent identity (2,1,1) = (3,0,0) + (−1,1,1). The recursion is a counting law; it does not upgrade the analogy to a combination law for F, U, S.
3. **The Klein plane of THM-4733.** The four Klein orbits of D(p^2 q r) under (C, Q) (`fus_shapes.py` part 4) are {p, pqr}, {p^2, qr}, {pq, pr} and {q, r, p^2 q, p^2 r}.
   * Each orbit has a unique representative in the closed quadrant f = a − b − c ≥ 0, g = b − c ≥ 0 of the (f, g) chart of [THM-4733](../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md).
   * The representatives are **p, p^2, pq, p^2 q**: the three asked shapes plus p. Equivalently they form p·Div(pq), the divisors of p^2 q divisible by p.
   * This is exact combinatorics; no arithmetic consequence is claimed.
4. **Atkin–Lehner reading** (CITED, classical: Ogg 1974). This reading applies to the cusp classes of X_0(N), which are indexed by d | N.
   * W_Q (Q ‖ N) acts on the exponent box by reflecting the Q-coordinates: a ↦ 2 − a for W_{p^2}, b ↦ 1 − b for W_q, c ↦ 1 − c for W_r.
   * THM-4733's complement C(d) = N/d is the Fricke involution W_N.
   * The prime exchange Q (q ↔ r) is **not** an Atkin–Lehner operator. The Klein group ⟨C, Q⟩ meets the Atkin–Lehner group ⟨W_{p^2}, W_q, W_r⟩ only in ⟨C⟩.
   * For squarefree N the cusps are exactly the divisors: F + 2 = τ(N) = #cusps (4 at pq). The Atkin–Lehner group has order 2^ω = S + 1 + [N squarefree] and rank ω = U + [N prime].
   * At p^2 and p^2 q the cusp counts are ∑_{d|N} φ(gcd(d, N/d)) = p + 1 and 2(p + 1). The middle divisor p carries φ(p) = p − 1 cusps, which the unweighted count F does not see.

## 4. Why F, U, S do not see supersingular counts, and where divisors and supersingular points meet

**Obstruction** (PROVED, one line).
* F, U, S are bounded functions of the exponent profile: F(p^2) = 1 for every p.
* Every supersingular count at these levels grows linearly in a prime: ss(p) ≥ (p − 1)/12, H(p, q) ≥ (p − 1)(q + 1)/12.
* So no identity of the form "count = expression in F, U, S" holds. The same applies to the F_p / F_{p^2} split, which depends on p through class numbers and residues (sections 2 and 6).

**The genuine meeting point** (CITED: Deligne–Rapoport, Katz–Mazur, Eichler; FINITE-EXACT: `dr_genus.gp`).
* For ℓ^a ‖ N, the special fibre of X_0(N) at ℓ consists of a + 1 copies of X_0(N/ℓ^a)_{F_ℓ}.
* The copies are indexed by the divisors 1, ℓ, …, ℓ^a, with multiplicities φ(ℓ^{min(i, a−i)}). These are the same weights that count cusps.
* All copies pass through the supersingular points of X_0(N/ℓ^a) mod ℓ, whose number is Eichler's class number
  ```text
  H(ℓ, M) = (ℓ−1) ψ(M)/12 + (ν_2(M)/4)(1 − (−4/ℓ)) + (ν_3(M)/3)(1 − (−3/ℓ)).
  ```

| shape | bad fibres | genus identity (tested) |
|---|---|---|
| p^2 | at p: 3 components (1, p, p^2; the middle one of multiplicity p − 1) through the ss(p) supersingular points | (no Deligne–Rapoport form; the fibre is non-reduced) |
| pq | at p: 2 × X_0(q) glued at H(p,q); at q: 2 × X_0(p) glued at H(q,p) | g(X_0(pq)) = 2g(X_0(q)) + H(p,q) − 1 = 2g(X_0(p)) + H(q,p) − 1: on 1332 ordered pairs < 160 |
| p^2 q | at q: 2 × X_0(p^2) glued at H(q, p^2); at p: 3 × X_0(q) through the H(p,q) points | g(X_0(p^2 q)) = 2g(X_0(p^2)) + H(q,p^2) − 1: on 408 pairs with p^2 q ≤ 40,000 |

The supersingular counts at the asked levels in closed form:
```text
H(p, q)   = (p−1)(q+1)/12 + (1/4)(1 − (−4/p))(1 + (−4/q)) + (1/3)(1 − (−3/p))(1 + (−3/q)),
H(q, p^2) = (q−1)p(p+1)/12 + (1/4)(1 − (−4/q)) ν_2(p^2) + (1/3)(1 − (−3/q)) ν_3(p^2),
ν_2(p^2) = 1 + (−4/p) (0 at p = 2),   ν_3(p^2) = 1 + (−3/p) (0 at p = 3).
```
* Symmetric form of the pq identity: 2 ss(q) + H(p, q) = 2 ss(p) + H(q, p). It links supersingular counts in two characteristics. It is checked on all pairs, and it is an algebraic consequence of the closed forms.
* Example rows (`dr_genus.out`):
  * N = 55: H(5,11) = 4, H(11,5) = 6, g = 5;
  * N = 275: H(11,25) = 26, g = 25;
  * N = 637: H(13,49) = 56, g = 57.

**The F_p / F_{p^2} split at these levels.** Every supersingular point of X_0(M) mod p is F_{p^2}-rational, as the case M = 1 generalises. The F_p-rational ones are counted in section 6.3.

## 5. The six identity families: three theta families and three Hauptmoduln

**Statements** are transcribed from McLaughlin–Sills–Zimmer, *Rogers–Ramanujan–Slater type identities*, Electron. J. Combin. DS15 (2nd ed. 2022), §§2.7, 2.8, 2.9, 2.14, 2.27. They are checked exactly to O(q^700) by `qseries_identities.py`: 17 identities, with the Rogers–Ramanujan pair as a positive control and two swapped products as hostile controls (both fail at q^1, as they should).
* Notation: f(a,b) is Ramanujan's theta function and θ_N(a; q) = f(−q^a, −q^{N−a}) = (q^a, q^{N−a}, q^N; q^N)_∞.

| family | DS15 / Slater | product side | level | identity index set |
|---|---|---|---|---|
| Rogers–Selberg | 2.7.1–3 (S.31–33) | θ_7(a; q)/f(−q^2) | 7 | a ∈ {1,2,3} = (Z/7)^×/±1 |
| Rogers mod 14 | 2.14.1–3 (S.59–61) | θ_14(2a; q)/f(−q) = θ_7(a; q^2)/f(−q) | 7 at 2τ | 2·{1,2,3} |
| Göllnitz–Gordon | 2.8.1, 2.8.4 (S.34, S.36) | θ_8(a; q)/ψ(−q) | 8 | a ∈ {1,3} = (Z/8)^×/±1 |
| Jackson–Slater | 2.8.9–10 (S.38, S.39) | f(q^a, q^{8−a})/f(−q^2) = θ_8(a; −q)/f(−q^2) | 8 at τ + 1/2 | {1,3} |
| Bailey mod 9 | 2.9.1–3 (S.40c–42c) | θ_9(a; q)/f(−q^3) | 9 | a ∈ {1,2,4} = (Z/9)^×/±1 |
| Dyson mod 27 | 2.27.1–4 (S.90–93) | θ_27(3b; q)/f(−q) = θ_9(b; q^3)/f(−q) | 9 at 3τ | 3·{1,2,4}, plus the non-unit b = 3 (θ_9(3; q) = (q^3;q^3)_∞) |

**Three theta families** (PROVED, definitional; the sum-side corollaries are verified to O(q^700)).
```text
RS_a-sum(q^2) f(−q^4) = R14_{2a}-sum(q) f(−q)          (a = 1, 2, 3)
B_a-sum(q^3)  f(−q^9) = D_{3a}-sum(q)  f(−q)           (a = 1, 2, 4)
GG_a-sum(−q)  ψ(q)    = JS_a-sum(q)    f(−q^2)         (a = 1, 3)
```
* The moduli 14 and 27 are the substitutions q ↦ q^2, q^3 (τ ↦ 2τ, 3τ) in the numerator theta function. They are not new geometric levels.
* Jackson–Slater is Göllnitz–Gordon at τ + 1/2.
* The concurrent S55 note ([signature_axis_cat_map_rr_mrb_20261022 §4](signature_axis_cat_map_rr_mrb_20261022.md)) found the same Dyson = Bailey(q^3) and Rogers mod 14 = M(2,7)(q^2) relations from the conformal side. This is independent agreement.

**Hauptmoduln from the product sides** (FINITE-EXACT: `level_hauptmoduln.gp`). Write f_g = θ_N(g; q)/(q^N; q^N)_∞ (Klein-form numerators). The common denominators of each family cancel in these ratios.
* For each coordinate below, the rational function j = A(x)/B(x) is found by linear algebra on O(q^220) and verified on O(q^520).
* No relation of the next lower degree exists. This is the hostile control: M7 has no degree-23 relation and M9 no degree-35 relation, so neither is a coordinate on a smaller curve.

| curve | coordinate (identity products) | deg j | cusps (denominator B) | numerator A |
|---|---|---|---|---|
| X_1(5) (control) | x = 1/M5 = q^{−1}(G/H)^5 = r^{−5} | 12 | x (x^2 − 11x − 1)^5 | (x^4 + 228x^3 + 494x^2 − 228x + 1)^3 (Klein's icosahedral form) |
| X_1(7) | M7 = q^{−1} RS_2^2 RS_3 / RS_1^3 = q^{−1} f_2^2 f_3/f_1^3 | 24 | x (x − 1)(x^3 − 8x^2 + 5x + 1)^7 | (quadratic)^3 (sextic)^3 |
| X_1(8) | x = 1/M8 = q^{−1}(GG_3/GG_1)^2 = v^{−2}, v the Göllnitz–Gordon continued fraction | 24 | x (x − 1)^2 (x + 1)^4 (x^2 − 6x + 1)^8 | (octic)^3, palindromic |
| X_1(9) | M9 = q^{−1} B_2 B_4 / B_1^2 = q^{−1} f_2 f_4/f_1^2 | 36 | x (x − 1)(x^2 − x + 1)^3 (x^3 − 6x^2 + 3x + 1)^9 | (cubic)^3 (nonic)^3 |

* The exponents in B are the cusp widths.
* The width-N cusps are defined over Q(ζ_N)^+, as they should be. Field discriminants: 49 for x^3 − 8x^2 + 5x + 1, 81 for x^3 − 6x^2 + 3x + 1, 8 for x^2 − 6x + 1.
* Norm maps to X_0(N):
  * t7 = (η(τ)/η(7τ))^4 = (M7^3 − 8M7^2 + 5M7 + 1)/(M7^2 − M7);
  * t9 = (η(τ)/η(9τ))^3 = (M9^3 − 6M9^2 + 3M9 + 1)/(M9^2 − M9);
  * t8 = 1/M8 + M8, with j = P8(t8)/Q8(t8) of degree 12.
* The classical relations j(τ) = (t^2+13t+49)(t^2+245t+2401)^3/t^7 and j(7τ) = (t^2+13t+49)(t^2+5t+1)^3/t (t = t7), and their level-9 analogues, are reconfirmed.
* The minimal-model characters of S55 have the same ratios, because (q;q)_∞ cancels. So M7 = q^{−1} χ_2^2 χ_3/χ_1^3 for the M(2,7) characters χ_a ∝ θ_7(a)/(q;q)_∞. The characters are coordinates on X_1(7). This is the bridge between the conformal lane and this one: the same q-series in a different carrier (META-PATTERNS: "the same representation is not the same carrier").

**Diamond operators** (FINITE-EXACT). The diamond ⟨d⟩ permutes the Klein forms, f_a ↦ f_{da}.
* On M7, ⟨2⟩ is x ↦ (x − 1)/x, which has order 3. It permutes the width-1 cusps 0, 1, ∞.
* On M9, ⟨2⟩ is the same Möbius map, x ↦ (x − 1)/x.
* On x = 1/M8, ⟨3⟩ is x ↦ 1/x. This is equivalent to the palindromic j-map.

## 6. Supersingular content at the identity levels

### 6.1 Counts: the mod-12 correction lives only at elliptic points (PROVED, classical; FINITE-EXACT)

**Mechanism.**
* The Eichler–Deuring mass is ∑_E 1/|Aut(E)/±1| = (p − 1)/12.
* A level structure of index d spreads a supersingular point with e = |Aut(E)/±1| into d/e points, except at elliptic points of the level group.
* Γ_1(N) is torsion-free for N ≥ 4, so #X_1(N)^ss = (p − 1) d_N/12 exactly, with d_N = [PSL_2(Z) : Γ̄_1(N)].
* For X_0(N) the correction terms are (ν_2/4)(1 − (−4/p)) and (ν_3/3)(1 − (−3/p)).

| curve (identity) | d | #ss mod p, p ∤ N, p ≥ 5 | correction | tested (`level_supersingular.gp`, T1/T2, all p < 400) |
|---|---|---|---|---|
| X(1) | 1 | ⌊p/12⌋ + {0,1,1,2} | mod 12 (j = 0 and 1728) | section 2 (p < 3000) |
| X_1(5) (Rogers–Ramanujan) | 12 | p − 1 | none | 75 primes |
| X_1(7) (RS; R14 at 2τ) | 24 | 2(p − 1) | none | 75 primes |
| X_1(8) (GG; JS at τ + 1/2) | 24 | 2(p − 1) | none | 76 primes |
| X_1(9) (Bailey; Dyson at 3τ) | 36 | 3(p − 1) | none | 76 primes |
| X_0(7) | 8 | (2(p − 1) + 4[p ≡ 2 mod 3])/3 | mod 3: ν_3(7) = 2 because (−3/7) = 1; ν_2(7) = 0 because (−4/7) = −1 | 75 primes |
| X_0(8), X_0(9) | 12 | p − 1 | none | 76 primes each |
| X_0(14), X_0(27) | 24, 36 | 2(p − 1), 3(p − 1) | none (ν_2 = ν_3 = 0) | Eichler's formula only; the identity products are pull-backs from levels 7 and 9 (section 5) |

This is the whole precise content of "the moduli against the supersingular count mod 12":
* the period 12 is lcm(4, 3), where 4 and 3 are the discriminants of Z[i] and Z[ω], the two orders with extra units;
* the correction survives at level N exactly when N has elliptic points of that type;
* among the identity levels only X_0(7) has them.

### 6.2 The F_p / F_{p^2} split: Frobenius is the diamond operator (PROVED; FINITE-EXACT)

**Proposition.** Let p ∤ N and N ≥ 4. On the supersingular locus of X_1(N) over F̄_p, the p^2-Frobenius acts as the diamond operator ⟨p⟩.

*Proof.*
1. Let E/F_{p^2} have Frobenius π_{p^2} = −p. Every subgroup of E is then π-stable, so every isogeny out of E is defined over F_{p^2}, and its target again has Frobenius −p.
2. The supersingular isogeny graph is connected, so every supersingular j has such a model.
3. Then σ(E, ±P) = (E, ±π P) = (E, ±pP) = ⟨p⟩(E, ±P). ∎

*Corollary.*
* Over F_{p^2}, the supersingular points of X_1(N) fall into Frobenius orbits of length k(p), the order of p in (Z/N)^×/±1.
* The only exceptions are over j = 0 or 1728, at points where an automorphism has an eigenvector modulo N.
* The supersingular points of X_0(N) are all F_{p^2}-rational, as is the case N = 1 (the brief's "F_p or F_{p^2}").

Tests (all primes 5 ≤ p < 400):
* **T3.** Every irreducible factor over F_{p^2} of the squarefree supersingular polynomial on the identity coordinate has degree k(p). The exceptional fibres are exactly where predicted:
  * X_1(7) over j = 0 for p ≡ 2 mod 3 and k = 3: degrees {1, 1, 3, 3}. The two linear factors lie over the two order-3 elliptic points of X_0(7).
  * X_1(5) over j = 1728 for p ≡ 3 mod 4 and k = 2: degrees {1, 1, 2, 2}.
  * X_1(8) and X_1(9): none. Neither i nor ρ has an eigenvector modulo 8 or 9.
* **T4.** x^{p^2} = μ^e(x) on every supersingular root, where ⟨p⟩ = ⟨g⟩^e for the diamond generator μ = ⟨g⟩ of section 5. The rule is consistent on 75–76 primes per level:

  | level | e = 0 | e = 1 | e = 2 |
  |---|---|---|---|
  | 7 (g = 2) | p ≡ ±1 | p ≡ ±2 | p ≡ ±3 |
  | 9 (g = 2) | p ≡ ±1 | p ≡ ±2 | p ≡ ±4 |
  | 8 (g = 3) | p ≡ ±1 | p ≡ ±3 | — |

**Reading for the identities.** The diamond operators permute the Klein forms, so Frobenius on the supersingular locus permutes each family's product sides through the class of p in (Z/N)^×/±1.
* The three Rogers–Selberg products are cycled for p ≢ ±1 (mod 7).
* The three Bailey products are cycled for p ≢ ±1 (mod 9).
* The two Göllnitz–Gordon products are exchanged for p ≡ ±3 (mod 8).

The identity index set is the Galois group of the F_{p^2}-structure. At level 1 that group is trivial, which is why every supersingular j lies in F_{p^2}.

### 6.3 Class-number content: F_p-rational supersingular points (PROVED modulo CITED inputs; FINITE-EXACT)

**Theorem.** Let p ≥ 5, p ∤ N, and let X be X_0(N) or X_1(N). Then
```text
#X^ss(F_p) = ( h(−4p) s_X(Z[√−p]) + [p ≡ 3 mod 4] h(−p) s_X(Z[(1+√−p)/2]) ) / 2,
```
where s_X(O) is the number of level structures of type X fixed by π on O/NO. For X_0(N) these are the π-stable cyclic subgroups of order N. For X_1(N) they are the ±-classes of points of exact order N with πP = ±P.

*Ingredients.*
* (a) **Delfs–Galbraith** (CITED). For a supersingular curve over F_p, End_{F_p}(E) is Z[√−p] (the floor) or, only for p ≡ 3 mod 4, Z[(1+√−p)/2] (the surface). There are h(O) F_p-classes with End_{F_p} = O, and each supersingular j ∈ F_p has two F_p-classes. The induced count n_p = (h(−4p) + [p ≡ 3 (4)]h(−p))/2 is checked in section 2.
* (b) **Tate's theorem** (CITED). T_ℓ E is free of rank one over O ⊗ Z_ℓ for ℓ ≠ p, since quadratic orders are Gorenstein. So π acts on E[N] as multiplication on O/NO.
* (c) **A twisted Burnside lemma** (elementary). The F_p-points of X over a supersingular j are the Frobenius-stable orbits of G = Aut(E)/±1 on level structures. Their number is (1/|G|) ∑_{g ∈ G} #Fix(g π). The two F_p-twists of j correspond to the classes of gπ, which gives the factor 1/2 over the two F_p-classes per j, including at j = 0 and 1728.

*Local numbers* (exact enumeration on (Z/N)^2, `level_fp_counts.gp`, `level_supersingular.gp`):

| X | s(floor) | s(surface) | #X^ss(F_p) |
|---|---|---|---|
| X_0(ℓ), X_0(ℓ^2), ℓ odd (here 7, 9) | 1 + (−p/ℓ) | 1 + (−p/ℓ) | (1 + (−p/ℓ)) n_p |
| X_1(7), X_1(9) | 6[p ≡ −1 mod N] | 6[p ≡ −1 mod N] | 6 [p ≡ −1 mod N] n_p |
| X_0(2) | 1 | 3 | (h(−4p) + 3[p ≡ 3 (4)] h(−p))/2 |
| X_0(4) | 0 (p ≡ 1 mod 4), 2 (p ≡ 3 mod 4) | 0 (p ≡ 3 mod 8), 4 (p ≡ 7 mod 8) | 3 [p ≡ 3 mod 4] h(−p) |
| X_0(8) | 4[p ≡ 7 mod 8] | 4[p ≡ 7 mod 8] | 4 [p ≡ 7 mod 8] h(−p) |
| X_1(8) | 4[p ≡ 7 mod 8] | 8[p ≡ 15 mod 16] | (2[p ≡ 7 mod 16] + 6[p ≡ 15 mod 16]) h(−p) |

Verification:
* **Fibre by fibre** on the Hauptmodul lines for X_0(7), X_0(9) and X_1(7): every supersingular j ∈ F_p for p < 400 matches, including the 77 fibres over j ∈ {0, 1728}.
* **Totals** for X_0(2), X_0(4), X_0(8), X_1(8) and X_1(9) on all 301 primes p < 2000.
  * For X_1(8) the obvious ansatz "a constant multiple of n_p for each p mod 8" fails at p ≡ 7 mod 8 (ratios 2 and 6 both occur).
  * The local computation explains the failure: at p ≡ 7 mod 8 the surface number depends on p mod 16. The mod-16 rule then holds on all 301 primes.
* **The local formula** on 1505 (curve, prime) pairs.
* **Independent check through Atkin–Lehner traces** (`al_genus_check.gp`).
  * The Deligne–Rapoport picture gives g(X_0(Np)/w_p) = g(X_0(N)) + (H(p,N) − #X_0(N)^ss(F_p))/2, the level-N form of g(X_0(p)^+) = (ss − ss(F_p))/2.
  * The left side is computed from PARI's Atkin–Lehner operator on S_2(Γ_0(Np)), independently of the Hauptmodul root counts.
  * Result: all 114 pairs agree (N ∈ {1, 2, 4, 7, 8, 9}, 5 ≤ p < 150, Np ≤ 300; 0 mismatches). Example: g(X_0(46)/w_23) = 0 with H(23, 2) = 6 = #X_0(2)^ss(F_23).

For ℓ = 7 and 9 the formula reads #X_0(ℓ^e)^ss(F_p) = (1 + (−p/ℓ)) n_p. For example, X_0(7) has no F_p-rational supersingular point when (p/7) = 1.

**Level-N Ogg lists** (p < 2000): the primes for which every supersingular point of the curve is F_p-rational.
* X_0(2): {5, 7, 11, 23}.
* X_0(4): {7}.
* X_0(7), X_0(8), X_0(9), X_1(7), X_1(8), X_1(9): empty.

By the genus relation above, this is the generalised Ogg condition g(X_0(Np)/w_p) = g(X_0(N)). For N ≥ 2 the list is finite, since n_p = O(√p log p) while H(p, N) is of order ψ(N) p/12. At the identity levels the split is F_{p^2}-heavy for every p.

## 7. Connection ledger

Columns: source → target, map, what is preserved and lost, cheapest decisive test, verdict.

| # | claimed connection | map | preserved / lost | decisive test | verdict |
|---|---|---|---|---|---|
| C1 | F, U, S at p^2, pq, p^2 q → supersingular counts | none possible | profile only / every prime | F(p^2) = 1 for all p versus ss(p) ≥ (p−1)/12 | REFUTED as an identity; dictionary via X_0(N) components and cusps (section 4) is CITED structure |
| C2 | F, U, S → F_p/F_{p^2} split | none | — | the split depends on (p mod N, class numbers); F is constant in p | REFUTED |
| C3 | six families → three theta families | q ↦ q^2, q^3, −q | theta numerators / denominators | 3 + 3 + 2 sum-side corollaries to O(q^700) | STRUCTURAL (exact) |
| C4 | product ratios → X_1(7), X_1(8), X_1(9) | Klein-form monomials M7, 1/M8, M9 | function field / fermionic sum sides | j = A(M)/B(M) of degree d_N to O(q^520); lower degrees fail | STRUCTURAL (FINITE-EXACT; classical in substance) |
| C5 | moduli → "supersingular count mod 12" | mass formula on X_1(N), X_0(N) | the elliptic-point correction | root counts on the coordinates, 75–76 primes each | PROVED + FINITE-EXACT: no correction on X_1(N); only X_0(7) keeps mod 3 |
| C6 | identity index (Z/N)^×/±1 → F_{p^2}-Frobenius | Frobenius = ⟨p⟩ | orbit lengths; exceptional fibres | factor degrees over F_{p^2}; x^{p^2} = μ^e(x) | PROVED + FINITE-EXACT |
| C7 | moduli → class numbers | (h(−4p) s_floor + h(−p) s_surface)/2 | local factors at N | fibre counts, totals, local formula, Atkin–Lehner traces | PROVED (mod CITED) + FINITE-EXACT |
| C8 | 14, 27 → the genus-one curves X_0(14), X_0(27) (27a: CM by discriminant −27) | none | — | C3: the products are θ_7(a; q^2), θ_9(a; q^3); no elliptic-curve structure of X_0(14) or X_0(27) is used | NUMEROLOGY |
| C9 | moduli 7, 9, 14, 27 → class-number-one discriminants | none | — | −9 and −14 are not discriminants; h(−36) = 2, h(−56) = 4 | NUMEROLOGY |
| C10 | shapes 9 = p^2, 14 = pq, 12 = p^2 q (the period of ss(p)); 7 = p, 8 and 27 = p^3, 60 = p^2 q r; \|PSL_2(F_l)\| = 6, 12, 60 = pq, p^2 q, p^2 q r (defects −2, −1, 0) | none for F, U, S | — | F(12) = φ(12) = 4 but F(18) = 4, φ(18) = 6; A_5 has no subgroups of order 15, 20, 30 although F(60) counts them | NUMEROLOGY; the genuine parts are 12 = \|D_{−4}\|·\|D_{−3}\| and #ss(X(l)) = (p−1)\|PSL_2(F_l)\|/12 |
| C11 | Ogg's fifteen primes → the identities | none | — | the level-7, 8, 9 Ogg lists are empty (section 6.3) | NUMEROLOGY for the list (agrees with S55 §7); the genuine generalisation is the level-N Ogg condition |
| C12 | S55's characters (conformal lane) → this lane | the same q-series | characters ↔ coordinates on X_1(N) | M7 = q^{−1} χ_2^2 χ_3/χ_1^3 identically | STRUCTURAL (same series, different carrier) |

## 8. Inheritance pass and board

* **Closest proved mechanism:** DB1/DB2 (divisor balance and its defect table, which already lists p^2, p^2 q at −1); [ogg_triangular_chowla_20260927](ogg_triangular_chowla_20260927.md) for Ogg's list and the class-number reading.
* **Canonical hostile:** p^2, where 1 ≠ 1 + 1. Also the literal exponent family p^q r.
* **Corrected near miss:** the "Chowla reading" of F = U + S as σ(n) − n − 1. It is a different object and is not used.
* **Least-used sidecar:** the cusp and Atkin–Lehner structure of X_0(N), indexed by divisors, which ties THM-4733's complement C to the Fricke involution.
* **Board:** exponent box / defect level sets / X_0(N) components and cusps / Klein-form coordinates / diamond = Frobenius / floor–surface class numbers.
* **META-PATTERNS cards used:**
  * "Search the statement before the method": the values were already in DB2.
  * "Type every analogy and every implication": the ledger in section 7.
  * "The same representation is not the same carrier": C12.
  * "Treat pulls and collisions as synthesis opportunities": S55 arrived mid-session.

## 9. What is not claimed, and next questions

**Not claimed.**
* No novelty. Everything proved here is classical in substance: Deuring, Eichler, Deligne–Rapoport, Katz–Mazur, Kubert, Delfs–Galbraith, Ogg's fixed-point counts for w_p.
* No consequence for LRC(14), Collatz or tournaments.
* No arithmetic consequence of F = U + S. THM-4733's "NOT PROVED: any arithmetic consequence for F = U + S" stands, and C1 shows that none of the identity type can exist.
* The KLS paper is recorded, not audited.
* No HYP was opened. Every formula that started as a fit (X_1(8) by p mod 16, the 2-power tables) is a theorem once the local numbers are computed.

**Next questions.**
1. Which supersingular paper did the user intend (section 1)?
2. Extend the dictionary to the other genus-zero levels with Rogers–Ramanujan-type families: Rogers mod 10 (DS15 §2.10) for X_1(10), and the mod 12 families (§2.12) for X_1(12).
3. Do the sum (fermionic) sides carry anything at supersingular points? For example, a congruence for the sum sides mod p analogous to E_{p−1} ≡ 1. None was tested.
4. Locate the literature form of the class-number formula of section 6.3. It should be the fixed-point count of w_p on X_0(pN) (Ogg, Kluit, Furumoto–Hasegawa); not checked here.

## 10. Reproduction and run log

Commands, run order and times are in the experiment [README](../../04-computation/experiments/supersingular_fus_qseries_20261022/README.md). All outputs end with `ALL CHECKS PASSED` unless stated below.
* `fus_shapes.out`, `supersingular_census.out`, `qseries_identities.out`, `level_hauptmoduln.out`, `level_supersingular.out`, `level_fp_counts.out`, `dr_genus.out`, `numerology_checks.out`: all pass.
* `al_genus_check.out`: 114 pairs, 0 mismatches, ALL CHECKS PASSED. The run took 270 s with about 1 GB peak memory; a 512 MB PARI stack overflowed when the level cap was 700.
* Audit: pending (an independent adversarial audit is run before this note is marked audited).
