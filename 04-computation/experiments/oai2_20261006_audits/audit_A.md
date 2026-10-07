# Audit A — THM-4558, THM-4559, THM-4561, THM-4562, MISTAKE-576, THM-418 update, results note §§3/4/7

Auditor: independent adversarial audit (read-only on worktree `math-wt-chessboard-20261006`), 2026-10-06/07.
All checks below were re-done with my own code (no repo script was imported or run):

| script (this folder) | what it checks | result |
|---|---|---|
| `auditA_gf.py`, `auditA_kappa.py`, `auditA_kappa_big.py` | own GF(p^n) arithmetic; κ(q)=χ(Cay(F_{q²},μ_{q+1})) by SAT (glucose4 / cadical195, different encoding and symmetry breaking from the repo) | κ = 4,3,4,4,4,4,3,5,4 for q = 2,3,4,5,7,8,9,11,16; κ(13) ∈ [5,6]; **κ(17) ∈ [5,6]** (a proper 6-colouring found, 18 s); k = 5 runs for q = 13, 17: see end |
| `auditA_hoffman.py` | λ_min of Cay(F_{p²},μ_{p+1}) (numpy + 50-digit mpmath at the minimiser), p < 140 | κ(p) ≥ 5 for all primes 23..139, ≥ 6 for 53..139 (≥ 7 at 71 and at every prime 89..139, ≥ 8 at 131, 139); margin at p = 31: 4.39·10⁻⁵ |
| `auditA_multiquad.py`, `auditA_ladder.py` (+ `.out`) | own decomposition/inertia groups of every prime ≤ 400 in each multiquadratic field, c ∈ I / c ∈ D∖I test | reproduces every upper bound in THM-4558's table and §4.4's list |
| `auditA_spindle.py` | generalized spindle, exact arithmetic in Q(√−3,√m), all exact unit edges, SAT | N = 3, 9, 21, 27: not 3-colourable |
| `auditA_tower.py`, `auditA_walsh_const.py`, `auditA_bent16.py` | THM-4561 closed form, transfer identity, F bounds, RS-switched sup, Walsh constant (exact recursion to k = 22), bent count | all claims reproduced except the F-range (below) |
| `auditA_jc_fibres.py` | THM-1300 fibre counts, brute force over GF(q) | q = 5,7,11,13,**25,49,121,125** all match |
| `auditA_misc.py` | note §7.1 mod-7 reduction; THM-4559 idempotent certificate for Pálvölgyi's heptagon with algebraic r = 3 | both pass |

External sources read (data only): Madore arXiv:1509.07023 (full text); Exoo–Ismailescu arXiv:1805.00157 (full text); Pálvölgyi arXiv:2609.23327 (full text incl. appendix); ACMVZ arXiv:2207.14179 (abstract); Dúcz arXiv:2606.12325; #172 sources in scratch (`hn/p172/01,02,07`) and `EuclideanRamsey.lean`; github.com/decalion89/chromatic-number-of-the-plane (README, `notes/local_colourings.md`, `notes/literature.md`; commit dates 2026-09-27..29); github.com/MildlyMeticulous/hn-2adic-obstruction (README, RESULTS.md; single commit 2026-07-30); A. Goucher, cp4space 2018-05-19 ("Royal Wedding and Polymath16").

---

## 1. THM-4558 (residue colourings of field planes) — **PASS WITH CORRECTIONS** (mathematics sound; one false claim; major missing attribution)

**Verified correct.**
* Lemma R (1)–(3): `v_𝔓(cu) = v_{c𝔓}(u) = v_𝔓(u)` ⇒ `v_𝔓(u)=0`; c̄ on k_𝔓; components = cosets of Γ ⊂ O_𝔓; edges change the residue by a nonzero element of U_𝔓. Split ⇒ not c-stable; ramified ⇒ c ∈ inertia ⇒ c̄ = id ⇒ U = {±1} (χ ≤ 2 in char 2, odd cycles of length p otherwise); inert ⇒ c̄ = Frobenius, U = μ_{q+1}. Correct.
* Every row of the values table, re-derived by hand and by `auditA_multiquad.py`: Q(√−d) (2 ramified iff d ≡ 1,2 mod 4; else an odd ramified prime); Q(ζ_8) (2 ramified over Q(√2)); Q(ζ_12) ((√3) inert, k = F_9); M (2 splits in Q(√33), D = Gal(M/Q(√33)) ∋ c, I = 1, k = F_4); M(i) (I = Gal(L/M), D = Gal(L/Q(√33)), c ∈ D∖I); L_N, N ≡ 2 mod 3 (prime over 3: e = f = 2, g = 1, inertia field Q(√−(4N−1)), c ∉ I); L_N, N = 3n, n odd (−(12n−1) ≡ 5 mod 8, 2 inert in both, splits in L⁺); Heegner compositum (all six d ≡ 3 mod 8, Frob_2 = c, the *only* c-stable prime below 400 is 2); Polymath field (11: inertia field Q(√−3,√5), decomposition field Q(√5), c ∈ D∖I, k = F_121). All correct.
* κ table values for q ≤ 11, 16 (third independent SAT computation); κ(3): Cay(F_9,μ_4) = K_3□K_3; Hoffman claims (23 ≤ p ≤ 79 ⇒ ≥ 5; 53 ≤ p ≤ 79 ⇒ ≥ 6; margin at 31 = 4.4·10⁻⁵).
* Knight torus: knight set = {x²+y² = 5} = (1+2i)μ_8 ⊂ F_49 (THM-4552 (ii)), so G_7 ≅ Cay(F_49, μ_8), χ = κ(7) = 4.
* ω_t: minimal polynomial t z² − (2t−1)z + t, discriminant −(4t−1) ⇒ Q(ω_t) = Q(√−(4t−1)); matches the Polymath16 definition (Goucher: "real part 1 − 1/(2t) and absolute value 1").
* Lower bounds: Moser spindle in M; generalized spindle built exactly for N = 3, 9, 21, 27 (n = 1, 3, 7, 9), each non-3-colourable. Polymath field: Goucher (cp4space 2018) states "Several of Marijn Heule's 5-chromatic graphs lie in Z[ω₁, ω₃, ω₄]" and that Z[ω₁,ω₃,ω₄] (and with ω₇) has homomorphic 5-colourings; hn-2adic confirms Heule/Parts record graphs are two Q(√3,√11) halves glued by ω₄. Note de Grey's original graph also uses ω₁₆ ∈ Q(√−7) (Exoo–Ismailescu §5).
* Corollary 2 (χ ≥ 4 ⇒ L/L⁺ unramified at all finite primes ⇒ L ⊂ narrow HCF, h⁺(L⁺) even) and "M is the genus field of Q(√33)": correct.

**Errors.**
1. Line 81: *"The 5-chromatic step needs the class-number-2 field `Q(√−15)`."* — **FALSE.** Exoo–Ismailescu (arXiv:1805.00157; DCG 64 (2020) 216–226), §5: "our 5-chromatic graph can be embedded in Q[√3,√11,√247] × Q[√3,√11,√247]" (rotation arccos(119/128), sin = 3√247/128). A quarter turn puts it in G(Q(√−3,√−11,√−247)) (decalion89 Thm 2 says exactly this), a field with no √−15; h(−247) = 6 (computed). Lemma R at 11 gives χ ≤ 5 there (verified), so that field is also exactly 5-chromatic. decalion89 further reports (unrefereed, DRAT-certified) 5-chromatic graphs in Q(√−3,√−7,√−11) and in Q(ζ_21), neither containing √−15. The repo's own results note §4.2 already cites the √−247 field — internal inconsistency.
   Corrected text: "The 5-chromatic step is not tied to one field: 5-chromatic graphs are known in Q(√−3,√−11,√−15) (Heule's ω₄; h(−15) = 2) and in Q(√−3,√−11,√−247) (Exoo–Ismailescu 2020; h(−247) = 6)."
2. Line 75 (Corollary 1): *"A 5-chromatic graph needs a rotation `ω_t` with `t` even, such as Heule's `ω₄`."* — overclaim: Exoo–Ismailescu's escaping rotation is (119+3√−247)/128 = ω_{64/9}, t not an integer. Corrected: "So no 5-chromatic graph has all its edge vectors in this field; the known ones leave it through ω₄ (de Grey, Heule), ω₁₆ (de Grey) or (119+3√−247)/128 (Exoo–Ismailescu)."
3. Line 54/50: *"17 | … | 5–7"* and *"(FINITE-EXACT; reader's SAT run and this session's CaDiCaL run agree)"* — the 13 and 17 entries come from the reader only (the session script tests only q ≤ 11, 16); κ(17) ≤ 6 (audit A found a proper 6-colouring; decalion89 also reports 5–6). decalion89 reports κ(13) = χ(G₁₃) = 6 by a completed computer proof (`notes/g13_chi.md`, 2026-09-29; not re-checked here). Corrected: "13 | 5–6 (6 per decalion89)", "17 | 5–6".
4. Line 70: *"the generalized spindle (two Eisenstein rhombi with tips at squared distance `N`, `3 \| N`)"* — rhombi only for N = 3. Corrected: "two copies of a rigid triangulated Eisenstein patch containing the pivot 0 and a point v ∈ √−3·Z[ω] with |v|² = N, the second rotated by ω_N".
5. Lines 68, 72: *"the prime over 2 is fixed by c"*, *"the prime over 11 (`e = f = 2`)"* — there are two such primes (g = 2) in M(i) and in the Polymath field. Say "each prime over 2/11".
6. Line 85/86: *"Every finite unit-distance graph lives in some number field"*, *"its compactness graph `H` has real-algebraic coordinates"* — only a realization (possibly with extra edges) does (Tarski–Seidenberg; Madore §5.4). Corrected: "has a realization in some conjugation-stable number field" / "can be realized with real-algebraic coordinates".
7. Line 81: *"What governs χ here is how 2, 3, 5, 7 and 11 decompose in `L/L⁺`"* — only upper bounds are governed by decomposition. Corrected: "What bounds χ from above here is …".
8. Line 79: *"the compositum of all six is 4-colourable"* — say "all six with d ≡ 3 (mod 8)". Adding the seventh rung √−7 (N = 2) gives fields with no small c-stable prime (Q(√−3,√−7,√−11): first c-stable prime 17), where decalion89 reports a 5-chromatic graph; the refutation of HYP-2278 (4) is unaffected (six rotations, χ = 4 ≠ 8).

**Missing attribution / overclaimed novelty (serious).**
* **decalion89 `notes/local_colourings.md` (commits 2026-09-27..29, i.e. before this THM; the THM cites the same repository v1.1.0 only for two values)** already contains: the exact setting (K ⊂ C stable under conjugation, L = K ∩ R, Γ(K), T(K) = {uū = 1}); **Proposition A = Lemma R** (non-split place: unramified ⇒ χ ≤ χ(G_q), G_q = Cay(F_{q²}, N₁) — i.e. κ(q); ramified ⇒ χ ≤ 3), credited there to Madore Prop. 3.2/¶6.6; **Theorem 1** χ(Γ(Q(√−3,√−11))) = 4 by the place over 2 (= the M row); **Theorem 2** χ(Γ(Q(√−3,√−11,√−247))) = 5; a χ(G_q) table (2:4, 3:3, 4:4, 5:4, 7:4, 11:5, 13:6, 17:5–6, 19:5, …); **Proposition B** χ(G_q) ≥ 6 for every prime q ≥ 53 (Hoffman + Weil) and ≥ 6 for q = 29, 37, 41, 43, 47 (three-point bound) — supersedes the THM's "κ(p) ≥ 6 for 53 ≤ p ≤ 79 (computed)"; and §5 the same field screen as results-note §4.4.
* **MildlyMeticulous, "A 2-adic obstruction to 5-chromatic unit-distance graphs", github.com/MildlyMeticulous/hn-2adic-obstruction (2026-07, unrefereed)**: Theorem A/A′ = Lemma R at σ-stable places over 2 in the z = x+iy formulation, with the norm-one refinement and "ramified ⇒ χ ≤ 2"; Corollary B ⇒ χ(Q(√3,√11)²) = 4 (M, M(i) rows); **Corollary B′** (K_max ∋ √d for all squarefree d ≡ 1,3 mod 8) contains both the Heegner-compositum bound and THM-4558 Corollary 1; **Theorem C/C1** ("spindling at Löschian r² escapes iff r² is even") is Corollary 1 for spindles. decalion89 itself credits this repository as earlier.
* The coset/additive-colouring mechanism predates Madore: Woodall 1973 (Q² mod 2), K. G. Fischer, Discrete Math. 82 (1990) 181–195 (Thm 1), Moorhouse 2010 draft (Lemma 4.2, 8.2); see also Axenovich–Choi–Lastrina–McKay–Smith–Stanton, Graphs Combin. 2014, Thm 2.3. So *"The reduction technique is Madore's (2015), generalized … to any conjugation-stable field L"* (line 4) misassigns both the technique and the generalization.
* Fischer 1994 is second-hand: decalion89 reports it from the zbMATH review ("we have not seen its full text"); title "A planar geometric graph of chromatic number four", Congr. Numer. 104 (1994) 73–79 (Q(√p,√q)², p ≡ 3, q ≡ 11 mod 16, pq ≡ 1 mod 32). The status line should say so.
* Line 17 *"Gibbs's homomorphic 4-colouring of the ring Z[omega_1, omega_3]"*: Goucher (2018) credits P. Gibbs with "a homomorphism to a four-element group"; decalion89 credits D. Speyer (Polymath16 thread 2, comments 4013/4206) with the reduction mod 2. Cite both.
* What remains new after these credits: the Heegner-rung values (χ = 3 for N ≡ 2 mod 3 via the prime over 3), the χ = 4 family N = 3n (n odd Löschian) with the generalized spindle (upper bound = hn-2adic C1), the explicit Heegner-compositum and Polymath-field statements, Corollary 2, and the refutation of the repo's HYPs.

**Typing.** Lemma R, κ(q) for q ≤ 11, M, M(i), the Q(√−3,√−11,√−247) analogue, Corollary 1 and the Heegner-compositum bound should be KNOWN (with the credits above), not new PROVED. The Hoffman bounds are floating-point (NUMERICAL with margin ≥ 4·10⁻⁵, or cite decalion89's interval-arithmetic version), not FINITE-EXACT.

---

## 2. THM-4559 (algebraic spherical sets are Ramsey, conditional on #172) — **PASS** (two wording fixes)

* Criterion transcription matches #172 Thm "Classification" and the Lean comparator `FieldCriterion` exactly: conditions only on the spatial block ("There are no separate conditions on the multiplied constant row or constant column of P"); F = Q(coordinates); s ≥ 2, injective, affine span = R^d.
* Deduction checked line by line: F number field ⇒ F⊗F ≅ F × (rest), m = projection onto F, ker m = (1−e)B so e·(x⊗1 − 1⊗x) = 0 and e(x⊗y) = e(xy⊗1). H ∈ Mat_{d+1}(F) with spatial block I_d exists: the centre solves 2(a_i−a_0)·c = |a_i|²−|a_0|² over F (unique by affine spanning; Cramer), r² = |a_0−c|² ∈ F; constant row/column of H are unconstrained by the criterion. P = e(H⊗1) satisfies both conditions. Correct.
* Representative independence: #172 Remark "alg:congruence" proves it directly. Algebraic squared distances ⇒ Gram matrix over R_alg (real closed) ⇒ congruent copy with real-algebraic coordinates. Correct; nothing hidden (centre and radius² lie in F).
* Independent exact check (`auditA_misc.py`): Pálvölgyi's heptagon with **algebraic** r = 3 (F = Q(√2,√5)): e idempotent, m(e) = 1, all 7 tensor evaluations vanish with e, 4 of 7 are nonzero without e. Consistent with "the heptagon needs a transcendental radius".
* No contradicting result found: #172 never states the corollary (its Prop. "quadratic independence" needs a rank hypothesis that the idempotent makes unnecessary); Pálvölgyi's Theorem 1 is for transcendental r > 2 (confirmed), and his unchecked appendix only says that derivation certificates cannot exclude circle configurations with at most three non-algebraic points after a Möbius map (Thm A.10), and that derivative-based colourings cannot exclude B ∪ {p*} with B algebraic (Thm A.17). It does not assert Ramseyness.
* Corrections: (a) *"At most 5 concyclic points are Ramsey (#172)."* reads as an upper limit (contradicting the theorem itself); write "Every set of at most 5 concyclic points is Ramsey (#172)". (b) *"Not found stated in #172 or in Palvolgyi's abstract"* can be strengthened to "nor anywhere in Pálvölgyi's paper (main text or appendix)". Optional: "transcendental" in "Every non-Ramsey spherical set is transcendental" should read "has a transcendental squared distance"; subsets of regular polygons (P_7 heptagon, hexagon) are Ramsey unconditionally (Kříž 1991).
* Typing CONDITIONAL / deduction PROVED: correct.

---

## 3. THM-4561 (Sierpinski skew-Hadamard tower) — **PASS WITH CORRECTIONS** (minor)

* Doubling: H^T = 2I − H gives c_j = 2z^j − r_j; lower rows (−c_j, c_j) = (r_j − 2z^j)(1 − w); columns c_j ∓ w r_j. Correct.
* Closed form: induction re-derived (l = k term gives (−1)^{y_k}[y′ = j]); checked k ≤ 9. Correct.
* Transfer identity re-derived: |y|² = A(2 + ε(z^m + z^{−m})), A = |x|²; constant coefficient of |y|⁴ = 6Σ C_u² + 8εS (the z^{±2m} terms vanish since deg A < m), coefficient of z^{−2m} = ‖x‖₄⁴ + 4εS. Hence R′ = (3/2+2εs)R, s′ = f(εs). CS: |S| ≤ Σ_{u≥1}C_u² = (‖x‖₄⁴ − m²)/2 ⇒ F′ ≤ 2F/(F+1). Walsh: (3/2+2t)(3/2−2f(t)) = (7+4t)/4 ≥ 3/2 on [−1/4,1/4]. All correct.
* Golay claim: r_i(1) = 0 for i ≠ 0 follows from the closed form (each spike sum is a nontrivial character sum); proof correct.
* RS bound re-derived: the l = 0 spike flips the Walsh part (±RS⊙w_i, Golay); for l ≥ 1 on y = a + 2^l t the RS form is const + a_{l−1}t_0 + Σ_j t_j t_{j+1} (path) and (−1)^{Σ_{m≥l} i_m y_m} is linear in t, so G_l is a Davis–Jedwab Golay sequence of length N/2^l, |G_l| ≤ √(2N)·2^{−l/2}; sum √(2N)(1 + 2(√2+1)) = (3+2√2)√(2N) ≈ 8.24√N; F ≥ 1/(2(3+2√2)² − 1) = 1/(33+24√2). Correct.
* FINITE-EXACT H_16: 0 (tower) vs 448 (Sylvester), best min |(Hd)_i| = 2: reproduced. NUMERICAL rates reproduced; best Walsh constant from the exact recursion over all sign patterns: 0.9377 (k=13) … 0.9074 (k=22), decreasing → ≈ 0.905; tower/Walsh ratio 1.47–1.49 (k = 11–13); columns reach F = 4 at k = 2, 3.
* Errors: (a) *"NUMERICAL (`k ≤ 12`): the worst row has `sup ≤ 3.0√N`, and `F ∈ [0.75, 3.0]` with mean 1.30."* — [0.75, 3.0] is the k = 12 range only ([0.7527, 3.0007]); over 3 ≤ k ≤ 12 the repo's own `.out` gives F ∈ [0.6975 (k=6), 3.2000 (k=4)]. Corrected: "F ∈ [0.70, 3.20] for 3 ≤ k ≤ 12 (at k = 12: [0.753, 3.001], mean 1.30)". (b) *"the RS-switched rows are √2-flat Golay mixtures"* — proven constant is (3+2√2)√2 ≈ 8.24 (observed ≈ 3.0); say "sums of √2-flat Golay pieces (sup ≤ 8.24√N)". (c) Item 1 *"so it does not flatten"* is not proved — type NUMERICAL (column F: 4, 4, 1.10, 0.68, 0.46, 0.30, 0.23 for k = 2..8). (d) The title states "max F → 0, sup/√N → ∞" without the "up to a perturbation sketch" qualifier carried in the status; add it (the sketch must also treat small k, where O(m^{−1/2}) is not small, by direct computation).

---

## 4. THM-4562 (slice lemma; THM-1300 fully non-slice) — **PASS WITH CORRECTIONS**

* **Error (lemma statement).** *"`V_j` is locally nilpotent **iff** `C[x] = R[F_j]` for some subring `R` containing every `F_i` (`i ≠ j`)."* — as written the right side always holds (take R = C[x]), so (⇐) would make every dual field locally nilpotent; the (⇐) proof uses ∂/∂F_j over R, which needs F_j transcendental over R. Corrected: "… iff `C[x] = R[F_j]` with `F_j` algebraically independent over `R` (i.e. `C[x] = R^{[1]}` in the variable `F_j`) for some subring `R` containing every `F_i`". Same wording in results note §6.3 and the reader's report.
* Everything else in the lemma is right: V_j F_i = δ_ij with V = (JF^T)^{−1}∂ (index convention checked); slice theorem; derivations determined on the separable algebraic extension C(x)/C(F); C^n ≅ X_j × A¹, G étale, injectivity equivalence; n = 3 via Fujita/Miyanishi–Sugie; n = 5: étale + injective ⇒ open immersion, O(X)^× = C^×, Hartogs ⇒ X ≅ A⁴. det JF = −2, weights (−2,−1,1), V_j polynomial — verified symbolically.
* Fibre counts re-derived for **every** q (not only primes): F_3 = x(2 − 3xy − x²z) (zero fibre {x=0} ∪ graph: q² + q(q−1); c ≠ 0 forces x ≠ 0: q(q−1)); F_1 is z-linear with coefficient u³ (u ≠ 0: q² − (q−1) points; u = 0 adds q(q−1) points only when c = 0); F_2 is z-linear with coefficient 3xu² (xu ≠ 0: (q−1)²; x = 0: F_2 = y, q points; u = 0: F_2 = −2y, q points iff c ≠ 0). F_1, F_3 hold for all q; F_2 needs char ∤ 6. Brute force over GF(q): q = 5, 7, 11, 13, 25, 49, 121, 125 all match. So the counts are PROVED for gcd(q,6) = 1, and the "Fibre geometry (reader; not re-derived here)" items (A² ∖ {xy = −1}, C^× × C) are re-derived (PROVED).
* Spreading out is valid (residue fields of a f.g. Z-algebra are F_{p^f}, which is why the all-q counts are needed — now covered). Stable-coordinate argument correct.
* *"By Dutta–Lahiri, no component is residual over a coordinate line either."* — needs the reference and definition: A. K. Dutta, A. Lahiri, On residual and stable coordinates, J. Pure Appl. Algebra 225 (2021) 106707, Thm 3.4 (residual coordinate over k[p] ⇒ stable coordinate); then the claim follows from "not a stable coordinate".
* Russell cylinder: CY = 16xy = 16z(z+1) = B(B+4); Y_1, Y_2 and THM-3561's surface c²e = b(b+4) match THM-3605 (3); #047's R = k[x,y,z,u]/(xy − z(z+1)) matches its §1. Correct.

---

## 5. MISTAKE-576 — **PASS WITH CORRECTIONS**

* Quotes of HYP-2277 (S699m), HYP-2278 (4) and HYP-2276 match `05-knowledge/hypotheses/INDEX-HISTORICAL-THROUGH-2026-07-21.md` (lines 3545–3548, mojibake-encoded). Refutations of the Heegner roadmap: correct.
* Facts: Falconer 1981 χ_m ≥ 5 — correct; ACMVZ arXiv:2207.14179 abstract: "upper bound of 0.2470", confirming m₁ < 1/4 ⇒ χ_m ≥ ⌈1/m₁⌉ = 5 — correct; Croft m₁ ≥ 0.2293 > 1/5 ⇒ density alone cannot give 6 — correct.
* **Error:** *"The 5-chromatic graphs need the class-number-2 field `Q(sqrt-15)` (Heule's `omega_4`)"* — false (Exoo–Ismailescu, Q(√−3,√−11,√−247), h = 6); replace as in §1 error 1.
* **Fairness of the HYP-2278 (1) quote:** the ellipsis removes the load-bearing step "χ_f(ℝ²)=1/m₁ ≤ 4.36 < 5 ⟹". State precisely: (i) the fractional statement is TRUE (χ_f(R²) ≤ 4.36 < 5; now 4 ≤ χ_f(R²) ≤ 4.36, arXiv:2311.10069); (ii) "χ_f(R²) = 1/m₁" is not a known identity (1/m₁ is the *measurable* fractional bound); (iii) the actual error is integrality: m₁ < 1/4 forces χ_m ≥ 5, and χ_m ≥ 5 was already Falconer 1981; (iv) "the χ ≥ 5 lower bound … is IRREDUCIBLY COMBINATORIAL … no analytic method can narrow it", read for χ(R²) (all colourings), is refuted only CONDITIONALLY (by #158's transfer theorem); unconditionally χ(R²) ≥ 5 is still known only combinatorially.
* **Incomplete:** the header lists HYP-2277 "claude S688", but its most directly refuted claim is not quoted: "the χ=4 junction field ranges over EXACTLY the class-number-one (UFD) imaginary-quadratic fields Q(√−7,−11,−19,−43,−67,−163)". THM-4558 shows Q(√−3,√−(4N−1)) is 3-chromatic for N ∈ {2,5,11,17,41}; these N are not even Löschian, contradicting the HYP's own requirement that N be an Eisenstein norm. Add it.

---

## 6. THM-418 update block — **PASS WITH CORRECTIONS** (wording)

* Correct: Q(√−d) ≤ 3, bipartite for d ≡ 1,2 (mod 4); Q(ζ_8) = 2 (Madore Prop. 3.9), Q(ζ_12) = 3 (Prop. 4.2); the two recursions are the ramified cases at (1+i) ((a+bi) mod (1+i) = (a+b) mod 2) and at (√−3) ((a+bω) mod √−3 = (a−b) mod 3, ω = e^{iπ/3}).
* *"(which contains both lattices at every scale and rotation)"* — false literally (irrational rotations/scales leave Q(ζ_12)). Corrected: "(which contains a similar copy of U(Z², D) and of U(Z[ω], D) for every norm D, via z ↦ z·ᾱ/D with N(α) = D)".
* Credit "Madore's technique, 2015; THM-4558's Lemma R" should add decalion89 Prop. A and hn-2adic Thm A′ (the d ≡ 1,2 mod 4 bipartite case also follows from hn-2adic Cor. B, since Q(√−d) ⊂ Q(√d)(i)).

---

## 7. Results note §§3, 4, 7 — **PASS WITH CORRECTIONS**

§3: as THM-4561 (F-range error; "√2-Golay in the Rudin–Shapiro signing" → "Golay-flat, sup ≤ 8.24√N"). Barker/Legendre claim correct (p = 5, 13 excluded by the count of minus signs, p > 13 by Turyn–Storer). Cyclotomic partition: failures in the `.out` all inside the bound; q ≥ 7 for s = 2 correct. "Turyn's conjecture" is #076's own terminology (citing Downarowicz–Lacroix). "In the plane, translation tiles admit periodic tilings (Bhattacharya; Beauquier–Nivat)": specify Z² (Bhattacharya) / polyominoes and discs (Beauquier–Nivat).

§4: all THM-4558/4559/MISTAKE-576 items apply. Additionally:
* §4.2 *"Our formulation works with an arbitrary conjugation-stable L instead of K(i). That is what lets primes over 2 (where x² + y² degenerates) and fields not containing i be handled uniformly."* — both are in decalion89 (Prop. A, conjugation-stable K ⊂ C; Theorems 1–2 for Q(√−3,√−11), Q(√−3,√−11,√−247)) and hn-2adic (Thm A at 2). Replace by the credits of §1.
* §4.2 *"Gibbs (Polymath16) 4-coloured the ring Z[ω₁, ω₃] by the same F_4 homomorphism."* — the F_4 reduction is credited to Speyer (decalion89); Gibbs per Goucher used "a four-element group". Rephrase.
* §4.3 *"So #158's sixth colour is necessarily topological"* — only density is excluded; write "cannot come from the density bound alone". *"J₀ enters #158 only through its decay"* — mark as the reader's (unverified) reading.
* §4.4 *"The compactness graph H has real-algebraic coordinates (Tarski)"* → "can be realized with real-algebraic coordinates". *"κ(17) ∈ [5, 7] and κ(29), κ(41) ≥ 5 (Hoffman)"* → "κ(17) ∈ [5, 6] (audit A; decalion89), κ(29), κ(41) ≥ 6 (decalion89, three-point bound)". First c-stable primes 17 / 29 / 41 for the three fields: verified. "Smallest" is right only in the sense of the smallest added t (over Z[ω₁,ω₃]: first non-excluded addition ω₂, next ω₁₀; over Z[ω₁,ω₃,ω₄]: ω₂, ω₅, then ω₉, ω₁₀, ω₁₁, ω₁₃); by inclusion Frac Z[ω₁,ω₂,ω₃,ω₄] ⊃ Frac Z[ω₁,ω₂,ω₃]. The whole paragraph reproduces decalion89 `local_colourings.md` §5 (same fields, same primes 17, 41, 83, 101 / 41, 101, 131), which also already exhibits a 5-chromatic graph in Q(√−3,√−7,√−11) — credit it.

§7.1: mathematics correct (checked: every primitive Pythagorean triple with m < 60 has 7 ∤ c and residue on the 8-point circle; (1+2i)·circle = knight set; all 8 residues occur, e.g. (3/5,4/5) ↦ (2,−2), so the map is onto). **Wording error:** *"So **`z ↦ (1 + 2i) z mod 7` maps every rational unit-distance pair to a knight move**"* — z mod 7 is undefined when 7 divides a denominator (z = 1/7). Corrected: "on each component z₀ + Γ (Γ = group generated by rational unit vectors ⊂ Z_(7)[i]), z ↦ (1+2i)(z − z₀) mod 7 maps every rational unit-distance pair to a knight move".
§7.2: untyped; "reduction modulo a conjugation-stable prime, a ring projection onto a finite field" is a residue map, not an idempotent — type the subsection ANALOGY.
§7.3: *"Isbell's 7-colouring uses the *split* prime `2 − ω` of norm 7"* (committed text) — with the note's ω = e^{iπ/3} (needed for "norm-7 directions `2 + ω, 3 − ω`", both checked: N = 7, 57 pairs at distance √7 reproduced), N(2 − ω) = 3. The uncommitted working copy (modified 00:09, during this audit) already reads "`2 + ω` of norm 7 (`ω = e^(iπ/3)`; `2 − ω` in the convention `ω = e^(2πi/3)`)", which is correct. Remaining: "Minkowski product" → "Minkowski sum (Erdős product)".

(Version note: THM-4558/4559/4561/4562 and THM-418 were audited at HEAD 9258bf9d52, unchanged since. MISTAKES.md and the results note had uncommitted edits at 00:08–00:09; the only edit touching the audited text is the §7.3 Isbell line. MISTAKE-576 and §§3, 4, 7.1, 7.2 are as quoted.)

---

## Required corrections

1. THM-4558 l.81, results note §4.3 bullet 3, MISTAKE-576: delete "The 5-chromatic step needs the class-number-2 field Q(√−15)" / "The 5-chromatic graphs need …"; replace with: 5-chromatic graphs exist in Q(√−3,√−11,√−15) (Heule; h = 2) and in Q(√−3,√−11,√−247) (Exoo–Ismailescu, DCG 2020, arXiv:1805.00157; h = 6; χ = 5 exactly by Lemma R at 11).
2. THM-4558 Corollary 1: replace "A 5-chromatic graph needs a rotation ω_t with t even, such as Heule's ω₄" by "no 5-chromatic graph has all edge vectors in this field; known ones leave it via ω₄, ω₁₆ or (119+3√−247)/128".
3. THM-4558 status/related, results note §4.2 Credit, THM-418 update: credit decalion89 `notes/local_colourings.md` (Prop. A = Lemma R; Thms 1–2; χ(G_q) table; Prop. B; §5 field screen; 2026-09-24..29) and MildlyMeticulous hn-2adic-obstruction (Thm A/A′, Cor. B/B′, Thm C/C1; 2026-07); replace "The reduction technique is Madore's (2015), generalized … to any conjugation-stable field L" by an attribution to Woodall 1973, Fischer 1990 (Thm 1), Moorhouse 2010 (Lemma 4.2), Madore 2015 (Prop. 3.2, ¶6.6), hn-2adic and decalion89; retype Lemma R, the M and M(i) rows, Corollary 1 and the Heegner-compositum bound as KNOWN; list what is new (Heegner rungs, N = 3n family, Polymath-field value, Corollary 2).
4. THM-4558 / note §4.2: state that Fischer 1994 ("A planar geometric graph of chromatic number four", Congr. Numer. 104, 73–79) is known only via the zbMATH review quoted by decalion89.
5. THM-4558 κ table: header provenance (q = 13, 17 not double-checked by the session; now checked by audit A); 17: "5–7" → "5–6"; note decalion89's χ(G₁₃) = 6; retype Hoffman bounds NUMERICAL (or cite decalion89 Prop. B: ≥ 6 for all q ≥ 53, and q = 29, 37, 41, 43, 47).
6. THM-4558 l.70: "two Eisenstein rhombi" → rigid Eisenstein patch through a point of √−3·Z[ω] at squared distance N (N = 3, 9, 21, 27 verified).
7. THM-4558 ll.68, 72: "the prime over 2/11" → "each prime over 2/11 (g = 2)"; l.72: "Heule's 5-chromatic graphs" → "several of Heule's 5-chromatic graphs (Goucher 2018; Polymath16)".
8. THM-4558 ll.85–86 and note §4.4: "lives in some number field"/"has real-algebraic coordinates" → "has a realization in"/"can be realized with".
9. THM-4558 l.17 and note §4.2: Moser-ring 4-colouring credit → Gibbs (four-element group, per Goucher 2018) and Speyer (reduction mod 2, Polymath16 thread 2); drop "by the same F_4 homomorphism" unless sourced.
10. THM-4559: "At most 5 concyclic points are Ramsey" → "Every set of at most 5 concyclic points is Ramsey"; extend the not-found statement to Pálvölgyi's full paper.
11. THM-4561 and note §3: F-range "∈ [0.75, 3.0] (k ≤ 12)" → "[0.70, 3.20] for 3 ≤ k ≤ 12; [0.753, 3.001] at k = 12"; "√2-flat Golay mixtures"/"√2-Golay" → "sums of √2-flat Golay pieces, sup ≤ 8.24√N"; type the column "does not flatten" remark NUMERICAL; add the sketch qualifier to the title.
12. THM-4562 Lemma 1 (and note §6.3, reader report): add "F_j algebraically independent over R (C[x] = R^{[1]})".
13. THM-4562: upgrade fibre counts to PROVED for all q with gcd(q,6) = 1 and the fibre geometry to PROVED; add the Dutta–Lahiri reference (JPAA 225 (2021) 106707, Thm 3.4) and definition.
14. MISTAKE-576: quote HYP-2278 (1) in full and say exactly what is wrong (integrality / χ_m vs χ_f; fractional part true; "no analytic method" refuted only CONDITIONALLY on #158); add the S688 HYP-2277 "exactly the class-number-one fields" quote and its refutation.
15. THM-418 update: "(contains both lattices at every scale and rotation)" → "contains a similar copy of U(Z², D) and U(Z[ω], D) for every norm D".
16. Note §4.3: "necessarily topological" → "cannot come from the density bound alone"; flag the J₀ remark as unverified.
17. Note §4.4: κ(17) → [5,6], κ(29), κ(41) → ≥ 6 (decalion89); qualify "smallest" (by added t); credit decalion89 §5 and its 5-chromatic graph in Q(√−3,√−7,√−11).
18. Note §7.1: per-coset base point in the bolded map statement.
19. Note §7.2: type ANALOGY; residue map is not an idempotent.
20. Note §7.3: "2 − ω of norm 7" → "2 + ω of norm 7 (conjugate 3 − ω)" (already applied in the uncommitted working copy); "Minkowski product" → "Minkowski sum".
21. THM-4558 l.79 (minor): "the compositum of all six" → "the compositum of the six with d ≡ 3 (mod 8)".

---

### Appendix: the k = 5 SAT runs (the "see end" of the header table)

* q = 13, k = 5: CaDiCaL 1.9.5, 1800 s, no answer (`auditA_kappa13_k5.out`). A parallel kissat run and a q = 17, k = 5 CaDiCaL run were stopped unfinished at about 1000 s each. This matches decalion89's report that CaDiCaL ran 25 min on G₁₃ without an answer; their χ(G₁₃) = 6 comes from a completed cube-and-conquer run, which audit A did not re-check.
* So audit A's own certified values are κ(13) ∈ [5,6] and κ(17) ∈ [5,6]: no 4-colouring (UNSAT), and an explicit, verified 6-colouring for both.
