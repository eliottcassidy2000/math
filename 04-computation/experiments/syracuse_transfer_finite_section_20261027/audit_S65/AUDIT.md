# Independent audit: finite-section transfer operators of the Syracuse map (opus-2026-10-27-S65)

**Audited (worktree `math-wt-s65-syracuse-transfer`, commit `cac59265c5`):**
- `05-knowledge/results/syracuse_transfer_finite_section_20261027.md` (the note);
- `01-canon/theorems/THM-4821-...md`, `01-canon/theorems/THM-4822-...md`;
- `04-computation/experiments/syracuse_transfer_finite_section_20261027.py` (+ `.out`);
- the related canon, for prior art: THM-4629, THM-4631, THM-4520, THM-4519, and the TCPC note §5.

**Method.** I did not use the author's code. All my code is in this folder:
- Python with python-flint 0.9.0 (exact integer matrices, exact minimal polynomials, modular ranks), numpy, scipy, mpmath, sympy.
- Every computed claim was re-derived by a different route.
- Negative controls (`neg_controls.py`) show that the checks fail on a wrong class map, a perturbed matrix, and a perturbed kernel.
- Nothing in the repository was edited, committed or pushed. `git status` in the worktree is clean.

**Overall verdict: SOUND WITH CORRECTIONS.**
- Both theorems are true as stated. Every proof step checks.
- Every number I re-computed agrees.
- There is one mathematical error, confined to the note: the displayed formula for the rows of Π in Theorem B(1) is mis-indexed and is *not* the law of Syrac. THM-4822 itself does not contain it.
- The other corrections are about prior art, citations and scope:
  - the J = 2, K = 3 and de Bruijn cases are already in the repository;
  - the forgetting statement B(1) is Tao's (1.21);
  - the universality in p is an isomorphism inherited from the classical parity-vector conjugacy;
  - THM-4629 2(c) should be 2(d) for the T-map Ulam matrix;
  - the journal year of Tao's paper is 2022;
  - the Jordan type of Theorem B counts transient states.

## Summary

| # | item | verdict | one-line reason |
|---|---|---|---|
| 1 | Definition of P_K | **SOUND** (one reproducibility correction) | Haar rows from lift counting equal the closed-form rows (K = 2..7, 7 multipliers); doubly stochastic; nnz = K 2^(K−1); a* row needs the Haar limit, as stated. The note's own lift-count check is not in the committed script. |
| 2 | Theorem A (THM-4821) 1–6 | **SOUND WITH CORRECTION** (typing/citation only) | Exact dense checks for K ≤ 11 and class-based checks for K ≤ 22 confirm everything, and all proof steps hold. Statement 6 is true and stronger: P_K^(p) is permutation-conjugate to P_K^(3) (classical parity-vector conjugacy), so "universality in p" is not new. |
| 3 | Theorem B (THM-4822) | **SOUND WITH CORRECTION** | Ranks, P^J = Π, Jordan type and Tao's mod-9 law are confirmed exactly; circulant eigenvalues and THM-4520's table (J = 2..7, from exact laws at 30 digits) also agree. The note's displayed law of Π is wrong. B(1) is Tao's (1.21). The Jordan type includes the transient states; on supp π it differs (index still J). |
| 4 | T-map contrasts | **SOUND WITH CORRECTION** (citation/wording) | (a) Characteristic polynomial exact for J ≤ 5 (and p = 5, 7, 17, 31). "THM-4631 at u = v = 1" should be "on the line u = v = z". (b) P^K = Π ≠ P^(K−1) exact for K ≤ 11. The citation should be THM-4629 2(d), the index K is classical de Bruijn, and the case is in-repo (procgen_repunit 2026-09-24 §2.1). |
| 5 | Floating-point claims | **SOUND** | numpy \|λ₂\| reproduced and shown to depend on the ordering. ℓ²(π) norms 0.674200, 0.971291, 0.996151, 0.999518, 0.999946 confirmed (J = 2, 3 to 15 digits with mpmath). The ε^(1/m) explanation is a fair order-of-magnitude heuristic. |
| 6 | Prior art and novelty | **SOUND WITH CORRECTION** | The "new here" items are new to the repository, but the note omits in-repo prior art for the small cases and the classical source of the p-universality. Tao's journal year is wrong. Lagarias is stated correctly in substance. |
| 7 | Overall answer | **SOUND** (scope clarifications suggested) | SLEM = 0 for every K, J; eigenvalue gap 1; singular gap ≈ 2^−(K−1); exact mixing time linear in K. No overclaim about Collatz. "Cannot see p = 3" is correct, and in the 2-adic case it is literally an isomorphism. |
| 8 | Reproduction | **SOUND** | The rerun is identical to the committed `.out` after normalising CRLF and the timing fields (74 s). The committed `.out` has mixed CRLF/LF line endings (cosmetic). |

## Details

### Item 1. The definition of P_K (`a1_definition.py`, `.out`)

**The Haar row from the definition.**
- Method: for x = a + 2^K t, I counted S_p(x) mod 2^K over t mod 2^M, M = K + 8. I kept only the lifts whose image mod 2^K is determined (v₂(px+1) ≤ M) and recorded the undetermined mass.
- Compared with the note's rule for K = 2..7 and p ∈ {3, 5, 7, 1, −1, 9, 15}:
  - every non-special row is equal exactly;
  - for a* the determined lifts are exactly uniform;
  - the undetermined mass is exactly 2^−(M−K+1), so it tends to 0.

**Integer representatives.**
- The note says C_K(a, b) = #{t mod 2^(K−1) : S(a + 2^K t) ≡ b}. This holds exactly on every row except a*.
- On a* with integer representatives the counts are not uniform. For K = 4, p = 3, a* = 5 the counts are {1: 3, 3: 1, 5: 1, 7: 1, 11: 1, 13: 1}.
- So the note's choice of the Haar limit for the a* row is necessary, and its warning about integer-representative matrices is right.

**Structure.**
- Doubly stochastic: S_p preserves Haar measure. The preimage of b + 2^K Z₂ has measure Σ_(v≥1) 2^(−v−K) = 2^−K.
- nnz = Σ_(v=1..K−1) 2^(K−1−v)·2^v + 2^(K−1) = K 2^(K−1). Checked for every case above, and in `a3` for K ≤ 18.
- The a* values agree with the `.out`: 85, 341, 341, 1365, ... at K = 8, 9, 10, 11, ...
- C_K is integral, with entries 2^(K−1−v) and 1.

**Reproducibility.** The note (§1 Checks) and the THM-4821 status cite "direct lift counts for K = 3..7 (scratch check)". This check is not in the committed script; `grep lift` finds only the de Bruijn docstring. It is therefore not reproducible from the commit. My a1 reproduces it independently.

### Item 2. Theorem A (THM-4821)

**Exact dense computations (`a2_dense_exact.py`, `.out`; no class machinery).**
- Setup: the integer C = 2^(K−1) P_K built from scratch, and J the all-ones matrix, so 2^(K−1)(P − Π) = C − J.
- p = 3, K = 2..11:
  - exact minimal polynomial minpoly(C) = (x − 2^(K−1)) x^(K−1) for K ≤ 10 (K = 10 took 36 s);
  - C^(K−1) = 2^((K−1)(K−2)) J, and 2^(K−1) C^(K−2) ≠ 2^((K−1)(K−2)) J;
  - (C − J)^(K−1) = 0 ≠ (C − J)^(K−2), computed directly without using PΠ = Π.
- p = 5, 7, 1, −1, K = 2..9: the same, all exact.
- Ranks:
  - rank P^n = 2^(K−1−n) and rank N^n = 2^(K−1−n) − 1, exactly for all n;
  - lower bound: rank modulo 2^61 − 1 and 2^62 − 57;
  - upper bound: an exact check that every row of C^n is constant on the classes mod 2^(K−n);
  - rows of N^n lie in rowspace(P^n) ∩ {sum = 0}.
- Jordan type of N on the whole space: the note's type on W plus one 1×1 block from the constants. For example K = 10: {1: 129, 2: 64, 3: 32, 4: 16, 5: 8, 6: 4, 7: 2, 8: 1, 9: 1}.
- Rows of C^n:
  - each row's support is a full residue class mod 2^m, with constant values;
  - each class of precision m ≥ 2 has a level-only multiplicity;
  - the number of distinct classes at level m is 2^(m−1).
- Hitting-time counts 2 C(K−2, h−1); mean K/2; worst TV 1 − 2^−(K−1−n).
- Dobrushin coefficients 1, …, 1, 0, by brute force over all pairs of rows for K ≤ 8.
- **Gram identity** (P^n)ᵀ P^n = Σ_m w_m(n) E_m, where E_m is the conditional expectation on classes mod 2^m. It holds exactly, entry by entry, for K ≤ 8 (p = 3) and K ≤ 7 (other p). Example K = 8: ‖(P − Π)^n‖² = 1, 63/64, 57/64, 21/32, 11/32, 7/64, 1/64, 0 = P(Bin(6, ½) ≥ n).

**Independent class-based code (`a3_class_method.py`, `.out`; my own Φ and v₂).**
- The class lemma, checked directly at every level:
  - the sum of the rows of C_K over a class C equals 2^(K−1−m+m′) 1_Φ(C), using scipy sparse aggregation;
  - p = 3 for K = 2..18; 13 multipliers (1, −1, 5, 7, 9, 11, 13, 15, 17, 33, 1023, −3, −5) for K ≤ 12.
- Trajectory checks: p = 3 for K = 2..22, and the 13 multipliers for K = 2..16.
  - R_n is exactly the set of all classes of precision 1..K−n.
  - The multiplicity of each class of precision m ≥ 2 in R_n is **binom(K − m − 1, n − 1)**. This is my closed form: choose the n − 1 intermediate precisions in (m, K). It makes the note's "level-only" statement explicit.
  - Hitting counts, mean K/2, worst TV, Dobrushin, and the ℓ² values Σ_(m≥2) w_m(n) = P(Bin(K−2, ½) ≥ n), all exact.
- K = 12 reproduced exactly:
  - hitting counts 2, 20, 90, 240, 420, 504, 420, 240, 90, 20, 2;
  - average TV profile 0.9995, 0.9941, 0.9673, 0.8867, 0.7256, 0.5, 0.2744, 0.1133, 0.0327, 0.0059, 0.0005, 0.

**Proof review (note §2).** Every step is sound.
- *(1)* x ↦ px + 1 is a Haar-preserving bijection of classes.
  - If v < m the valuation is constant on the class and division by 2^v is uniform.
  - If v ≥ m, the odd part of a Haar element of 2^m Z₂ is Haar on Z₂^×.
  - The author's merge-identity check (the "same | kids" test) is a valid induction up the tree, given the finest level.
- *(2)* c ↦ pc + 1 is a bijection from the odd to the even residues mod 2^m. The residues of valuation exactly v (1 ≤ v ≤ m−1) number 2^(m−1−v) and divide onto all odd residues mod 2^(m−v). Hence Φ hits each class of precision 2..m−1 once and the root twice.
  - The induction for multiplicities is mult_(n+1)(m′) = Σ_(m=m′+1..K−n) mult_n(m) for m′ ≥ 2. It is level-only, with solution binom(K−m−1, n−1).
  - The span claim is immediate.
- *(3)* ΠP = PΠ = Π (doubly stochastic), so P^n = Π + N^n and rank N^n|_W = rank P^n − 1. The block counts follow from consecutive rank differences. The dimension check Σ s 2^(K−2−s) + (K−1) = 2^(K−1) − 1 holds.
- *(4)* The subset count is right: each set T ⊆ [2, K−1] of intermediate precisions comes from exactly two residues (backward: unique preimage at each level, two at the root).
  - In parity-vector coordinates (below), h − 1 is the number of 1s among the parities ε₁..ε_(K−2), so h − 1 ~ Bin(K−2, ½) directly.
  - TV(μ_C, π) = 1 − 2^(1−m). Rows are disjoint iff two incomparable classes occur, iff K − n ≥ 2.
- *(5)* The conditional-expectation argument is right, because ‖P^n f‖² = Σ_m w_m(n)‖E_m f‖² holds for every f. I verified the stronger exact operator identity (P^n)ᵀP^n = Σ_m w_m(n) E_m.
  - The maximiser f = ±1 on 1/3 mod 4 attains every ‖E_m f‖ = ‖f‖.
  - Because the E_m are nested projections, the identity gives the **full singular spectrum**: 1 once, and √P(Bin(K − m, ½) ≥ n) with multiplicity 2^(m−2) for m = 2..K.
  - Checked against a float SVD of P^n for K = 3..11 and all n, with maximum deviation 3.9e−14 (`e3_full_svd_powers.out`).
  - n = 1 gives σ₂ = √(1 − 2^−(K−2)) as claimed (`e1_float.out` agrees to 1e−10 for K ≤ 12).
- *(6) Universality, which is stronger than stated and not new* (`a4_conjugacy.py`).
  - Let Q_p be the parity-vector map of T_p: x ↦ x/2, (px+1)/2. By Terras's argument it is a bijection Z/2^K → {0,1}^K for every odd p, compatible across K.
  - Then ψ = Q_3^(−1)∘Q_p is a tree automorphism of Z/2^K, and P_K^(p)(a, b) = P_K^(3)(ψa, ψb) exactly.
  - Reason: S_p is the first return to the odd residues of T_p, and in parity coordinates it is "shift past the next 1" for every p.
  - Verified for K = 2..13 and p ∈ {5, 7, 1, −1, 9, 15, −3, 1023}. The T-map de Bruijn matrices are conjugated by the same ψ.
  - So the finite sections for different p are *the same Markov chain up to relabelling*. This is the finite-section shadow of the classical fact (Terras 1976; Lagarias 1985 for p = 3, with the same proof for every odd p) that all T_p are conjugate to one Bernoulli shift under Haar measure.

**Table of §2.** Every entry checks:
- nnz 1,024 / 24,576 / 524,288 / 10,485,760;
- index K − 1;
- block counts 32, 16, 8, 4 / … / 131072, 65536, 32768, 16384;
- 63/64 … 262143/262144;
- singular gaps 7.84e−3, 4.88e−4, 3.05e−5, 1.91e−6;
- numpy 0.0029 and 0.0179.

**A precision note on "binomial cutoff at n = K/2".**
- Exactly, the average TV is P(B > n) + ½ P(B = n) with B ~ Bin(K−2, ½). I checked this exactly for K = 2..18 (`a5_avgtv_circ.out`).
- It crosses ½ at n = (K − 2)/2, not K/2. The mean of h is K/2.
- The difference is O(1) against a √K window, so this is wording, not an error.

### Item 3. Theorem B (THM-4822) (`b1_tao_chain.py`, `b2_fourier.py`, `b3_tao_J6_exact.py`, `b4_recurrent_scope.py` and their `.out` files)

**Exact matrices.** M = (2^L − 1) P was built by summing the weights 2^(L−a), a = 1..L, along y = (px+1)2^−a, so each entry is 2^(L−a₀)/(2^L − 1) as the note says. Cases (all exact with flint):
- p = 3, J = 1..5, plus J = 6 with all 729 starts (37 s);
- p = 5, J ≤ 3;
- p = 7, J ≤ 3 (2 not primitive, cosets of size 3, 21, 147);
- p = 9, 11, 13, 15, 17, 21, 23, 31 at J = 2 (index of ⟨2⟩ in the units mod p^2: 1, 1, 1, 2, 2, 6, 2, 6).

In every case:
- rows of P^n depend only on x mod p^(J−n);
- P^J = Π and P^(J−1) ≠ Π;
- rank P^n = p^(J−n), from a modular lower bound equal to the exact distinct-row count;
- rank N^n = p^(J−n) − 1;
- the Jordan type on Z/p^J is the note's type plus one 1×1 block, e.g. p = 3, J = 5: {1: 109, 2: 36, 3: 12, 4: 4, 5: 2};
- minpoly(M) = (x − (2^L − 1)) x^J exactly for p^J ≤ 125.

**Stationary law.**
- For p = 3, J ≤ 4, three routes agree exactly: the row of P^J, the nullspace of Pᵀ − I, and Tao's sum Σ_m 3^(m−1) 2^−(a_1+…+a_m) by exact DP. The same holds for p = 5, 7 at J = 2.
- The mod-9 law is (0, 8, 16, 0, 11, 4, 0, 2, 22)/63 exactly.
- By the recursion Syrac_J = (3 Syrac_(J−1) + 1)2^−a (Tao's Lemma 1.12) for J = 1..7:
  - consistency mod 3^(J−1) holds exactly;
  - the max primitive |Fourier coefficient| at J = 2..7 is 0.37792363, 0.25223684, 0.17699888, 0.12927357, 0.096106391, 0.075870038 (mpmath), attained at ξ = 4, 8, 65, 211, 64, 256;
  - these match THM-4520 (5)'s table at all six J, including J = 7 (0.07587), which the note did not test;
  - TV(Syrac, uniform on units) = 0.16667, 0.27778, 0.34601, 0.37054, 0.38887, 0.39855, 0.40686, as in the `.out` up to J = 6.

**The error in the note.**
- Theorem B(1) says "The rows of Π are the law of Σ_(j=1..J) p^(j−1) 2^-(a_j + ... + a_J)". That is wrong for J ≥ 2.
- The closed form displayed one line earlier is correct. It gives x_J = Σ_j p^(j−1) 2^−(a_(J−j+1)+…+a_J), whose law (reversing the i.i.d. a's) is Σ_j p^(j−1) 2^−(a_1+…+a_j). That is Tao's (1.29).
- The note's sum puts p^(j−1) on the last J−j+1 variables instead of the last j. At p = 3, J = 2:
  - its mod-3 marginal is (0, 5/9, 4/9), against Syrac's (0, 1/3, 2/3);
  - its mod-9 law is (0, 188, 208, 0, 395, 262, 0, 152, 118)/1323.
- It also fails at J = 3, 4 and at p = 5, 7. THM-4822 does not contain this formula. The script computes the law by iterating the chain, so the computed laws are right.

**Proof review.**
- (1) The closed form x_n = p^n x₀ 2^−A_n + Σ_(j=1..n) p^(j−1) 2^−(a_(n−j+1)+…+a_n) is correct (induction). It is Tao's (1.7)/(1.21)/(1.24); see item 6.
- (2) The depth-filtration argument is correct:
  - U F_d ⊆ F_(d−1);
  - row x̄ lives on the coset (px+1)⟨2⟩ of the units mod p^d, since px + 1 ≡ 1 mod p is a unit;
  - rows on distinct cosets are disjoint, and on one coset they are distinct cyclic shifts, since x̄ ↦ px + 1 is injective;
  - any set of distinct rows of an invertible circulant is independent.
  - Circulant convention: with weight w(j) at cyclic offset j from px+1 (the note's convention), the eigenvalues are Σ_j w(j) ζ^j = ζ/(2 − ζ) = 1/(2ζ^(−1) − 1). As ζ runs over the L-th roots of unity this is {1/(2η − 1)}, on |z − 1/3| = 2/3, with min |z| = 1/3.
  - Checked numerically for L = 2, 6, 18, 54. |det Circ| = (2^L − 1)^(L−1) ≠ 0 exactly for L = 2, 3, 4, 6, 18, 21, 54.
  - Putting w(1) at offset 0 would give 1/(2 − ζ), the circle |z − 2/3| = 1/3. That is also nonzero, so invertibility does not depend on the convention, and the note's convention agrees with THM-4520 (2).
  - Primality of p is not used. Composite p = 9, 15, 21 pass.
- (3) The ranks give the Jordan type. Its dimension total is p^J − 1, and the index is J because rank N^(J−1) = p − 1 > 0.
- *Every start needs exactly J steps*: the three rows of P^(J−1) are independent and π = ⅓ row₁ + ⅔ row₂, so no row equals π.

**Scope (not an error).**
- The Jordan type is for the chain on all of Z/p^J, which includes p^J − r p^(J−1) transient states (r = |⟨2 mod p⟩|).
- On the recurrent class R = supp π = {y : y mod p ∈ ⟨2⟩} (for p = 3, the units), the Jordan type on the π-mean-zero space is:
  - r − 1 blocks of size J;
  - r(p − 2) + 1 blocks of size J − 1;
  - r(p − 1)² p^(J−2−s) blocks of size s ≤ J − 2.
- Verified exactly for 12 cases, including p = 7, 17, 31 and p = 9, 15.
- For Tao's chain (p = 3) this is 1, 3 and 8·3^(J−2−s), against 2 and 4·3^(J−1−s) on Z/3^J. The index is J in both.
- The TCPC note's "rank 2" at mod 9 is this units restriction, and is consistent.

**Checks of the note.**
- "Equal rows of P^J from every start for J ≤ 4, five starts for J = 5, 6" is accurate. The note's Status bullet ("verified exactly … J = 1..6") overstates it slightly.
- The audit closes the gap: J = 5 and J = 6 are exact from every start.

### Item 4. T-map contrasts (`c1_tmap.py`, `e2_debruijn_svd.py`)

**(a) T-map random chain on Z/p^J.**
- I computed the exact characteristic polynomial of the integer matrix 2P with flint. It equals (μ − 2) Π_d (μ^(L_d) − 1)^(c_d), with c_d = φ(p^d)/L_d.
  - p = 3, J = 1..5 (the note checked J ≤ 3);
  - also p = 5 (J ≤ 3), 7 (J ≤ 3, c_d = 2), 17 (J ≤ 2), 31 (J ≤ 2, c_d = 6).
- With μ = 2λ this is the note's (λ − 1) Π_d (λ^(L_d) − 2^−(L_d)). So λ = 1 is simple, every other eigenvalue has modulus ½, and the chain is never exactly mixed. SLEM = ½ exactly for every J.
- The proof sketch is right: halving permutes the primitive characters of depth d (one cycle of length L_d for p = 3), and the odd branch lowers depth.
- Identity with THM-4631: at j = 0 and u = v = 1 the operator M is the Koopman operator U, but det(I − M) at u = v = 1 is just 0. The characteristic polynomial comes from the line u = v = z: det(I − zU) = (1 − z) Π_d (1 − (z/2)^(L_d))^(c_d), hence det(λ − U) = λ^(p^J) det(I − U/λ). THM-4631's identity is polynomial in independent (u, v), so this specialisation is legitimate. "At u = v = 1" is loose wording.

**(b) T-map Ulam (de Bruijn) matrix on Z/2^K, K = 2..11, exact.**
- (2P)^K = J, i.e. P^K = Π, and P^(K−1) ≠ Π.
- rank P^n = 2^(K−n), with an exact upper bound from class-constancy.
- Jordan type on the mean-zero space: one block of size K and 2^(K−1−s) blocks of size s.
- σ₂ = 1 exactly for every K ≥ 2: P Pᵀ is a sum of 2×2 averaging blocks, so the singular values are 1 and 0, each with multiplicity 2^(K−1). The note's table leaves this cell as "—"; it can be filled.
- Attribution: THM-4629 2(c) is the *halving-branch* averaging E_j (single nonzero eigenvalue ½, rank E_j^q = 2^(j−q)), not this matrix.
  - The full T-map Ulam matrix is M|V₀ at u = v = 1, whose nonzero spectrum {(u+v)/2} = {1} is THM-4629 proof 2(d) (also THM-4631's V₀ paragraph).
  - The exact index K (P^K = Π) is not stated in THM-4629. It is classical: the binary de Bruijn graph has A^K = J, with exactly one walk of length K between any two words (Terras's parity bijection).
  - It is already in the repository: procgen_repunit_20260924 §2.1, "mod 8: the window chain of 3 i.i.d. bits, mixed exactly in 3 steps", and "for T the state mod 2^k is the window of the next k i.i.d. parity bits".
- Script docstring slip (line 395): "(2P)^K = 2^(K-1) J?" should be "(2P)^K = J". The code does not use it.

### Item 5. Floating point (`e1_float.py`, `.out`)

- numpy `eigvals` |λ₂| in the natural ordering: 0.0029, 0.0046, 0.0089, 0.0125, 0.0179 (K = 8..12), the same as the note.
  - Under three random relabellings: 0.0027/0.0023/0.0030 (K = 8) … 0.0215/0.0156/0.0224 (K = 12).
  - All 2^(K−1) − 1 zero eigenvalues come out with |λ| > 1e−6. This confirms they are pseudospectral artefacts.
- (1.1e−16)^(1/(K−1)) = 0.0053 … 0.0354. The observed values are about half of this, so the ε^(1/m) explanation is right in order of magnitude and growth.
- Tao's chain:
  - numpy |λ₂| = 0.0000, 0.0000, 0.0001, 0.0002, 0.0013 (J = 2..6);
  - ℓ²(π) norm of P − Π on the units = 0.674200, 0.971291, 0.996151, 0.999518, 0.999946;
  - mpmath at 40 digits: 0.674199862463242 (J = 2) and 0.971290824687078 (J = 3).
- The note's 0.674, 0.971, 0.996, 0.9995 and 0.99995 are correct roundings. The "not proved" label is appropriate.

### Item 6. Prior art and novelty

**In-repository prior art not cited by the note or the canon files.**
- `05-knowledge/results/procgen_repunit_20260924_digits_rotation_repunits.md` §2.1, Theorem H (PROVED), the exact Haar Markov chains of the Syracuse residues:
  - mod 9: the matrix, Tao's stationary law (8, 16, 11, 4, 2, 22)/63 on (1, 2, 4, 5, 7, 8), and "P^2 has identical rows, so the chain forgets its state in two steps". This is the J = 2 case of Theorem B.
  - mod 8: "rows 1 and 5 are uniform; 3 → {1,5} and 7 → {3,7}; P^2 is uniform". This is the K = 3 case of Theorem A: P_3^(K−1) = Π.
  - T mod 2^k: "the window of the next k i.i.d. parity bits … mixed exactly in 3 steps" at k = 3. This is the de Bruijn contrast of §4.
- The TCPC note §5.2. The note cites it only for Tao's law and the IFS, but it also states:
  - the mod-9 transition matrix (on the units) has rank 2 and P² = 1π exactly ("exactly mixed after two steps");
  - "A_(n+1) mod 3^k is a function of the last k valuations", which is Theorem B(1).
- `graph_theory_as_common_carrier_20261013.md` §5 and HYP-9261 (5): the T-map class graph mod 2^j is B(2, j); E_j is rank one plus nilpotent of index j; M|V₀ has nonzero spectrum (u+v)/2.
- THM-4631 already notes that p does not enter the 2-adic (de Bruijn) part: "the only place p enters is the odd branch … the reduction goes through for every odd p". That is the T-map analogue of the "universality in p".

**What remains new in the repository.** I found none of the following stated before:
- the exact Jordan type, hitting-time law, TV and ℓ² profiles of the forward Syracuse Ulam operator for general K;
- the general-J ranks and Jordan type of Tao's chain via Circ(w).

The note hedges "as far as the repository and these sources go", which is honest. The universality item should be credited to the classical conjugacy.

**Literature and citations.**
- *Tao.* I fetched arXiv:1909.03562 (v7) and read §1.3.
  - Journal: Forum Math. Pi **10 (2022)**, e12, so "Tao 2019, Forum Math. Pi 10, e12" mixes the arXiv year with the journal volume.
  - Lemma 1.12 is "Recursive formula for Syracuse random variables". Remark 1.13 identifies the law of Syrac(Z₃) as "the unique stationary measure for the discrete Markov process on Z₃ that maps each x to (3x+1)/2^a … with transition probability 2^-a". Both match the note's numbering.
  - The mod-9 law (0, 8/63, 16/63, 0, 11/63, 4/63, 0, 2/63, 22/63) is stated right after Lemma 1.12. Correct.
  - Tao's (1.21), Syr^n(N) = F_n(a^(n)(N)) mod 3^k for 0 ≤ k ≤ n, together with (1.22) and (1.24), *is* the forgetting statement B(1) and P^J = Π. It should be cited for B(1).
  - Tao's (1.29) is the correct reversed form of the law; compare the error above.
- *Lagarias 1985* (Amer. Math. Monthly 92, 3–23). The T-map on Z₂ is conjugate to the one-sided Bernoulli shift by the parity-vector map (measure-preserving). The note states this correctly in substance, without a theorem number, and types it as CITED from memory.
  - The same parity-vector argument (Terras 1976) works for every odd multiplier. That is the source of THM-4821's statement 6 (see item 2).
- *Ulam's method.* P_K is exactly Ulam's matrix P(a, b) = m(A_a ∩ S^(−1)A_b)/m(A_a). The usage is correct.
- *Not cited; possibly relevant, from memory and not re-read:*
  - Matthews and Watts (Acta Arith. 1984, 1985) on Markov chains of generalised Syracuse maps modulo m under the d-adic measure;
  - for the de Bruijn A^K = J, the classical theory of 0–1 matrices with A^k = J (Hoffman 1967; Lam and van Lint 1978).
  - My web search found no source for the exact Jordan types. Novelty beyond the repository remains unestablished, as the note says.

### Item 7. The overall answer

- Correct for the 2-adic carrier:
  - SLEM = 0 for every K (certified by (P − Π)^(K−1) = 0);
  - the eigenvalue gap is exactly 1, uniformly;
  - the singular gap 1 − √(1 − 2^−(K−2)) ≈ 2^−(K−1) decays;
  - the exact time to equilibrium is linear: mean K/2, and worst case K − 1 for every ε < ½, since the worst-start TV stays at least ½ until the last step.
- For Tao's chain: SLEM = 0, and the exact time is J from *every* start.
- The two readings — the eigenvalue gap overstates mixing, the singular gap understates it — are right. Fill's bound with σ₂ would suggest a number of steps of order K 2^K.
- No overclaim about Collatz: "None of this bears on the Collatz conjecture itself" and "NOT CLAIMED: any statement about integer orbits".
- The universality interpretation is correct and, in the 2-adic case, literally true: P_K^(3) and P_K^(5) are the same Markov chain up to a relabelling of residues, so no finite-section quantity at all can separate 3x + 1 from 5x + 1, spectral or not.
- Suggested scoping:
  - Say why. Under 2-adic Haar measure every T_p (p odd) is the same Bernoulli shift. The 3x+1 / 5x+1 difference enters the integer problem through the archimedean drift E log(p 2^−v) = log p − 2 log 2, negative for p = 3 and positive for p = 5, which no 2-adic Haar section can see.
  - For the p-adic chain, "behaves the same way" means "the same form, with p in the counts". The state spaces differ, and the stationary law (where the arithmetic lives) depends on p.

### Item 8. Reproduction

- `python syracuse_transfer_finite_section_20261027.py` finished with exit 0 in 74 s, printing ALL CHECKS PASSED.
- The output, after stripping CR characters and the `[..s]` and `elapsed .. s` timing fields, is byte-identical to the committed `.out` (104 lines; `committed_norm.txt` vs `rerun_norm.txt`).
- The float lines agree too, on this machine and LAPACK build.
- The committed `.out` has mixed CRLF/LF line terminators (cosmetic).

## Exact corrections (sentence → replacement)

**Errors.**

1. Note §3, Theorem B (1).
   - Current: "The rows of Π are the law of Σ_(j=1..J) p^(j−1) 2^-(a_j + ... + a_J): for p = 3, Syrac(Z/3^J)."
   - Replace with: "The rows of Π are the law of Σ_(j=1..J) p^(j−1) 2^-(a_(J−j+1) + ... + a_J), equivalently (reversing the i.i.d. a's) of Σ_(j=1..J) p^(j−1) 2^-(a_1 + ... + a_j): for p = 3, Syrac(Z/3^J) in the form of Tao's (1.29)."
   - Reason: the stated sum has mod-3 marginal (0, 5/9, 4/9) at J = 2, not (0, 1/3, 2/3).

2. Note §1.
   - Current: "(Tao 2019, Lemma 1.12 and Remark 1.13; CITED via the TCPC note)"
   - Replace with: "(Tao, Forum Math. Pi 10 (2022), e12, arXiv:1909.03562 (2019), Lemma 1.12 and Remark 1.13; CITED, numbering checked against arXiv v7)".

3. Note §6.
   - Current: "The Syracuse valuations are i.i.d. geometric (Tao 2019, Forum Math. Pi 10, e12, §1)."
   - Replace with: "The Syracuse valuations are i.i.d. geometric (Tao, Forum Math. Pi 10 (2022), e12 [arXiv 2019], §1.3, Proposition 1.9)."

4. THM-4822 status.
   - Current: "CITED: the chain's stationary law is Syrac(Z/3^J) (Tao 2019, Lemma 1.12 and Remark 1.13, as recorded in the TCPC note)."
   - Replace with: "CITED: the chain's stationary law is Syrac(Z/3^J) (Tao, Forum Math. Pi 10 (2022), e12, Lemma 1.12 and Remark 1.13); statement 1 is Tao's (1.21)/(1.24), and its J = 2 case is in the TCPC note §5.2 and procgen_repunit_20260924 §2.1."

5. Note §3 proof.
   - Current: "*Proof.* (1) is induction on n."
   - Replace with: "*Proof.* (1) is induction on n; it is Tao's (1.21)/(1.24) (Syr^n(N) mod 3^k depends only on the valuations for k ≤ n), and the case J = 2 (P² = Π mod 9) is in the TCPC note §5.2 and procgen_repunit_20260924 §2.1."

**Citation and attribution.**

6. Note §4 table.
   - Current: "| T-map Ulam (de Bruijn) | Z/2^K | 0 | yes | K | — (THM-4629 2(c)) |"
   - Replace with: "| T-map Ulam (de Bruijn) | Z/2^K | 0 | yes | K | 1 exactly (P Pᵀ is a sum of 2 × 2 averaging blocks); spectrum THM-4629 proof 2(d) at u = v = 1; index K classical (de Bruijn A^K = J) and procgen_repunit_20260924 §2.1 |".

7. THM-4821 status.
   - Current: "Canon THM-4629 proof 2(c) and THM-4631 prove the corresponding de Bruijn nilpotency for the inverse-tree (T-map) operator."
   - Replace with: "Canon THM-4629 proof 2(c)–(d) and THM-4631 prove the corresponding de Bruijn structure for the inverse-tree (T-map) operator (2(c): the halving-branch averaging E_j, rank E_j^q = 2^(j−q); 2(d): the full de Bruijn walk has the single nonzero eigenvalue (u + v)/2); the exact index K of the T-map Ulam matrix is classical (de Bruijn A^K = J) and is in procgen_repunit_20260924 §2.1."

8. Note §6, Canon bullet.
   - Current: "THM-4629 2(c) / THM-4631 already prove the rank-one-plus-nilpotent structure of the 2-adic de Bruijn averaging in the inverse-tree (T-map) operator, and THM-4631 finds the doubling-cycle modes of the 3-adic levels. The T-map rows of §4 are those results restated (checked again here)."
   - Replace with: "THM-4629 2(c)–(d) / THM-4631 already prove the rank-one-plus-nilpotent structure of the 2-adic de Bruijn averaging in the inverse-tree (T-map) operator, and THM-4631 finds the doubling-cycle modes of the 3-adic levels; the T-map rows of §4 are those results restated (checked again here). Small cases of both theorems are already in the repository: procgen_repunit_20260924_digits_rotation_repunits.md §2.1 (Theorem H) has the mod-8 Syracuse chain with P² uniform (K = 3 of Theorem A), the mod-9 chain with P² of identical rows and Tao's law (J = 2 of Theorem B), and the T-map window chain mod 2^k mixed exactly in k steps; the TCPC note §5.2 has the mod-9 rank 2 (on the units), P² = 1π and 'A_(n+1) mod 3^k is a function of the last k valuations' (Theorem B(1))."

9. "THM-4631 at u = v = 1". This phrase appears in the note's Status ("the doubling cycles of THM-4631 at u = v = 1"), the §4 table ("— (THM-4631 at u = v = 1)"), the THM-4822 status, and script section G.
   - Replace each with: "THM-4631 with j = 0 on the line u = v = z (det(I − zU) = (1 − z) Π_d (1 − (z/2)^(L_d)))".
   - Reason: at u = v = 1 the determinant is 0.

**Novelty and universality.**

10. THM-4821 statement 6 and note §2 item 6.
    - Current: "Statements 1–5 do not depend on p." / "None of 1–5 depends on p."
    - Replace with: "For every odd p, P_K^(p) = Q P_K^(3) Q^(−1) for the permutation Q of the odd residues induced by ψ = Q_3^(−1)∘Q_p, where Q_p is the parity-vector map of x ↦ x/2, (px + 1)/2 (a tree automorphism of Z/2^K; Terras 1976, Lagarias 1985); hence statements 1–5, and every relabelling-invariant property, do not depend on p."
    - Optionally replace the proof outline line with: "*(6)* S_p is the first return of T_p to the odd residues, and in parity coordinates it is the same shift for every odd p."

11. Note §6, "New here", third item.
    - Current: "the universality in p."
    - Replace with: "the universality in p, which is a corollary of the classical parity-vector conjugacy (all x ↦ x/2, (px+1)/2 with p odd are conjugate to one Bernoulli shift under Haar measure): the finite sections for different p are the same Markov chain up to relabelling."

12. INDEX.md entry.
    - Current: "Universality: identical for (px+1)/2^v and (px+1) 2^-a mod p^J for every odd p (checked p = 1, 5, 7, -1), so these untwisted finite sections cannot see p = 3."
    - Replace with: "Universality: for (px+1)/2^v the Ulam matrices of all odd p are permutation-conjugate (classical parity-vector conjugacy), and for (px+1) 2^-a mod p^J the statements have the same form with p in the counts, so these untwisted finite sections cannot see p = 3."

**Typing and reproducibility.**

13. Note Status.
    - Current: "It is verified exactly for K = 8..20 and J = 1..6, and modulo a prime for J = 7, 8."
    - Replace with: "It is verified exactly for K = 8..20 and for J = 1..4 from every start (J = 5, 6 exactly from five starts and from every start modulo a prime), and modulo a prime for J = 7, 8."
    - Alternatively keep "J = 1..6" and add "(J = 5, 6 from every start in the independent audit)".

14. Note §1 Checks.
    - Current: "The script checks the rows against direct lift counts for K = 3..7 (scratch check)."
    - Replace with: "The rows were checked against direct lift counts for K = 3..7 in an uncommitted scratch script (reproduced independently by the audit for K = 2..7 and seven multipliers)."
    - Matching edit in THM-4821 status: "the rows were checked against direct lift counts for K = 3..7." → "the rows were checked against direct lift counts for K = 3..7 (uncommitted scratch check; reproduced by the independent audit)."
    - Or add the lift-count check to the committed script.

**Scope clarifications (recommended additions, not errors).**

15. THM-4822 statement 3 and note §3 Theorem B (3). After "(p − 1)² p^(J−1−s) blocks of size s (1 ≤ s ≤ J − 1)", append:
    "This counts the chain on all of Z/p^J, including its transient states. On the recurrent class R = supp π = {y : y mod p ∈ ⟨2⟩} (r = |⟨2 mod p⟩|; for p = 3 the units), the Jordan type of N on the π-mean-zero space is r − 1 blocks of size J, r(p − 2) + 1 blocks of size J − 1 and r(p − 1)² p^(J−2−s) blocks of size s ≤ J − 2 (p = 3: 1, 3 and 8·3^(J−2−s)); the index is J in both."

16. Note §4, "The finite sections are blind to the multiplier". After "...behaves the same way for every odd p.", append:
    "In the 2-adic case the finite sections for different p are the same Markov chain up to relabelling: the finite-section shadow of the classical fact that every x ↦ x/2, (px + 1)/2 is the same Bernoulli shift under Haar measure. The multiplier enters the integer problem through the archimedean drift log p − 2 log 2 (negative for p = 3, positive for p = 5), which no 2-adic Haar section sees."

17. Note §2, after the K = 12 profile.
    - Current: "a binomial cutoff at n = K/2 with a window of order √K."
    - Replace with: "a binomial cutoff: the average total variation is exactly P(B > n) + P(B = n)/2 with B ~ Bin(K − 2, 1/2), equal to 1/2 at n = (K − 2)/2, with a window of order √K."

18. Note §5. After "returns |λ₂| = 0.0029, 0.0046, 0.0089, 0.0125, 0.0179", append:
    "(natural ordering; the values depend on the ordering and the LAPACK build: under random relabellings the audit saw 0.0023–0.0030 at K = 8 and 0.016–0.022 at K = 12)".

19. Script docstring of `section_G_tmap_2adic`.
    - Current: "Exact integers: (2P)^K = 2^(K-1) J? -> P^K = Pi, P^(K-1) != Pi (THM-4629 2(c))."
    - Replace with: "Exact integers: (2P)^K = J (all ones), i.e. P^K = Pi, and P^(K-1) != Pi (classical de Bruijn; spectrum THM-4629 proof 2(d))."
    - Matching edit in the print at line 597: "(THM-4629 2(c))" → "(THM-4629 2(d); de Bruijn A^K = J)".

**Optional strengthenings the note may add (verified).**
- The multiplicity of each class of precision m ≥ 2 among the rows of P^n is binom(K − m − 1, n − 1).
- The full singular spectrum of P_K^n: 1 once, and √P(Bin(K − m, 1/2) ≥ n) with multiplicity 2^(m−2), m = 2..K.
- σ₂ = 1 exactly for the de Bruijn T-map Ulam matrix.
- THM-4520 (5)'s J = 7 value 0.07587 is reproduced too.

## Auditor's files (all in this folder)

- `a1_definition.py` / `.out`: item 1, Haar lift counting vs closed-form rows (7 multipliers, K = 2..7).
- `a2_dense_exact.py` / `.out`: item 2, exact dense flint (minpoly, ranks, Jordan type, Gram identity, TV, Dobrushin), p ∈ {3, 5, 7, 1, −1}.
- `a3_class_method.py` / `.out`: item 2, class lemma by sparse aggregation and trajectories (K ≤ 22, 14 multipliers).
- `a4_conjugacy.py` / `.out`: items 2(6) and 7, permutation conjugacy P_K^(p) ≅ P_K^(3) (K ≤ 13, 8 multipliers).
- `a5_avgtv_circ.py` / `.out`: the average-TV closed form, and the circulant eigenvalue convention.
- `b1_tao_chain.py` / `.out`: item 3, exact chain for 20 (p, J) cases; three routes to the law; the wrong formula exhibited; circulant determinants. The "units-only" lines there restrict to x ≢ 0 mod p. That is the recurrent class only for prime p with 2 a primitive root mod p (p = 3, 5, 11, 13). The correct recurrent-class computation is `b4`.
- `b2_fourier.py` / `.out`: exact Syrac laws for J ≤ 7 and Fourier maxima vs THM-4520.
- `b3_tao_J6_exact.py` / `.out`: J = 6, all 729 starts exact.
- `b4_recurrent_scope.py` / `.out`: Jordan type on supp π (12 cases).
- `c1_tmap.py` / `.out`: item 4, exact characteristic polynomials (p = 3, 5, 7, 17, 31) and the de Bruijn checks (K ≤ 11).
- `e1_float.py`, `e2_debruijn_svd.py`, `e3_full_svd_powers.py` (+ `.out`): item 5 and singular spectra.
- `neg_controls.py` / `.out`: the checks fail on wrong inputs.
- `author_rerun.out`, `committed_norm.txt`, `rerun_norm.txt`: item 8.

Slips in my own tooling:
- `a2`'s first version used "distinct rows" as a rank upper bound, which overcounts when nested classes occur (K = 3: rows root, 1 mod 4, 3 mod 4 have rank 2). It was replaced by an exact class-constancy bound before the final run.
- `b4`'s first version double-subtracted the constants block. Fixed before the final run.
- `b1.out`'s printed "index of <2> in units" uses (p^J − p^(J−1))/L. That is right for prime p but wrong for composite p: the true indices are p = 15: 2 and p = 21: 6, not 3 and 10. Only that label is affected, not any check.
