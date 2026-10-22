<!-- Independent audit of the S55 pre-audit commit e212bd63c2, returned by the auditor agent on 2026-10-22 and saved verbatim by opus-2026-10-22-S55. The scripts it names are in this folder. -->

# Audit of opus-2026-10-22-S55: THM-4762, THM-4763, HYP-9372, HYP-9373 and the signature-axis note

The mathematics holds up: every numerical claim I recomputed with my own code reproduces. The two serious problems are both prior art. The decomposition at the heart of THM-4763 is already a theorem in the FGST paper it cites. HYP-9373's one concrete prediction is already a published theorem in FGST's companion paper. Beyond that there are four MEDIUM errors and a list of LOW fixes. I did not commit, push or edit any document. All scripts and outputs are in `C:\Users\Eliott\AppData\Local\Temp\math-wt-s55-signature-mtc\04-computation\experiments\catmap_modular_20261022\audit_S55\` (run with `python3` from that folder). I returned the report here rather than writing a REPORT.md.

## Verdict table

| Document | Verdict | Worst severity |
|---|---|---|
| THM-4762 | Parts (1)–(4) **CONFIRMED**. All 135 census rows reproduce, using my own builders and full modular-relation checks. **CORRECTED**: A_1(−1) is wrongly listed as having identical C; the Dong–Lin–Ng journal is wrong; π(N)/2 is not an integer for N ≤ 2 (5 of the 88 "saturating" cases); the minimal-model set is misdescribed; the rapidity is off by a factor 2. | MEDIUM |
| THM-4763 | Mathematics **CONFIRMED** analytically for p = 2..6. **PRIOR ART**: the decomposition is FGST 2006, eq. (1.2) and Theorem 2.3. The session's own Lucas-trace check was tautological. "R_1 is trivial" (p = 2) is wrong. | HIGH (attribution) |
| HYP-9372 | (a) **CONFIRMED**, and can be sharpened. (b) as titled ("exactly when") is **REFUTED**; the one-directional form is **CONFIRMED** in 128 of 132 cases. (c) **CONFIRMED** and essentially provable for Ising-type factors. | MEDIUM |
| HYP-9373 | Logic is sound and the general claim is rightly OPEN. **PRIOR ART**: the "predicted next case" (W(p,p′) with Sym² C² = C³) is FGST, Nucl. Phys. B 757 (2006), Thm 1.2(3). The HYP mis-cites it as "J. Math. Phys. 2006" and calls it UNVERIFIED. "PROVED for K = 0" is overstated. | HIGH (prior art) |
| Synthesis note | The 18 identities, c_eff values, q → 1 ratios, all claimed relations, the mock-theta identity, the Ising ring and the LTE claim are **CONFIRMED**. **CORRECTED**: the coefficient-growth claim is false for Rogers–Selberg and Bailey mod 9 (Bailey "no growth" is false), plus several LOW fixes. | MEDIUM |

## A. THM-4762 (1): CONFIRMED (`a_pisano.py`, `a_pisano.out`)

- **The identity.** Exact sympy computation gives t s t⁻¹ s⁻¹ = Q² = [[2,1],[1,1]]. The element st has trace 1 (order 6) and s has order 4.
- **ord(A mod N) = π(N)/2 for N ≥ 3.** Proof: det Qᵐ = (−1)ᵐ, so Qᵐ ≡ I forces m even when N ≥ 3. Checked for all N ≤ 2000 with no failure; π(N) is even for every 3 ≤ N ≤ 2000. For N = 2, ord(A mod 2) = 3 = π(2).
- **The specific periods.** No N ≤ 2000 has π(N) = 15. π(4) = 6, π(6) = 24, π(30) = 120, and 30 is the least N with π(N) = 120 (the next are 45, 60, 66, 80).
- **Extra result used in E.** For every 3 ≤ N ≤ 2000, the number of central powers of A in ⟨A⟩ is 1 or 2, never 4. Also π/α = 4 holds exactly when α is odd, with no violations. The counts of π/α are 1: 234, 2: 1516, 4: 248.
- **Rapidity (LOW, `x_extra.out`).** The adjoint action of A has eigenvalues φ^{±4} and 1. So its image in SO(2,1) is a boost of rapidity 4 log φ = 1.9248. The "2 log φ" in the text is the spinor half-rapidity.

## B. THM-4762 (2)–(3): CONFIRMED, citation venue corrected

**Normalisation independence: correct.**
- Any two linear lifts differ by a character of SL2(Z), and A is a commutator, so ρ(A) = C for every lift.
- Swapping S for S⁻¹ = S̄ also changes nothing, because s⁻¹ = −s and S² commutes with T.
- The pointed-lattice claim is a tautology: S and the normalised T are built from the discriminant form (G, q) alone.

**Ng–Schauenburg, CMP 300 (2010).** Thm 6.8, as restated in Dong–Lin–Ng: the kernel of the projective representation is a congruence subgroup of level N = ord(t̃), where t̃ is the twist matrix normalised so that θ₀ = 1. This equals the order of [T] in PGL, so the level must be computed with T₀₀ = 1. `catmap_mtc.py` does this (h_unit = 0).

**Dong–Lin–Ng.** The paper is in *Algebra & Number Theory* 9 (2015) 2121–2166, not Adv. Math. Its Theorem II says that for every linear modular representation ρ of any MTC, ker ρ is congruence of level n = ord ρ(t), with N | n | 12N. Theorem 3.10 covers rational C2-cofinite VOAs with T = t̃·e^{−2πic/24}. So the attribution of the linear statement is right; only the venue is wrong.

**"C^{π(N)/2} is scalar": correct.**
- A^{π(N)/2} ≡ I mod N exactly. If it were only ≡ −I, it would give the charge conjugation S², which need not be scalar.
- For N ≤ 2, π(N)/2 is not an integer. The statement should use ord_N(A), which is 1 for N = 1 and 3 for N = 2. The same applies to N′.

**"ord C | π(N′)/2": correct.** It follows from DLN Thm II applied to ρ(s) = S, ρ(t) = T̄·e^{−2πic/24}. This is a genuine representation because (ST̄)³ = e^{2πic/8}·S², which I verified for all 135 data.

## C. Census: CONFIRMED for all 135 rows (`c_modular_data.py`, `.out`, `c_modular_data_all.out`, `c_census_audit.json`)

**My builders.**
- SU(3)_k from the Kac–Peterson formula in orthonormal R³ coordinates.
- D(S_3) from the Coste–Gannon–Ruelle group formula, with no hard-coded S.
- Lattices from my own enumeration of L*/L.
- Minimal models from my own Kac table.

**Checks on every datum.**
- S is symmetric and unitary, and S² is a permutation matrix.
- (ST)³ = Θ·S² with Θ equal to both e^{2πic/8} and the Gauss sum over D.
- Verlinde coefficients are non-negative integers (rank ≤ 100).

**How the orders were found.** I took eigenvalues at 50 digits for rank ≤ 16, then confirmed C^o = I, C^{po} = λI, and minimality at each prime divisor.

**Required list.** Every entry equals the census:

| Model | N | ord_N(A) | proj ord C | ord C | N′ |
|---|---|---|---|---|---|
| Fibonacci | 5 | 10 | 5 | 10 | 60 |
| Ising M(3,4) | 16 | 12 | 3 | 3 | 48 |
| Lee–Yang M(2,5) | 5 | 10 | 5 | 10 | 60 |
| M(2,7) | 7 | 8 | 4 | 4 | 42 |
| M(4,5) | 80 | 60 | 15 | 30 | 240 |
| SU(3)_1 | 3 | 4 | 4 | 4 | |
| SU(3)_2 | 15 | 20 | 20 | 20 | |
| D(S_3) | 6 | 12 | 6 | 6 | |
| A1, A1+U, A1(−1) | 4 | 3 | 3 | 6 | |
| A2, A2+U | 3 | 4 | 4 | 4 | |
| golden lattice | 5 | 10 | 10 | 10 | |

- **SU(2)_1..12**, (proj ord, ord): (3,6), (3,3), (15,15), (12,12), (12,12), (12,12), (6,12), (30,30), (15,30), (12,12), (21,21), (24,24).
- **D(Z_N), N = 2..10**: proj ord = ord = 3, 4, 3, 10, 12, 8, 6, 12, 30.

**Full rebuild.** All 135 rows agree in all seven fields. Projective order divides ord_N(A), and ord C divides ord_{N′}(A), in every case. Deficiency counts are r = 1: 88, r = 2: 43, r = 4: 4 (M(3,4), M(4,5), M(4,11), SU(2)_2).

**D(S_3).** The session's hard-coded S-matrix and twists match my group-theoretic data exactly, with the same labelling.

**Saturation count (LOW).** 88 cases saturate ord_N(A), but only 83 have N ≥ 3. The other five have N ≤ 2, where π(N)/2 is not an integer: M(2,3), E8 and U (N = 1), and D4 and D(Z_2) (N = 2).

**Minimal-model description (LOW).** The script actually uses p < p′ ≤ 12 with pp′ ≤ 66, giving 23 models. There are 46 minimal models with pp′ ≤ 66.

**Definite vs indefinite lattices (MEDIUM).**
- C(A1) = C(A1+U) and C(A2) = C(A2+U) exactly (difference 0.0).
- A1(−1) has discriminant form (Z/2, q = 3/4), which is not isomorphic to (Z/2, 1/4). Its C is the complex conjugate of C(A1), and max|C(A1) − C(A1(−1))| = 1.0. It is a wrong example of "identical C", though the orders agree.

**D(Z_N) is a theorem, not evidence (`x_extra.out`).** After a partial Fourier transform in the charge label, C is a 0/1 permutation matrix of a trace-3 element of SL2(Z/N) conjugate to A (checked for N = 3, 4, 5, 6, 8). So ord C = projective ord C = ord_N(A).

**Methodology (LOW).** `catmap_mtc.py` never checks (ST)³ ∝ S² or Verlinde integrality. All data turned out to be correct anyway.

## D. THM-4763: mathematics CONFIRMED, attribution HIGH (`d_wp_analytic.py`, `.out`)

**Construction.** By Poisson summation:
- E_g(−1/τ) = (2p)^{−1/2} Σ_d e(gd/2p) E_d.
- O(−1/τ) = τ·M_s·O, with M_s = −i√(2/p)·[sin(πgd/p)].
- Both pick up the phase e(g²/4p − 1/24) under τ → τ + 1.

So in the basis (E, τO, O), ρ(s) = R(s) ⊕ s⊗M_s and ρ(t) = R(t) ⊕ t⊗M_t.

**Checks, with no fitting.**
- F(γτ) = ρ(γ)F(τ) at random τ for γ = s, t, A, A² and a random word. Maximum relative error is 1e−40 to 2e−40 for p = 2..6, using my own η (pentagonal series).
- The relations ρ(s)² = I, (ρ(s)ρ(t))³ = ρ(s)² and M_s² = −I hold.
- ρ(A) = R(A) ⊕ A_def ⊗ M(A) exactly.
- The eigenvalue moduli are {φ⁻², 1, φ²}.
- The trace law holds for n = 1..12, with M(A) computed independently.
- For random hyperbolic words the spectral radius equals the dilatation, and ρ(t) has Jordan blocks.

**Traces reproduce `log_catmap.out` exactly.**
- p = 2: 3, 7, 21, 47, 123, 325, 843, 2207.
- p = 3: 3, −5, −36, −49, 123, 640, 843, −2209.
- p = 4: 1, −1, 55, −1, 1, 971, 1, −1. Here the Ising-type odd block R_3(A) has order 3.

**p = 2 (LOW).** R_1 is not trivial. It is the η² multiplier character, with M(s) = −i and M(t) = e^{2πi/12}. It is trivial on A, so the trace formula is unaffected.

**The session's trace-law check was tautological (MEDIUM).** `log_catmap.py` sets trR = tr(Clogⁿ)/L_{2n} and then tests Tr Cⁿ = Tr CEⁿ + L_{2n}·trR. That only tests block-diagonality. The law itself is true (shown above), but the status line overstates what was verified.

**FGST citation, CMP 265 (2006) 47–93.**
- The abstract does state the equivalence of the centre representation with the extended characters, and the product of two commuting representations.
- More importantly, eq. (1.2) and Theorem 2.3 state Z_cft = R_{p+1} ⊕ C² ⊗ R_{p−1} explicitly. Here R_{p+1} is the lattice-VOA characters, R_{p−1} is the sl(2)_{p−2} characters, and "C2 is the standard two-dimensional SL(2, Z)-representation".
- §5 proves the same structure on the quantum-group centre, and §2.2 calls the (3p−1)-dimensional span the space of torus conformal blocks.
- So the decomposition is KNOWN. THM-4763's new content is only the cat-map corollary (Lucas traces, spectral radius), which is one line from FGST.

**Identification of the span: accurate.** FGST's 2p characters span exactly ⟨E, O⟩ (rank 2p for p = 2..6). The span ⟨E, O, τO⟩ is invariant, has dimension 3p − 1 and contains them, so it equals FGST's Z_cft.

**LOW wording.** "W(p) ≃ Ū_q sl(2)" should be "Kazhdan–Lusztig-dual". "Affine VOAs at positive integer level ≃ quantum groups" should say "semisimplified tilting categories (Finkelberg)". The rapidity "log λ(γ)" is the spinor convention.

## E. HYP-9372 (`e_center.py`, `e_center.out`)

**Spot checks** (my data, both directions). All agree with `center_check.out`.

| Case | N | Central powers of A, j (A^j) | Scalar powers of C | proj ord | Faithful (C^j scalar ⟹ A^j central) |
|---|---|---|---|---|---|
| Ising, SU(2)_2 | 16 | 6 (9I), 12 | 3, 6, 9, 12 | 3 | no |
| M(4,5) | 80 | 30 (41I), 60 | 15, 30, 45, 60 | 15 | no |
| M(4,11) | 176 | 30 (89I), 60 | 15, 30, 45, 60 | 15 | no |
| M(4,7) | 112 | 12 (97I), 24 | 12, 24 | 12 | yes |
| M(4,9) | 144 | 6 (89I), 12 | 6, 12 | 6 | yes |
| D(Z_5) | 5 | 5 (4I), 10 | 10 | 10 | yes |
| SU(3)_2 | 15 | 10 (11I), 20 | 20 | 20 | yes |

**Cases center_check.py skipped.** The 9 census entries with N ≥ 2 that it omitted (8 lattices and Fib ⊠ D(Z_2) ⊠ D(Z_3)) are all faithful in the one-directional sense. So that form holds in **128 of 132**, not "119 of 123".

**The converse fails (MEDIUM).** For these cases A^j is central but C^j is not scalar:
- D(Z_3), D(Z_5), SU(3)_1, SU(3)_2, A2, A2+U, A4, the golden lattice and the Deligne product. Here −I acts as a nontrivial charge conjugation.
- The self-dual SU(2)_4: A⁶ = 17I mod 24, but C⁶ is not scalar.

So the title's "exactly when" is false. `center_check.py` only tests the subset direction, which is the form actually supported.

**Sharpening of (a).** |⟨A⟩ ∩ Z| is 1 or 2 for every N ≥ 3. Proof: Q^α = λI with λ² = (−1)^α, so π/α = 4 exactly when α is odd, and then only Q^{2α} = −I lies in ⟨Q²⟩. Hence faithful cases have r ∈ {1, 2}, and r = 4 forces non-faithfulness. The HYP's {1, 2, 4} is true but loose.

**Mechanism (c): confirmed and essentially proved.**
- Use the basis (χ₀+χ_ψ, χ₀−χ_ψ, √2·χ_σ), which is proportional to √(θ₃/η), √(θ₄/η), √(θ₂/η).
- In this basis S is the permutation of the last two vectors, and T is monomial over the first two. These are exactly s and t mod 2 acting on the three even spin structures.
- C is monomial over the 3-cycle. Its entries are e^{−iπ/8}, 1, e^{iπ/8} for Ising and e^{∓3iπ/8}, 1 for SU(2)_2. The product is 1, so C³ = I.
- So C does not "kill the phases"; their product around the cycle is 1.

**Deligne factorisation.** M(4,p′) ≅ Ising_ν ⊠ X holds exactly, where X = {(r,1) : r odd} is a valid modular datum with S_X = 2·S restricted to X:
- p′ = 5: ν = 7, and X is the conjugate Fibonacci.
- p′ = 7, 9, 11: ν = 13, 3, 9.
- The projective order is lcm(3, po X) = 15, 12, 6, 15, matching the census.

**Masking text (LOW).** It uses ord_{p′}(A) = 8, 12 where the relevant numbers are po(X) = 4, 6. M(4,9) still has r = 2.

**π/α ∈ {1, 2, 4}** is classical (Vinson 1963; Wall 1960), and the counts above agree.

## F. HYP-9373

**Miyamoto.** The citation is correct: Duke Math. J. 122 (2004) 51–91. It proves SL2(Z)-invariance of the span of pseudo-traces, which are polynomials in log q with q-series coefficients. That supports the τ-polynomial structure. It does not give the Sym^k ⊗ R_k grading or finiteness of R_k, and the HYP rightly lists that as the open input.

**Prior art (HIGH).** FGST, Nucl. Phys. B 757 (2006) 303–343, Theorem 1.2(3):
- The SL2(Z) representation generated by the W_{p₊,p₋} characters has dimension (3p₊−1)(3p₋−1)/2.
- It decomposes as R_min ⊕ R_proj ⊕ C²⊗(R ⊕ R′) ⊕ C³⊗R_min, with C³ = S²(C²).
- Their footnote 2 ties C^{n+1} to the top power n of τ (1 for (1,p), 2 for (p₊,p₋)). That is exactly the HYP's "Why" argument.
- Their Conjecture 1.3 says this space is the space of torus amplitudes.

So the K = 2 case is known, and φ⁴ is a corollary. The J. Math. Phys. paper is vol. 48 (2007) 032303, on the KL-dual quantum group.

**"PROVED for K = 0" (MEDIUM).** This needs rationality, via the Dong–Lin–Ng congruence theorem. A non-rational C2-cofinite VOA with K = 0 is not covered.

**Otherwise.** OPEN status is appropriate, and the tensor-product consistency claim is correct.

## G. Synthesis note (`g_rr.py`, `g_growth.py`, `g_mrb_collatz.py` and their `.out` files)

**Identities.** All 18 verify exactly to q^250 with my own series code:
- Rogers–Ramanujan G, H; Rogers–Selberg A, B, C; Rogers mod 14 (Slater 59–61); Göllnitz–Gordon 1, 2.
- Bailey mod 9 (Slater 40–42, including B with its (1 − q^{3n+2}) factor); Dyson mod 27 (Slater 90–93).
- Jackson–Slater, with both the session's product and Slater (39), (−q³, −q⁵, q⁸; q⁸)/(q²; q²).

**Relations: all CONFIRMED.**
- Jackson–Slater sum = the Rocha-Caridi Ising character χ_{1,1}, i.e. q^{1/48}χ₀. It is also the even part of Σ q^{n²/2}/(q;q)_n.
- Rogers–Selberg A, B, C = the M(2,7) characters with r = 3, 2, 1, times (q;q²)_∞.
- Rogers mod 14 = the same M(2,7) characters at q², times (−q;q)_∞.
- Dyson A, C, D = Bailey A, B, C at q³, times (q⁹;q⁹)/(q;q). This holds even at the level of the sum sides.
- The Dyson products are the M(2,27) characters with r = 12, 9, 6, 3.
- The fifth-order mock theta identity (q;q)_∞·f₀ = the Hecke-type sum holds.

**c_eff and q → 1 ratios.**
- c_eff at q → 1 is 2/5, 1/14, 11/14, 3/8, 0, 8/9, 1/2 respectively. All CONFIRMED.
- The q → 1 ratios match 1/φ, the heptagon and nonagon sine ratios, and tan(π/8) to ≤ 1.4e−10.

**Coefficient growth: CORRECTED (MEDIUM).** I computed the growth exponent at every cusp h/k with k ≤ 36.
- For the positive families (RR, Rogers mod 14, Göllnitz–Gordon, Dyson, Jackson–Slater), q → 1 dominates and the claim holds.
- For Rogers–Selberg, q → −1 dominates with c* = 11/28 (not 1/14). Exact coefficients to n = 6000 fit c = 0.3927, with alternating signs and log|a₆₀₀₀| ≈ 116.
- For Bailey mod 9, q → e^{±2πi/3} dominates with c* = 8/27. Coefficients fit c = 0.2960, with sign pattern + − − and log|a₆₀₀₀| ≈ 100.
- So "coefficients grow like exp(2π√(c_eff·n/6))" is false for these two families, and "no growth" for Bailey mod 9 is false.

**"Each sum side has a sign-twisted mock-modular companion" (LOW).** Only Rogers–Ramanujan is verified; the mod-7 partner is cited. The general statement is an overgeneralisation.

**Ising ring.** Exact in Q(√5) for n ≤ 20, with K = −ln5/4 (antiferromagnetic):
- e^{−nK}·Z_periodic = 2ⁿLₙ and e^{−nK}·Z_antiperiodic = 2ⁿ√5·Fₙ.
- The eigenvalue ratio is −1/φ².
- The unnormalised Z carries an extra factor 5^{−n/4} (LOW).

**Collatz.**
- v₂(3^k − 1) = 1 for odd k and 2 + v₂(k) for even k, for all k ≤ 2000.
- 2^k − 1 reaches 3^k − 1 after k steps on odd values (checked k ≤ 79).
- v₂ = k only for k = 1 and k = 4, so "unique" needs k ≥ 2 (LOW).

**MRB.**
- S₅ = 0.006447893 > 1/200, and the minimum even partial sum is 0.22387392 at N = 58.
- But S₃ = −0.028036, so "the minimum odd partial sum for N ≤ 59 is S₅" holds only for 5 ≤ N ≤ 59 (LOW).
- The trace argument holds: the minimal polynomials of n^{1/n} for n ≤ 30 have zero trace coefficient.
- I recomputed D_p by complete-period block summation. It agrees with the session to ≤ 6e−5 and sits at about 0.1 of the stated bound for all six primes.
- The source PDFs and abstract are not in the repo, so fidelity to the source documents is UNVERIFIABLE.

**§7 list (LOW).** The census has 12 entries with projective order 15: M(3,5), M(3,11), M(4,5), M(4,11), M(5,11), SU(2)_3, SU(2)_9, SU(2)_20, SU(2)_29, SU(2)_59, SU(2)_60 and D(Z_22). The note lists five.

## H. Prior art

**Cat map as [T,S] with Pisano bounds.** I found no explicit statement, which agrees with the session. Related work it should cite:
- Jeffrey, CMP 147 (1992): torus bundles with hyperbolic SL2(Z) monodromy in TQFT.
- Charles, arXiv:1005.3447: modular-group elements act as Fourier integral operators.
- The abelian case is the quantum cat map: Hannay–Berry, and Dyson–Falk, Amer. Math. Monthly 99 (1992) 603–614.
- Chandra (arXiv:2607.24857) does define π(N) as the order of A in SL(2, Z/NZ), so the session's citation is accurate.
- The Pisano bound itself is a one-line corollary of Ng–Schauenburg and Dong–Lin–Ng.

**Non-semisimple side.**
- The decomposition is KNOWN (FGST CMP 2006, Thm 2.3), and so is the K = 2 structure (FGST NPB 2006, Thm 1.2).
- Infinite-order mapping-class actions in non-semisimple TQFTs: Blanchet–Costantino–Geer–Patureau-Mirand, arXiv:1404.7289.
- I found no explicit Lucas-trace or spectral-radius statement, but it is an immediate corollary of FGST.

**Rogers–Ramanujan dictionary.**
- Jackson–Slater as the Ising vacuum character: Andrews–van Ekeren–Heluani, arXiv:2005.10769; Kedem–Klassen–McCoy–Melzer, hep-th/9211102.
- Bailey mod 9 and Dyson mod 27: Bailey, Proc. LMS 49 (1947).
- Dyson B is Dyson's "favourite" identity: Andrews, Ann. Comb. 23 (2019).
- The product-level relations are one-line rewritings via (q;q) = (q;q²)(q²;q²) and θ_{3b,27}(q) = θ_{b,9}(q³). They carry no novelty.

## Required edits, ranked

1. **[HIGH] THM-4763 status.** Current text: "The decomposition and its consequences are elementary from the classical modularity of unary theta series … is CITED: Feigin-Gainutdinov-Semikhatov-Tipunin, Comm. Math. Phys. 265 (2006), hep-th/0504093."
   Replace with: "The decomposition ρ = R_(p+1) ⊕ C² ⊗ R_(p−1), with C² the standard representation, is KNOWN: FGST, CMP 265 (2006) 47–93, eq. (1.2) and Thm 2.3. FGST also identify the span with the torus conformal blocks and prove the equivalence with the centre (Thm 1.1, §5). It is re-derived here. NEW: the cat-map corollaries (Lucas traces, spectral radius φ²)."
   Also line 43: replace "Their statement that the ribbon element's Jordan decomposition "splits this representation…" is this tensor structure." with "FGST state this decomposition explicitly (eq. (1.2), Thm 2.3); the product of two commuting representations in their abstract is a reformulation of it."
2. **[HIGH] HYP-9373 status.** Current text: "The W(p,p') case is untested here. Its SL2(Z) structure (Feigin-Gainutdinov-Semikhatov-Tipunin, J. Math. Phys. 2006) is cited as the likely source of a C^3 = Sym^2 C^2 factor (UNVERIFIED)."
   Replace with: "KNOWN: FGST, Nucl. Phys. B 757 (2006) 303–343, Thm 1.2(3) gives R_min ⊕ R_proj ⊕ C² ⊗ (R ⊕ R′) ⊕ C³ ⊗ R_min with C³ = S²(C²); their footnote 2 ties C^(n+1) to the top power n of τ. Hence φ⁴ on that space. (J. Math. Phys. 48 (2007) 032303 is the KL-dual quantum group paper.)"
   Title: replace "The predicted next case is K = 2 … rank-3 Jordan blocks." with "K = 2 is realised by W(p,p′), whose C³ = Sym² C² piece is known (FGST, NPB 757, Thm 1.2), giving φ⁴." Credit FGST footnote 2 in the "Why" section, and drop the W(p,p′) question in note §8.
3. **[MEDIUM] HYP-9372 title (b), and the `center_check.py` docstring.** Replace "C^j is scalar exactly when A^j is a central matrix lambda I" with "C^j scalar implies A^j = lambda I central". Add a note that the converse fails for D(Z_3), SU(3)_1, A_2 and SU(2)_4.
4. **[MEDIUM] THM-4763 status.** Replace "the Lucas trace law holds for n = 1..10" with "the Lucas trace law holds for n = 1..12 (audit_S55/d_wp_analytic.out; the original log_catmap.py test was tautological)". Alternatively, fix `log_catmap.py` to take M(A) from the (O,O) block of C and test against that.
5. **[MEDIUM] Note line 57.** Replace "The coefficients grow like `exp(2π√(c_eff·n/6))`." with: "c_eff is the q → 1 exponent. It governs coefficient growth only for the positive families. Rogers–Selberg grows with c* = 11/28 from q → −1 (alternating signs); Bailey mod 9 grows with c* = 8/27 from q → e^{±2πi/3}."
   Line 67: replace "weight-0 modular units (no growth; signed coefficients)" with "c_eff = 0 at q → 1, but the coefficients grow like exp(2π√(4n/81)), signed with period 3".
6. **[MEDIUM] THM-4762 line 56 and note line 38.** Delete `A_1(−1)` from the identical-C examples. Add: "A_1(−1) has the conjugate discriminant form (Z/2, 3/4), so its C is the complex conjugate of C(A_1), with the same orders."
7. **[MEDIUM] HYP-9372 line 28.** Replace "`|⟨A⟩ ∩ Z| ∈ {1, 2, 4}`" with "|⟨A⟩ ∩ Z| ∈ {1, 2} (π/α = 4 exactly when α is odd; then only −I lies in ⟨A⟩)".
   Line 30: replace "This explains (a) in the faithful cases." with "Faithful cases therefore have r ∈ {1, 2}; r = 4 forces non-faithfulness."
8. **[MEDIUM] HYP-9373 line 6.** Replace "PROVED for K = 0 (THM-4762 with the congruence theorems)" with "PROVED for rational C2-cofinite VOAs (Dong–Lin–Ng)".
9. **[LOW] HYP-9372 status and body.**
   - Line 7: replace "119 of the 123" with "128 of the 132 (the audit adds the 9 lattice and Deligne entries)".
   - Line 10: replace "(c) is a proposed mechanism" with "(c) is proved for Ising-type factors (monomial 3-cycle, C³ = I); M(4,p′) ≅ Ising_ν ⊠ X is verified for p′ = 5, 7, 9, 11".
   - Line 46: replace "kills those phases and keeps the permutation" with "is monomial over the 3-cycle A mod 2, with phase product 1".
   - Line 40: use po(X) = 4, 6 in the masking explanation, and note that M(4,9) has r = 2.
10. **[LOW] THM-4762 line 6.** Replace "Dong-Lin-Ng (Adv. Math. 2015)" with "Dong-Lin-Ng (Algebra Number Theory 9 (2015) 2121-2166, Thm II)".
11. **[LOW] THM-4762, N ≤ 2 cases.**
    - Title (3): replace "C^(pi(N)/2)" and "pi(N')/2" with "ord_N(A)" and "ord_N′(A)", which equal π(N)/2 and π(N′)/2 when N, N′ ≥ 3.
    - Census claim: replace "equals pi(N)/2 in 88 cases" with "equals ord(A mod N) in 88 cases (83 with N ≥ 3)".
12. **[LOW] THM-4762 line 64.** Replace "with `pp' ≤ 66`" with "with `p < p' ≤ 12` and `pp' ≤ 66` (46 minimal models have pp′ ≤ 66)".
13. **[LOW] Rapidity.** In THM-4762 line 37 and note line 18, "rapidity `2 log φ`" becomes "SO(2,1) rapidity 4 log φ (spinor eigenvalues φ^{±2})". In THM-4763 line 66, "log λ(γ)" becomes "2 log λ(γ)". In HYP-9373 line 41, "2K log φ" becomes "4K log φ".
14. **[LOW] THM-4763 line 53.** Replace "`R_(p−1)` is trivial" with "R_1 is the η² multiplier character (order 12), trivial on A".
15. **[LOW] THM-4763 lines 63–64.** Replace "≃" with "Kazhdan–Lusztig-dual" in line 64, and with "≃ semisimplified tilting categories (Finkelberg)" in line 63.
16. **[LOW] Note line 90.** Replace "for `N ≤ 59`" with "for 5 ≤ N ≤ 59 (S₃ = −0.0280)".
17. **[LOW] Note line 99.** Replace "the ring partition functions become" with "the normalised e^{−nK}·Z (K = −ln5/4) become".
18. **[LOW] Note line 109.** Replace "the unique coincidence `v_2 = k`" with "the unique coincidence with k ≥ 2 (k = 1 gives 1 → 2)".
19. **[LOW] Note line 115.** Replace the list with the full set of 12 entries from section G.
20. **[LOW] Note line 76.** Restrict the companion claim to Rogers–Ramanujan (verified) and the mod-7 family (Andrews 1986).
21. **[LOW] Prior art.** Add Jeffrey 1992 and Charles 2010 to THM-4762, and state that the D(Z_N) case is a theorem (permutation representation). Add Andrews–van Ekeren–Heluani, Kedem–Klassen–McCoy–Melzer, Bailey 1947 and Andrews 2019 to note §4. Add modular-relation checks to `catmap_mtc.py`.

## Files

All in `C:\Users\Eliott\AppData\Local\Temp\math-wt-s55-signature-mtc\04-computation\experiments\catmap_modular_20261022\audit_S55\`:

| Script | Output | Covers |
|---|---|---|
| `a_pisano.py` | `a_pisano.out` | A |
| `c_modular_data.py` | `c_modular_data.out`, `c_modular_data_all.out`, `c_census_audit.json` | B, C |
| `d_wp_analytic.py` | `d_wp_analytic.out` | D |
| `e_center.py` | `e_center.out` | E |
| `g_rr.py` | `g_rr.out` | G (identities, relations) |
| `g_growth.py` | `g_growth.out` | G (c_eff, growth, ratios) |
| `g_mrb_collatz.py` | `g_mrb_collatz.out` | G (Ising ring, LTE, MRB) |
| `x_extra.py` | `x_extra.out` | C (D(Z_N)), A (rapidity) |

Sources:
- [FGST, CMP 265 (2006), hep-th/0504093](https://arxiv.org/abs/hep-th/0504093)
- [FGST, Nucl. Phys. B 757 (2006), hep-th/0606196](https://arxiv.org/abs/hep-th/0606196)
- [FGST, J. Math. Phys. 48 (2007), math/0606506](https://arxiv.org/abs/math/0606506)
- [Ng–Schauenburg, arXiv:0806.2493](https://arxiv.org/abs/0806.2493)
- [Dong–Lin–Ng, arXiv:1201.6644](https://arxiv.org/abs/1201.6644) and [ANT 9 (2015)](https://msp.org/ant/2015/9-9/ant-v9-n9-p03-s.pdf)
- [Miyamoto, Duke 122 (2004)](https://projecteuclid.org/journals/duke-mathematical-journal/volume-122/issue-1/Modular-invariance-of-vertex-operator-algebras-satisfying-C-2/10.1215/S0012-7094-04-12212-2.full)
- [Chandra, arXiv:2607.24857](https://arxiv.org/abs/2607.24857)
- [Dyson–Falk reference (OEIS A281230)](https://oeis.org/A281230)
- [Jeffrey, CMP 147 (1992)](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-147/issue-3/Chern-Simons-Witten-invariants-of-lens-spaces-and-torus-bundles/cmp/1104250751.pdf)
- [Charles, arXiv:1005.3447](https://arxiv.org/abs/1005.3447)
- [Blanchet–Costantino–Geer–Patureau-Mirand, arXiv:1404.7289](https://arxiv.org/abs/1404.7289)
- [Andrews–van Ekeren–Heluani, arXiv:2005.10769](https://arxiv.org/abs/2005.10769)
- [Kedem–Klassen–McCoy–Melzer, hep-th/9211102](https://arxiv.org/abs/hep-th/9211102)
- [Bailey mod 9 (MathWorld)](https://mathworld.wolfram.com/BaileyMod9Identities.html), [Dyson mod 27 (MathWorld)](https://mathworld.wolfram.com/DysonMod27Identities.html)
- [Andrews, Dyson's "favorite" identity](https://pure.psu.edu/en/publications/dysons-favorite-identity-and-chebyshev-polynomials-of-the-third-a/)
- [Knopp–Mason, arXiv:0910.3976](https://arxiv.org/abs/0910.3976)
- [Creutzig–Gannon, arXiv:1605.04630](https://arxiv.org/abs/1605.04630)
