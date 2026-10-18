# Referee report: kind-pasteur-2026-10-18-S2 (THM-4708, HYP-9330, HYP-9331, results note)

Saved verbatim by the session from the subagent referee's final message (2026-10-18). The subagent environment
returned the report as text; HTML entities from the transport (&amp;lt; &amp;gt;) were converted back to < and >.

Audited worktree: `math-wt-kps-s2-bridges-20261018` at `5958495dfb`. The referee's code is in `04-computation/cm_detector_sporadic/audit_20261018/` (absolute path `C:\Users\Eliott\AppData\Local\Temp\math-wt-kps-s2-bridges-20261018\04-computation\cm_detector_sporadic\audit_20261018\`). Nothing was committed and no audited file was edited.

**Verdict: MAJOR 4, MINOR 6, NOTE 10**

**Summary.**
- **The mathematics holds up.** I reproduced every theorem-level statement by independent methods, and no number in the four documents failed.
  - THM-4708 parts 1–4 are correct, including both ten-element lists and all twenty D labels.
  - The finite data behind HYP-9330 and HYP-9331 reproduce exactly. They survive extensions of 10–100× in range, and HYP-9331 passed two out-of-sample tests the author did not run.
- **The MAJOR findings are about attribution and reporting, not mathematics:**
  - M1: the THM-4708 lists and the level-5 parametrization are already in the literature (Beukers 2024 App. C; Z.-W. Sun 2011), uncredited.
  - M2: half of HYP-9330(C) is in Caruso–Fürnsinn–Vargas-Montoya–Zudilin (2025).
  - M3: the stated inert-prime valuation law is contradicted by the author's own output.
  - M4: a pre-registered prediction (P2) is mis-scored as PASSED.

## 1. Method and independence

- **Code written from scratch.** All code was written from the definitions in the audit brief. The author's modules were never imported or run. The author's `.out`/`.json` files were read only after my runs, for comparison. They agree with mine numerically everywhere, which is also how M3 and M4 were traced to the author's own files.
- **Methods chosen to differ from the author's:**
  - **Characters:** computed from the Rocha–Caridi formula, not from Andrews–Gordon products.
  - **m_5, m_7:** compared with Monster McKay–Thompson coefficients in the OEIS b-files A007251 (5A) and A007264 (7A).
  - **CM classification (no Hilbert class polynomials, no j-maps, no resultants):**
    - For each of the 42 orders with h(D) ≤ 2, recounted by enumerating reduced forms, I took every one of the (N+1)·h(D) points of Γ0(N)\H lying over its CM j-values.
    - At each point I evaluated m_N directly from η-products at 300 digits.
    - From these values I rebuilt the Galois-stable polynomial ∏(X − m_N(P)), rounded it to Z, factored it with sympy, and kept the linear factors.
  - **(B1') detector:** truth tables indexed by m mod p², built for each prime with numpy int64 Horner evaluation.
  - **F_p polynomial algebra:** hand-written Euclid, radical and Yun square-free decomposition (no sympy or flint).
  - **Hulek–Verrill numbers hv_n:** computed from the multinomial definition as S_r(n) = Σ_k C(n,k)² S_{r−1}(n−k), and checked against the literal 5-fold sum for n ≤ 8.
  - **Valuations:** exact integer arithmetic mod p⁸.
- **Literature (primary-source bodies):**
  - Read myself: Beukers arXiv:2403.03301 v3 (and v1, v2 compared) for App. B and C, Rem. 1.8, 1.14 and 1.17, Thm 1.4, 1.13 and Conj. 1.27–1.28; Z.-W. Sun arXiv:1103.4325v3 (TeX source); CFVZ arXiv:2510.23298 (Thm 2–4, Table 2); Straub arXiv:2301.12248 (Thm 1.3); BEKK arXiv:2204.06506 (Table 1); Campbell–Cooper–Ye arXiv:2602.09352 (§4.5, level 5); OEIS A274786, A282877, A007251, A007264.
  - Gathered by two sub-agents from primary sources: Chen–Yui, Chan–Cooper 2012, A183204, and the remaining Beukers items. Every item used in a MAJOR finding I re-read myself.
- **Pre-registration provenance:** checked with `git log`, the remote-tracking push reflog, an append-only diff of PREREGISTRATION.md, and worktree mtimes.

## 2. The ten requested checks

| # | Item | Verdict | Evidence |
|---|---|---|---|
| 1 | THM-4708 part 1 | **CONFIRMED** | ∏_s χ_{1,s}^{(2,2K+1)} = (η((2K+1)τ)/η(τ))^{K−1}, exactly to O(q³⁰⁰) for K = 2..8, with leading power K(K−1)/12. The Andrews–Gordon form of each χ_{1,s} also equals Rocha–Caridi. Proof re-derived: Σ_{s≤K}(2K+1−2s)² = K(2K−1)(2K+1)/3. (`a1`) |
| 2 | Part 2 | **CONFIRMED** (attribution: M1) | m_5 − T_5A ≡ 16 and m_7 − T_7A ≡ 9, exactly through q¹⁴⁹⁹. Both generating-function identities hold exactly through q¹⁶¹. Both recurrences hold for 1 ≤ n < 699 against the binomial-sum definitions. OEIS's own A274786 formula Σ(−1)^j C(2n,j)C(j,n)³ equals C(2n,n)b_n for n < 300. The four level-5/7 j-maps hold as exact q-series to O(q¹⁴²). (`a2`, `a2b`, `a9`) |
| 3 | Part 3 | **CONFIRMED** by a different method (attribution: M1; rigor: n1) | Both lists come out exactly, all values integers. Every one of the 20 values sits on a Heegner pair with D(τ) = D(Nτ), so all D labels are right. Worst distance of any rebuilt coefficient to Z: 2.1e−234 (N=5), 1.9e−234 (N=7). Orders with N dividing the conductor give irrational values: m_5 = −108 ± 48√5 (D(5τ) = −75) and m_7 = −112 ± 24√21 (D(7τ) = −147). (`a3`, `a3c`) |
| 4 | Part 4 / HYP-9330(E) | **CONFIRMED** | (B1') and (B1'_ord) survivors for all 15 families over 0 < \|m\| ≤ 10⁵ are identical to the note's lists. They are unchanged up to 10⁷. All 18 nonzero rational CM values pass (B1') at all 165 primes 7 ≤ p ≤ 997. (`a4`, `a4c`, `a4d`) |
| 5 | HYP-9330(A) | **CONFIRMED** | The proof is correct; it needs only the Lucas congruence. Numerically, (B1') ⟺ V_p ≡ u₁q_pF_p (mod p) at 6644/6644 (m,p) pairs with p ≤ 257. (`a7`) |
| 6 | HYP-9330(B),(C) | **CONFIRMED as finite data** (novelty: M2) | No counterexample anywhere: all 11 families at all 75 primes 7–397, plus 11 further families at 59 primes 7–293 (n4). Simple roots, fibre by fibre, sit exactly at primes inert in the fibre's CM field: Apéry D=−24; Domb Q(√−3); s7 Q(√−7); s10 m=−4 ↔ D=−20 and m=16 ↔ D=−40; A274786 D=−20. (`a6`, `a6b`, `a6c`) |
| 7 | HYP-9330(D),(F) | (D) **CONFIRMED** (wording: m5); (F) **ISSUE** (M3) | Slope = 1 at every ordinary prime ≤ 257 where it is defined, at all CM points and cusps (22–29 primes per weight-2 point, 25 per weight-1 CM point). Valuation law: see M3. (`a7`, `a7b`) |
| 8 | HYP-9331 | **CONFIRMED** (wording: m3) | Survivors {−7,1,9,25}, {−5832,729}, {256}, {3125}, unchanged to 10⁷. Valuations as stated. Two out-of-sample confirmations (n5). (`a4c`, `a7`, `a10`, `a10b`) |
| 9 | Literature | **ISSUES** (M1, M2, m2, m4) | Details in §3–§5. |
| 10 | Typing / overclaim | **ISSUES** (M3, M4, m1, m6) | Statuses are otherwise appropriate. |

## 3. MAJOR findings

### M1. THM-4708 does not credit prior tabulations of its values and of its level-5 parametrization, including in a paper it cites as "checked"

**What the literature already contains:**
- **Beukers arXiv:2403.03301 v3, App. C §9.11 "Cooper s7"** gives:
  - the mirror map t = h₇/(1+13h₇+49h₇²) on Γ0(7)+⟨7⟩;
  - the rational CM values −1, 1/2, 1/8, −1/8, 1/24, 1/27, −1/64, 1/125, −1/10648, each with its τ.

  These are exactly the nine nonzero m_7 values of part 3.
- **Beukers App. C §9.14 "Apéry2 ⊙ central"**, g_n = C(2n,n)Σ_k C(n,k)²C(n+k,k), gives:
  - the operator (1−44t−16t²)θ³ + …, which is equivalent to the part-2 recurrence;
  - the mirror map t = h₅/(1+22h₅+125h₅²);
  - −4Θ = E_{2,5}, which is Σ U_n m_5^{−n} = (5E₂(5τ) − E₂(τ))/4;
  - nine rational CM values with τ: 1/4, −1/3, −1/28, 1/36, 1/18, 1/72, **1/588**, −1/828, −1/15228.
- **I evaluated m_N at all 18 printed τ** with my own η code (`a3b`). Seventeen agree to < 1e−55. At τ = √(−3/5) the value is **147**, so v3's "1/588" is a misprint. v1/v2 print 4/147 in their C(2n,n)/4ⁿ normalization, which is correct. **The THM's 147 is right.**
- **Z.-W. Sun, arXiv:1103.4325v3 (March 2011), §3**, defines a_n = Σ C(n,k)²C(n+k,k) and lists for Σ_{k<p} C(2k,k)a_k/m^k the pairs (m,d) = (−3,15), (4,11), (18,1), (−28,35), (36,19), (72,10), (147,15), (−828,115), (−15228,235). These are all nine nonzero level-5 values, 13 years earlier. His Conj. 3.12 also includes "≡ 0 mod p² if (−115/p) = −1".
- **Further sources:**
  - The A274786 recurrence is in OEIS (V. Kotesovec, Dec 2017).
  - The m_7 η-expression is OEIS A282877 (M. Somos, 2017).
  - Read only by a sub-agent: the level-5 family appears in Chan–Cooper 2012 (Thm 3.1/3.5, Table 7), Cooper's 2017 book and Campbell–Cooper–Ye 2026. I confirmed the CCY §4.5 discussion myself: the family C(2n,n)·b_n, G = 1 − 44x − 16x², and Sato's series at x = −1/15228.
  - **Huber–Schultz–Ye 2023** (Acta Arith. 207), cited on the OEIS A274786 page, says in its abstract that "A complete list of rational and quadratic series corresponding to singular values of the parameters is derived", in terms of Thompson series. Its body was not accessible to me or to the sub-agent.

**What the THM says instead:**
- Part 2 for A274786 is "a q-series identity checked exactly … found as the unique one-dimensional nullspace".
- The CM values are "expected to appear in the Ramanujan-Sato literature".
- The results note §5 says "presumably".

**What remains new:**
- Beukers presents App. C as "listings of the rational CM-values that we found" (Rem. 1.17). He claims completeness only for App. B (Rem. 1.8(3)).
- So the proof of completeness (over all CM τ), the inclusion of m = 0, and the link to the detector survivors are the THM's contribution. Completeness is routine, though (n3).

**Repair:**
- Part 2 (N=5): mark as CITED, giving Beukers v3 App. C §9.14, Chan–Cooper 2012, Cooper 2017, the OEIS recurrence (Kotesovec 2017), the product formula (Bala 2024) and A282877.
- Part 3: state that the nine nonzero values for each level are tabulated, without a completeness claim, in Beukers App. C §9.11/§9.14, and the level-5 ones in Z.-W. Sun 2011 §3. Note the v3 misprint 1/588 → 1/147.
- Read HSY 2023 and state precisely how its completeness differs from the THM's (all CM τ versus the Ramanujan–Sato singular points), or drop "no statement of completeness was located".
- Make the same corrections in results note §5 and §7.

### M2. HYP-9330(C), second bullet, together with P7's "post hoc" finding and the stated mechanism, is prior art: CFVZ, arXiv:2510.23298 (27 Oct 2025), which OEIS A274786 links

**What CFVZ contains:**
- **Thm 2 (proved):** for the Apéry numbers, A_p = B_p² if p ≡ 1, 5, 7, 11 (mod 24), and A_p = (t² − 34t + 1)·B_p² if p ≡ 13, 17, 19, 23 (mod 24). Here A_p is the note's F_p. The second residue set is exactly the primes inert in Q(√−6).
- **Thm 3 (Domb) and Thm 4 (Almkvist–Zudilin):** the same structure, proved.
- **Table 2 (computational):**
  - A274786: extra factor 1 − 44t − 16t² exactly when (−5/p) = −1;
  - A183204 = s7 (level 7, (−7/p));
  - A005260 = s10 (level 10, four cases by (−5/p) and (−10/p));
  - s18 and others.
- **Mechanism:** f(t(x)) = φ(x)h(x)² for a rational parametrization t(x). This is the note's heuristic ("F_p behaves like f², simple roots from the ramified 2:1 Fricke quotient").

My per-fibre data agree with CFVZ in every case.

**Related prior art on V_p:** Straub, arXiv:2301.12248, Thm 1.3, covers all 15 sporadic sequences: A(pn+k) ≡ A(k)A(n) + p·n·A′(k)A(n) (mod p²). Hence V_p ≡ u₁·Σ_{k<p} A′(k)x^k (mod p). So the "first-carry defect" is u₁ times a truncated formal-derivative series. This is a published version of the HYP's proposed "proof route for (B)". It covers Zagier A–F, s7 and s10, but not A274786.

**What the HYP claims:** "No statement of (B), (C) or the weight dichotomy was found (Beukers 2024; KLMSY 2016; OEIS A274786)". Results note §7 lists "(B) and (C)" as new.

**What still appears new** (none of these is in CFVZ, which treats only F_p mod p, or in Straub):
- rad(F_p) | V_p in weight 2;
- gcd(F_p, V_p) = 1 in weight 1 (square-freeness of F_p in weight 1 is Igusa-type);
- the resulting dichotomy for the mod-p² block test.

**Repair:**
- Cite CFVZ.
- Mark (C)'s second bullet as PROVED (CFVZ Thm 2–4) for Apéry, Domb and AZ, and OBSERVED (CFVZ Table 2) for s7, s10 and A274786.
- Cite Straub Thm 1.3 for the meaning of V_p.
- Restrict the novelty claim to rad(F_p)|V_p, gcd(F_p,V_p)=1 and the detector dichotomy.
- Append to PREREGISTRATION.md that the P7 post hoc pattern is in CFVZ.

### M3. "At inert primes the defect valuation equals the weight" is contradicted by the author's own output (Cooper s10, m = −2)

**Where the claim appears:**
- the HYP-9330 title;
- HYP-9330(F): "= 2 apart from sporadic coincidences of frequency about 1/p";
- results note §1 and §3: "Inert primes: v_p(Δ) = 1 in weight 1 and 2 in weight 2".

**What the data show:**
- The author's `defect_valuation_probe.out`, line 49, reads `w2 Cooper s10 m= -2 … supersingular v-histogram {3: 22, 4: 1}`. That is v_p ≥ 3 at every non-ordinary prime p ≤ 199; the chance expectation Σ1/p is about 0.6.
- My reproduction (`a7b`): v_p ≥ 3 at all 39 non-ordinary primes ≤ 400, and v_p = 3 at all ordinary primes as well.
- The point has CM by Z[i]: its non-ordinary set is the inert set of D = −4, and τ = (2+i)/5 (Beukers §9.12, t = −1/2).
- **The other 41 weight-2 CM points behave as claimed:** 27 have no non-ordinary prime with v ≥ 3, 10 have one, and 4 have two, which is chance level.
- The weight-1 half is correct: v = 1 at 22/22 non-ordinary primes at each of the 8 CM points.

**Repair:**
- Restate (F): v_p ≥ w at non-ordinary primes (forced by (B)/(C)), with equality generically, and a systematic exception at s10, m = −2.
- Delete the clause from the HYP-9330 title and from results note §1 and §3.
- Add s10 m = −2 (non-ordinary primes) to results open problem 4.

### M4. Pre-registered P2 is reported "PASSED in weight 2" although its control clause failed in the author's own output

- **The prediction** (P2, pushed in `8bf465cb8c`): "Non-survivor controls (m in {2, 3, 5, 7, 11, 100, −100, 999}) have density < 0.10."
- **The author's output** (`sporadic_cm_detector_M100000.out`): **Apéry m = 100: 0.13** and **Domb m = 999: 0.11**. I reproduced both exactly (`a4b`: 0.130 and 0.110).
- **The reports:** the PREREGISTRATION outcome ("P1 and P2 PASSED for all zeta(3)-type families … density 0.45-0.54") and the results-note table both score only the survivor clause.
- **No mathematical consequence:** survivor densities are 0.446–0.560, and the largest non-survivor density over all |m| ≤ 10⁵ is 0.264. But a pre-registered prediction must be scored as written.
- **Repair:** append a correction to PREREGISTRATION.md (the file has been append-only, see n6): "P2 in weight 2: survivor clause PASSED; control clause FAILED for Apéry (m = 100: 0.130) and Domb (m = 999: 0.110); overall PARTLY FAILED". Fix the P2 row of the results-note table and the "pre-registered P1, P2" wording in HYP-9330's status.

## 4. MINOR findings

**m1. Prime count.** "All 41 primes 7 ≤ p ≤ 199" appears twice in results note §1 and once in the P7 outcome ("of the 41 primes"). There are **43** such primes. The author's code (`primes_between(7, 199)`) ran all 43, and its simple-root counts (21/22/21/35) equal mine over 43 primes. Replace 41 with 43. (The "37 of 41" in HYP-9331 counts ordinary primes and is correct.)

**m2. Class-number citation.**
- "Arno" is not a source for h = 2 or for non-maximal orders: Arno 1992 is about class number 4, and Arno–Robinson–Wheeler 1998 about small odd class numbers.
- Correct sources:
  - fundamental discriminants with h = 2: Baker 1971, Stark 1975 and Montgomery–Weinberger 1974 (or Watkins 2004);
  - orders: the formula h(f²D_K) = h(D_K)·f·∏_{ℓ|f}(1 − (D_K/ℓ)/ℓ)/[O_K^×:O^×], which forces φ(f) ≤ 6 (e.g. Cox, Thm 7.24).
- My recount: 27 fundamental discriminants up to |D| ≤ 20000 (largest 427), and 13 + 29 orders.

**m3. HYP-9331 wording versus the author's own `weight3_attractor_scan.out`.**
- "Generic m = 5: v_p = 1 at ordinary primes, in every family" is false as written: v = 2 occurs at 2 ordinary primes for HV (p = 13, 67), at 1 for the ₄F₃ family, and at 3 for C(2k,k)⁴. Add "apart from sporadic coincidences".
- For C(2k,k)⁴ at m = −256, "there [p ≡ 3 mod 4] v_p = 1" fails at p = 127, where v = 3 and V(1/m) ≡ 0. The non-ordinary set equals {p ≡ 3 mod 4} exactly (23 primes ≤ 199), and v = 2 at all 20 ordinary primes.

**m4. Beukers citation precision.**
- (i) The numbering matches **v3 only** (16 Jun 2025). In v2 the lift theorem is 1.3, the inert-prime remark is 1.8(3), and the conjectures are 1.28–1.29. Cite "v3".
- (ii) Beukers does not call Thm 1.4 a "Dwork congruence"; he reserves that name for Thm 1.29.
- (iii) Rem. 1.14 lists as exceptions the 1728t case, s7, s10, **s18**, and **both** Hadamard products: Apéry2⊙central = A274786, and Franel⊙central. "s7, s10 and C(2n,n)b_n are exactly Beukers' Remark 1.14 exceptions" should read "are the three of our five weight-2 families on Beukers' exception list".
- (iv) "Match Beukers' appendix tables" holds for A, C, D, E and F. For B, v3 prints t = +1/9, +1/3 in its twisted variable (x = −t), i.e. m = −9, −3. That is a sign slip; v1/v2 are correct. The HYP's B = {3, 9} is the true set, and it should say so.
- (v) Thm 1.7 and 1.15 are stated for split p via α, c, d. "t₀^σ = t₀" is Beukers' introductory heuristic, not a hypothesis of the theorems.

**m5. HYP-9330(D) and results §3: "Truncated below 2p at a fixed point … these give our split-prime slope law".**
- Write t^σ = t^p(1 + pD(t)). Beukers' congruence mod p² then gives V_p ≡ u₁·(F_p·D mod t^p) (mod p).
- So V_p(t₀) = u₁F_p(t₀)D(t₀) − u₁t₀^pR(t₀), where R is the degree-≥p part of F_p·D_{<p}.
- The slope law therefore needs both D(t₀) ≡ q_p(m) and R(t₀) ≡ 0. Say "suggested by", not "expected from".
- Separately, at m = ±1 (the A and C cusps; the Apéry, Domb and s10 CM points) q_p(m) = 0, so the "slope" is undefined. (D) should say "where defined".

**m6. Physics framing in THM-4708.**
- The "Reading" says "these are the attractor / rational-CFT points of the torus E_τ". A CM τ is necessary but not sufficient for an RCFT (Gukov–Vafa also need ρ CM in the same field), as S43's own HYP-9311 says.
- The AD ↔ Virasoro (2,2K+1) identification needs a citation (e.g. Cordova–Shao 2015; Beem et al. 2013).
- "The three odd genus-0 levels 5, 7, 9 … are exactly the levels 2K+1" is numerology: every odd N is 2K+1, and the operative condition is (N−1) | 24.
- Suggested title: "Rational CM values of the Fricke Hauptmoduln T_5A+16 and T_7A+9", with the Argyres–Douglas reading kept as an identification/ANALOGY remark.

## 5. NOTES

**n1. Status of part 3 ("PROVED, computer-assisted").** Justified. The author's only floating-point step (H_D by 120-digit rounding, margin 6.9e−88) is bypassed entirely by my method (300 digits, margin 2e−234). For full rigor, either check H_D against a published table (all have degree ≤ 2) or use interval arithmetic.

**n2. Degree argument.** Correct.
- D may be read at τ or at Nτ: both j(τ) and j(Nτ) lie in Q(t), so both orders have h ≤ 2.
- A resultant over J = j(Nτ) alone is complete, because the second root b/t of the quadratic is the Fricke image τ′, where j(Nτ′) = j(τ).
- Orders with N dividing the conductor (−75, −100, −147) are among the 42 and are handled; they give irrational values (table, check 3).
- The sentence "it is the same at both Fricke-conjugate points" is true of the 20 listed points (all Heegner), not of CM points in general. Reword it.

**n3. Structural cross-check and the novelty of completeness.**
- The 20 discriminants are exactly:
  - {h(D) = 1 with N split or ramified}, together with
  - {h(D) = 2, N ramified, N ∤ conductor, prime above N non-principal}.
- D = −20 drops out at N = 5 because (√−5) is principal, which gives the two W₅-fixed points 22 ± 10√5.
- This is the standard description of rational CM points on X0⁺(N) (cf. Galbraith 1999; arXiv:1607.04558 §3 for h = 1). Completeness is therefore routine and deserves one sentence and a citation rather than a headline.

**n4. Extensions with no counterexample found.**
- Detector survivors for all 15 families are unchanged up to |m| ≤ 10⁷.
- The 18 CM values pass at every prime p ≤ 997.
- (B)/(C) hold at all 75 primes 7–397. They also hold at all 59 primes 7–293 for 11 families the author did not test:
  - weight 1: C(2k,k)², (3k)!/k!³, (4k)!/((2k)!k!²), (6k)!/((3k)!(2k)!k!);
  - weight 2: C(2k,k)³, C(2k,k)²C(3k,k), C(2k,k)²C(4k,2k), C(2k,k)C(3k,k)C(6k,3k), Almkvist–Zudilin, Cooper s18, Franel⊙central.
- This supports HYP-9330's general dichotomy.

**n5. HYP-9331 out-of-sample confirmations** (families chosen from BEKK Table 1, predictions stated before running):
- **AESZ 11**, ₄F₃(1/4,1/3,2/3,3/4) with u_k = (4k)!(3k)!/((2k)!k!⁵): survivors in 0 < |m| ≤ 10⁶ are exactly **{−432, 1728}**, i.e. BEKK's z₀ = −1/(2⁴3³) plus the conifold. v = 3 at 39/43 primes and 4 at 4.
- **AESZ 36**, A_n = C(2n,n)²Σ_k C(n,k)C(2k,k)C(2n−2k,n−k) (verified against the AESZ operator): survivors exactly **{−64, 64, 128}**, i.e. z₀ = −1/2⁶ plus the two conifolds. v = 2 at all 43 primes.
- BEKK Table 1 has no entry for AESZ 1 (quintic) or AESZ 3 (C(2k,k)⁴), consistent with the HYP's "none known" controls.
- Recommend adding these cases. The remaining 21 points of Table 1 are the natural next test.

**n6. Pre-registration provenance verified.**
- Push reflog (prediction pushed / outcome committed): P1–P3 at 13:06:23 / 13:11:22; P3′–P5 at 13:13:01 / 13:23:17; P6–P7 at 13:26:14 / 13:29:15; P8 at 13:31:36 / 13:33:00; P9 at 13:43:29 / 13:45:05.
- PREREGISTRATION.md is strictly append-only across its 10 commits (0 deleted lines).
- Worktree mtimes are consistent for P8 (script 13:32:04, outputs 13:32:10) and P9 (13:43:41, 13:43:47). Earlier mtimes were reset by rebase checkouts.
- Apart from M4 and m1, the outcomes are reported faithfully, including the failures of P1, P2 (weight 1), P3, P3′ and P7.

**n7. Mirror maps and j-maps.** Verified exactly to O(q⁶⁰): t_B(q) = −P₉(−q) for Zagier B, and t_D = q∏(1−qⁿ)^{5(n/5)} for Zagier D. The four j-maps hold far beyond the valence bound (agreement through q⁰ already proves them).

**n8. THM-4708 part 1 is folklore-level.** K = 2 is the classical G(q)H(q) = (q⁵;q⁵)∞/(q;q)∞, and general K is one line from Andrews–Gordon. No explicit general statement was found; the closest (Milas–Mortenson–Ono) concerns the Wronskian of the characters. It is fine as a lemma, but should not be advertised as new.

**n9. Integrality and inert-prime literature.**
- Integrality of the values can be cited from Chen–Yui, Thm I / Thm 3.4: T_g(τ) is an algebraic integer for fundamental Thompson series at CM points with (a,N) = 1, which covers 5A and 7A. They give no table of rational values.
- For context on the inert-prime picture: Z.-W. Sun's conjectures that F_p(1/m) ≡ 0 (mod p²) at inert p (e.g. his Conj. 3.12 at m = −828) are a related but different statement from rad(F_p) | V_p.

**n10. Beukers v3 erratum list,** worth recording since the repository relies on these tables:
- §9.14: t(√(−3/5)) = 1/147, not 1/588.
- §8.6 (B): the signs of 1/9 and 1/3 should be negative in the twisted variable.
- Thm 1.5 table: F's cusp values are printed in the untwisted sign.
- §8.6: the upper limit of B's sum should be ⌊n/3⌋, not ⌊n/2⌋.
- The first is mine; the last three are sub-agent findings, of which I re-verified the B sign.

## 6. Reproduction table

Run from `04-computation/cm_detector_sporadic/audit_20261018/` as `PYTHONIOENCODING=utf-8 python3 <script> <args>`. Requires Python 3, numpy, sympy 1.14 and mpmath 1.3. Outputs are the `.out` and `.json` files of the same name.

| Script and arguments | Checks | Range covered | Time |
|---|---|---|---|
| `audit_seqs.py` | Self-test: Zagier A–F against closed forms; hv_n against the literal multinomial sum | n < 30; n ≤ 8 | < 1 s |
| `a1_characters.py 300` | Check 1 (Rocha–Caridi products) | K = 2..8, O(q³⁰⁰) | 1 s |
| `a2_hauptmoduln.py 1500 160 700` | Check 2: m_N against OEIS b-files (`b007251.txt`, `b007264.txt`, trimmed to 1600 terms); GF identities; recurrences | q¹⁴⁹⁹; q¹⁶¹; n < 699 | 10 s |
| `a2b_a274786_definition.py 300` | OEIS A274786 definition = C(2n,n)b_n | n < 300 | 2 s |
| `a3_cm_classification.py` | Check 3: independent CM classification over all 42 orders | N = 5, 7; 300 digits | 12 s |
| `a3b_beukers_points.py` | m_N at Beukers' 18 printed CM points (the 1/588 misprint) | — | 2 s |
| `a3c_nonheegner.py` | Non-Heegner orbits are irrational | D = −3 with N ∣ f | 1 s |
| `a4_detector_scan.py 100000` | Checks 4 and 8: (B1'), (B1'_ord) and densities, 15 families | 0 < \|m\| ≤ 10⁵, 7 ≤ p ≤ 199 | 45 s |
| `a4b_p2_controls.py` | M4: P2 control densities | 7 ≤ p ≤ 500 | 10 s |
| `a4c_large_scan.py 10000000` | Survivors at large height, 15 families | 0 < \|m\| ≤ 10⁷ | 2.5 min |
| `a4d_cm_values_large_p.py 1000` | The 18 CM values pass (B1') | 7 ≤ p ≤ 997 | 3 s |
| `a6_polys_modp.py 7 199` and `a6_polys_modp.py 200 400` | Check 6: (B), (C), Lucas congruence, multiplicities | 75 primes | 4 s + 1 min |
| `a6b_singular_fibres.py 7 199` and `… 200 400` | Simple roots per singular fibre against inert primes | 75 primes | 4 s + 1 min |
| `a6c_more_families.py 7 199` and `… 200 300` | (B)/(C) on 11 untested families | 59 primes | 3 s + 30 s |
| `a7_slopes_valuations.py 257` (also run with 199) | Checks 5, 7 and 8: lemma (A), slopes, valuations | p ≤ 257 | 12 s |
| `a7b_inert_valuation_probe.py 400` | M3: inert-prime valuations at all 42 weight-2 CM points | p ≤ 400 | 40 s |
| `a8_mirror_maps.py 60` | Mirror maps of Zagier B and D | O(q⁶⁰) | 20 s |
| `a9_jmaps_exact.py 150` | j-maps for levels 5 and 7 | O(q¹⁴²) | < 1 s |
| `a10_aesz11_attractor.py 1000000` and `a10b_aesz36_attractor.py 1000000` | n5: out-of-sample HYP-9331 tests | 0 < |m| ≤ 10⁶ | 3 s each |

Session note: the two largest JSON dumps (`a7_slopes_valuations_p199.json`, `a7_slopes_valuations_p257.json`, about 1 MB together) are not committed; they are regenerated by the `a7` command above.
