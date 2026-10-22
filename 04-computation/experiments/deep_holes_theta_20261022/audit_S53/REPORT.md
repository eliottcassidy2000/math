# S53 audit: THM-4758, THM-4759, HYP-9368 and the note, plus THM-4760

Independent audit agent, session opus-2026-10-22-S53, 2026-10-22. The subagent harness blocked report-file writes, so
the main session saved this text verbatim from the agent's final message. The checking scripts and their `.out` files
are in this directory, and the audited text is snapshotted in `snapshot/`. The audit modified no session files.

**Versions.** Audited HEAD fc8da89dd0 plus the uncommitted edits as of 06:42, and read the registered version
c31e9a6a3f. Every sentence fix below was re-checked against HEAD 92a4c0a0af. The session read the interim outputs
during the audit and had already applied three fixes (19e71373ee, effcde315f, e23fbb78f4); these are marked
**FIXED at HEAD**.

## Summary table

| # | Item | Verdict |
|---|---|---|
| 1 | THM-4758 swap theorem; "PROVED" | CONFIRMED WITH CAVEAT (4 text fixes) |
| 2 | THM-4759 (1)–(2) | CONFIRMED; classical (add references); (2) is a one-line corollary |
| 3 | THM-4759 (3), weight-1 traces 24 and −k | CONFIRMED |
| 4 | THM-4759 (4) against CDH | CONFIRMED |
| 5 | THM-4759 (5), (5′) | CONFIRMED WITH CAVEAT (one missing line in the proof; component rule is a Frame-shape statement, FIXED at HEAD) |
| 6 | THM-4759 (7), (10) | (7) CONFIRMED. (10) had a real rationality gap as registered, now CLOSED. Lattice part of (10) is KNOWN. Two overclaims remain |
| 7 | THM-4759 (8), Borel–de Siebenthal descent | CONFIRMED; p = 2 is KNOWN |
| 8 | M24 point | CONFIRMED, and classical, so it can be PROVED (CITED) |
| 9 | Carnahan | CONFIRMED WITH CAVEAT (theorem numbering; class names now verified) |
| 10 | HYP-9368 | WRONG as registered: a counterexample to (A) at order 5. Definition FIXED at HEAD; one sentence remains |
| T1 | THM-4760 (1)–(2), cusp criterion and rationality certificate | CONFIRMED; condition (ii) is redundant |
| T2 | THM-4760 spot checks (d6⁴ and e8³, shell 1) | CONFIRMED by re-run and by an independent implementation |
| T3 | Carrier table (E7 in D10E7², E8 in D16E8) | CONFIRMED |
| T4 | Ligozat orders and class names | CONFIRMED |
| T5 | Is THM-4759 (10) complete now? | YES, given the exact finite enumerations |

## 1. THM-4758: CONFIRMED WITH CAVEAT

**What checks out.**
- **Sturm lemma.** Correct: each factor F_j|g is a combination of components, so its order is > k/12, so the norm
  form N_j has order > K/12 and the valence formula gives 0.
- **Weil representation.** Signature 24 means no phase. O(D) commutes with ρ_M. For a real vector, invariance under ρ
  and under ρ* are equivalent.
- **Certificate.** Re-run: 23/23 pass.
  - The S-invariance test is exact (polynomial identity mod Φ_h).
  - The float-pruned enumeration is safe: every vector is re-normed exactly, and the PARI qfminim count shows none is
    missed.
  - The decomposition into the indicators 1_{S_d} is unique.
- **Lines A and B.** Re-derived by hand: |π(x)|² = |x|² − 2(x·z)/h.
- **Numbers.** 7275 ✓, 64746 ✓, A1²⁴ ✓, and the classes 2B–13B ✓.

**Fixes needed.**
- (a) The status line cites "Duke 1998, Thm 4.1". It should be Borcherds, *Invent. Math.* 132 (1998), Thm 4.1.
- (b) **Error in (3′).** "coincides with that of a fixed-point-free Conway element of Frame shape h^k/1^k" is false
  for h = 25: Co₀ has no element of order 25. Add "for h ≠ 25".
- (c) The remark "For prime h … W(R) realizes the torus on D" is unjustified. W(R) does not fix ρ, so it does not act
  on M, and the argument really needs rationality. Mark it SKETCH.
- (d) The note's §1 divisor-lattice sentences are wrong.
  - "E8³ gives 2, 3, 5, 6, 10": the actual Coxeter numbers are h_d = 2, 3, 5, **5, 9**, 14 for d = 2, 3, 5, 6, 10, 15.
  - "At odd d … drop": E6⁴ drops at d = 4.

## 2. THM-4759 (1)–(2): CONFIRMED

- **Derivations.** Re-derived the Kostant/V_Q trace argument and the Poisson computation, checked anti-self-duality
  for all five families by hand, and checked the glue lemma (alcove maps permute the affine simple roots, so they fix
  ρ/h).
- **Independent enumeration** (`per_component_check.py`): D5, D6, E6, E7 and E8 pass in both the shifted and the
  twined form.
- **Classical references to cite:**
  - Lepowsky, *Ann. ÉNS* 12 (1979), §17 and Table 4.
  - Kac–Peterson, "112 constructions" (1985), eqs. (11.6)–(11.7).
  - Kac's book, Prop. 10.9.
- **Typing of (2).** It follows in one line from Kac–Peterson, Kostant and the glue lemma. It is new only as the
  identification with CDH's T^X, which is not found in print. Lam–Miyamoto (*IMRN* 2006, Art. 35967) was not checked.

## 3. THM-4759 (3): CONFIRMED

The trace of the principal element on a simple Lie algebra is Tr(c|𝔥) = −1. Proof: dim g_(j) = ℓ + mult(j) for
j ≢ 0. So the weight-1 trace on V_N is 24 + Σ(−1 − ℓ_R) = −k, and on V_Λ it is 24.

## 4. THM-4759 (4): CONFIRMED

Checked against CDH:
- Eq. (2.23): T^X = q⁻¹ − d^X + O(q), where d^X is the number of components.
- Prop. 2.1 is proved "by inspection".
- Remark 2.2 asks for a conceptual proof.
- All 23 lambencies match Table 2, and its Monster classes match the THM.

## 5. THM-4759 (5), (5′): CONFIRMED WITH CAVEAT

- **Exponent law.** Every hand formula holds for A_n and D_n with n ≤ 300 and for E. The equality cases are exactly
  as stated.
- **The other steps are correct:** the twisted-module bound, the norm quantization in (2d/h)Z, and attainment from
  the pole at 1/d.
- **Uniqueness of the minimal vector for D_n (n even)** can be PROVED, not just verified: round each coordinate
  j/(n−1) of ρ/(n−1) to the nearest integer.
- **Component rule.**
  - It is correct as a Frame-shape statement.
  - A single component's eta quotient picks up e^{−2πiℓ/24} under W_e; measured e^{−iπ/3} for D4.
  - The full T^X is invariant exactly on S_X for all 23 types (`we_niemeier_check.py`). FIXED at HEAD.
- **Gap in the Leech-side step** (line 536). The proof shows N + dρ/h and Λ + d′c agree below q¹, then concludes about
  Λ + dc. Add one line: Θ^{[c]} is rational, so θ_{Λ+d′c} = θ_{Λ+dc}.

## 6. THM-4759 (7) and (10)

**(7): CONFIRMED.**
- Γ₀(e)-invariance holds because rationality follows from Θ^{[c]}.
- The cusp analysis is right, F − G is W_e-odd (hence 0), and n₂ = 24 follows.

**(10): a real gap in the registered version, now CLOSED.**
- The gap: the proof needed rationality of the twining and never established it.
- Closed by THM-4760 (2)–(3). Also closed independently:
  - `units_hplus1.py` shows every unit multiple u·c″, for all 23 types, has profile {2 − 2/e: 1, 2: 24}.
  - A Γ₁(e) principal-part argument then forces Z_v to equal the rational series θ_{Λ+v}(eτ)/Δ(eτ).

**Monster names.** All coefficients and all cited OEIS A-numbers are verified (6B A007255, 8A A007265, 9A A007266,
10D A058100, 14C A058504, 15C A058510, 26B A058597, …). The groups 6+6, 8+, 9+, 10+10, 14+14, 15+15 and 26+26 are
monstrous: 6B, 8A, 9A, 10D, 14C, 15C, 26B.

**Supersingular remark.** Sound as an implication: (10) only re-derives that X₀(p)/W_p has genus zero; Ogg's theorem
does the rest.

**(h+1)-construction, lattice side: KNOWN.** Chenevier–Lannes Thm 3.4.2.10 proves vois_{h+1} = s_ρ(vois_h) ≅ Leech,
with counts; Thompson did D24. HEAD now cites this. The Fricke-twining statement itself was not found in the
literature.

**Overclaims remaining at HEAD.**
- Line 204: "These are exactly the types with a Fricke point" contradicts (10).
- Line 373 (and note line 271): "Each is the (h+1)-point of the corresponding Fricke descent" is proved only at order 3.

## 7. THM-4759 (8): CONFIRMED

- **Subsystems.** Checked on the extended-diagram labels and by Kostant height counts 2Σ⌊m/p⌋: E8 gives 112, 72, 40;
  E7 gives 56; E6 gives 18; D_n gives n(n−2).
- **Neighbours and divisor table.** All eleven neighbours match, and the divisor table is consistent.
- **Simplification.** L_e = vois_p(N; ρ_N), so the inclusion R(L_e) ⊇ X_p is immediate.
- **Prior art for p = 2:** Chenevier–Lannes §III.4.4 (Prop. 4.4.2: h(vois₂) = h/2 − ι/8 + 2; ψ(R) = R/2 iff
  ι(R) = 24) and Borcherds' thesis Prop. 4.5.3. HEAD now cites both.
- Line 265 should tag "the Collatz map 3n + 1" as ANALOGY.

## 8. M24 point: CONFIRMED

- **Distribution.** Derived by hand and by a Golay-code count: 5: 24, 4: 759, 3: 6072, 2: 21528, 1: 43056, 0: 53682.
- **Orbit sizes.** 33,965,568,000 + 452,088,000 = 34,417,656,000, exactly.
- **Classical.** Conway, SPLAG Ch. 10, Theorem 29 (transitive on vectors of type 6₂₂ and of type 6₃₂) and Table 10.4
  (·6₃₂ = M24, ·6₂₂ = PSU₆(2).2). The status line 53 should change from VERIFIED to PROVED (CITED).

## 9. Carnahan: CONFIRMED WITH CAVEAT

- **Hypotheses.** The main result is Thm 4.5, with hypotheses fixed-point-free and E₀ > 0; Cor. 4.6 is the class
  bijection. These match the THM.
- **E₀.** Carnahan's printed formula lacks the factor 1/24; the THM uses the correct one. For every Coxeter Frame
  shape E₀ = 1/h.
- **The appendix** lists 38 anomaly-free classes, including the 22, and the names agree with CDH Table 2.
- **Fixes.** Delete the "summarising tool" caveat (line 444) and cite Thm 4.5.

## 10. HYP-9368: WRONG as registered

**The problem.** The definition of a strict Fricke point ignored the unit multiples uv. The censuses were at orders 3,
4 and 6, the only orders with units {±1}, so they could not see this.

**Counterexample at order 5.** `order5_strict_census.py` sampled 231,699 norm-40 vectors, giving 45,890 strict
points in two classes. `verify_order5_example.py` confirms the example exactly with PARI.
- **39,322 points:** these have the invariant of c″(A3⁸).
- **6,568 points:** these have n₂ = 23 and Λ + 2v = {6/5: 1, …}, so the twining is irrational and they are
  inequivalent to c″(A3⁸). Their back-nudge y/4 is not a deep hole. The session has since identified this class as
  the A2¹² shell-1 point.

**Consequences.**
- This refutes (A) as stated.
- It refutes "THM-4759 (10) proves that then Z_v − 24 is the Fricke Hauptmodul".
- It refutes step 1 of the prime-order sketch, which was labelled PROVED.

**Fixed at HEAD (19e71373ee):** the definition now requires the Hauptmodul property.

**Remaining fixes.**
- Line 94 is now circular and should be dropped.
- Order 3 can be upgraded to PROVED (CITED).

**Upgrade available.** At prime levels (B′) is PROVED: a Hauptmodul twining is rational, so dim V′₁ = 24p forces
h = p − 1. Hence no translation twines by T_{p+} for p = 2, 29, 41, 59, 71.

## THM-4760 (scope addition)

**T1. Criterion: CONFIRMED.**
- Re-derived: F(u,w)(τ+1) = e^{−πiu²}F(u,u+w), the Poisson law, and F(0,v)|γ = εF(cv,dv).
- Γ₀(e)-invariance: trivial character (det e² is a square) and no phase. A Jacobi-form check agrees.
- The Φ_c | P reduction is correct.
- The "only if" is correct; a simpler route is that Γ₀(e)-invariance gives Z_v = σ_d(Z_v) directly.
- The α_d certificate and the Sturm lemma via the valence formula on Γ(N) are correct.
- **Note:** condition (ii) is redundant, because G − Z_v is a constant that is odd under W_e.

**T2. Spot checks: CONFIRMED.** The re-run of `fricke_cusp_certificate.py` and the independent
`cusp_criterion_indep.py` agree; the PARI count matched the float enumeration in all 12 runs.
- d6⁴ shell 1: Hauptmodul, q⁻¹ + 24 + 14q + 36q² (12H).
- e8³ shell 1: irrational, units ±1 mod 8; Λ + 3v has minimum 15/8.
- d24 shell 1: irrational, units ±1 mod 12.
- Shell 0 of d24, a24 and e8³: Hauptmodul.

**T3. Carriers: CONFIRMED.**
- The E7 carriers in D10E7² are Hauptmoduln (20F).
- The D10 carrier is irrational (units ±1 mod 10).
- The E8 carrier in D16E8 is a Hauptmodul; the D16 carrier is irrational.
- Shell-1 completeness holds, because the glue group acts transitively on each carrier's |N/Q| vectors.

**T4. Ligozat orders.** All seven, plus 6B, 10D, 15C, 14C and 26B, are confirmed independently (`ligozat_check.py`).

**T5. Is THM-4759 (10) complete?** YES.

**Typing.**
- THM-4760 (3′) is labelled FINITE-EXACT but its large shells were sampled; it should read VERIFIED (sampled).
- Add PARI count cross-checks to `leech_h_plus_one`, `fricke_h_plus_one` and `fricke_cusp_certificate`.

## Corrections still needed at HEAD (92a4c0a0af)

1. **THM-4758, status line:** Duke 1998 → *Invent. Math.* 132 (1998), Thm 4.1.
2. **THM-4758 (3′):** add "for h ≠ 25".
3. **THM-4758, W(R)-torus remark:** mark SKETCH.
4. **THM-4759, line 53:** M24 remark → PROVED (CITED SPLAG Thm 29 / Table 10.4).
5. **THM-4759, line 204:** "…types whose deep hole has a Fricke multiple".
6. **THM-4759, line 373 and note, line 271:** "has the same neighbour and theta data as" (equivalence proved only at
   order 3).
7. **THM-4759, line 536:** add the rationality line for θ_{Λ+d′c} = θ_{Λ+dc}.
8. **THM-4759, line 444:** delete the caveat; "38 classes, including these 22"; cite Thm 4.5.
9. **THM-4759, line 265:** tag ANALOGY.
10. **Note, lines 97–98:** use the rule h_d = d − 1 iff h/d is a common bad prime ∥ h.
11. **HYP-9368:** drop line 94; set order 3 to PROVED; add (B′) at prime levels as PROVED.
12. **THM-4760 (3′):** VERIFIED (sampled).

## Prior art

- **Chenevier–Lannes**, *Automorphic Forms and Even Unimodular Lattices*, Ergebnisse 69 (2019), arXiv:1409.7616:
  Thm 3.4.2.10, Prop. 3.4.2.16, §3.4.4, Prop. 3.4.1.1. https://arxiv.org/abs/1409.7616
- **Borcherds**, thesis, arXiv:math/9911195: Lemma 2.2.7, Prop. 4.5.3.
- **Lepowsky**, *Ann. ÉNS* 12 (1979), §§16–17. http://www.numdam.org/item/ASENS_1979_4_12_2_169_0/
- **Kac–Peterson** (1985), eqs. (11.6)–(11.7). https://math.mit.edu/~kac/not-easily-available/112.pdf
- **Kac**, LNM 848 (1981), Thm 1.
- **Conway**, SPLAG Ch. 10, Thm 29 and Table 10.4.
- **Möller–Scheithauer**, arXiv:1910.04947, Cor. 5.12–5.13: confirms the root-count citation in (10).
- **Unchecked:** SPLAG Ch. 4 eq. (144), the A1²⁴ hole theta series, possibly already in closed form (equal to
  θ₂(τ/2)²⁴/4096 + 48Δ); and Lam–Miyamoto 2006.

## Recommended MISTAKE-662 addition

**Audit (a): a Fricke-point definition blind to Galois conjugates.**
- HYP-9368 (A) constrained only v and its non-unit multiples.
- Its censuses were at orders 3, 4 and 6, where the only units are ±1.
- At order 5, 6,568 of 45,890 sampled strict points have n₂ = 23 and Λ + 2v of minimum 6/5. Their twining is
  irrational and they are inequivalent to c″(A3⁸).
- This refuted (A) as stated, the claim that THM-4759 (10) applies to every strict point, and step 1 of the
  prime-order sketch.
- Lesson: for order n with φ(n) > 2, test every generator of ⟨v⟩, or the rationality of the twining.

**Audit (b): smaller canon errors.**
- THM-4758 (3′) invokes a Conway element of Frame shape 25/1, which does not exist.
- The Borcherds theorem is cited to the wrong journal ("Duke 1998").
- The note's divisor-lattice sentences ("odd d"; "E8³ gives … 6, 10") are stale.
- "Exactly the types with a Fricke point" contradicts (10).
- The (h+1)-point identification of the shallow multiples is stated as fact but proved only at order 3.
- The M24 orbit is left as VERIFIED although Conway's theorem makes it citable.
- The component rule was stated as invariance of functions; per component it holds only up to e^{−2πiℓ/24} (fixed in
  effcde315f).
