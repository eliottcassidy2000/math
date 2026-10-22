# Independent audit of opus-2026-10-21-S52 (2026-10-22)

The auditor was an independent subagent. It used its own code: the five scripts in this directory, each ending in
ALL CHECKS PASSED. It wrote no files outside this directory. Subagents cannot write report files, so the parent
session saved this report from the auditor's final message without changing its content.

## Verdict

**The mathematics is correct.** Every numerical claim was reproduced with independent code, usually by a different
method. Two checks go beyond the session's own:
- **Leech roots.** A complete enumeration by height shows that ν = (½²⁵ | 5/2) lies in the closed Leech chamber and
  that exactly 25 Leech roots are orthogonal to it. This does not use Conway's pairwise argument.
- **The quotient.** ν^⊥/ν = N(A24) is computed directly, without the Niemeier–Venkov step.

**The presentation needs correcting before close-out.**
- THM-4752 (b) calls U(T) "the Lorentzian distance" to the nearest cusp. That is false.
- Statement (c) is already published in Borcherds' table of norm-0 vectors of II_{25,1}. It lists exactly 2ν for
  A24, with height 25 and 25 orthogonal simple roots. The thesis gives z·ρ = −h and ρ = (ρ_N, h, h+1). Only the
  tournament reading is new.
- The note overclaims about THM-869/872 (whose question is about n = 24) and about THM-484 Remark 1.
- The note misreports the Nebe-48 PDF, gets the Smale sentence wrong, and misstates THM-4748 and HYP-9359 in §6.

**Recommended typing for THM-4752.** The dictionary is PROVED and elementary.
- The lattice side is KNOWN: Conway 1983; Borcherds' thesis §§2.2–2.5; Borcherds' norm-0 table.
- The tournament side uses only elementary or classical inputs: Kendall–Babington Smith; the rearrangement
  inequality (Kano–Sakamoto); Landau; Watson.

Apart from M5 and R1, the PDF referee verdicts stand. The e8.pdf "trace lock" refutation is confirmed by an argument
that needs no SEA.

## Mandatory corrections

**M1. THM-4752 (b) is false as stated.** It says "So U(T) is the Lorentzian distance from v_T to the nearest of the
25! Leech cusps σ(ρ)".
- The difference v_T − ρ_π = (s − π | 0) is spacelike, with |v_T − ρ_π|² = 2(BL_T(π) − c3).
- So the squared distance to the nearest cusp is |s↑ − ι|² = 2(U − c3), not U.
- Counterexample (A4): a random p = ½ tournament has U = 956, but its squared distance is 710 (c3 = 601). The two
  agree only when U = 2c3, as for the regular tournament.
- What is true: U is the pairing −v_T·ρ_s, a horospherical height. The "nearest cusp" identification still holds,
  because both quantities are monotone in BL.
- Fix: U = min_σ(−v_T·σ(ρ)). The minimising cusp is also the nearest one, at squared distance |s↑ − ι|² = 2(U − c3).

**M2. Statement (c) is already published.**
- Borcherds' norm-0 table, https://math.berkeley.edu/~reb/lattices/norm0, gives for each vector its height, then the
  numbers of simple roots at inner products 0, 1, 2, then coordinates.
  - The Leech line is 2·(70 | 24, …, 0).
  - The A24 line is `25 25 5 5 a24 2 5 1 … 1`: that is 2ν, with height 25 and 25 / 5 / 5 simple roots at inner
    products 0 / −1 / −2.
- The thesis (arXiv math/9911195) contains:
  - Lemma 2.4.2: ((ρ, h, h+1), z) = −h;
  - Lemma 2.5.1: h = −(z, w);
  - Lemma 2.2.4 (Venkov): ρ_N² = 2h(h+1).
- So (d)'s transitive class, (12, …, −12) = ρ_A24 with 1300 = 2·25·26, is Borcherds' decomposition read in
  tournament coordinates.
- Fix:
  - type (c) as KNOWN and reproduced;
  - requalify "No prior instance";
  - replace "written in tournament coordinates" with "in Conway's standard coordinates, read as tournament
    coordinates".

**M3. The THM-869/872 overclaim.** Their question concerns n = 24. At n = 25, centred scores in A24 are already
THM-872's odd-n case. The only new step is realising N(A24) canonically as ν^⊥/ν; this is canonical because ν is the
unique primitive future null vector fixed by S₂₅. Say this, and say that the n = 24 question remains open.

**M4. THM-484 Remark 1 is misread.** That remark concerns the THM-472 spectral identity Σμ_j² = C(n,2), not score
vectors. Delete "settles" and rewrite the related-list entry.

**M5. The Nebe-48 PDF is misreported.** The PDF claims "at least 6,648,079,969,440" facets, and its Theorem 4.1 says
this is a lower bound; norm-12 relevance is left open. Fix: "at least …; correct as a lower bound".

**M6. The Smale sentence is wrong.** Θ_n *is* the group of h-cobordism classes. Smale's theorem identifies Θ_n with
S_n, the diffeomorphism classes of homotopy spheres, for n ≥ 5 (Milnor 2011, footnote 2, p. 805).

**M7. §6 misstates the merged results.**
- (a) "2^k − 1 needs three back arcs" holds only for k ≥ 5: 15 = 2⁴ − 1 has a realization with feedback arc number 2.
- (b) The uniqueness claim should read "the fan and its converse … among strong classes with H < 2100", typed
  FINITE-EXACT.
- (c) The arithmetic support should be typed HYP-9359 (NUMERICAL) and say "almost all".
- (d) N(n,h) is only eventually a polynomial, for n ≥ ⌊3 log₅ h⌋.

**M8. Housekeeping.**
- The hypotheses INDEX still lists THM-4752 / THM-4753 / HYP-9363 as reserved stubs.
- The results INDEX has no entry for the note.

## Recommended improvements

- **R1. The Cohn–Elkies value.** Interpolating Table 1 of Afkhami-Jeddi–Cohn–Hartman–de Laat–Tajdini
  (arXiv:2006.02560) gives an LP bound at d = 48 of about 4.5e−7, about 20× the extremal density. The PDF's 3.01e−8
  and its "23 % gap" are inconsistent with this. (NUMERICAL; the same interpolation reproduces the Leech value at
  d = 24 to within 2e−4.)
- **R2. An SEA-free refutation.** [p+1∓5]P ≠ O for a random point P proves a_p ≠ ±5 for all ten curves.
  - The q = 5/4, 15-digit order is certified: N = 37824·c, with c a probable prime greater than 4√p.
  - [N]P = O holds for all ten session values.
- **R3. The span.** The roots span ν^⊥ over R, but over Z they span an index-5 sublattice modulo ν. Say that ν is
  future-pointing. Optionally, compute directly: A24 plus glue [10], with glue minima 4, 6, 6, 4.
- **R4. Canonicity.** Every S₂₅-invariant w with w·ν = −25 gives the class s − 12. The score set is the level −25
  slice intersected with x₂₅ = 70.
- **R5. Watson.** At n = 1 the vector is zero; say "nonzero null vector".
- **R6. Citations and typing.**
  - Cite Kano–Sakamoto precisely: Networks 13 (1983) 143–151, doi 10.1002/net.3230130109.
  - Type U = c3 + ½|s↑ − ι|² as elementary, possibly folklore.
  - Monsuur (EJOR 2005) is unread and paywalled.
- **R7. Context for (e).** Π(x_i + x_j) = s_δ, so score vectors are the weights of V_δ, and tournament counts are
  Kostka numbers (Gessel 1979; repo camion_busch §6.1).
- **R8. The Leech ζ wording.**
  - Say L(Δ,13) "is present (non-critical value)" rather than that it "cannot be removed".
  - The shell-sum check is circular, since it uses the same σ₁₁ − τ decomposition. The independent θ-continuation and
    E4³ − 720Δ routes agree to 1e−40.
- **R9. The im J wording.**
  - That the im J orders are "exactly the numerators" is automatic (the coefficient is the reciprocal of B₂ₛ/(−4s)).
  - "240 three ways" is two ways: the roots and the kissing configuration are the same vectors.
  - The Atiyah attribution is unverified.
- **R10. HYP-9273.** "Holds through q = 31" is its existing status from S33/S34; credit those sessions.
- **R11. The projective-plane test has little power.** The factor-2 band contains 98.6 % of all odd d.
  - A regression of log S(d) on a cubic in log d plus the indicators [3|d], [9|d], [5|d] gives coefficients
    0.233 / 0.071 / 0.050, with residual s.d. 0.107.
  - The worst PG size is 133, at z = +2.74: an excess, not a depletion. The 0.67 ratio at 31 is the 3-divisibility
    effect.
- **R12. The measure caveat (Session27).** The degree law concerns N(n,h) as n → ∞ and says nothing about strong-class
  counts such as S(31)/S(33).
- **R13. MISTAKE entry.** Record M1–M7, especially the missed Borcherds table.
- **R14. The Nebe-48 PDF's reference list looks garbled.** Not checked further.

## Confirmed items (independent code; all pass)

**THM-4752**
- All vector identities: ρ, ν, r*, and the Ã24 Gram matrix.
- On twelve 25-vertex tournaments and 720 random rankings:
  - v_T² = −2c3;
  - BL = −v_T·ρ_π;
  - U equals the Landau slack sum, which equals c3 + ½|s↑ − ι|²;
  - |s − 12|² = 1300 − 2c3.
- Exhaustive over all tournaments and all rankings for n ≤ 6; 1500 tournaments × 5040 rankings at n = 7.
- The identity holds at orders 3–40.
- The Leech-root enumeration, with −70 ≤ t ≤ 70 (33,678 roots), is complete for r·ν ≥ −22:
  - no root has r·ν > 0, and only the 25 cycle roots have r·ν = 0;
  - all pairwise products are ≤ 0 for t ≤ 30;
  - **observation:** #{r : r·ν = −k} = 5·p(k−1), exactly for 1 ≤ k ≤ 22.
- ν^⊥/ν = N(A24), computed directly.
- Watson's uniqueness checked to 10⁷.

**lattice_cusps**
- Vinberg rerun with the opposite sign convention: 10 and 19 roots, the cusps, ρ² = −1240, h = 30.
- Siegel mass checked with sensitivity: doubling any |G₁| or |G₂| breaks the identity, and M₁₂ is forced for A₂¹².
- The E8 holes, and Milnor's pair through norm 8.

**pdf_referee**
- The extremal rank-48 theta series; the norm-12 shell is a₆ = 437,824,977,408,000.
- The density, and the PDF's covering-radius bound.
- L-values:
  - Λ(k) agrees to 45 digits across two split points;
  - the Manin/Zagier period ratios hold;
  - L(Δ,5) = 0.666709188434003643826…
- ζ_Λ(10) = −1.975523970070365355620842, by θ-continuation alone.
- ζ_Λ(26) by three routes.
- McLaughlin: the reductions are right, and a constructed Leech lattice gives 552 = 2·276.

**bitdepth**
- a_p ≠ ±5 is proved for all ten curves.
- One curve's order is certified.
- The committed .out (b8ac6c0b13) is complete and matches the scratch run.

**numerology**
- H(P_q) for q ≤ 27 by an independent DP; H(P_31) comes from the repo (audited in S34).
- Aut orders 3, 21, 55, 171, 253, 1053, 465.
- The cube pattern holds only at q = 3, 7, 11. H/|Aff| ≡ 1 (mod 8) at all seven q, while H/|Aut| ≡ 3 at q = 27.
- The bP orders match Milnor 2011 Theorem 5, and |Θ₂₃| = 48·|bP₂₄|.
- Also confirmed: Adams' im J orders, Mp₂(Z)^ab = Z/24, j, and the Heegner values.
- The strong inventory has 723,544 classes. The logic is sound: Moon's strong-deletion argument plus insertion
  monotonicity make it complete.

**Concurrent work.** Session27 is credited accurately and is not claimed as new.

## Prior-art findings

1. Borcherds, norm-0 vectors of II_{25,1}: https://math.berkeley.edu/~reb/lattices/norm0 (format at
   https://math.berkeley.edu/~reb/lattices/index.html).
2. Borcherds, *The Leech lattice and other lattices* (thesis), arXiv math/9911195: Lemmas 2.2.4, 2.4.2, 2.5.1.
3. No web source on tournament score vectors in II_{25,1}. In the repo:
   - `ogg_triangular_chowla_20260927.md:283–284` names (0, …, 24; 70) "the Leech vector in II_(25,1)" (NUMEROLOGY);
   - `strong_tournament_balance_golden_rigidity_20261014.md:156`.
   - For classical tournament/Weyl context, see Gessel 1979 and Hausel's slides
     (https://hausel.pages.ist.ac.at/wp-content/uploads/sites/229/2023/12/ViennaCombinatorics23.pdf).
4. Backward length:
   - Kano–Sakamoto, Networks 13 (1983) 143–151;
   - Monsuur, EJOR 164 (2005) 778–784 (paywalled; priority for the closed form unresolved).
5. Cohn–Elkies LP numerics: arXiv:2006.02560, Table 1.
6. Milnor 2011: https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/milnor2011.pdf (Theorem 5 p. 806;
   footnote 2 p. 805).

## Scripts (each has a .out and ends ALL CHECKS PASSED)

| Script | Checks | Runtime |
|---|---|---|
| `a_lorentzian_audit.py` | THM-4752 (a)–(f), the Leech-root enumeration, ν^⊥/ν, canonicity, the M1 counterexample | 13 s |
| `b_vinberg_mass_audit.py` | Vinberg, Siegel mass with sensitivity, E8 holes, Milnor's pair | 4 s |
| `c_pdf_numerics_audit.py` | Rank-48 theta series, LP interpolation, L-values and period relations, Epstein zeta values, McLaughlin | 5 s |
| `d_bitdepth_audit.py` | Its own elliptic-curve arithmetic: a_p ≠ ±5, an order certificate, a PARI cross-check | 2 s |
| `e_numerology_audit.py` | Paley H and Aut, the bP / Adams / η / j / Heegner values, projective-plane fairness | 150 s |

`e_numerology_audit.py` regenerates the strong spectrum with the S50 audit tools and writes
`e_strong_spectrum_S.json`.
