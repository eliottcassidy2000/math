# Audit H: THM-4611, THM-4609 (5)–(6), the HYP-9244 update and the results note (3.4–6b)

Independent adversarial audit, run 2026-10-08 as a session subagent. The harness refused a report file, so the main session saved the subagent's final message here verbatim.

All my scripts and outputs are in this folder. None of them imports the author's modules. I edited no audited file, committed nothing and downloaded nothing.

Versions audited (file times): THM-4611 and the results note at 11:49, HYP-9244 at 11:50, THM-4609 at 10:15. They grew during the audit, so earlier wording problems I found and that are now fixed are noted as such.

**Bottom line.** The mathematics holds, and every exact certificate I re-ran passed:
- the 5 fixed-length certificates;
- the 6 adaptive certificates for the {±1} maps, including (1,6,4,31,1);
- the adaptive Z_4 certificate;
- the 3 new certificates for audit F's maps;
- all 18 census certificates.

The claims that no fixed-length certificate exists were backed only by a float optimizer. I proved them exactly with dual certificates, for (1,2,3,7,1) at k = 2, 3 and for all six {±1} maps at every k ≤ 4. There are no major findings. The fixes are typing, wording and stale counts, plus one false sentence: that the twelve rank-3 census maps *need* adaptive lengths.

## 1. Verdicts

| File | Verdict |
|---|---|
| THM-4611 (statements 1, 2, 2′, 3, 3′, 4, 5, 6) | **CORRECT WITH FIXES** (minor) |
| THM-4609: statement 5, the non-translation rows of statement 6, the "translation-only = all m_i congruent mod d" setting | **CORRECT** (nits) |
| HYP-9244: status and "Update 2026-10-08" | **CORRECT WITH FIXES** |
| Results note, sections 3.4, 3.5, 4, 5, 6b | **CORRECT WITH FIXES** |

## 2. Findings

### THM-4611

**1. (minor) Statement 3′ says "None has a fixed-length certificate for k ≤ 4" and types it FINITE-EXACT, but the author's evidence is float.**
- The evidence is the optimizer's negative "best margin" in `block_census_fixed_k4.out`.
- The claim is **true**; I now prove it exactly:
  - k = 1, all six maps: some block has exact rank 1.
  - k = 2 for (1,6,11,11,4): some block has rank 2.
  - Every other pair (map, k ≤ 4), plus (1,2,3,7,1) at k = 2, 3: a dual certificate. It consists of rational z_l and weights c_l ≥ 0 such that Σ c_l [A_l − 2(A_l z_l)(A_l z_l)^T/(z_l^T A_l z_l)] is negative definite.
  - Why that suffices: if a form P = Q⁻¹ balances A_l, then ⟨B_l, P⟩ > 0. A negative definite sum therefore rules out every form. Six cuts were enough in each case.
- Evidence: `f_dual_part1.out`, `f_dual_part2.out`, `f_dual_sixth.out`. The certificates are stored in `dual_certificates.json`, and `l_verify_certs.py` re-checks them from scratch: it recomputes the block set, then confirms membership and negative definiteness.
- Fix: cite these certificates, and type the quoted optimal margins (−0.045 … −0.107) as NUMERICAL.
- I also confirmed the Reading's optima for (1,2,3,7,1) by bracketing; the upper end of each bracket is proved exactly by a dual certificate.

  | k | Optimum lies in |
  |---|---|
  | 2 | (−0.27465, −0.27462] |
  | 3 | (−0.09107, −0.09106] |
  | 4 | (+0.01354, +0.01438] |

**2. (minor) Statement 5 says "The twelve rank-3 maps need adaptive lengths 3–5". This is false for (1,1,28,11,7).**
- That map has a fixed k = 4 certificate, Q = [[23,−4,−8],[−4,20,−4],[−8,−4,24]], found by census v1.
- I verified it exactly on all 156,556 nonzero blocks; its margin is 0.0053 (`c_census_v1_fixed.out`).
- Census v2 tries only fixed k = 2 before going adaptive, so fixed k ≥ 3 was never tried.
- Fix: write "were certified with adaptive lengths 3–5 (fixed k ≥ 3 not tried)".

**3. (minor) Statement 6 and results note 6b say that for Z_7 (1,2,3,5,1,1,1) "lengths 2–4 do not suffice, leaf optimum −0.036".**
- `auditF_maps_cert_Z7.out` contains only the pilot line; the −0.036 is not recorded.
- One pilot-driven partition failing does not show that no adaptive 2–4 certificate exists. Partitions are not monotone: a lift's block can be worse than its parent's.
- Fix: write "the pilot construction with lengths 2–4 did not certify it" and save that run, or drop the sentence.

**4. (minor) The Reading says multiplier-only certificates "fail for every map tested". HYP-9244 item 3 and results note 6b repeat this.**
- `adversarial_cert.py` tests one form per map (the statement-3 form) for k ≤ 8.
- A positive adversarial value for one Q does not exclude every Q.
- Fix: restate it as HEURISTIC/NUMERICAL ("with the statement-3 forms the adversarial bound stays positive for k ≤ 8"), or optimize over Q.

**5. (minor) The status claims property (i) "for every certified map (`property_i_check.py`)", but that script covers 13 of the 27 maps.**
- Missing from it:
  - census maps with dependent multipliers, where (i) is not automatic: (1,6,8,21,3), (1,3,9,7,11), (1,2,3,14,8);
  - census map (1,11,1,21,2);
  - the six rank-4 census maps;
  - Z_4 (1,3,5,7).
- `block_census.py` never checks (i).
- I checked (i) directly for all 18 census maps, for Z_4 and for the statement-6 maps; it holds in every case.
- Fix: extend `property_i_check.py` and add the check to the census.

**6. (nit; already fixed in the current text) The 11:21 version's escape-word reason was false.**
- It read "(M, e) ≠ (1, 0) loses one digit of agreement per identity step".
- Counterexample on Z_5 (1,2,3,7,1) with standard constants: the state (2/7, 0) is reached from (1, 1) by the digits 0, 3, 0, 2. Digit 0 then fixes it, so the identity coupling persists forever along 000… with constant agreement (`m_identity_persist.out`).
- The current argument ("otherwise Mv + e = v on a cylinder", via the Terras property) is correct; keep it.

**7. (nit) "9–20% of the level-3 states needed longer blocks" (THM-4611 Reading; "need" also in note 6b).**
- These counts come from the pilot refinement rule, not from necessity.
- With the final forms, an exact greedy refinement refines only 304, 393, 682, 677, 360 and 632 of the 6,250 level-3 states (4.9–10.9%), and 1.0–9.3% of their lifts.
- Fix: write "were refined".

**8. (nit) Wording: "whatever happens at the hidden digits" and "no independence enter".**
- The proof does use i.i.d. uniform fresh digits (the Terras property). What it does not use is mixing of the hidden state.
- M is a function of the debt, so the only free hidden datum is e.
- Suggest: "uniformly in the hidden state at block starts".

**9. (nit) The proof of 3 describes only the rank-3 exact test.** For rank 4 the code uses exact integer 3×3 and 4×4 minors; say so.

**10. (nit) "semi-decidable question: does every map have a certificate?"**
- Existence of a certificate is semi-decidable *per map*: each finite check is a semialgebraic feasibility problem.
- The universal question is not semi-decidable. The same phrasing appears in HYP-9244.

**11. (nit) Status line.**
- The every-k Z_4 claim is PROVED by the identity-run computation; I confirmed it to k = 8. The status lists only "k ≤ 5" as FINITE-EXACT.
- The adaptive Z_4 certificate is not typed in the status.
- "Not yet independently audited" can now be updated.

### HYP-9244

**12. (minor) "Rank-4 maps with coupling group the squares mod 7 or 11 have explicit forms (8I − J for d = 7)" overgeneralizes.**
- Forms exist only for these unit-set orbits:
  - d = 7: {0,1,2};
  - d = 11: {0..6} and {0..5,7}.
- For d = 7, the other rank-4 orbit {0,1,3} has a squares-group coupling whose covariance has rank 2, so no one-step form exists there (`g_hierarchy.out`).
- For d = 11, the orbits {0..5,8} and {0..4,6,7} are undecided.
- Fix: name the unit sets.

**13. (nit) "sharp for the standard form" (here and in THM-4609 (5)).**
- For the odd-order row, sharpness at rank 4 needs a non-trivial odd-order subgroup, i.e. d − 1 not a power of 2.
- For Fermat primes (5, 17, 257, 65537) that class is just the translations, and Q = I already works from rank 4.
- Fix: scope the claim to the checked d.

**14. (nit) Item 3:** see findings 4 and 10. "The sampled maps of that class have no fixed-length certificate for k ≤ 4" is now exactly true; cite audit H.

### THM-4609

**15. (nit) "Checked exactly for every unit set with d = 7, 11, 13".**
- `rank_hierarchy.py` checks ranks 3–8 only, and only unit sets containing 0. The second restriction is a valid shortcut: the coupling class is invariant under the affine group.
- I checked every unit set and every rank 3..d−1 for d = 7 and 11.
- Fix: write "ranks ≤ 8, unit sets up to translation".

**16. (nit) THM-4609's Reading is stale.**
- Its "still open" list (Z_5 (1,2,3,7,1) and d = 4 (1,3,5,7)) is now covered by THM-4611.
- "odd-order groups need rank 5" holds only for the standard form; the squares-group forms work at rank 4.
- The status still says statement 5 is unaudited.

### Results note

**17. (minor) Section 6b and surrounding text are stale.**
- 6b says "all five sampled {±1} maps", and its table has five rows. Six were certified; (1,6,4,31,1) is missing.
- §6 open items 3 and 4 do not reflect THM-4611.
- The header omits 3.5 (the non-translation rows) and 6b from the sections added after audit G.

**18. (nit) §5 "The sketch route … works for every μ".** It is a sketch; write "should work".

**19. (nit) 6b "a sticky reflection, often of rank one or two".**
- On Z_5 the reflection covariance rank is always ≤ 2, and ≥ 1 under (i). So "always".
- "Every subset is reflection-symmetric" is special to Z_5: 28 of 128 subsets fail on Z_7, and 1,364 of 2,048 on Z_11.

**20. (nit) §3.4 repeats "checks every unit set of each rank"** (see finding 15).

### Prior art

**21. (minor) THM-4611 should credit the classical tools behind statements 1, 2 and 2′.**
- State-dependent (k-step) Lyapunov drift: Malyshev–Men'shikov (Trudy MMO 39, 1979); Fayolle–Malyshev–Menshikov (CUP 1995); Meyn–Tweedie (Ann. Appl. Probab. 4, 1994).
- The |x|^(−α) many-dimensional method: Lamperti (1960); Menshikov–Popov–Wade (CUP 2016); Georgiou–Menshikov–Mijatović–Wade (Adv. Appl. Probab. 48A, 2016).
- The trace criterion: Peres–Popov–Sousi (2013), Thm 1.3. Already cited.
- Walks with internal states: Comets–Menshikov–Popov (Ann. Probab. 26, 1998); Georgiou–Wade (SPA 124, 2014, internal states not assumed Markov); Krámli–Szász (1983); Markov additive processes.

The new content is arithmetic: the exact hidden-residue block recursion, the certificate computations, adaptive refinement, sticky reflections, and the Z_4 identity-run obstruction. I found no prior source for that part.

## 3. What I verified, and how

**A. Statement 1** (`a_recursion_check.py`, 112 s; `dp_crosscheck.py`).
- I simulated the pair chain with exact fractions over all d^k digit words.
- Σ_w ηη^T, the definition Σ_s d^(k−1−s) Σ_w D(π_s), and the recursion's A_k(h) agree exactly, and Σ_w η = 0. The debt change was re-derived by factorizing M.

| Map | k checked |
|---|---|
| Z_5 (1,2,3,7,1), Z_5 (1,1,2,3,7), Z_4 (1,3,5,7) | 2, 3, 4 |
| Z_7 (1,1,1,2,3,5,11) | 2, 3 |

- Starts tested per (map, k):
  - 30 random exact states;
  - 30 perturbed twins, (M·g with g ≡ 1 mod d^k, e + d^k·u). They give identical coupling sequences, confirming that e′ mod d^(s−1) depends only on (M, e) mod d^s;
  - special states (identity runs; M ≡ 1 mod d^k with M ≠ 1 and e = 0), up to k = 5 on Z_5 and k = 6 on Z_4.
- Result: 0 mismatches. My vectorized DP also agrees with the memoized recursion on 6,600 random (k, h).

**B. Statements 2 and 2′ (proof review).** I checked:
- the Taylor bound, with a remainder uniform under |η| ≤ kB or KB;
- the finite set of blocks and zero blocks;
- the sampled supermartingale, its stopping index, bounded V, and the overshoot terms ((r₁+kB)/r)^α and KB;
- that the τ_n are stopping times, and the strong Markov property;
- the escape word in its current text.

No gap remains. None of statements 1, 2 or 2′ uses primality of d: π is a permutation for any unit M̄, d divides N, and the Terras bijection mod d^n holds for composite d. This answers the d = 4 question; the ratio group mod 4^s is all units, of size 2^(2s−1).

**C. Statement 3** (`c_certificates.py`). I used my own chunked DP, exact deduplication, and an exact positive-definiteness test of S′. The test uses principal-minor sums in Python integers, not Sylvester's criterion.
- All counts match: 250; 156,557; 3,907,806 (7,812,500 states, 50 s, 0.43 GB); 2,471,816 (31 s, 0.45 GB); 1,048.
- Only the state (1, 0) has a zero block; the minimum exact block rank equals the rank; all five forms are exactly balanced.
- Margins of the integer forms are 0.0587, 0.0030, 0.0307, 0.0298 and 0.0517, matching the new column. α_max is 0.125, 0.0061, 0.063, 0.061 and 0.109.
- Minimum one-step coupling ranks are 2, 1, 1, 1, 2, all attained at reflections. Λ, independence, G and (i) all check out.
- Negative controls agree with float checks away from 0: changing one entry of the (1,2,3,7,1) form makes 14 blocks fail (`neg_control.out`).

**C′. Statement 3′ and the Z_4 claim** (`e_adaptive_check.py`, `h_greedy.py`).
- I used a greedy exact refinement: refine a node exactly when its block is nonzero and not exactly balanced. For a fixed Q this succeeds if and only if some refinement-tree partition does.
- Lifts are built as a canonical lift times the kernel of H_{s+1} → H_s, together with e + d^s t. Fibre sizes are equal and coverage is counted exactly.
- With the author's forms, all six {±1} maps (levels 3..5) and Z_4 (levels 1..6) pass.

  | Map | Nodes refined per level | Distinct leaf blocks |
  |---|---|---|
  | Z_4 (1,3,5,7) | 7, 50, 400, 1,200, 1,819 | 26,004 |
  | Z_5 (1,6,4,31,1) | 632, 1,463 | 28,428 |

- The only zero-block leaf is the state (1, 0).
- Level-K blocks match the recursion on 400 random nodes (`h_sanity.out`).
- I reviewed `block_adaptive.py` itself. The lifts are complete, `next_level` implements the recursion correctly, final-level nodes all become leaves, zero blocks are correctly excluded, and every leaf goes through the exact check. The `exact_check_vec` change does not affect passing forms.

**D. Statement 4** (`j_identity_runs.py`).
- A_k(1, b d^(k−1)) = Σ_paths D(j ↦ j + bΠm̄) for every b. I checked this to k = 8 on Z_4, k = 6 on Z_5 (1,2,3,7,1) and (1,4,1,11,34), and k = 4 on Z_7 (1,1,1,1,2,3,5).
- On Z_4, A_k = 4^(k−1) D(j ↦ j+2), which has rank 2, for every k ≤ 8.
- Full Z_4 tables up to k = 6 have minimum block rank 2.

**E. Sticky reflections** (`n_sticky.py`).
- All 32 subsets of Z/5 are reflection-invariant.
- From M ≡ −1 and a sticky e ≡ b, every digit keeps M′ ≡ −1 (200 exact random states per map).
- Reflection covariance ranks are all ≤ 2.

**F. THM-4609** (`g_hierarchy.py`, 11 s).
- I checked every unit set and every rank from 3 to d−1 for d = 7 and 11. Ties were decided with exact Bareiss elimination.
- Q = I works for translations from rank 4, for the largest odd-order group from rank 5, and for all units from rank 6. In each class it fails at the rank just below.
- The structural facts never fail: at most one fixed point; tr D = 2·#(moved non-units); λ_max ≤ 4, with equality iff an even cycle lies inside the non-units.
- The three squares-group forms and all six statement-6 rows check out exactly: residues, coupling group, rank, independence, Πm < d^d, Λ, r_i, and that the form balances every coupling with a in G.
- The generalization holds: if all m_i are congruent mod d, every ratio is ≡ 1, so M ≡ 1 mod d at all times.

**G. Census and statement 6** (`i_census_check.py`, `h_greedy_part3.out`, `h_greedy_part4.out`).
- Class sizes {(3,2): 2,652, (3,4): 44,208, (4,4): 3,408} recomputed independently.
- All 18 census certificates and property (i) verified.
- The three statement-6 certificates verified: Z_5 (1,8,3,7,12) with r = (0,−3,4,4,2), Z_5 (1,2,3,7,6), and Z_7 (1,2,3,5,1,1,1). The Z_7 run took 4 s and 1.29 GB and ran alone.
- Audit F's rank ≥ 3 list in `partC_maps.out` is exactly these three plus (1,3,5,7), (1,6,11,16,1), (1,2,3,7,11) and the HYP map. So "all of audit F's maps" holds.

**H. Other claims in the note.**
- The A4-metric facts hold: under 5I − J, translations are strictly balanced, and every coupling with a ≠ 1 sits exactly at equality.
- The lift counts 2,950 and 3,000 of 3,125 reproduce for (1,2,3,7,11) and (1,6,11,7,2) (`o_a4_lifts.out`).
- The twisted-obstruction lemma is correct, and Z_3 (1,1,5) has no constant-c obstruction (a two-line argument).
- The numbers in §4 and §5 match `twisted_check.out` and `rank1_two_valued.out`.

## 4. Prior art

Searches were brief: web search plus the arXiv abstract pages for 1203.3459, 1402.2558 and 1506.08541. The references are in finding 21. I could not give theorem numbers for the Fayolle–Malyshev–Menshikov or Menshikov–Popov–Wade books, because that would need PDF downloads.

## 5. Not checked

- The author's specific partitions and leaf counts (578/1,823 …, and 2,529,032 for Z_4). They depend on unrecorded pilot forms; I verified existence through the greedy partition instead.
- The float optimal margins at k = 4 for the six {±1} maps. Only their sign is proved, and exactly.
- THM-4610, which HYP-9244 cites for rank one; it is outside this audit.

## Scripts and outputs

| Script | Output | Covers |
|---|---|---|
| `hcore.py`, `hdp.py` | — | core routines |
| `a_recursion_check.py`, `dp_crosscheck.py` | `.out` | statement 1 |
| `c_certificates.py` | `c_certificates_light.out`, `c_certificates_Z5_11237.out`, `c_certificates_Z7_1111235.out` | statement 3 |
| `c_census_v1_fixed.py` | `.out` | finding 2 |
| `neg_control.py` | `.out` | exact-test controls |
| `e_adaptive_check.py`, `h_greedy.py`, `h_sanity.py` | `e_adaptive_check.out`, `h_greedy_part1–4.out`, `h_sanity.out` | 3′, Z_4, statement 6 |
| `f_dual_certificates.py`, `f_dual_sixth.py`, `k_dump_certs.py`, `l_verify_certs.py` | `f_dual_*.out`, `dual_certificates.json`, `l_verify_certs.out` | finding 1 |
| `f_optimum.py` | `f_optimum_reading.out`, `f_optimum_k4.out` | optimum brackets |
| `g_hierarchy.py` | `.out` | THM-4609 |
| `i_census_check.py` | `.out` | census |
| `j_identity_runs.py` | `.out` | statement 4 |
| `m_identity_persist.py` | `.out` | finding 6 |
| `n_sticky.py` | `.out` | sticky reflections |
| `o_a4_lifts.py` | `.out` | A4 metric and lifts |

I ran heavy jobs one at a time after checking `ps`. One run (`j_identity_runs.py`) built a full Z_4 k = 6 table at roughly 1–1.5 GB while two lighter jobs were running, which bent the one-heavy-job rule.

Sources: [arXiv:1203.3459](https://arxiv.org/abs/1203.3459), [arXiv:1402.2558](https://arxiv.org/abs/1402.2558), [arXiv:1506.08541](https://arxiv.org/abs/1506.08541), [Meyn–Tweedie 1994](https://projecteuclid.org/euclid.aoap/1177005204), [Comets–Menshikov–Popov 1998](https://projecteuclid.org/journals/annals-of-probability/volume-26/issue-4/Lyapunov-functions-for-random-walks-and-strings-in-randomenvironment/10.1214/aop/1022855869.full), [Malyshev–Men'shikov 1979](https://mathnet.ru/eng/mmo369), [MPW 2016 (Cambridge)](https://resolve.cambridge.org/core/books/abs/nonhomogeneous-random-walks/references/0D79EF27E4305AEF5C713B846C6575FE), [FMM 1995 (catalogue)](https://catalogue.i2m.univ-amu.fr/bib/8982), [Georgiou–Wade abstract](https://www.maths.dur.ac.uk/users/andrew.wade/abstracts/gw1.html)
