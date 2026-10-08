# Run-compressed pair chains, the orphan law, and the reducible locus

**Session.** mac-mini-2026-10-07-twoanchor (continuation).

**Owner prompt.** "keep pursuing possible next steps and explore around and synthesize new ones that you pursue as well". The previous part is [two anchors](twoanchor_reset2_friezes_20261007.md) (THM-4600–4602).

**Concurrent work.** The codex sessions build directly on THM-4601:
* [general head phases](collatz_general_head_phases_20261007b.md);
* [two-anchor head decoder](collatz_twoanchor_head_decoder_20261007b.md);
* [signed ladders](collatz_complement_ladders_20261007c.md);
* [two-twos grounding](collatz_two_twos_grounding_20261007c.md).

They ground specific Mersenne children; this note does not duplicate that. Its ladder lemma completes their decoder (§1).

COLLATZ IS STILL OPEN.

**Audits.** Two independent audits ran on 2026-10-07 (audit_C, audit_D in `04-computation/experiments/runcompress_20261007/`). The core mathematics and all data were confirmed. Their corrections are applied below; MISTAKE-587 records them.

## 0. Results and types

| Item | Statement | Type |
|---|---|---|
| [THM-4603](../../01-canon/theorems/THM-4603-run-compressed-pair-chains-drift-lemma-ladder-completeness-and-long-run-transitions.md) (1) | **Drift lemma.** If `u − c = 3^k(v − c′)` with c, c′ T-periodic points, the orbits run their patterns in lockstep and the debt drifts by the difference of odd densities per Terras step. THM-4600 is the zero-drift case. | PROVED; 720 + 2,748 exact tests |
| THM-4603 (2) | **Ladder completeness.** In an admissible state, every equal-time merge in which both orbits have an odd value before the merge is a ladder collision, `z′ = 4^i z + (4^i−1)/3`. So the ladder compiler of THM-4601 (iv) (all i, child and source ladders) is complete, and the concurrent signed-gap decoder becomes complete for all common futures. | PROVED (elementary, on the classical ladder); 1651 + 2,948 merges |
| THM-4603 (3) | **Long-run transitions.** Conditional on eventual periodicity of the rational limit orbit (open; contains Collatz on negative integers), a long run of any periodic pattern leads to absorption, re-anchoring, zero drift, or linear drift. Universal state: a child ones-run re-anchors at −1 with debt −5; a source ones-run drifts at +1/2. | conditional classification; FINITE-EXACT table (550 periodic limit pairs) |
| [THM-4604](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md) | **Collatz words lie on the reducible locus.** Carries form a twisted 1-cocycle, trivialised on ⟨w⟩ by c_w, and w ↦ G_w is faithful. Traces depend only on (\|w\|, Σw) and equal 2cosh(δ/2). Word pairs lie on the Cayley cubic tr[A,B] = 2; Markov triples live at −2; no G_w is ±I. Frieze minors of the carry configuration are carries (carry exchange relation). | PROVED (elementary; KNOWN in substance); the frieze reading is DICTIONARY |
| [HYP-9242](../hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md) | **Orphan law.** Deletion orphans tend to 0 like a negative power of log n. The upper bound (log n)^(−1/2) up to logs holds at sketch level for random t. Effective exponents are 0.42 (K = 2, 3), 0.54–0.63 (K = 9–33) and 0.59 on the Mersenne line. Orphans' post-run orbits are +1 SD long against random orbits, and 40% are rescued 3-adically. | NUMERICAL; upper bound at sketch level |
| [HYP-9241](../hypotheses/HYP-9241-pair-independence-product-law-for-collatz-translation-partners.md) (update) | Direct product test: q_2/(q_01 q_12 q_02) = 0.95 → 1.00; q_3/(all six pairs) = 0.97–1.13 for T ≤ 256. | NUMERICAL |

## 1. Ladder completeness (THM-4603 (2))

**Why it holds.**
* Just before two orbits first coincide, their last odd values z_u ≠ z_v both lead to the first common odd value Z.
* For integers they are two of the classical odd U-preimages `(2^a Z − 1)/3`, a running over a progression of step 2. So they form the ladder `z ↦ 4z + 1`.
* For 2-adic orbits the same conclusion needs admissibility: the equal-time identity forces `2^(t_u) ≡ 2^(t_v)` mod 3, so the letter gap is even.
* Merges are child ladders `z_v = 4^i z_u + (4^i−1)/3` (letters c, c+2i) or source ladders. The heads satisfy `|v| = |u| + k_0` and `Σu − Σv = ±2i`, plus the anchored identity. The merge time is `max(Σu, Σv) + 1`.

**Scope.** Both hypotheses are needed:
* the reset state (1,1) merges at s = 1 with no odd value before the merge;
* inadmissible starts give half-ladders with odd gaps (audit C).

**Consequences.**
* MISTAKE-586's "the i = 1 grammar is incomplete" is sharpened: the full ladder grammar is complete.
* The concurrent signed-gap decoder proves completeness of `F_v(Y) = S^r(F_u(X))` only at a fixed marked state, "not for arbitrary common-future diagrams". THM-4603 (2) shows every common future has that form.

**Frequencies.** Among the 1651 merges from the universal state:

| ladder | i = 1 | i = 2 | i = 3 | i = 4 |
|---|---|---|---|---|
| child | 870 | 50 | 18 | |
| source | 672 | 33 | 7 | 1 |

Source ladders are almost as common as child ladders. Audit C replicates the shares: i = 1 0.936, source 0.431.

## 2. Drift and long-run transitions (THM-4603 (1), (3))

**The drift lemma.**
* Debts change at the rate of the difference of the odd densities of the periodic points the orbits are near.
* This explains the collapse of the reset children in THM-4601 (v): the source sits on +1 (density 1/2) and the child's limit on the halving point 0 (density 0), giving +1/2. It also explains re-anchoring (zero drift).

**From the universal state (3, 1−27).**
* The realisable transitions are listed in THM-4603 (3); the exact table is `transitions.out`.
* Runs beginning with the letter 2 cannot follow a two-run. In particular the session's first absorption list included four non-realisable entries (audit C).
* Key cases, checked on actual sources:
  * child ones-run: re-anchor at −1 with debt −5 after 12 steps (limit pair (−53, −1));
  * source ones-run: drift +1/2 after an 8-step transient;
  * child letter-3 run: debt 4; letter-4 run: debt 9; letter-5 run: drift +7/15;
  * child (1,1,4): absorbed at 23; source (1,2,5): absorbed at 20.
* The label "SHIFT" in `transitions.out` covers both a same cycle out of phase and different cycles of equal odd density. 224 of its 280 cases are the latter.
* The same correction was made independently by codex-complement ([run-classification scope](collatz_run_classification_scope_20261007e.md)). Its exact control is the child run on (2,4) from the universal state: distinct cycles of periods 12 and 6, both of density 1/3. Without a periodicity witness the outcome is UNRESOLVED.

**What this gives.**
* Positive-drift runs explain why no bounded post-run depth suffices: a source ones-run of length L, a new Mersenne-like excursion, raises the debt by about L/2.
* They do not cause the s^(−1/2) tail, which is the diffusive return of the zero-drift debt walk (THM-4581). Survival to s = 3000 is 0.316 overall and 0.312 among paths without a source ones-run ≥ 14 (audit C).

## 3. The orphan law (HYP-9242)

**Positive integers.**
* A source n is certified by h_D exactly when they merge at equal time before 1. For orbits reaching 1 this is equivalent to `o(n) = o(h_D)` together with `σ_T(n) = σ_T(h_D) + D`. This is elementary: aligned orbits that first reach 1 together agree from 8 on.
* An **orphan** has no such D.

**Mersenne line** (orbits of 2^K − 1 computed exactly for K ≤ 12800 and reproduced by audit C):

| K window | odd K | orphans | fraction | fraction × √K_mid |
|---|---|---|---|---|
| [200, 400) | 100 | 6 | 0.0600 | 1.01 |
| [400, 800) | 200 | 13 | 0.0650 | 1.55 |
| [800, 1600) | 400 | 13 | 0.0325 | 1.09 |
| [1600, 3200) | 800 | 17 | 0.0213 | 1.01 |
| [3200, 6400) | 1600 | 27 | 0.0169 | 1.14 |
| [6400, 12800] | 3200 | 30 | 0.0094 | 0.89 |

* The fitted exponent is 0.59 (95% CI [0.42, 0.75]). The data do not distinguish 1/2 from 0.6–0.7 or from log-corrected laws.

**General residual sources** `n = 2^K t − 1`:
* Fitted exponents (audit C, 2000 per row, B ≤ 3200): about 0.42 for K = 2, 3; 0.54 for K = 9; 0.605 for K = 17; 0.62 for K = 33.
* The effective exponent grows with the number of children, and small sizes bias it downward. K = 2 measures 0.42, although its single-chain exponent is 1/2.
* There is an upper bound at sketch level for random t: `P(orphan) ≤ C L^(−1/2)(log L)^2`, from the child D = 1 alone. Absorption before log2 h_1 steps is a certificate, and the driving bits are exactly uniform.

**Which K are orphans.**
* Partners share σ_T − K exactly. So comparing orbit lengths within the Mersenne line is biased: tercile rates of 0.31% / 0.92% / 3.0% are reproduced by a null model (audit C).
* The genuine signal: against random orbits of the same size, orphans' post-run lengths are +1.0 SD long (median +0.94, 86% positive, n = 83), and non-orphans' −0.01. That is a positive aggregate debt drift.
* It is not a single long ones-run (Mann–Whitney z = −0.38), and not a leading-digit effect (1.9–2.3% in each quarter of frac(K log2 3) for K in [600, 6400]).

**Full-orbit success per shell** (m = v2(K−1), odd K in [600, 6400]): the reset pair succeeds 82–91% and h_3/h_4 82–88%, both flat in m.

**Rescue by 3-adic predecessors.**
* Within backward depth 36 a smaller node certifies n for 20 of the 53 orphans in [1000, 6400], against 0.403 of all odd K (p = 0.41).
* Depth 3 occurs exactly when K ≡ 5 mod 6, via (8n−5)/9.
* Doubly uncovered: 33 of the 2700 odd K (1.2%).

## 3b. The depth spectrum of the rule on the Mersenne line

`mersenne_depth_spectrum.py` (in `depthspec_20261007/`) covers odd K in [1001, 3999], children D = 1 and D = 3, with depth counted in Terras steps after the run end.

| quantity | value |
|---|---|
| merges with D = 1 / D = 3 / either | 86.9% / 81.3% / 90.5% |
| median best deletion depth | 204 steps (≈ 0.09 K) |
| median descent depth | 2.81 K |
| deletion certificate shorter than descent | 87.2% |

Survival P(best depth > s), with "no merge" counted as infinite:

| s | 10 | 30 | 100 | 300 | 1000 | 3000 | 10000 |
|---|---|---|---|---|---|---|---|
| P | 0.95 | 0.83 | 0.65 | 0.49 | 0.31 | 0.17 | 0.11 |
| P·√s | 3.0 | 4.5 | 6.5 | 8.5 | 9.7 | 9.4 | 11.1 |

* The variable depth of the rule is heavy-tailed: P·√s is roughly flat (9.4–11) for s between 10³ and 10⁴, then cut off at the orbit length.
* Deletion certificates are typically far shorter than descent.

## 4. The reducible locus (THM-4604): what traces, minors and friezes see

**Collatz words lie on the reducible locus.**
* The steps are upper-triangular matrices, `G_w = [[3^|w|, B_w],[0, 2^(Σw)]]`, composed chronologically.
* The carry is a twisted 1-cocycle, positive, odd and prime to 3.
* The cycle point c_w trivialises it on ⟨w⟩.
* The centraliser of G_w is the torus fixing c_w and ∞; the anchored states are its elements with multiplier 3^k.
* w ↦ G_w is faithful: the carry decodes the word.

**Traces see only the semisimplification.**
* `t(w) = 2cosh(δ_w/2)` with `δ_w/ln 2 = |w| log2 3 − Σw`. When Σw = round(|w| log2 3) this is HYP-9230's θ_f with f = |w|.
* Every pair of words lies on the Cayley cubic `x² + y² + z² − xyz = 4` (tr[A,B] = 2). Example: the odd/even pair gives `(5/√6, 3/√2, 7/√12)`.
* Markov triples live at tr[A,B] = −2, with integer points 3·Markov.
* Conway–Coxeter friezes are quiddities with `M(a_1)⋯M(a_n) = −I`, where `tr[M(a), M(b)] = 2 + (a−b)²`.
* No Collatz transfer matrix is ±I. Read as quiddities, however, some valuation words do close friezes: (1,1,1) and (1,2,2,2,2,2,2,1,7) are the U-words of 7183 and 2583211, which merge at 24245.

**Minors see the cocycle.**
* The frieze minors of the carry configuration are carries: `det(v_i, v_j) = 3^i 2^(A_i) B_(w[i+1..j])`.
* The Plücker relation becomes the **carry exchange relation** `B_xy B_yz = B_y B_xyz + 3^|y| 2^(Σy) B_x B_z`, verified on 3000 random triples.
* This is the cocycle-side type-A structure (a positive frieze with 6-smooth coefficients). It holds for every triple, so on its own it does not detect merges.

**Reading (DICTIONARY).**
* Merges are incidences `F_u(n) = F_v(h)` at integers, decided by carries and ternary labels. Traces cannot detect them. Example: 483 and 469 have U-words (1,7) and (7,1), equal traces, and merge at 17.
* The frieze program yields coordinates because its minors are cocycle values.
* It yields no rewrites because ear moves are relations of a ↦ M(a), while a ↦ G_a has none. The ear move changes (|w|, Σw) by (+1, +3) and flips the endpoint label mod 3.
* On the Christoffel tree of O/E words the Markov-type mutations walk the Farey tree of slopes with values 2cosh(δ/2). This is the linear GL(2,Z)-action on (δ_u, δ_v) (Goldman 2003).

## 5. HYP-9240 and digits

* The +1 barrier's landing N_D and entry debt are functions of D and of the first s_0 binary digits of 3^(−D), where s_0 ≈ 1.8–2.2 D.
* Real-size bounds give only `N_D + 1 < (3/2)^D`; observed landings are far smaller (≤ 880 for D ≤ 3000).
* An all-D proof needs 2-adic control of a nonlinear function of those digits. The zeroless-power problems are an analogy, not a reduction. D ≤ 8000 is settled exactly.

## 6. Next directions

1. A debt-payment mechanism for orphans: negative-drift runs, 3-adic predecessors, or something new for the doubly uncovered 1.2%.
2. Make the orphan upper bound rigorous, by bringing THM-4581 (4) beyond sketch level, and decide whether the many-children exponent exceeds 1/2.
3. Complete ladder grammars at −1 (the tiling compiler) and at the short-J states. Then give the exact Haar coverage per depth for the whole residual branch.
4. Use the carry exchange relation as an organising principle for collisions. Do ladder collisions of two words satisfy a mixed exchange relation?

## 7. Reproduction

**`04-computation/experiments/runcompress_20261007/`:**
* `ladder_and_drift.py`, `drift_check_actual.py`, `fricke_check.py`: ALL CHECKS PASSED.
* `transitions.py` (+ `.out`).
* `mersenne_sigma*.py` (data `mersenne_sigma_12800.txt`).
* `orphan_law_mersenne.py`, `orphan_law_general.py`, `mersenne_partners_by_shell.py`, `orphan_structure.py`, `orphan_time.py`, `orphan_rescue.py`, `pair_product_R3.py` (+ `.out`).
* `audit_C/`, `audit_D/` (scripts, outputs, REPORT.md). These include audit C's null model (`c4_length_selection.py`) and its exponent fits (`c4_general_exponent.py`).

**`04-computation/experiments/depthspec_20261007/`:** `mersenne_depth_spectrum.py`, `analyze_depth_spectrum.py` (+ outputs).
