# Independent audit: kind-pasteur-2026-11-09-S18 (THM-5070, THM-5071, HYP-9640)

**Verdict: PASS WITH CORRECTIONS.** 2 MAJOR findings (one credit, one over-claim) and 15 MINOR findings. No wrong
number and no wrong proof was found.

**Scope.**
- Commit b60c7cdea8: THM-5070 (1)-(5), THM-5071 (1)-(5), HYP-9640 (A)-(D), and the packet (README, PROOF,
  PREREGISTRATION, scripts, `out/`).
- Coordinator addendum, commit 140c2996b4: THM-5070 (6) and HYP-9640 (E), with `bases_gaps.py`, `hunt.c`,
  `verify_witnesses.py` and the base 4-8 outputs.

**Method.** Every number below was recomputed with the audit's own code (`audit/a*.py`, `audit/a*.c`). The packet's
scripts were read but not reused, with one exception: `chain_blocks.py` was run once, only to compare its matrix set
with the audit's.

**What holds.**
- All mathematics is correct.
- Every finite claim reproduces exactly, including the packet's BFS state counts.
- Several claims were checked further than the packet checked them (see "Independent checks").

The corrections concern credit, status typing and over-reaching glosses.

---

## Verdict by item

| item | packet status | audit |
|---|---|---|
| THM-5070 (1) transfer | PROVED | **correct**. The proof was re-derived. A black-box linear solve on every block class of order <= 5 in bases 2-8 gives a nonnegative integer matrix (0 failures). 1,300 random chains of orders 8-16 match brute force in h, E and O, bases 2-8 (0 failures). Credit: m3 |
| THM-5070 (2) apex bridges | PROVED | **correct**: 2,800 checks, 0 failures; the a(W) recurrence matches Session191 one_apex/PROOF.md (1). Wording: m11 |
| THM-5070 (3) tower | PROVED | **correct**: char. poly, (t-2)(t^2+t+1), eigenvector (3,1,2), H_3 = (3*2^n + c_n)/7 for n <= 14 |
| THM-5070 (4) coverage | FINITE-EXACT | **correct, every number and state count reproduced**. Both order-5 rows were also confirmed by direct runs to 2*10^5 (the packet used containment). The 15-generator family was extended to 3*10^6. Proof text: m4 |
| THM-5070 (5) tower directions | FINITE-EXACT | **finite statement correct**: also only tower powers at length exactly 4. **The gloss over-claims: M2** |
| THM-5070 (6) bases 4-8 to 30,000 | FINITE-EXACT | **correct**, conditional on THM-4940 (6) for m <= 80. All 61 witnesses were re-verified independently. The attainment side of the m <= 80 table was independently confirmed. Wording and robustness: m5, m6 |
| THM-5071 (1)-(2) weighted formula | PROVED | **correct**: exhaustive on all 1,052,741 loopless digraphs of order <= 5, plus 320,000 random of orders 6-8, by three routes, 0 failures. Typing and citation: m1 |
| THM-5071 (3) parity | PROVED | **correct**: same 1,372,741 digraphs. Credit to Berge: m2 |
| THM-5071 (4a) sink back arcs | PROVED | **correct**: every S at orders 2-10 |
| THM-5071 (4b) one reciprocal pair | PROVED | **mathematically correct**, ingredients rechecked. **Credit: M1** (the 7 and 1 cases and the method are Session234's) |
| THM-5071 (4c) census | FINITE-EXACT | **correct**. Credit: m8 |
| THM-5071 (5) cores | FINITE-EXACT | **correct**: 3 + 16 classes, all 19 decompositions agree. Title vs status: m9 |
| HYP-9640 (A)-(E) | OPEN | statements and implications sound. (A) now holds to 3*10^6. Stale pointers in HYP-9530: m7 |

---

## MAJOR findings

### M1. THM-5071 (4b) presents known results as new, without credit (credit and status typing)

**Claim.** (4b) says that a semicomplete digraph with exactly one reciprocal pair never has H in {1, 4, 7}, at any
order. It is typed PROVED in this session, with a "self-contained four-cell lemma".

**Session234 (read first by S18, per README "Concurrent work") already proves most of this:**
- `session234_digraphs_20261109/PROOFS.md`, **DG3**: H(D) = prod H(C_i) over the ordered strong components.
- **DG4**: a strong semicomplete digraph of order >= 3 has a Hamiltonian cycle. Keeping its arcs gives a **strong
  spanning tournament, for any number of reciprocal pairs**, so H(D) >= H(T) >= f(n).
- Section 3: Moon plus path insertion gives a nondecreasing f, and f(3..7) = 3, 5, 9, 15, 25.
- **DG5 refinement**: "a semicomplete graph with H=7 has **exactly two** reciprocal pairs." The Session234 README
  repeats it: "Seven paths force an absent pair or exactly two reciprocal pairs".
- `GENERAL_GRAMMAR.md` section 1: every non-singleton strong component has H >= 2.

**What this means for (4b).**
- The cases **H != 7** and **H != 1** were already proved in the repository.
- The method is Session234's DG3 + DG4 + Moon floor.
- The four-cell lemma is the b = 1 special case of DG4. The audit checked DG4's general form on every strong
  semicomplete digraph of order 3-5 (51,969 at order 5; `a7_defects.out` (ii)).
- The only new case is **H != 4**.

**Where the credit is missing.** The packet even quotes Session234's "7: a missing pair or exactly two reciprocal
pairs", yet these places present (4b) as this session's theorem:
- THM-5071 status;
- PROOF.md (4b);
- README section 2 ("This is proved for all orders");
- README "What survives" ("never for 7 (4b)");
- HYP-9640 (C) status;
- both INDEX rows.

**The mathematics is fine.** The audit re-verified:
- f(3..7) = 3, 5, 9, 15, 25;
- Moon's deletion property at orders 4-7;
- the one-pair lemma at orders 3-7 (892,208 strong digraphs at order 7);
- no strong one-reciprocal-pair digraph of order <= 4 has H in {1, 4, 7}.

**Fix.**
- Type (4b) as: "PROVED. The cases 1 and 7 are KNOWN (Session234 DG3/DG4 and the DG5 reciprocal-count refinement).
  New: the exclusion of 4, by the same route."
- Cite DG3 in Step 1 and DG4 in Step 3. The four-cell lemma may stay, labelled as the b = 1 case of DG4.
- In README "What survives", attribute "never for 7" to Session234.

### M2. THM-5070 (5) is stated beyond its finite check

**What is checked, and correct.** Among products of length <= 3 of the 57 base-3 blocks, the matrices that preserve
the direction (3, 1, 2) are exactly M_T^1..M_T^12. The audit reproduced this, finding 119,516 distinct products. At
length exactly 4, only M_T^4..M_T^16 occur (`a5b_tower_len4.out`).

**Glosses that assert more:**
- THM-5070 (5): "So these blocks give **no** finite-state digit system (affine numeration) along the tower."
- HYP-9640 (A) Route: "THM-5070 (5) **excludes** an affine digit system along the tower."
- `05-knowledge/results/INDEX.md`: "**no** tower-direction digit system **exists** among order-4 blocks."
- README (5): "This is **also the reason** the 2-dimensional 'sink-normalised' subsystem and the apex-bridge family
  alone leave **infinitely many** gaps (exploratory runs, not claimed)."

**Why these are over-claims.**
- The first three are universal statements about products of every length, and "digit system" is never defined.
- The fourth infers "infinitely many gaps" and a causal link from finite exploratory runs. THM-5070's own Remarks say
  only "gaps that keep growing".
- THM-4940's audit treated "stated beyond its proof" as MAJOR, so the same grade applies here.

**Fix.** Qualify each with "among products of length <= 3" (the audit adds: and of length 4). In the README, replace
"infinitely many" and "the reason" with the canon's wording.

---

## MINOR findings

**m1. THM-5071 (1): citation and typing.**
- (1) is the specialisation at (1, 0, 0, ...) of Grinberg-Stanley, arXiv:2307.05569v1, **Theorem 1.30** (Section
  1.8.2). There phi(sigma) sums (length - 1) over the cycles of sigma that are D-cycles.
- Irving-Omar, arXiv:2412.10572v3, Corollary 14 restates it; their citation bracket reads "[12, Theorem 1.3]".
- The packet cites the paper without a theorem number. The specialisation is one line (F_S(1, 0, ...) = [S empty]),
  so (1) should be typed **KNOWN/CITED**, with the theorem number.
- (2) is correctly PROVED as a regrouping. Credit that the pairing of a D-cycle with its reversed complement cycle is
  GS's own device (their proof of the tournament p-expansion, Section 3).
- The sign convention agrees with GS. Pages read: arxiv.org/abs/2307.05569, the arXiv and ar5iv HTML of 2307.05569,
  and arxiv.org/html/2412.10572v3.

**m2. THM-5071 (3): Berge not credited.**
- Berge's theorem says ham(D) = ham(complement of D) (mod 2) for every digraph. GS recover it, and it appears in
  Grinberg-Stanley's IPAC 2023 slides. arXiv:2402.07606 derives it from the reciprocity
  u_X(-m) = (-1)^|V| u_Xbar(m); the audit saw this only in a search summary and did not check a theorem number.
- (3) refines Berge. The defect cycles of D and of its complement are the same set: the reciprocal-pair cycles of one
  are the missing-pair complement cycles of the other. So (3) implies Berge.
- Replace "generalising Redei" with "refining Berge's parity theorem (Redei is the tournament case)".
- Whether GS Section 6 already states (3) could not be read (the HTML is truncated there). Check before claiming
  novelty.

**m3. THM-5070 (1): missing precedent.**
- THM-4952 (A) already gives a finite-state transfer for I(Omega(T), x) in **every base** on span-w ribbon words.
- A port chain with blocks of order <= k is a ribbon word of span <= 2k - 1 in its natural vertex order. Internal
  back arcs span <= k - 1, and the glue arc w -> u spans <= 2k - 1.
- The 3-dimensional nonnegative port state is a genuinely new compression, but THM-4952 (A) should be cited as the
  general-base precedent.
- The S86 packet asked for exactly this engine (`packing_bases_repunits_20261108/README.md`, "A ternary coverage
  engine ... a base-3 analogue of THM-4875 (D) is needed"). Credit it.

**m4. THM-5070 Proof (4) / PROOF.md (4): the exactness argument is incomplete as written.**
- Pruning at h > LIM needs **h' >= eta h >= h**, which holds because eta = I(Omega(B), x) >= 1. This is not stated.
- "min(., CAP) commutes with every row" should read cap(row . s) = cap(row . cap(s)) for nonnegative integer rows,
  with h never capped.
- The argument is sound once these are added. Empirically, the audit got identical miss sets with the packet cap, a
  larger cap and no cap at all, in bases 2, 3, 4 and 8 (`a4_cover_satcheck.out`).

**m5. `bases_gaps.py`: latent unsoundness.**
- The script parses only the printed misses, but `chain_cover.c` prints at most 200 (`if (nmiss <= 200)`). A run
  with more than 200 misses would silently count the unprinted ones as attained.
- This did not happen here: the largest list is 46, in base 8.
- Fix: parse the `count` line and assert that it equals the length of the parsed list.

**m6. THM-5070 (6) / README wording.**
- (i) The "order-join closure" source settled **nothing** in (80, 30000] in any base. After closure the open sets
  equal the chain misses above 80 (`a9_bases.out`). Say so, or drop "three sources".
- (ii) "61 independently verified witnesses": they were verified by a second script of the same session. This audit
  is the first independent verification, and all 61 pass (`a9_bases.out`).
- (iii) The m <= 80 part rests wholly on THM-4940 (6). Its attainment side for b = 4-8 is now independently confirmed.
  These non-gaps <= 80 are missed by order-5 chains:

  | b | non-gap | where the audit found it |
  |---|---|---|
  | 5 | 62 | order <= 7 census |
  | 6 | 47 | order <= 7 census |
  | 6 | 59 | new witness, order 8 |
  | 8 | 43 | new witness, order 10 |
  | 8 | 47 | new witness, order 8 |
  | 8 | 80 | order <= 7 census |

  The new witnesses are in `a9_small_hunt.out`, rechecked by `a1_brute`. The gap side remains THM-4940 (6)
  (audited). It is consistent with the audit's chains and its census of order <= 7.

**m7. Stale pointers in HYP-9530.**
- HYP-9530 (A) status still says "complete for alpha-sums <= 9678 (values <= 29,035)" and "Every m <= 9678 is already
  settled".
- (B) does not point to THM-5070 (6) or HYP-9640 (E).
- Add a status line pointing to THM-5070 (4) and (6) (values to 10^6; bases 4-8 to alpha-sum 30,000).

**m8. P4 and (4c): partial priority.**
- The H = 7 half of P4 (order 5, one missing pair, no reciprocal pair) was already in Session234's census before S18's
  pre-registration commit 189865f292: `census/general_census.csv` row `5,7,1,0,1,600`. The audit recounted the 600
  labelled digraphs (`a7_z1count.out`).
- Credit Session234 in (4c) and in the P4 outcome. The 21-at-order-6 half and the full 0..129 census are new.

**m9. THM-5071 (5): title and status disagree.**
- The title says "Every strong semicomplete core with H = 7 has no odd tournament cycle", which is an all-order claim.
  The status says FINITE-EXACT at orders <= 5.
- The all-order form follows because H = 7 cores have order 4 (Session234 DG4/DG5). Cite that, and type it accordingly.

**m10. Wording in THM-5070 Remarks.** "The same 271 / 264 order-5 matrices realise both Chan-Pak spectra": it is the
same **blocks** (387 classes). The matrices depend on the base: 264 distinct at x = 2, 271 at x = 3, and 270 or 271
at x = 4..8.

**m11. Wording in THM-5070 (2).**
- "Every odd cycle of T passes through u" is equivalent to T - u being transitive, since any cycle in a tournament
  implies a 3-cycle. So the hypothesis means exactly one-apex sides.
- The general display is therefore no more general than the one-apex formula. State the equivalence.

**m12. Citations in the README.**
- "7 would need alpha = (3, 0) ... impossible (THM-1370)": the alpha = (3, 0) exclusion is THM-4940 (3) / THM-343.
  THM-1370 excludes H = 7 by strong floors.
- "The diagonal never enters" is Session234 DG1. Credit it in the README and in the THM-5071 Remarks.

**m13. README section 2.** "The values excluded by one 11 pair are exactly {1, 4, 7}. This is proved for all orders."
"Exactly" is HYP-9640 (C), which is OPEN above 197. Reword: "1, 4 and 7 are excluded at all orders (proved); that
nothing else is excluded is HYP-9640 (C)".

**m14. README section 3.** "The weighted formula shows it would need weight 3 on every odd cycle" is a heuristic,
since nothing forces a construction to act cycle by cycle. Write "suggests", consistent with the ANALOGY typing.

**m15. PREREGISTRATION P5 wording.** It promised "10^5 random digraphs of orders 6 and 7"; the run used 10^5 at
order 6 and 10^4 at order 7. This is harmless, since P5 is a gate and the audit's check is far larger. Say so in the
outcome line.

---

## Observations (not findings)

**The SEAM remark made explicit.** At x = 2, SEAM.md's coordinates are a fixed linear image of THM-5070's:
(H, S, R) = (h, O - E, E + O).
- On the converse, H counts Hamiltonian paths, S counts those ending at the port, and R sums H(T - V(P)) over paths
  P from the port.
- Checked on 60 random port chains of orders <= 8 (`a3_seam_coords.out`).
- This makes "KNOWN at x = 2 through OCF, in other coordinates" concrete. The identity S = O - E is observed, not
  proved here.

**HYP-9640 (A) evidence extended.** The 15 base-3 generators reach every v = 1 (mod 3) up to **3*10^6** except the
same ten values.
- The run used 293,669,502 states.
- All 999,990 reached values are certified by exact, unsaturated word products (`a4_cover_gen15_3e6.out`).

**Untracked files.** At the start of the audit, `bases_gaps.py` and `out/bases_4to8_30000.out` were untracked in this
worktree. They were committed in 140c2996b4 and are audited above under the addendum.

---

## Independent checks run (all outputs in `audit/out/`)

**Transfer matrices (`a1_blocks.py`).**
- Method: black-box linear solve, with no case analysis. Each block class (B, w, u') of order <= 5 (387 classes) is
  glued after 8 test tournaments. States are computed by brute force: DFS cycle enumeration, memoised packing
  polynomial and DFS path enumeration. N is solved from 3 independent states and validated on 5 more, for
  x = 2..8.
- Result: 0 failures of linearity, integrality or nonnegativity.
- Distinct matrices:

  | block order | x = 2 | x = 3 | x = 4..8 |
  |---|---|---|---|
  | <= 4 | 56 | 57 | 57 |
  | <= 5 | 264 | 271 | 270 / 271 |

**Comparison with the packet (`a1_compare.py`).**
- The audit's sets equal `chain_blocks.py`'s sets for x = 2, 3 and k <= 4, 5.
- Every line of `generators_base3_15.txt` (15) and `generators_base2_13.txt` (13) was rebuilt from its block. Each
  block is a tournament of order <= 4 and gives the printed matrix.

**Random chains (`a1_brute.c` + `a1_chains.py`).**
- Brute force in C: cycle DP, packing DP and a path-from-port DP, giving (h, E, O) for x = 2..8, plus the direct
  Hamiltonian path count.
- 1,300 random chains of orders 8-16, built from every block class: 0 failures in h, E and O, all bases, and OCF.

**Apex bridges and tower (`a2_apex_tower.py`).**
- Apex bridges: 2,800 checks (k, l <= 6, x = 2..8) of the formula, of E, O, q_e, q_o and of a(W). 0 failures.
- Tower:
  - det(tI - M_T) on a 9 x 10 grid;
  - the factorisation at x = 3;
  - M_T (3, 1, 2) = 2 (3, 1, 2);
  - brute-force (h, E, O) of the tower for n = 1..14 against M_T^(n-1) (1, 0, 1) and (3*2^n + c_n)/7.

  0 failures.

**Coverage (`a4_cover.c`).** The audit's own BFS uses Fibonacci hashing, parent pointers, and an exact 128-bit
recomputation of the certificate word for every reached value. All state counts equal the packet's.

| family | base | bound | states | misses |
|---|---|---|---|---|
| order <= 4 (57) | 3 | 2*10^5 | 22,540,354 | the ten |
| 15 generators | 3 | 10^6 | 64,420,757 | the ten |
| 15 generators (extension) | 3 | 3*10^6 | 293,669,502 | the ten |
| order <= 4 (56) | 2 | 2*10^5 | 55,566,866 | 7, 21, 47, 67, 127, 173, 491 |
| 13 generators | 2 | 2*10^5 | 14,389,973 | the 12 listed |
| order <= 5 (271), direct | 3 | 2*10^5 | 92,242,587 | 10, 31 |
| order <= 5 (264), direct | 2 | 2*10^5 | 322,705,694 | 7, 21 |
| order <= 5, bases 4-8 | 4-8 | 1 + 30000b | — | identical to `out/bases_4to8_30000.out` |

- Saturation exactness: identical miss sets with the packet cap, cap = LIM and no cap, in bases 2, 3, 4, 8.
- Certificate failures: 0 in every run.

**Tower directions (`a5_tower_direction.py`, `a5b_tower_len4.py`).** Products of length <= 3 that preserve (3, 1, 2)
are exactly M_T^1..12. Products of length 4 that preserve it are exactly M_T^4..16.

**Weighted cycle formula (`a6_gs.c`).**
- Coverage: every loopless digraph of order <= 5 (1,052,741), plus 200,000 / 100,000 / 20,000 random digraphs of
  orders 6 / 7 / 8. The random samples mix uniform, semicomplete, oriented and near-tournament digraphs.
- Three routes, all compared with H:
  - (1) through vertex-set cycle weights and a set-partition DP;
  - (1) by explicit enumeration of the n! permutations, at orders <= 6;
  - (2) by explicit DFS of the cycles of D and of its complement, classified by pair type.
- Results: 0 formula failures and 0 parity failures. (4a) holds for every S at orders 2-10.

**Defects (`a7_defects.c`).**
- f(3..7) = 3, 5, 9, 15, 25. Moon's deletion property holds at orders 4-7.
- Strong-resolution lemma: holds for b = 1 at orders 3-7, and for any b at orders 3-5.
- One reciprocal pair, orders <= 7: the values below 198 are exactly {2..197} minus {4, 7}; 21 first appears at
  order 6.
- One missing pair, orders <= 7: every value 0..129; 7 first at order 5, 21 first at order 6.
- `a7_witnesses.py` checks the explicit witness matrices. `a7_z1count.py` recounts Session234's 600.

**Cores (`a8_cores.py`).**
- 3 strong H = 7 classes, all of order 4 with b = 2, none with an odd tournament-only cycle.
- 16 strong H = 21 classes at order 5: 6 with b = 2 and 10 with b = 3, as in Session234.
- All 19 packing decompositions equal the packet's lines.

**Bases 4-8 (`a9_small.c`, `a9_bases.py`).**
- A complete census of I(Omega, b), b = 2..8, over all tournaments of order <= 7.
- Hill-climbed witnesses for b = 6, m = 59 and for b = 8, m = 43 and 47.
- Closure from the audit's own chain misses: the open lists equal the packet's.
- All 61 packet witnesses verified by `a1_brute`: tournament property, I(Omega, b) = target, and OCF.

**Literature.** Fetched pages:
- the arXiv abstract and HTML of 2307.05569 (v1 numbering: Theorem 1.30 general, Section 3 tournament case);
- arXiv HTML 2412.10572v3 (Corollaries 14, 19, 20);
- a web search for Berge's theorem.

## Reproduction

Binaries go in a temporary directory, never in this folder. Let `S` be a scratch directory.

- **Block matrices.** `python3 a1_blocks.py S 5` writes `S/mats_x{2..8}.txt` and `S/blocks_all.txt`.
- **Compile.** `gcc -O2 -o S/a1_brute a1_brute.c`; the same for `a4_cover.c`, `a6_gs.c`, `a7_defects.c` and
  `a9_small.c`.
- **Plain matrix files for a4_cover.**
  - One family: `awk -F'|' '$1+0 <= K {print $2}' S/mats_xX.txt > S/m.txt`.
  - Generator files: `grep -v '^#' ../out/generators_base3_15.txt | cut -d'|' -f1 > S/g.txt`.
- **Coverage.** `S/a4_cover S/m.txt X LIM CAPMODE`, where CAPMODE is 0 (packet cap), c (cap c) or -1 (exact).
- **The rest.**
  - `python3 a1_chains.py S/blocks_all.txt S/a1_brute TRIALS NMAX KMAX SEED`
  - `python3 a9_bases.py out/a4_cover_bases4to8.out ../out/bases_4to8_witnesses.txt S/a1_brute`
  - The other scripts take the arguments shown in their docstrings.

The largest run (3*10^6) needs about 8 GB and 1 minute; the others need seconds to 2 minutes.
