# Independent audit: THM-4863, HYP-9455 and the S72 path-weighted deletion tower note

**Audited:**
- `05-knowledge/results/primegame_deletion_tower_path_weighted_decks_20261104.md` (the note);
- `01-canon/theorems/THM-4863-...md`;
- `05-knowledge/hypotheses/HYP-9455-...md`, with its pre-registration commits `ce33106a23` and `e58a09c18b`;
- the five session programs and their outputs in `04-computation/experiments/`.

**Auditor:** adversarial audit for opus-2026-11-04-S72, 2026-11-04, worktree `C:/tmp/wt-s72` at `cd2bc6a5e1`.

## Verdict: PASS WITH CORRECTIONS

The census reproduces exactly with independently written code, in full through order 10 and with exact comparisons (no hashing). The proofs of Lemma W, Theorem H, Lemma L, Lemma J and the Phi properties are correct; Lemma L has one small gap (m1).

The required corrections concern prior art and typing:
- (B), the line-class count, is already canon (THM-851).
- Most of the census structure the note calls "New here" is in the literature.
- Theorem H's literature status must be filled in. Its key step is Pouzet's cycle-cover reconstruction lemma, and its determinant analogue is known.
- Several repository precedents are uncited.

No numerical claim was found false.

## MAJOR issues

### M1. (B), the line-class count (2^m + 2^h)/4, is existing canon but is typed as a new PROVED result

**Locations:**
- THM-4863: the title clause "(B) Line classes number (2^m + 2^h)/4", the status line "PROVED: ... and (B)", and section (B);
- the note: the Status bullet "the Burnside count of line classes" (PROVED), the §4.5 heading "(PROVED Burnside ...)", and the §6 row "line classes (2^m + 2^h)/4 | PROVED".

**Prior art:**
- **THM-851** (codex-2026-07-15-S15, `01-canon/theorems/THM-851-node-coloured-defect-coherent-closure.md`):
  - eq. (13) tabulates "line orbits 1, 3, 20, 272, 8320" for n = 3..7;
  - eq. (14) proves |L_n/<sigma>| = (2^(m-1) + 2^(h-1))/2 = 2^(h-2)(2^r + 1) by Burnside, with the same fixed-point argument: a black line cannot be fixed by reflection.
- **kind-pasteur-2026-07-01-S15:** `07-reflections/the-quarter-tiling-model-and-the-reconstruction-wall-at-n7-kps.md`, line 29, has the Klein-four quotient <sigma, phi> with total (2^m + 2^D)/4 and values 3, 20, 272, 8320, 525312 for n = 4..8.
- **Session130:** `colored_incidence/REPORT.md` and `primegame_lines/THEORY.md` §1 split n = 5 as 8 blue lines plus 12 black pairs, which is 20.
- **THM-796 (4.5)-(4.7):** the formula is a one-line corollary, b_n + k_n/2.

**Required:**
- Retype (B) as CITED (THM-851 (14); THM-796 §4).
- Remove (B) from THM-4863's PROVED list and title, or mark it CITED there.
- Keep the following, which I reproduced:
  - the DICTIONARY reading of the owner's "1" and "3";
  - the FINITE-EXACT remark that for n <= 5 line classes correspond one-to-one with node adjacencies (1, 3, 20 classes on 1, 3, 20 adjacencies). It fails at n = 6: 272 classes on 183 adjacencies.

### M2. The "New here" census structure is largely in the literature

**Locations:**
- note §3.2: "New here. In this census, every even-order pair consists of mutual converses. Every odd-order pair consists of two self-converse tournaments whose Hamiltonian cycle counts differ by exactly one.";
- the §6 row for McKay's census;
- THM-4863: title clause (P) and section (P);
- the results-INDEX entry;
- HYP-9455 (C): "visible at order 4 and confirmed at orders 8 and 10".

**What is already known:**
- **(a) Odd-order pairs are self-converse.** McKay 2022 (Australas. J. Combin. 83(3), 448-457), p. 453, item (f), states: "Stockmeyer's tournaments of order 2^t + 1 are self-complementary."
  - For tournaments, self-complementary means self-converse.
  - The order-5 and order-9 pairs are the only pairs at those orders, so they are Stockmeyer's. I checked this by construction (below).
  - The order-3 case is trivial.
  - The note quotes McKay (f) one paragraph earlier, so this half of "New here" is not new.
- **(b) Even-order pairs are converse pairs.** Boudabbous–Ille, "The {-2,-1}-selfdual and decomposable tournaments", Discuss. Math. Graph Theory 38 (2018) 743-789, p. 749.
  - They attribute to Stockmeyer 1977 a tournament tau of order 2^n + 2 (n >= 3) with three properties: every card is self-dual, tau is prime, and tau is not self-dual.
  - Its same-deck partner is therefore its dual. This covers the order-10 pair.
  - The order-4 pair is the Harary–Palmer pair 1 -> C3 / C3 -> 1, visibly a converse pair (Stockmeyer 2021, Figure 1, as the note says).
- **(c) Hamiltonian cycle counts differ by exactly one.** Bankoussou-mabiala, Boussaïri, Chaïchaâ, Chergui, Lakhlifi, "On the spectral reconstruction problem for digraphs", arXiv:1910.13914 (2019).
  - Lemma 6 (§3) is stated there as a particular case of Pouzet, "Note sur le problème de Ulam", J. Combin. Theory Ser. B 27 (1979) 231-236, Equality (A). It says that for same-deck digraphs, det G - det H = (-1)^(n+1)(C(G) - C(H)), where C counts Hamiltonian cycles.
  - p. 3: Dumont checked, as reported by Pouzet 1979, that |det B_n - det C_n| = 1 for n <= 6. Here (B_n, C_n) is Stockmeyer's pair of order 2^n + 1.
  - Proposition 1: the two determinants have opposite parity for all n >= 3, that is, at every order 2^n + 1 >= 9.
  - Hence "hc differs by exactly one" at orders 3, 5, 9 (and 17, 33, 65) is in the literature in determinant form.
  - I verified the determinant law on every collision in the census (determinants 0/1, 2/3, 11/12, 54/55; numbers in the reproduction section). The law also holds for the converse-merged collision at order 7, by the same argument, since det A^T = det A.
  - Pouzet 1979 was not accessible (paywall). It is cited through arXiv:1910.13914.

**What remains new, as far as found:**
- that the six pairs at orders 6 and 8 are converse pairs. I found this stated nowhere. One order-6 pair (H = 43) is Stockmeyer's 1981 D_2,1, and one order-8 pair (H = 633) is his G_8; for neither family did I find the converse relation stated, and the 1981 paper was not accessed;
- the Hamiltonian-path form of the statement, which is Theorem H (see M3).

**Supporting computation (outside the census):**
- **Stockmeyer's 1981 tournament family D_m,n:** I built it from the description in Egrot–Hirsch, arXiv:2008.01327, Definition 6.1 and Figure 9, family D.
  - It gives same-deck converse pairs at orders 6, 10, 12 and 18.
  - It gives same-deck self-converse pairs at orders 3, 5, 9 and 17.
  - The order-17 pair has hc 405,661,087 / 405,661,086 and H differing by 2.
- **The power-of-two family G_p (Stockmeyer 2021):** it gives converse pairs at p = 4, 8 and 16.
- So HYP-9455 (C)'s required consequence for Stockmeyer's even families holds at every family order computed (6, 8, 10, 12, 16, 18). Part of this is already stated in the literature.

**Required:**
- Delete "New here", or restrict it to orders 6 and 8.
- Type the census structure as FINITE-EXACT reproduction plus CITED: McKay (f); Boudabbous–Ille 2018; Pouzet 1979 Equality (A) and Dumont, via arXiv:1910.13914, Lemma 6, p. 3 and Proposition 1.
- Add these sources to THM-4863's CITED list.

### M3. Theorem H: literature status and missed precedents

**Locations:**
- note §3.3: "We did not find Theorem H ... Its novelty is left to the audit";
- the §6 row for Theorem H;
- THM-4863 (H).

**Proof:** correct; see "Checks of the proofs" below.

**Literature:** I found no statement of Theorem H anywhere, in either form:
- H reconstructible at even order and H - 2hc at odd order;
- the converse-merged version.

But it is an immediate corollary of published ingredients:
- **(i) The odd-cycle formula:** Grinberg–Stanley, arXiv:2307.05569, Theorem 1.39 and Lemma 6.5 (THM-002).
- **(ii) The classical cycle-cover reconstruction lemma:**
  - **Pouzet 1977** ("Quelques remarques sur les résultats de Tutte concernant le problème de Ulam", Publ. Dép. Math. Lyon 14(2), 1-8, numdam PDML_1977__14_2_1_0), I-2.2(a): for an arbitrary square matrix, the sum over permutations with q >= 2 cycles of prescribed lengths is reconstructible. I-2.2 says explicitly that the argument does not depend on the graphs being oriented.
  - This is exactly steps 3-4 of the proof. Tutte 1979 and Kocay 1981 are the other standard sources.
- **(iii) The determinant analogue:** the same argument with the sign character in place of the odd-cycle weight 2^(#cycles) is Pouzet 1979 Equality (A), as quoted in arXiv:1910.13914 Lemma 6.

So Theorem H is "new as stated" but a short corollary. The note credits only Kelly's lemma and Kocay's lemma.

**Repository precedents the note missed:**
- **The codex reflection** `07-reflections/tournament-subset-decks-higher-order-faces-and-kneser-middle-duality-codex-20260821.md`, §3, eqs. (3.6)-(3.12). It shows that at n = 7 neither H nor the Hamiltonian 7-cycle count is an affine function of the 6-deck, and the defects are paired:
  - one example has H totals 144 vs 146 and hc totals 11 vs 12;
  - a 5-deck hostile has (H, hc) = (109, 8) vs (111, 9).

  These are instances of the Theorem H mechanism. At n = 7 no spanning odd-cycle collection has two or more cycles, so H - 2hc is an affine function of the (merged) deck.
  - I recomputed from the reflection's bitstrings (`codex0821_check.py`): the H - 2hc sums are 31 + 91 = 122 and 29 + 93 = 122, and the 5-deck pair has H - 2hc = 93 = 93.
  - The reflection's sentence "There is no affine functional of that deck recovering even the Hamilton-path count" is true for H, but must be reconciled with the fact that H - 2hc is affine there.
- `05-knowledge/results/iso_class_info_geometry_s20z.out`, line 475, remarks that Stockmeyer's same-deck pairs may differ in H.

**Required:**
- Type Theorem H as: "PROVED. Statement not found in the literature; a short corollary of the odd-cycle formula and the Pouzet/Tutte/Kocay cycle-cover reconstruction lemma. Determinant analogue: Pouzet 1979 Equality (A)."
- Cite Pouzet 1977 I-2.2(a) for steps 3-4.
- Cite and reconcile the codex-20260821 reflection.
- State that the even-order half has no non-trivial instance through order 10, so its "checks" there are vacuous:
  - merged decks are injective at orders 4, 6, 8 and 10;
  - every even-order ordinary same-deck pair is a converse pair.

## MINOR issues

**m1. Lemma L: the "distinct nodes" argument is incomplete** (note §3.4; THM-4863 (L)).
- The two strong-component sequences ({1}, C_1, ..., C_k) and (C_1, ..., C_k, {1}) coincide when every C_i is a singleton, that is, when S is transitive.
- Add one line: S is not transitive, because a transitive tournament is self-converse and S is not isomorphic to S^op.

**m2. "This explains the parity law through order 10" overclaims** (note §3.4). HYP-9455 makes the same claim: "the note, THM-4863 (P), (L) and (H), explains it as follows".
- Lemma L accounts for the odd-order collisions without Hamiltonian cycles.
- Nothing explains the even-order injectivity; it is a census fact through order 10.
- The sentence "the lift of a self-converse pair {A, B} has cards [A] and [B], which differ" does not show that the two lifted decks differ, because [A] could occur among the other cards of 1 -> B. In the census it holds only because merged decks are injective at orders 4, 6 and 10.
- Reword as: accounts for the odd-order collisions; even-order injectivity is observed, not explained.

**m3. Pre-registration record** (note §7; HYP-9455).
- **The first pre-registration is clean:**
  - commit `ce33106a23` (09:51) contains P1-P5 and the census program;
  - the outputs first appear in `e58a09c18b` (09:58);
  - the census source is unchanged since;
  - my rerun of the committed program reproduces the committed output line for line, apart from CRLF line endings.
- **The second pre-registration is weaker:**
  - `e58a09c18b` names `deletion_tower_order10_20261104.c` and its hash protocol;
  - but that program, together with its output, was first committed in `cd2bc6a5e1` (10:21);
  - so the predictions were committed before the results, but the program was not frozen with them.
- Record this in §7.

**m4. "McKay's complete census through order 13"** (note §3.2): wording only.
- McKay tabulates counts only through order 8 (p. 454).
- From order 9 he states that only Stockmeyer's families are known, and that for tournaments this is the full set through order 13.
- The counts 1, 1 at orders 9 and 10 are therefore inferred from the family orders, with one tournament pair per order 2^s + 2^t.
- McKay's §4 is titled "Directed graphs", and his own searches (e)-(f) use reduced decks.

**m5. Uncited repository precedents.** None anticipates Lemma W, the tower injectivity or the parity law, but each should be cited where it is used:
- **THM-801 §4**, text after its table: the face-row entries "count fixed-Hamiltonian-path presentations whose high endpoint, gap, or low endpoint chart lands in a specified lower converse class". This is the reading that Lemma W makes explicit as a vertex-deck formula.
- **Session130**, `colored_incidence/REPORT.md` §2, "Hamiltonian-path meaning": |F_n(u)| = s(T) H/|Aut T| with Aut acting freely on paths, and Theorem 4 (the endpoint-deletion profile r_u = b_u + deg K_u).
- **Session129**, `tournaments/REPORT.md` §3, the "self-converse sandwich": T = 1 -> S -> 1 has an axis profile because every Hamiltonian path runs from source to sink. This is a special case of Lemma J.
- **THM-3324 (5), THM-3369, THM-3372 and THM-3377:** ordered-join compilers in which deletion respects the join factor. Lemma J is the same mechanism for W.
- **The October-19 Hasse/deletion theorem**, recovered in `session105_metagraph_recursion_20261029/README.md` §2: the unweighted, unmerged deletion map D_n has finite cokernel, so it has full rank over Q at every n. This is the unweighted analogue of HYP-9455 (B) and is useful context.
- **THM-796 §7:** Phi is a signed marginal of the coupling tensor Xi_n.
- **The 2026-10-31 PRIMEGAME lanes.** §0's "Two lanes answered it then" undercounts them:
  - S71 lists four lanes;
  - session130 `primegame_lines/THEORY.md` reads the X0 two-cycles 11M <-> 13M as blue lines and has the n = 5 line counts;
  - session132 has a FRACTRAN program for the literal line counts.

**m6. Uncited external context.**
- **(a) Harary–Palmer 1967** (Monatsh. Math. 71, 14-23): non-strong tournaments on at least 5 vertices are reconstructible (as cited in Achour–Boudabbous–Boussaïri, arXiv:1204.2513).
  - A tournament without a Hamiltonian cycle is exactly a non-strong one (Camion).
  - So every Hamiltonian-cycle-free merged collision is necessarily a merged-only collision, which fits Lemma L.
- **(b) The converse-merged deck problem is {-1}-half-reconstruction** ("demi-reconstruction", reconstruction up to duality) in the Lopez–Pouzet school:
  - Hagendorf–Lopez, C. R. Acad. Sci. Paris Sér. I 317 (1993) 7-12;
  - Boudabbous–Dammak, C. R. Acad. Sci. Paris Sér. I 326 (1998) 1037-1040, who introduced {-k}-half-reconstruction (zbMATH and Dammak–Salem arXiv:1311.1765, p. 5; the paper itself was not accessed);
  - Dammak, Proyecciones 22 (2003).

  Stockmeyer's self-dual odd pairs already give {-1}-half-reconstruction failures at every order 2^t + 1. I found no statement of the even/odd parity pattern in this literature, so the FINITE-EXACT parity law through order 10 appears new.

**m7. §4.1 wording.**
- "They are separated only by the involution kappa": TT3 and C3 are separated by H, scores and Aut. What is meant is that deletion data do not separate them.
- "is accurate in this form" overstates a re-indexing: the owner's single node at n = 3 becomes the single node at n = 2. Say "can be read as". The DICTIONARY typing is right.

**m8. §4.2 and the §6 row for the order-4 line fraction.** The note calls the reading "exact, labelling-independent".
- The identity itself is PROVED: Phi(TT4) = e_C3 - e_TT3, Phi(pair) = 0, Phi(S4) = e_TT3 - e_C3.
- Naming TT3 and C3 by any two primes p and q gives {p/q, q/p, 1}.
- So the match with {11/13, 13/11} is session129's reciprocal-pair dictionary. It is independent only of which of 11 and 13 is assigned to which node.
- The §6 row, the results-INDEX entry ("= PRIMEGAME's first block with 1") and the THM-4863 title should say: identity PROVED, PRIMEGAME reading DICTIONARY.

**m9. Title and INDEX wording.**
- "path-weighted decks determine every merged tournament node through order 10" is false at order 3 (TT3 and C3). Say "at orders 4-10".
- "fails exactly at the odd orders" should read "at the odd orders 3-9 and at no even order 4-10".

**m10. §4.4, "Ten is the number of PRIMEGAME's registers 2, ..., 29".** S71, following Conway, calls only the data primes 2, 3, 5, 7 registers; 11, ..., 29 are control flags. Say "ten primes".

**m11. HYP-9455 (C) and (D).**
- In (C), "Equivalently: every same-deck pair of even order is a converse pair, and ...": the first clause is implied by (C), not equivalent to it; the second clause alone restates (C).
- (D) is already established at many odd orders:
  - at orders 2^t + 1, by Stockmeyer's self-converse pairs (McKay (f));
  - through Lemma L, at every order m + 1 where order m has a converse same-deck pair. My Stockmeyer runs give such pairs at orders 6, 10, 12 and 18, hence collisions at 7, 11, 13 and 19.
  - It is open, for example, at order 15.

**m12. Small inaccuracies.**
- §1: "as THM-796 (6.2) states" — the equality of the two faces is stated in the sentence after (6.3), not in (6.2).
- §4.2, "for which it is automatic": a node with a single tiling has Phi = 0 when [kappa d t] = [d t]. I checked that this never happens at orders 6-8, but it is not automatic.

**m13. Update the bookkeeping per M1-M3.**
- THM-4863's status and CITED list need to add:
  - canon: THM-851 and THM-801;
  - repository packets: session130;
  - literature: Pouzet 1977 and 1979, arXiv:1910.13914, Boudabbous–Ille 2018, Harary–Palmer 1967.
- Move (B) to CITED.
- The results-INDEX line "Theorem H: ... (THM-002 + Kelly)" should add the Kocay/Pouzet partition step.

**Concurrent lane, after the S72 checkpoint, not prior art:** opus-2026-11-04-S73's THM-4866 (one-arc joins: H = x_h y_h + x_R y_R - x_e y_e and hc = s_a e_b) uses endpoint counts much as Lemma J does. Cross-reference the two when both land.

**Incidental, outside this packet:** `05-knowledge/results/reduction_dag_recursion_modes_boxeph_S196.out`, near line 31, says the deck "determines T's char_A" through char_A' = sum of the char_(T-v).
- The derivative identity determines char_A only up to its constant term, which is (-1)^n det A.
- Same-deck tournaments can differ in that term: TT3 has characteristic polynomial x^3 and C3 has x^3 - 1, and the odd same-deck pairs have determinants 2/3 (order 5) and 54/55 (order 9).
- This deserves a separate correction.

## Checks of the proofs

**Lemma W.** Correct.
- Aut acts freely on directed Hamiltonian paths: an automorphism mapping a path to itself fixes its first vertex, hence the whole path.
- So the fibre of a class T is HP(T)/Aut, and the number of tilings deleting into v is (1/|Aut|) times the sum of h_e(x) over x with [T - x] = v.
- For a converse pair, the T^op sheet contributes the start-weighted deck, because h_e^(T^op) = h_s^T and [T^op - x] = [T - x].
- For a self-converse T, an anti-automorphism phi gives h_e(x) = h_s(phi x) and [T - phi x] = [T - x], so the start- and end-weighted decks agree.
- W itself is converse-invariant.
- Checked against all tilings at orders 3-8, through both faces, together with the fibre sizes.

**Theorem H.** Correct.
- **Step 1:** H = sum over S of g(S), where g(S) sums 2^(#cycles) over spanning odd-cycle collections of T[S], and g is converse-invariant.
- **Non-spanning terms:**
  - Kelly's identity s([S], T) = (1/(n-k)) times the sum over x of s([S], T - x) holds for the merged class [S].
  - s([S], X) depends only on the merged class of X.
- **Spanning terms with c >= 2 cycles:**
  - marking one cycle counts each collection c times;
  - N([A],[B]) = s([A],T) s([B],T) minus the number of overlapping pairs;
  - overlapping pairs live in proper subsets, are grouped by the class of their union, and p_A,B is converse-invariant.
- **The only remaining term** is 2hc(T), and it exists only for odd n.
- **The (1/c) bookkeeping is right.** When [A] = [B], N counts ordered pairs, consistently with s([A],T)^2.
- **Attempts to break it:**
  - On every merged-deck collision and every ordinary same-deck pair through order 9, I computed the full cycle-cover decomposition by subset dynamic programming. The odd-cycle formula holds, the spanning covers with c >= 2 agree, and the non-spanning parts agree. Example at order 9: c = 3 covers 10/10, non-spanning part 2417/2417, hc 210/211.
  - Outside the census, the order-17 Stockmeyer pair also satisfies H - H' = 2 = 2(hc - hc').
  - No counterexample is possible, given the proof.

**Lemma L.** Correct, apart from the gap in m1.

**Lemma J.** Correct.
- No arc goes from B to A, so every Hamiltonian path traverses A completely and then B.
- Hence h_s^T = H(B) h_s^A on A, h_e^T = H(A) h_e^B on B, and the cards are (A - x) => B and A => (B - y).

**Phi.** Correct.
- The sum over all nodes is 0, because kappa is a bijection of X_n and d kappa = kappa d.
- The two faces agree, because sigma permutes each fibre, d_L sigma = sigma d_H, sigma commutes with kappa, and [sigma z] = [z].
- At order 4, kappa swaps the two order-3 nodes, which gives the identity.

**Burnside count.** Correct, but it is THM-851's (M1).

## Independent reproduction

All code is in this folder and was written from scratch. No session code is imported.

**My canonical form differs from the session's:**
- individualization-refinement, starting from (score, number of 3-cycles through v, sum of out-neighbour scores) and refining on out-neighbour colour counts;
- full search tree;
- maximum code, with bit = 1 if the later vertex beats the earlier one;
- |Aut| = the number of maximal leaves.

**Validation of the canonical form:**
- invariance under 3 random relabellings, on all classes at orders <= 8 and on 20,000 classes at each of orders 9 and 10;
- brute-force |Aut| over all n! permutations at orders <= 7;
- class and self-converse counts equal OEIS A000568 and A002785 through order 10;
- the Python port reproduces the C codes.

**Commands:**

```
gcc -O2 -fopenmp -o audit_tower audit_tower.c && mkdir -p data
./audit_tower 9 data tiles 8 > run9.out     # 31 s on 12 threads (orders 1-9; Lemma W, Phi and line classes from tilings to order 8)
./audit_tower 10 data > run10.out           # 572 s (orders 1-10, exact sorting, no hashing)
python3 audit_analysis.py data > analysis.out     # 3 min; needs data/ from the order-9 run (data/ deleted after the audit, ~32 MB)
python3 support_check.py; python3 codex0821_check.py
cd stockmeyer; gcc -O2 -o hpaths hpaths.c; python3 stockmeyer_1981.py; python3 stock_gp_powers_of_two.py 4 8 16
gcc -O2 -o dtc ../deletion_tower_census_20261104.c; ./dtc 9    # session program rerun: identical to the committed .out (CRLF aside)
```

(`stockmeyer_1981.py` calls `./hpaths.exe`, the name the binary gets on Windows.)

Data dumps and binaries were deleted after the audit; all are regenerable with the commands above.

**Results, all equal to the note's claims unless marked:**

| n | classes / self-converse / merged | primitive W distinct | rank R_n (mod p, R^T D R, random D) | merged-deck collision classes | ordinary same-deck pairs |
|---:|---|---|---|---:|---|
| 3 | 2 / 2 / 2 | 1 of 2 | 1 of 1 | 1 | 1: H 1/3, hc 0/1, both self-converse |
| 4 | 4 / 2 / 3 | 3 | 2 | 0 | 1 converse pair, H 3 |
| 5 | 12 / 8 / 10 | 10 | 3 | 2 | 1: H 13/15, hc 2/3, both self-converse |
| 6 | 56 / 12 / 34 | 34 | 10 | 0 | 4 converse pairs: (H, hc) = (43,5), (31,3), (29,3), (23,1) |
| 7 | 456 / 88 / 272 | 272 | 34 | 5 | 0 |
| 8 | 6,880 / 176 / 3,528 | 3,528 | 272 | 0 | 2 converse pairs: (633,48), (645,49) |
| 9 | 191,536 / 2,752 / 97,144 | 97,144 | 3,528 | 3 | 1: H 2917/2919, hc 210/211, both self-converse |
| 10 | 9,733,056 / 8,784 / 4,870,920 | 4,870,920 | not tested | 0 | 1 converse pair: H 15,483, hc 939 |

- **Lemma W against every fixed-path tiling, both faces and fibre sizes:** 0 mismatches at n = 3..8, including all 2,097,152 tilings at n = 8. W/(2|Aut|) or W/|Aut| is integral at every node through n = 9. Node tilings sum to 2^binom(n-1,2) through n = 9.
- **Merged-deck collisions** (H, hc and type):

  | n | collision | H | hc | type |
  |---:|---|---|---|---|
  | 7 | strong | 141 / 143 | 14 / 15 | self-converse vs converse pair |
  | 7 | lifts (4) | 43, 31, 29, 23 | 0 | lifts |
  | 9 | Stockmeyer pair | 2917 / 2919 | — | self-converse |
  | 9 | lifts (2) | 633, 645 | 0 | lifts |
  | 5 | same-deck pair | 13 / 15 | — | self-converse |
  | 5 | lift | 3 / 3 | — | lift |

- **Lemma L and the classification:**
  - Every even-order converse same-deck pair (at orders 4, 6, 8) lifts to two distinct merged nodes with equal merged decks.
  - The lifts are exactly the 1, 4 and 2 hc-free collisions at orders 5, 7 and 9.
  - The remaining classes are {TT3, C3}, the order-5 pair, the order-7 strong pair (not an ordinary same-deck pair; both members strong) and the order-9 pair.
- **Theorem H decomposition:** checks listed above. Determinants are 0/1 (order 3), 2/3 (order 5), 11/12 (order 7 strong pair), 54/55 (order 9), and the determinant law holds on every pair.
- **Lemma J:** all 1,295 joins of ordinary classes with |A| + |B| <= 8, 0 failures (the note checked 55).
- **Fibonacci rows:**
  - consecutive-Fibonacci rows: 1/2, 3/7, 1/29, 2/261, 3/3515 at orders 4-8;
  - their nodes' strong-component structures are TT2 => S4, S4 => C3, TT3 => S4, S4 => 1 => C3, 1 => C3 => S4 and TT4 => S4. Each is S4 => TT_k, S4 => C3, S4 => C3 => 1 or S4 => 1 => C3 up to converse;
  - source extensions turn (5, 8) into (2, 3, 8, 13) and (2, 3) into (1, 1, 3, 5);
  - order 7: the pure-blue node with 3 tilings is 1 -> T5 -> 1 (H 15, |Aut| 5), and its Phi is composite.
- **Phi:**
  - order 4: Phi(TT4) = e_C3 - e_TT3, Phi(pair) = 0, Phi(S4) = e_TT3 - e_C3;
  - order 5: 9 values on 10 nodes; only (1,1,C3) and (1,S4) coincide;
  - injective at orders 6, 7, 8;
  - single ratios occur exactly at TT4 and S4 (order 4), at the three pure-blue nodes (order 5), and at the pure-blue single-tiling nodes (2, 3, 3 of them at orders 6-8);
  - the two faces agree and the column sums are 0 at every order 4-8.
- **Line classes:**
  - 1, 3, 20, 272, 8320, 525312 by orbit enumeration at n = 3..8;
  - they match node adjacencies one-to-one for n <= 5;
  - 183 adjacencies at n = 6.
- **Support sets:** merged-deck supports number 1, 2, 6, 34, 264, 3528 and 97141 at n = 3..9, as HYP-9455 states.
- **PRIMEGAME calibration:**
  - codes TT6 = 2, TT3 => C3 = 10, C3 => C3 = 5, 1 -> S4 -> 1 = 7^5;
  - 34 distinct codes;
  - canonical assignment: 1 incidence;
  - my own 20,000 random assignments (seed 4242): min 0, median 1, mean 0.714, max 3, P = 0.579;
  - with primitive rows: P = 0.81;
  - NUMEROLOGY confirmed.
- **Stockmeyer families outside the census:**
  - the 1981 family D_m,n gives converse pairs at orders 6, 10, 12, 18 and self-converse pairs at 3, 5, 9, 17. The order-17 pair has H 12,363,896,559 / 12,363,896,557 and hc 405,661,087 / 405,661,086;
  - the 2021 family G_p gives converse pairs at p = 4, 8, 16. The G_16 pair has H 1,498,610,505 and hc 52,638,402. This script was written by my literature subagent from Stockmeyer's definition and rerun by me.

## Literature consulted

**Accessed** (local copies in the audit scratchpad):
- McKay 2022 (AJC and arXiv versions);
- Stockmeyer, arXiv:2110.06158;
- Bankoussou-mabiala et al., arXiv:1910.13914;
- Pouzet 1977 (numdam);
- Boudabbous–Ille 2018;
- Egrot–Hirsch, arXiv:2008.01327;
- Achour–Boudabbous–Boussaïri, arXiv:1204.2513;
- Dammak–Salem, arXiv:1311.1765.

**Not accessed** (paywall or captcha); judged only from abstracts, zbMATH or citing papers:
- Stockmeyer 1977, 1981 and 1988;
- Kocay 1985;
- Pouzet 1979;
- Boudabbous–Dammak 1998.

**Searches that found no statement of Theorem H or of the parity pattern:**
- deck reconstruction of Hamiltonian path counts for tournaments;
- the Rédei–Berge literature (Grinberg–Stanley, Irving–Omar, Mitrović–Stojadinović);
- half-reconstruction summaries.
