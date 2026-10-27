# Independent adversarial audit of opus-2026-10-27-S63 (THM-4810 / HYP-9420 / note `h_spectrum_atlas_complement_budgets_paley_clocks_20261027.md`)

Audited tree: worktree `/tmp/math-wt-s63` at HEAD `e9b023a7ff` (no file of the worktree was modified). All auditor code, outputs and this report live in `C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a\scratchpad\s63_audit\` (files `my_*.cpp`, `my_*.py`, `my_*.out`, `verify_s73_witnesses.py`, `compare_sweep.py`, `chanpak/chanpak.txt`). Every number below was recomputed with the auditor's own programs (none of the author's binaries or scripts were executed; the author's `.out` files were only read for comparison).

## Verdict: PASS WITH CORRECTIONS

- **Every computed number in the package reproduces exactly** (all 436 table rows of the atlas at orders 3-8, the augmentation results at orders 9-10, all 4,407 certificate lines, the complement-move census at orders 3-7 including every level count and histogram, f(P_q) for q = 7, 11, 19, 23, 31, b(P_7), b(P_11), G1 at order 21, the P_11 depth-one and depth-two censuses, the aliquot table). The proofs of Lemma 1, Theorem P, Theorem C2, Proposition B and Theorem L are correct (with the wording repairs listed under MINOR).
- **But the package misses prior art that is already in this repository's canon and in the literature**, and in consequence its "what is new" framing, HYP-9420 (1) and (3), the "unresolved top of S_9" remark, and the "verified to 6,000" headline are wrong or regressive (MAJOR M1, M2, M5). Two further statements are false as written (M3: the value 661 is not a single-class value; M4: "the regular tournament Z_5{1,4}" is neither regular nor a valid circulant and the named obstruction is labelled-only). The mathematics survives all of these; the text, the typing of novelty, and the index/ledger rows must be corrected.

---

## MAJOR findings

### M1. Sections 1-2, THM-4810 (A)(B), HYP-9420 (1)(3): the strong fibres, the exact S_9, the low ends of S_9 / S_10 and the strong witnesses to 6,000 are already in canon, uncited; the "unresolved" top of S_9 is resolved in canon; the headline range regresses canon

Quotes (note, section "What is new here"): *"(1) The complete atlas of the Hamiltonian-path count on tournaments of order 3-8 ... the complete strong fibres F_s(v) for every odd v <= 74 and the first-instance orders; the exact low end of the strong spectra at orders 9 and 10 by augmentation (section 2), which decides the strong-untouchable odd values below 200; strong witnesses with certificates to 2,881 at order 9 and beyond at order 10."* Section 2: *"The hill-climbing phase found none of the missing values, so the top region of S_9 is not resolved here."* and *"This is THM-1370's completeness conjecture, in its stronger strong-tournament form, verified to 6,000."* HYP-9420 (3): *"Complete for v <= 74; test: the fibres of 75 <= v <= 200 from the augmentation lists."*

Evidence (all files in the audited worktree):

1. **THM-4744** (`01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, opus-2026-10-21-S50, audited): "(3) FINITE-EXACT, complete in the count: every strong tournament with fewer than 625 Hamiltonian paths has order at most 12, and the strong raw-count polynomials P_d are determined for every odd d < 625; for odd d in [3, 623], P_d = 0 exactly for d = 7, 21." Its data (`04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_375.out` line 145 and `strong_ears_625.out` lines 116-117) give `P_d(1)` = the number of strong classes of all orders with H = d, i.e. exactly the S63 fibre F_s(d), for every odd d < 625: `odd d < 625 with no strong tournament of any order: [1, 7, 21]`. I checked all 36 entries of the S63 table against this list: **identical** (`S63 F_s(v), v<=73, equals THM-4744 P_d(1): True`). The same file shows `P_d(1) = 1` exactly for d = 3, 5, 11, 13 among all odd d <= 373, and every P_d(1) for 375 <= d <= 623 is at least 56; so **HYP-9420 (3) is already verified for every odd v <= 623 by existing canon**, not merely "complete for v <= 74", and its proposed "test" (75 <= v <= 200) was done in S50. The session24 synthesis (`05-knowledge/results/session24_research_synthesis_20261021.md` line 77) also records "a globally complete strong-atom inventory for every odd H <= 73" (the same range as "v <= 74").
2. **THM-4097** (`01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`, codex 2026-08-25, independently audited): "An exact subset-DP census of all 1,526,032 such ears from the 6,008 strong order-eight classes independently reproduces the full 1,482-value order-nine strong H-spectrum. It contains every odd value from 85 through 2881." The data file `05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out` lists exactly 1,482 odd values (min 75, max 3357). I compared it with the S63 order-9 certificate: `S63 certificate values: 1482; cert - canon: []; canon - cert: []`. **So the S63 search found all of S_9, and every value it could not find (2883-2977 minus twelve, 2997, 3077-3357 minus twenty-one) is provably absent from S_9 by canon.** The sentence "the top region of S_9 is not resolved here" and pre-registration P2's "the search cannot show absence" must be replaced by the citation. The S50 census `h_multiplicities.out` (`n = 9: 191536 classes, 1520 distinct H values, max H = 3357 ... strong classes 178133`; `odd holes below the max: 159 (first [7, 21, 2883, 2887, 2889, 2893, 2895, 2897])`) also contradicts the note's "where exhaustion over all labelled tournaments is out of reach" as a justification for sampling: the exhaustive class census at order 9 exists in canon.
3. **HYP-9029, THM-4102, THM-4104** (`05-knowledge/hypotheses/HYP-9029-strong-interval-tiling-law.md`, `01-canon/theorems/THM-4102-...`, `THM-4104-...`): "strong(7) >= [65,105], strong(8) >= [69,609], strong(9) >= [85,2881] (islands: 75, 81), strong(10) >= [249,14649] (selected bank), strong(11) >= [429,80265] (selected bank) ... Thus every odd in [65,80265] has a STRONG witness; multiplication bridges the remaining allowed values through 80405." Together with THM-4744 for [3, 623], canon already knows that **the strong-untouchable odd values below 80,266 are exactly 1, 7, 21**. The S63 headline "below 6,001" (note section 2, THM-4810 (B), HYP-9420 (1) status, both INDEX rows, the ledger bullet) is a regression of canon by a factor 13, and "in its stronger strong-tournament form" is not new: HYP-9029 is the strong-tournament form. THM-4810 (A)'s corollary "every odd 3 <= v <= 609 except 7 and 21 is attained by a strong tournament on at most 8 vertices" is implicit in HYP-9029's strong(7)/strong(8) intervals plus THM-1370 (4) and THM-4744.
4. The low end of S_10: THM-4744's order-10 inventory (`strong_ears_225.out` line 6: `order 10: 239 classes; least H = 125; first values [125, 135, 145, 147, 153, 155, 159, 161, 165, 167]`; 8,769 strong 10-classes with H <= 625 in `strong_ears_625.out`) already contains S_10 ∩ [1, 625] exactly; S63's B = 200 computation is contained in it (and agrees on the printed values).
5. The "augmentation principle" of section 2 is exactly the "strong ear" method of THM-4097 and THM-4744 (Moon 1966 + the insertion lemma), used there to orders 9-13.

Suggested correction: cite THM-4744, THM-4097, THM-4102, THM-4104, HYP-9029, the S50 note `h_spectrum_multiplicities_20261021.md` and the session24 inventory in section 0 and in THM-4810's status; retype "(1)" of "What is new here" as *an independent labelled-sweep / Burnside reproduction of the class and strong-class counts per H (new: the labelled counts, the first-instance table) that agrees with THM-4744's inventories and THM-4097's exact S_9*; delete "the top region of S_9 is not resolved here" and the "cannot show absence" outcome of P2 in favour of "S_9 is exact in THM-4097: the 1,482 witnessed values are all of S_9"; restate HYP-9420 (1) with the canon evidence (strong form to 80,265; all-tournament form to 80,555 in the literature, see M2) and HYP-9420 (3) as verified to 623 by THM-4744; fix the two INDEX rows and the ledger bullet accordingly.

### M2. The completeness conjecture is Chan-Pak's published Conjecture 6.2 (origin: MathOverflow 232751, bof and Gordon Royle, 2016; Royle: all odd values below 80,557 by order <= 11); the note attributes it to THM-1370 and never cites Chan-Pak or Royle

Quote (note, section 2): *"This is THM-1370's completeness conjecture, in its stronger strong-tournament form, verified to 6,000."*; section 0: *"the completeness conjecture spectrum(infinity) = odd \ {7, 21}, OPEN"*; HYP-9420 (1): *"THM-1370's completeness conjecture ... is the all-tournament version"*.

Evidence: the author's own cited THM-115 (`01-canon/theorems/THM-115-h21-permanent-gap.md`, banner of 2026-10-20) says "Only the assertion that 7 and 21 are the *only* missing positive odd raw counts remains open; it is Chan--Pak's published Conjecture 6.2." I extracted the paper (`chanpak/chanpak.txt`, from https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf): S. H. Chan, I. Pak, *Computational complexity of counting coincidences*, Theoretical Computer Science 1015 (2024) 114776, arXiv:2308.10214, section 6.5: "Conjecture 6.2. Function h takes all odd values except 7 and 21." and "This conjecture was made by the MathOverflow user bof and Gordon Royle (March 2016). Royle reported that the conjecture holds for all odd integers up to 80557 (ibid)." The repository notes `05-knowledge/results/spectrum_basins_20261020e/spectrum.md` and `rotation_owners_spectrum_20261020e/spectrum.md` already record this (with the inclusive endpoint 80,555 and Royle's use of tournaments on at most eleven vertices) and THM-1370's own "Priority correction, 2026-10-20" retracts novelty for the omissions.

Suggested correction: cite Chan-Pak (Conjecture 6.2) and Royle's MO answer wherever "completeness conjecture" appears (note sections 0, 2, 3; THM-4810 (B) and OPEN line; HYP-9420 (1)); state the known ranges: all-tournament form verified to 80,555 (Royle, orders <= 11; CITED), strong form to 80,265 in canon (HYP-9029 / THM-4104).

### M3. False statement: "661 = P_8-like maximum among them" (single-class values at order 8)

Quote (note, section 1): *"only three of the 320 values of order 8 have a single class (661 = P_8-like maximum among them)."*

Evidence: the author's own table `h_spectrum_sweep_20261027_n8.out` has the row `661	241920	6	6	241920` (six classes). My sweep (`my_sweep_n8.out`): `values with exactly one class: 1 629 637 (3 values)`. There is also no Paley tournament on 8 vertices ("P_8-like" is meaningless). Correction: "the three values of order 8 with a single class are 1, 629 and 637; the maximum 661 is attained by six classes".

### M4. Theorem P's example "the regular tournament Z_5{1, 4}" is not regular, not a circulant, and the obstruction is a labelled-tournament statement only (note section 4; THM-4810 (C))

Quotes: note: *"at n = 5 the regular tournament Z_5{1, 4} is unreachable from the transitive one"*; THM-4810 (C): *"At n = 5 the regular tournament Z_5{1, 4} is in the other class from the transitive tournament and is unreachable from it even with the converse."*; program output `complement_moves_bfs_20261027_n5.out`: *"regular tournament on 5 vertices (Z_5 with {1,4}): f = 2, distance from transitive = 255"*.

Evidence: {1, 4} = {1, -1} is not a tournament connection set (S and -S must be disjoint). The author's `paleyCode()` sets bit k iff (j - i) mod 5 is a square for the pair i < j; the resulting labelled tournament (code 665) has scores (2, 1, 2, 3, 2) and f = 2 (my `my_bfs_n5.out`: `author's 'Z_5{1,4}' code ...: scores 2 1 2 3 2, forward arcs 5, f = 2, L = 2 -- regular? NO`); a regular 5-tournament has f >= L = 3 by the note's own Proposition B, so "f = 2" already contradicts "regular". The intended object, the circulant Z_5{1, 2} (code 947), is regular with 7 forward arcs (odd) and is indeed outside Tr's component (`regular circulant Z_5{1,2}: ... forward arcs 7 (parity odd; Tr has 10), f = 3, L = 3, dist from Tr = -1`). However, of the 24 labelled regular 5-tournaments, `12 with even forward parity, 12 with odd`: a transposition of two labels flips the forward-arc parity (it toggles the pair {i, j} and twice each pair {i, k}, {j, k} with k between i and j), so every isomorphism class has copies in both halves. Correction: replace by "the labelled circulant Z_5{1, 2} (7 forward arcs) is unreachable from Tr; the parity obstruction separates labelled tournaments only - every isomorphism type has labelled copies in both halves". Also fix the program's label.

### M5. Prior art for f(P_7) = 7 and f(P_11) = 20 (section 5 table, THM-4810 (E))

Quote: *"(3) The exception budgets of the Paley clocks P_q beyond P_7: f(P_q) = 7, 20, 64, 92, 180 at q = 7, 11, 19, 23, 31"* presented under "What is new here", typed FINITE-EXACT.

Evidence: OEIS A003141 ("Minimal number of arcs whose reversal yields a transitive tournament", maximised over n-vertex tournaments; Bermond 1972; Sanchez-Flores, Neumann-Lara and Bermond 1994) reads 0, 0, 0, 1, 1, 3, 4, 7, 8, 12, 15, 20, 22, 28 (n = 0..13): the maximum minimum-feedback-arc count is 7 at n = 7 and 20 at n = 11, which P_7 and P_11 attain. My census confirms max f = 7 over all 7-tournaments (`my_bfs_n7.out`: `max f = 7`). Session73's own ports report (cited by the note) already states "Subset DP at 3,7,11 gives 1,7,20". The values f(P_19) = 64, f(P_23) = 92, f(P_31) = 180 I could not find in the literature (search only). Correction: cite A003141 / Bermond for 7 and 20 (and note that P_7 and P_11 attain the known extremal values; by A003141 the maximum at n = 19, 23, 31 is not known to me, so the new values stand as FINITE-EXACT for the Paley tournaments, without a novelty claim for q = 7, 11).

---

## MINOR findings

1. **Miscounts in the order-9 witness report** (note section 2 and section 7). Quote: *"2883-2977 (39 values), 2997"* and *"40 values in [2883, 2997] were not found"*. My recount of the committed certificate (`my_cert_n9.out`): 36 values of [2883, 2977] are missing and twelve are witnessed (2885, 2891, 2903, 2917, 2919, 2931, 2959, 2961, 2963, 2967, 2971, 2973); with 2997 that is 37 missing in [2883, 2997]; 120 missing in [3077, 3357] (21 witnessed, as listed); the total 218 (<= 3400) is right. By M1 all of these are provably absent from S_9 (THM-4097).
2. **Theorem C2 mixes d(T) with the distance to the transitive class.** Quote: *"d(T) <= 2 iff some total order r has backward-pair set D_r(T) = E(P_1) xor E(P_2)"* while d(T) is defined as the distance from the fixed labelled Tr (order 0 < ... < n-1). T xor D_r(T) is Tr_r, the transitive tournament of order r, so the criterion characterises d(T, transitive class) <= 2 (as the P_11 sentence correctly says). Write "d(T, transitive class)" in C2 (note and THM-4810 (C)).
3. **Garbled sentences with a literal "...":** Lemma 1: *"k moves at most k(n - 1) pairs when k is even and at most ... the complement of a path when k is odd"* (for odd k the toggled set is the complement of E(P_1) xor ... xor E(P_k), hence at least n(n-1)/2 - k(n-1) pairs); section 3 table: *"(every proper divisor is at most n / 2 ... the classical height bound)"*; header of `complement_moves_bfs_20261027.cpp`: *"one move toggles n-1 ... exactly the M-(n-1) pairs off the path"*.
4. **Proposition B's "clock-linear bound" is weaker than the cited session73 G3 tariff.** The note proves f(tau_P T) <= f(T) + n - 1 and d >= ceil(f/(n-1)); session73's ports report (section 5, cited by the note) proves |f(tau_P T) - f(T)| <= n - 3 for n >= 4, hence d(T, transitive class) >= ceil(f/(n-3)) (3 for P_11 instead of 2) and its D2 bound ceil((f - (n-1))/(n-3)) for forest frames. My census checks the sharper bound too: `session73 tariff: violations of d >= ceil(f/(n-3)): 0; tight: 2` at n = 4, 5, 6, 7. Cite G3 and use n - 3 in the table column and in Proposition B (Theorem L supersedes both anyway).
5. **"as soon as q >= 17" is not the threshold.** (q^2 - 1)/8 > 2q - 4 iff q^2 - 16q + 31 > 0 iff q >= 14; for odd q it holds from q = 15 (28 > 26) and fails at q = 13 (21 < 22). The statements "for q >= 17" (note, THM-4810 (C)) are true but not sharp; replace by q >= 15.
6. **Output quirk:** `strong_low_spectrum_augment_20261027_B200.out` prints `order 10: 0 labelled strong candidates with H <= 200` because `list10` is filled only with the `11` flag; the true count is 3,752 (`my_augment_B200.out`). Harmless but misleading.
7. **Feedback density of P_23:** 92/253 = 4/11, the same density as P_11; the note writes 92/253 unreduced while stressing that C_3 and P_7 share 1/3.
8. **Depth-two census cost.** The score multiset of tau_Q(T') depends only on T' and the two endpoints of Q (internal vertices get n-1-d, the start n-d, the end n-2-d), so the full 2L histogram over the 8,349,303,325 path pairs is a per-image count of Hamiltonian paths by (start, end) and takes 18 s (`my_p11.out`) instead of "about one hour on 12 threads". Not an error; worth recording since HYP-9420 (4) proposes deeper searches. Also at depth one every one of the 95,095 images has L = 14 exactly (`2L histogram: 28:95095`), so the minimum f = 14 (attained by 55 images) is tight against L.
9. **Pre-registration (section 7).** Git confirms the order: `strong_spectrum_witness_prereg_20261027.md` and the witness program were committed in `604ec71e1b` (2026-10-26 21:21:14 -0600), the witness outputs and certificates in `e9b023a7ff` (21:50:44). P1-P3 outcomes are reported honestly. But P1 (nine values 611-655 at order 9) and P2 were decidable from THM-4097's exact S_9, and P3 from THM-4102/4104; the pre-registration tested nothing that canon did not already contain (see M1).
10. **Provenance nit:** `strong_spectrum_witness_20261027_n9.out` says "certificate ssw9_cert.txt written" while the committed file is `..._n9_certificate.txt` (renamed; contents verified line by line).
11. **THM-1370's "near-top holes {611,...,655}"** is a loose range in THM-1370 itself (the actual holes are the nine values 611, 615, 617, 619, 623, 625, 635, 647, 655); the note's "THM-1370 recorded ... the near-top holes at order 8" is fine but could say "the hole range".
12. **A_n = multiplicative closure of the S_m:** canon has this as THM-1862 / THM-1936 (cited by HYP-9029); the note derives it in a sentence without citation. Cite.
13. **L is an integer:** sum_i (s_(i) - i) = 0, so sum_i |s_(i) - i| is even; stating this makes "|L(tau_P T) - L(T)| <= 1" literally correct (the proof gives a change of at most 2 in 2L).
14. The typing "strong witnesses at orders 9 and 10: NUMERICAL (each witness is exact)" is conservative but acceptable; after M1 the better type is "FINITE-EXACT membership certificates, reproducing THM-4097's exact S_9 and a subset of THM-4102's selected order-10 image".
15. The note's section 3 is carefully typed DICTIONARY / ANALOGY and does not overclaim; its cited facts (s(8) = 7, s^-1(21) = {18, 51, 91}, 28 hermit, 4 = s(9), 13 = s(27), 52 untouchable, the sixteen single-source values below 45, the (m-1)^2 cutoff) all check (`my_aliquot.out`). THM-1370 does state the all-tournament coverage of [1, 609] ("by n = 8 it covers every odd value in [1, 609] except 7 and 21"), as the note says; the strong version is HYP-9029's (M1), not THM-1370's.
16. All 29 relative links in the note and the one in THM-4810 resolve; HYP-9420 has no relative links. Placeholders [[AUDIT]], [[AUDIT-THM]], [[F11]], [[F11-THM]], [[SLA300]], [[AUDIT-HIDX]], [[AUDIT-RIDX]], [[AUDIT-LED]] are open slots, not errors.

---

## Numbers reproduced with the auditor's own code (exact outputs)

Machine: 12 cores, g++ 15.2 `-O2 -std=c++17 -fopenmp`, Python 3.13. Sources: `my_sweep.cpp`, `my_augment.cpp`, `my_cert.cpp`, `my_bfs.cpp`, `my_paley.cpp`, `my_p11.cpp`, `my_aliquot.py`, `verify_s73_witnesses.py`, `compare_sweep.py` (all in the audit directory).

### 1. Atlas, orders 3-8 (`my_sweep.cpp`: labelled sweep; Burnside over all 11,024 non-identity all-odd permutations at n = 8, each with its own pair-orbit enumeration; DP self-checked against brute-force permutation counting on every tournament for n <= 6: `0 mismatches`)

`compare_sweep.py` against the author's `h_spectrum_sweep_20261027_n{3..8}.out`: every row (H, labelled, classes, strong classes, strong labelled) identical at every order - `n=3 ... differing rows: 0` through `n=8: author rows 320, mine 320, differing rows: 0`; `ALL TABLES IDENTICAL`. Totals: (8, 2, 1), (64, 4, 1), (1024, 12, 6), (32768, 56, 35), (2097152, 456, 353), (268435456, 6880, 6008). Spectra: `f=25 max=189 gaps=[63, 107, 119, 149, 161, 163, 165, 167, 169, 173, 177, 179, 181, 183, 185, 187]` at n = 7; `f=45 max=661 gaps=[47, 55, 67, 611, 615, 617, 619, 623, 625, 635, 647, 655]` at n = 8; |A_n| = 2, 3, 7, 19, 77, 320; |S_n| = 1, 1, 4, 13, 67, 297. Single-class values: n = 7: `1 61 85 87 125 127 147 155 157 171 175 189 (12 values)`; n = 8: `1 629 637 (3 values)`. `strong classes with H <= 200: 951` at n = 8. The fibre table F_s(v), v <= 73, summed from my tables equals the note's table and THM-4744's P_d(1) (36/36 entries). Run time n = 8: 4 m 38 s.

### 2. Augmentation (`my_augment.cpp`; canonical form = minimum code over score-sorted relabellings; two routes at order 7)

```
order 6: 56 classes, 35 strong (expect 56, 35)
order 7 via Moon chain from strong 6-classes: 353 strong classes (expect 353)
order 7 exhaustive: 456 classes, 353 strong (expect 456, 353); Moon-chain set identical: yes
order 8: 6008 strong classes (expect 6008)
order 8: strong classes with H <= 200: 951
order 9: 2919 labelled strong augmentations with H <= 200; f(9) = 75; attained: 75 81 85 87 ... 199
order 9: odd gaps in [f, 200]: 77 79 83
order 9: 532 distinct strong classes with H <= 200
order 10: 3752 labelled strong augmentations with H <= 200; f(10) = 125; attained: 125 135 145 147 153 155 159 161 165 167 169 171 175 177 ... 199
order 10: odd gaps in [f, 200]: 127 129 131 133 137 139 141 143 149 151 157 163 173
```
(20 s.) Agrees with the author's B200 output and with Busch's f(9) = 75 = 3 x 25, f(10) = 125 = 5^3. The logic (Moon 1966: a strong n-tournament, n >= 4, has a vertex whose deletion is strong; insertion lemma H(T - v) <= H(T); hence S_n ∩ [1, B] is exactly the set of strong augmentations of strong (n-1)-tournaments with H <= B) is sound; augmenting one representative per class in all 2^(n-1) ways is complete up to isomorphism.

### 3. Certificates (`my_cert.cpp`: bit k = 1 iff i -> j for the k-th pair (i, j), i < j, lexicographic; strong by two closures; H by subset DP)

```
n = 9: 1482 certificate lines, 0 bad, 1482 distinct verified values, max 3357
verified odd values as maximal runs: 75 81 [85,2881] 2885 2891 2903 [2917,2919] 2931 [2959,2963] 2967 [2971,2973] [2979,2995] [2999,3075] [3079,3085] 3095 3101 3105 3109 3159 3243 3249 3255 3267 3275 [3279,3281] [3297,3299] 3303 3333 3357
n = 10: 2925 certificate lines, 0 bad, 2925 distinct verified values, max 5999
verified odd values as maximal runs: 125 135 [145,147] [153,155] [159,161] [165,171] [175,5999]
```
So the order-9 certificate covers every odd value of [85, 2881] and 75, 81; the order-10 certificate covers every odd value of [175, 5999] and the twelve listed values; with the exhaustive [3, 609] every odd 3..6000 except 7, 21 is strong-attained at order <= 10 (and 7, 21 are absent by THM-1370; 1 only on transitive tournaments). The 1,482 order-9 values equal THM-4097's exact S_9 (see M1).

### 4. Complement-move census (`my_bfs.cpp`: explicit neighbour lists, components by repeated BFS, f by subset DP, L from sorted scores, Lipschitz check on every directed edge, diameters for n <= 6 by BFS from every state)

- n = 3: `components: 2; sizes: 4 4`; `BFS from Tr: reached 4, eccentricity 2; levels: 1 1 2`; diameter 2; `Paley P_3: ... f = 1, L = 1, dist from Tr = 1; BFS from P_3: reached 4, eccentricity 1; levels: 1 3`.
- n = 4: `components: 1`; `eccentricity 7; levels: 1 1 4 11 19 18 8 2`; `max f = 1; ... tight: 2; histogram: 0:2 1:4 2:9 3:16 4:15 5:10 6:6 7:2`; diameter 7.
- n = 5: `components: 2; sizes: 512 512`; `states whose component contains both parities: 0`; `eccentricity 6; levels: 1 1 8 62 221 197 6:22`; `(d - L) histogram: 0:10 1:31 2:128 3:141 4:152 5:30 6:20`; diameter 6; the Z_5 facts quoted in M4.
- n = 6: `components: 2; sizes: 16384 16384`; `eccentricity 7; levels: 1 1 16 314 3288 8869 3695 200`; `max f = 4; tight: 2; histogram: 0:2 1:16 2:314 3:3273 4:8839 5:3563 6:255 7:122`; `max L = 3; (d - L) histogram: 0:82 1:1405 2:4090 3:5150 4:3249 5:2031 6:255 7:122`; diameter 7.
- n = 7 (35 s): `2097152 states, 165150720 directed move edges, relation asymmetries: 0, Lipschitz violations |L'-L|>1: 0`; `components: 1`; `eccentricity 13; levels: 1 1 30 1364 35263 196755 415382 521492 465574 292614 126214 35870 6112 480`; `max f = 7; violations of d >= ceil(f/(n-1)): 0; tight: 2; histogram: 0:2 1:30 2:1364 3:35245 4:197667 5:415081 6:520133 7:465386 8:292580 9:126328 10:36170 11:6380 12:660 13:126`; `max L = 6; violations of d >= L: 0; (d - L) histogram: 0:15771 1:88123 2:227920 3:359210 4:416707 5:376509 6:284218 7:177064 8:90306 9:42656 10:13072 11:4888 12:582 13:126`; `Paley P_7: code 1956555, scores 3 3 3 3 3 3 3, H = 189, f = 7, L = 6, dist from Tr = 7; BFS from P_7: reached 2097152, eccentricity 7; levels: 1 189 20076 502365 905649 532582 122850 13440`.
All of these equal the author's `complement_moves_bfs_20261027_n{3..7}.out`. At every order 3..7 the Lipschitz violations are 0 and d >= L has 0 violations, as Theorem L requires; d >= ceil(f/(n-3)) also has 0 violations (n = 4..7).

### 5. Paley budgets (`my_paley.cpp`)

```
P_7: regular 1, f = 7  [(q^2-1)/8 = 6, q(q-1)/4 - f = 3.5, ceil(f/(q-1)) = 2, ceil(f/(q-3)) = 2, 2q-4 = 10]
P_7: H = 189
P_7: cycle debt b = 1 attained by 21 orders; f-optimal orders (E = 7): 21, of which connected backward graph: 21; min E over all orders 7
P_11: regular 1, f = 20  [(q^2-1)/8 = 15, ... ceil(f/(q-1)) = 2, ceil(f/(q-3)) = 3, 2q-4 = 18]
P_11: H = 95095
P_11: cycle debt b = 10 attained by 2420 orders; f-optimal orders (E = 20): 2420, of which connected backward graph: 2420; min E over all orders 20
P_19: regular 1, f = 64  [(q^2-1)/8 = 45, q(q-1)/4 - f = 21.5, ... 2q-4 = 34] (0.1 s)
P_23: regular 1, f = 92  [(q^2-1)/8 = 66, q(q-1)/4 - f = 34.5, ... 2q-4 = 42] (1.2 s)
P_31: regular 1, f = 180  [(q^2-1)/8 = 120, q(q-1)/4 - f = 52.5, ... 2q-4 = 58] (192.2 s)
G1: f(C3) = 1, f(P7) = 7, f(P7[C3]) = 70 (G1 predicts 7*1 + 9*7 = 70), f(C3[P7]) = 70 (G1 predicts 3*7 + 49*1 = 70)
```
The brute force over all 11! orders confirms the note's deduction "every f-optimal backward graph of P_11 is connected" directly (2,420 f-optimal orders, all connected), which the author's program did not count separately. The deficits 3.5, 7.5, 21.5, 34.5, 52.5 and the ratios to q^(3/2) (0.19, 0.21, 0.26, 0.31, 0.30) follow.

### 6. P_11 depth one and two (`my_p11.cpp`, 18 s)

```
P_11: f = 20, 2L = 30, H = 95095
depth one: 95095 images, 95095 distinct labelled tournaments
depth one: exact minimum f = 14; f histogram: 14:55 15:4235 16:28325 17:58740 18:3740
depth one: minimum 2L = 28; 2L histogram: 28:95095
endpoint reduction validated on first 25 images by explicit path enumeration: 0 mismatches
depth two: 8349303325 path pairs (= sum of H over the depth-one images); minimum 2L = 26 (L = 13); 2L histogram: 26:6621199255 28:1618277760 30:109826310
consequences: every depth-two state has f >= L >= 13; forest frames need f <= 10: none at depth two
```
Identical to the author's `paley11_two_moves_20261027.out` (8,349,303,325 pairs; 26: 6,621,199,255; 28: 1,618,277,760; 30: 109,826,310; depth-one minimum 14). Hence f >= 13 at depth two, no state with f <= 10, no forest frame at depth two, L >= 12 at depth three and L >= 11 at depth four by Theorem L, so at least five moves to a forest frame and fifteen to the transitive class.

Session73's witnesses (`verify_s73_witnesses.py` on the committed receipts): `P is a Hamiltonian path of P11: True; tau_P(P11) == receipt intermediate_rows: True; Q is a Hamiltonian path of tau_P(P11): True; tau_Q(tau_P(P11)) == receipt final_rows: True; f(intermediate) = 15 (receipt says 15); f(final) = 13 (receipt says 13); final state sorted scores [4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6] 2L = 26 L = 13; R is a Hamiltonian path of the depth-two state: True; tau_R(depth-two state) == receipt rows: True; f(depth-three state) = 12 (receipt says 12); ... L = 12`. So the minimum f at depths 1, 2, 3 is exactly 14, 13, 12 (= 15 - k), with L tight at each witness, as the note and HYP-9420 (4) state.

### 7. Aliquot fibres (`my_aliquot.py`, sieve to 2700, completeness by n <= (m-1)^2)

`untouchable m <= 45: [2, 5]`; `single-source m <= 45: [3, 4, 7, 9, 10, 11, 12, 18, 24, 26, 28, 30, 34, 36, 38, 39]`; `s^-1(7) = [8]; s^-1(21) = [18, 51, 91]; s^-1(28) = [28]; s^-1(3) = [4]; s^-1(13) = [27, 35]; s^-1(4) = [9]`; `52 untouchable: True`. The fibre lists for m = 2..45 are identical to the author's `aliquot_fibres_to_45_20261027.out`.

### 8. Proof checks (by hand)

- Lemma 1: (all \ E(P)) xor (all \ E(Q)) = E(P) xor E(Q); the same P twice is the identity. Correct.
- Theorem P: (n-1)(n-2)/2 toggled pairs per move; even iff n = 1, 2 (mod 4); converse toggles n(n-1)/2, odd iff n = 2, 3 (mod 4). Correct (and confirmed by the census: parity classes = components at n = 5, 6; connected at n = 4, 7). The obstruction is labelled-only (M4).
- Theorem C2: correct for the transitive class (minor 2); the disjointness remark when f = 2(n-1) is correct.
- Proposition B: f(tau_P T) <= f(T) + n - 1 correct (relative to the reversed best order the backward set changes by E(P)); the score bound: for any order, s(v_i) = (n-1-i) - (backward arcs into v_i) + (backward arcs out of v_i), so sum_i |s(v_i) - (n-1-i)| <= 2|D_r|, and the sorted matching minimises the left side over bijections; f >= L follows; regular q: sum |(q-1)/2 - i| = (q^2-1)/4. Correct; the threshold is q >= 15 (minor 5).
- Theorem L: scores of tau_P(T) are the reflected scores with +1 at the start and -1 at the end of P (internal vertices lose one out-arc and gain one); reflection preserves sum |s_(i) - i|; the sorted L1 matching is 1-Lipschitz per coordinate. Correct; confirmed exhaustively at n = 3..7 (0 violations on 165,150,720 directed edges at n = 7).

---

## Prior art (with URLs)

Literature:
- S. H. Chan, I. Pak, *Computational complexity of counting coincidences*, Theoret. Comput. Sci. 1015 (2024) 114776; author PDF https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf ; arXiv https://arxiv.org/abs/2308.10214 ; DOI https://doi.org/10.1016/j.tcs.2024.114776 . Section 6.5, Conjecture 6.2: "Function h takes all odd values except 7 and 21"; credited to MathOverflow user bof and Gordon Royle (March 2016); "Royle reported that the conjecture holds for all odd integers up to 80557".
- MathOverflow 232751, *The number of Hamiltonian paths in a tournament* (2016; Royle's answer updated January 2022): https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament (not fetchable from this environment; content taken from Chan-Pak and from the repository's `spectrum_basins_20261020e/spectrum.md` and `rotation_owners_spectrum_20261020e/spectrum.md`, which report coverage of all allowed odd values strictly below 80,557 by tournaments on at most eleven vertices, and the strong-order argument for the omission of 21).
- A. H. Busch, *A note on the number of Hamiltonian paths in strong tournaments*, Electron. J. Combin. 13 (2006) N3: https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3 (f(n) = 3·5^((n-3)/3), 5^((n-1)/3), 9·5^((n-5)/3) for n = 0, 1, 2 mod 3; so f(9) = 75, f(10) = 125, confirming section 2's floors). J. W. Moon, *The minimum number of spanning paths in a strong tournament*, Publ. Math. Debrecen 19 (1972) 101-104; J. W. Moon, *On subtournaments of a tournament*, Canad. Math. Bull. 9 (1966) 297-301 (Theorem 2, s(n, n-1) = 2).
- OEIS A003141 (maximum over n-vertex tournaments of the minimum number of arcs to reverse): https://oeis.org/A003141 - 0, 0, 0, 1, 1, 3, 4, 7, 8, 12, 15, 20, 22, 28; references J.-C. Bermond, *Ordres à distance minimum d'un tournoi et graphes partiels sans circuits maximaux*, Math. Sci. Hum. 37 (1972) 5-25; Sanchez-Flores, Neumann-Lara, Bermond, Graphs Combin. 10 (1994); I. Charon, O. Hudry, *An updated survey on the linear ordering problem for weighted or unweighted tournaments*, Ann. Oper. Res. 175 (2010) 107-158 (score-based lower bounds for the Slater index are classical in this survey's area; the exact form f >= (1/2) sum |s_(i) - i| was not located online during this audit).
- OEIS A000568 (tournament classes 1, 1, 2, 4, 12, 56, 456, 6880, 191536) and A051337 (strong classes 1, 0, 1, 1, 6, 35, 353, 6008, 178133): the "known sequences" the note's totals match.

Repository (all uncited by the package): THM-4744 (complete strong inventories P_d for odd d < 625, orders <= 12; strong-untouchables below 625 = {1, 7, 21}; P_d(1) = 1 only at 3, 5, 11, 13 up to 623), THM-4097 (the exact 1,482-value S_9, every odd of [85, 2881]), THM-4102 / THM-4104 / HYP-9029 (strong witnesses for every odd in [65, 80265]; all allowed odd values through 80,405), `05-knowledge/results/h_spectrum_multiplicities_20261021.md` and `h_multiplicities.out` (exhaustive class census by H through order 9; the twelve single-class values at order 7), `session24_research_synthesis_20261021.md` (strong-atom inventory for H <= 73), THM-115's banner and the two 2026-10-20 spectrum notes (Chan-Pak / Royle), THM-1862 / THM-1936 (spectrum = multiplicative closure of strong values), session73 ports report section 5 (G3: |f(tau_P T) - f(T)| <= n - 3). For the complement grammar itself I found no earlier reachability or component census in the repository (the closest objects, `strong_atom_continuation_vitali_20261026/README.md` and the session50 metagraph packets, define the same complement on free-bit masks but study responses, not reachability), and no earlier forward-arc parity invariant or score-monotonicity statement; those parts of THM-4810 (C)(D) appear to be new within the repository, with the caveat that the score bound f >= L is elementary and probably classical.

---

## What I could not check

- The MathOverflow page itself (domain not reachable from this environment); Royle's 80,557 statement is quoted via Chan-Pak and the repository's two spectrum notes.
- Whether f(P_19) = 64, f(P_23) = 92, f(P_31) = 180 appear in the literature (web search only; nothing found), and whether the score-sequence bound f >= (1/2) sum |s_(i) - i| is explicitly in Charon-Hudry or Slater (not located).
- The author's depth-two census was reproduced exactly through the endpoint reduction (validated on 25 images by explicit enumeration), not by re-enumerating the 8.3 x 10^9 pairs; the f-values of individual depth-two states were not recomputed (unnecessary: L >= 13 gives f >= 13, and session73's f = 13 witness was verified).
- The diameter of the move graph at order 7 (the note only claims >= 13; my BFS from Tr gives eccentricity 13, no all-pairs run at n = 7).
- Beyond the first ten printed values, the identity of THM-4744's order-10 strong value list below 225 with S63's exact S_10 ∩ [1, 200] (the ten printed values agree; the full list is not printed in `strong_ears_225.out`).
- The literature on Paley tournaments attaining A003141's maxima at n = 7 and 11 specifically (the OEIS entry does not name the extremal tournaments).
