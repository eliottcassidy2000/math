# Independent audit: THM-4904 and HYP-9494 (kind-pasteur-2026-11-07-S14)

Audited checkpoint 58a5a8a3d1. The packet is `05-knowledge/results/root_change_carmichael_chamber_20261107/`. The canon files are
`01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md` and
`05-knowledge/hypotheses/HYP-9494-root-change-carmichael-hierarchy-two-root-change-coverage-c1-depth-three-adic-factor.md`.

I wrote every program in this folder from the definitions. The packet's C files and `literal_witnesses.py` were read
only to compare definitions and outputs; no code was reused. No file outside this folder was edited.

## Verdict: PASS WITH CORRECTIONS

The mathematics of THM-4904 (1)-(5) is correct, and every number the packet reports was reproduced except one. One
conjecture, HYP-9494 (C), is false as stated (M4); P3, which concerns K = 2, is unaffected. My own letter-by-letter
enumerator, which uses no closed form for the final block, reproduces:

- the full two-root-change census to 10^9, so P1-P4 are reproduced and P1 and P3 are independently confirmed;
- the depth table, with my own uncapped search;
- the literal witnesses for 1877 and 1897;
- the chamber, parity and singleton data;
- the exploratory 3-adic means.

The pre-registration is honest. The pre-registration commit was pushed at 11:00:57 and the raw 10^9 output was written
at 11:26:28. `rc_census.c` is byte-identical to the registered file.

There are four majors:

- **M1:** an unproved "sharp for K = 0" sits inside a PROVED theorem.
- **M2:** prior work is not credited, and Lemma 4 is superseded by a stronger result already in the repository.
- **M3:** one FINITE-EXACT count is wrong: 3,983 should be 4,083.
- **M4:** HYP-9494 (C), stated for every K >= 2, is false. M_7(1897) = 1, and the packet's own `depth_K7_1e6.txt`
  shows it.

## MAJOR findings

**M1. "Sharp for K = 0" is stated as part of PROVED THM-4904 (5), but only data supports it.**

- Theorem 5 proves an upper bound.
- "Sharp" asks for #{1 < N <= X : M_0(N) = 1} >> sqrt(X). That is HYP-9494 (A), which is OPEN.
- What can be proved is that the K = 0 chamber itself has Theta(sqrt X) words. U at most doubles per letter, so
  F <= 2U <= 2^(|w|+2). Every A in {O,Z}^a with 2a + 2 <= log2 X therefore gives a chamber word of value <= X, which is
  at least sqrt(X)/2 - 1 words. At most 2 sqrt(X) such words exist (Appendix B).
- So no better bound can come from counting chamber words. Nothing proves that the singletons fill a positive fraction
  of the chamber.

Exact fixes:

- **THM-4904 title.** Replace "(5) singleton values number O_K(sqrt(X) (log X)^floor(K/2)), sharp for K = 0 - density
  zero" with "(5) singleton values number O_K(sqrt(X) (log X)^floor(K/2)) - density zero (the K = 0 chamber has
  Theta(sqrt X) words; that the singletons are also of order sqrt X is HYP-9494 (A))".
- **THM-4904, Reading.** Replace "The bound is sharp for K = 0: about 0.49 sqrt(X) singletons, about half of the 1,023
  = 2^10 - 1 chamber words below 2^20." with "For K = 0 the chamber itself has Theta(sqrt X) words. Observed singletons
  are 500 of the 1,023 chamber words below 2^20 and 4,083 of the 8,191 below 2^26 (0.498 sqrt X); see HYP-9494 (A)."
- **PROOF section 4, heading.** Replace "with a sharp square-root bound" with "with a square-root chamber bound".
- **PROOF section 4, closing text.** Replace "The bound is sharp for K = 0. Below 2^20 ... about 0.49 sqrt(X)." with:
  "For K = 0 the chamber bound cannot be improved by counting chamber words. U at most doubles per letter, so every A in
  {O,Z}^a with 2a + 2 <= log2 X gives a chamber word of value <= X. Whether the singletons are themselves of order
  sqrt(X) is HYP-9494 (A). Below 2^20 the 0-chamber has exactly 1,023 words and 500 singletons. Below 2^26 it has 8,191
  words and 4,083 singletons, 0.498 sqrt(X)."
- **README item 2.** Change "sharp for K = 0" to "for K = 0 the chamber has Theta(sqrt X) words (singletons: HYP-9494
  (A))".
- **CARMICHAEL_MAP, dictionary row.** Change "Sharp for K = 0 (about 0.49 sqrt(X) observed)" to "the K = 0 chamber has
  Theta(sqrt X) words; about 0.50 sqrt(X) singletons observed (HYP-9494 (A))".

**M2. Priority. Parts of THM-4904 restate, or are weaker than, results already in the repository, and the packet cites
none of them.** THM-4904 is typed "PROVED (new): (1)-(5)".

- **(a) The K = 0 case of (2), (3) and (5) is Session 206.** The source is
  `session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md`, sections 1-2 (audited, 2026-11-06), recorded in THM-4888
  (B), the S13 packet.
  - Its converse involution T(2^(n-1)+a) = 2^n - 1 - rev_(n-1)(a) is exactly the dagger on {O,Z} prefixes. The tail
    word b_0 ... b_(n-2) Z goes to dagger(b_0 ... b_(n-2)) Z.
  - H(T(m)) = H(m) is proved there by reversing the tournament.
  - Its fixed points, "the n-1 lower bits equal their reverse complement", number 2^((n-1)/2). These are exactly the
    K = 0 chamber A A^dagger.
  - Its orbit count gives the O(sqrt X) fixed-parameter bound.
  - The K = 0 parity law, the dagger-fixedness of K = 0 singletons and their O(sqrt X) count follow at once.
  - What is new is the extension to K >= 1 through Session199's adjoint.
- **(b) Theorem 3 (3) reuses Session200 Theorem 3.** The "neither begins nor ends with J" step is the argument of
  `session200/carmichael-fibers/PROOF.md`, Theorem 3: a terminal J or O gives a twin value, and the dagger moves the
  start to the end. Here it is run inside a root-change budget. The PROOF cites only Session200's Theorem 1, the twin.
- **(c) Lemma 4 is strictly weaker than Session200 root-weighted-height.** The source is
  `session200/root-weighted-height/PROOF.md`, section 1 (audited, 2026-11-05).
  - W = U - d satisfies O: W' = 2W + d, Z: W' = 2W, J: W' = W + d, with W(t) = 1.
  - Hence W >= 2^(n(w)) and F >= 2^(n(w)) + 2 for every word, with no factor lost per root change.
  - This replaces c0^(#J+1). Theorem 5 then holds with C_K = 0 and n(A) < (1/2) log2 X.
  - On my 76,333,010 states (and 916,133,475 at 2^26) the ratio U / (2 c0^(#J+1) 2^n) never falls below 1.1298. The
    ratio W / 2^n reaches 1 exactly, at the words Z^n.

Exact fix for the THM-4904 status: replace "PROVED (new): (1)-(5), given ..." with:

> PROVED: (1)-(5), given ... New for K >= 1. The K = 0 case of (2), (3) and (5) is Session 206's converse involution T
> (ROOT_OBSTRUCTION sections 1-2; THM-4888 (B)): T is the dagger on {O,Z} prefixes, its fixed points form the K = 0
> chamber, and its orbit count gives O(sqrt X). The J-exclusion in (3) is the argument of Session200 carmichael-fibers
> Theorem 3, run under a root-change budget. (4) is superseded by the KNOWN bound W = U - d >= 2^(#O/Z) (Session200
> root-weighted-height section 1), with which Theorem 5 holds with C_K = 0.

Further fixes:

- Add the same three sentences at the ends of PROOF sections 2, 3 and 4.
- Add `session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md` and `session200/root-weighted-height/PROOF.md` to
  `related:`.
- Optional: replace Lemma 4 by the W-bound. Appendix B gives the two-line proof.

**M3. A wrong FINITE-EXACT count, and the constant derived from it.**

- PROOF section 4 says "Below 2^26 there are 3,983 singletons, about 0.49 sqrt(X)".
- The packet's own window list sums to 4,084 including N = 1, that is 4,083 above 1. That list is PROOF section 6 and
  `explore/rc_K0_2p26_windows.txt`.
- My independent count is **4,083** (`chamber_audit_K01_2p26.txt`), from 8,191 = 2^13 - 1 chamber words.
- So c = 4,083/8,192 = 0.498. The singleton share of the chamber is 4,083/8,191 = 0.4985, which matches the root's
  Poisson miss rate 0.4985 better than the text claims.

Fixes:

- PROOF section 4: change "3,983" to "4,083" and "about 0.49 sqrt(X)" to "about 0.50 sqrt(X) (0.498 at 2^26)".
- PROOF section 7 (A): change "c about 0.49" to "c about 0.50".
- HYP-9494 title: change "c about 0.49 (half the chamber)" to "c about 0.50 (4,083 of the 8,191 chamber words below
  2^26)".
- HYP-9494 (A) body: change "the total is about 0.49 sqrt(X)" to "the total is 4,083 = 0.498 sqrt(X) below 2^26".
- THM-4904 Reading: fixed by M1.
- README item 5: change "about 0.49 sqrt(X)" to "about 0.50 sqrt(X)".
- Record a MISTAKES entry.

**M4. HYP-9494 (C) is refuted as stated.**

- The HYP-9494 title says "(C) for K >= 2 no value N > 1 has exactly one word". PROOF section 7 (C) says "K >= 2
  grammars have no singleton above N = 1".
- But **M_7(1897) = 1**. The only counted word is ZZJJJJJJJOOZ: its prefix w = ZZJJJJJJJOO is the one prefix with
  exactly 7 letters J and F = 1897, and 1897 has depth 7, so no word with fewer J's exists.
- w is dagger-fixed, of the form A J A^dagger with A = ZZJJJ, exactly as THM-4904 (3) predicts.
- The packet's own `depth_K7_1e6.txt` shows it: window [1024,2048), c0..8 = 3,**1**,0,0,2,0,0,0,2. My
  `census_K7_1e6.txt` and `k_singletons_2048.out` agree.
- The mechanism is general. If depth(N) = K, then M_K(N) = B_K(N), with no doubled 2A_(K-1) term. A depth-K value
  whose layer of exactly-K-J prefixes is a single prefix, necessarily dagger-fixed, is therefore a singleton of the
  K-grammar. Nothing in the packet's heuristic excludes this for K >= 3. For K > depth(N) the value is never a
  singleton, because A_(K-1)(N) >= 1.
- For comparison, M_7(1877) = 4 (no dagger-fixed word).

Below 10^9 this is the only word-singleton above 1 with K >= 2:

- K = 2: none (P3).
- K = 3: none. Its singletons would need depth 3, and my K = 3 census to 4*10^6 shows none.
- K = 4, 5, 6: none (census to 10^6).
- K >= 8: none. No value below 10^9 has depth >= 8, so M_K >= 2 on F_<=7.

Fixes:

- **HYP-9494 title.** Change "(C) for K >= 2 no value N > 1 has exactly one word" to "(C) for K = 2 no value N > 1 has
  exactly one word (for larger K singletons can sit at values of exact depth K: M_7(1897) = 1, the only word-singleton
  above 1 with K >= 2 below 10^9)".
- **HYP-9494 (C) section and PROOF section 7 (C).** Make the same change, and say why (C) is only for K = 2: M_K = B_K on
  F_<=K minus F_<=K-1.
- **CARMICHAEL_MAP, dictionary row.** Change "K = 2 has none through 10^9 (pre-registered P3)" to "K = 2 has none
  through 10^9 (pre-registered P3); K = 7 has the singleton 1897".
- **MISTAKES.** Record a MISTAKES entry: the counterexample sat in the packet's own output.

## MINOR findings

**m1.** Theorem 5 and THM-4904 (5) are off by one at N = 1. M_K(1) = 1 for every K, from the empty word, but 1 is not
a chamber value. So "#{N <= X : M_K(N) = 1} <= #{chamber words of value <= X}" fails for small X, for example K = 2,
X = 1.

- Write "#{1 < N <= X : M_K(N) = 1}".
- Say "chamber words w with F(t.w) <= X". "Value" otherwise reads as D(t.w).

**m2. A missed strengthening (optional).** The identity JZ = OJ (Session 205; THM-4888 (D)1) pairs words with the same
number of J's, and it gives two consequences:

- A singleton prefix contains neither of the factors JZ and OJ.
- So for K = 1, A ends in Z. For K = 2, B ends in Z, and C is empty or begins with O.

Checked: no singleton prefix contains JZ or OJ, among 520 singletons below 2^20 (K <= 2) and 4,122 below 2^26
(K <= 1). This is worth a remark in section 3 and in HYP-9494 (C)'s proof routes.

**m3. Corollary 5.1 attribution.** Two results credited elsewhere belong to Session198 c1-monotone
(`session198/c1-monotone/PROOF.md`, 2026-11-05 12:07, frozen as Session199's `C1_SOURCE_READOUT.md`):

- The readout "C1 spectrum = F(S)" is its section 1. Session199 section 1 calls it inherited.
- The 29-gap list through 50,001 is its section 3. It predates THM-4884 (2b), 15:50 the same day.

Fixes:

- Write "This extends Session197 (to 5,001) and Session198 c1-monotone section 3 (to 50,001; replayed in THM-4884
  (2b))".
- Cite Session198 c1-monotone section 1 beside Session199 Theorem 1.
- For THM-4884's owner: its claim "(2b) ... the identification of the two is what is new" is also pre-empted by
  Session198.

**m4. Root-change count of the depth-7 witnesses.** The literal witnesses for 1877 and 1897 are closed by J, so each
construction has eight letters J. Say "source words with at most seven root changes (the closing J is not counted,
Session199's convention)" in README item 4 and PROOF section 5.1.

**m5. The two-run chamber claim in PROOF section 6 is false.** It says "Session200's diagonal family is exactly the
K = 0 chamber of two-run words".

- The two-run K = 0 chamber is {O^a Z^a} together with {Z^a O^a}.
- The second family has F(t.Z^a O^a) = (2^a + 1)^2 = 9, 25, 81, 289, 1089, ..., and 25 and 1089 are K = 0 singletons.
- Write: "Session200's diagonal O^a Z^a is the half of the two-run K = 0 chamber that begins with O. The other half,
  Z^a O^a, has values (2^a + 1)^2; 25 and 1089 are singletons."

**m6. A mislabelled column.** VERIFY_CHAMBER_2p20.txt and `verify_chamber.c` print "with word w (value of wZ)" but list
the counted words wZ. For example, "3:Z" means w is empty. Relabel the column "counted word wZ".

**m7. README overstates what `verify_chamber.c` checks.** It claims an independent check of THM-4904 (1), but the file
computes M_K by the formula [N=1] + 2A + B. (1) is tested only through the agreement of its window histograms with
`rc_census`. Say so.

**m8. Wording of the finite checks, and a proof of (phys).**

- "Re-verified exhaustively here: 7,174,452 words ..., 76,333,010 states" in the THM status should read "re-checked on
  ...". These are finite checks.
- (phys) has a five-line self-contained induction (Appendix A), so THM-4904 need not rest on a citation for it.

**m9. Wording in the exploratory parts.**

- HYP-9494 (E) says "This explains why the misses and near-misses lie in N = 1 (mod 3)". Write "This is consistent
  with ...".
- In PROOF section 5.3, change "This accounts for" to "This is consistent with".
- The class means differ by 6% (199.96 against 211.73). That does not by itself explain a 104:26 split. (E) stays
  EXPLORATORY.
- The working tree also holds an uncommitted update to HYP-9480, made at 11:48 after the audited checkpoint. It says the
  bias "comes from" the 3-adic factor; apply the same fix there. That update and `explore/diag_root.*` (also 11:48,
  uncommitted) were not audited, and I did not touch them.

**m10. CARMICHAEL_MAP precision. No misquote was found.** Every bound, year, theorem number and the Carmichael-Klee
lemma agree with Ford's text. Four points:

- **(i) The text read is the 2013 arXiv v2 revision of Ford 1998.** It contains post-1998 material: Ford-Konyagin 1999,
  Ford 1999, the corrected Theorem 3 and Baker-Harman-Pintz 2001. Write "Ford 1998, read in its arXiv v2 (2013)
  revision", here and in README item 1.
- **(ii) Pomerance 1974 is mis-phrased.** "Then n is a counterexample" should read "then phi(x) = phi(n) has the unique
  solution x = n (phi(n) is a counterexample)". Neither this result nor the "very improbable" remark is in Ford, so
  cite Pomerance, Proc. AMS 43 (1974) 297-298, directly.
- **(iii)** Carmichael's 10^37 is from 1922 (Bull. AMS 28, 109-110).
- **(iv)** The 2p rule also appears in Session218's README, as the fiber of 2q.

The 2025 NNTDM reference is correct (31(3), 504-534): authors Shannon, Shiue, He and Saito, cases 2^t, 2p^s, 4p and
2pq.

**m11. A prior depth table.** Session199 hidden-structure PROOF section 6 already tabulated minimum J cost, which is
root-change depth, for the balanced sub-readout with at most three J's, to 1,000,001:

- cost 0/1/2/3 counts 213,873 / 266,697 / 19,137 / 220;
- 73 gaps;
- 293 misses with at most two J's, the largest 879,001.

The packet cites GAPS_THREE.json only. Credit the depth table there.

**m12. Make the pre-registration evidence durable.** The proof that no output existed when P1-P4 were pushed currently
lives only in a scratch directory and the local reflog. Its evidence:

- the reflog push of 6d0e88f0e8 at 11:00:57;
- the raw `rc_1e9.txt` and `.err` written at 11:26:28;
- `real 30m54.574s`, so the run started at about 10:55:34;
- all stdout is printed after the enumeration.

Add these four facts to PREREG_OUTCOMES.txt or the README.

## Recomputed numbers

The code is in this folder. COMPARISON.txt lists 56 automated comparisons: 55 match and 1 mismatches (M3). The packet's 10^9 run took 1,845.7 s. Mine
ran segmented on 12 threads, with run times in the err files: (2^26, 10^9] took 875.5 s and 10^8 took 65 s.

| Quantity | Packet | Auditor (own code) | Match |
| --- | --- | --- | --- |
| Session199 identities O^T Q = QZ, Z^T Q = QO, J^T Q = QJ | cited | exact integer check | yes |
| (adj) F(t.w) = F(t.w^dagger) | 7,174,452 words, length <= 14 | 64,570,080 words, length <= 16, via row-vector transfer (eq. 4) | yes |
| prefixes with #J <= 2, F <= 2^20 | 76,333,010 | 76,333,010 | yes |
| (phys) failures | 0 | 0 on 76.3 M and 916.1 M states (also E >= 0, d >= E, D-E+d <= U, D and d odd, E = U mod 2) | yes |
| Lemma 4 failures; min ratio | 0 | 0; 1.1298. W >= 2^n: 0 failures, min ratio 1 | yes |
| c0 | 0.44257 | 0.442573109396 | yes |
| (1) M_K = [N=1] + 2A_(K-1) + B_K, K = 0, 1, 2 | by construction | 0 mismatches, words counted from the child D coordinate (2^20; K <= 1 to 2^26) | yes |
| parity law | 0 failures at 2^20 | 0 failures: K <= 2 at 2^20, K <= 1 at 2^26, K = 2 at 2^24; Fix computed two ways | yes |
| singletons above 1 below 2^20, K = 0/1/2 | 500 / 20 / 0 | 500 / 20 / 0, identical words | yes |
| chamber words below 2^20, K = 0/1/2 | 1,023 / 901 / 6,344 | 1,023 / 901 / 6,344 | yes |
| Theorem 3 on every singleton | 0 failures | 0 failures (520 at 2^20, 4,122 at 2^26); no JZ or OJ factor | yes |
| window histograms at 2^20 (60 rows) | VERIFY_CHAMBER | identical | yes |
| K = 0 singletons below 2^26 | **3,983 (text)**; window sum 4,083 | **4,083**; chamber 8,191 | **no (M3)** |
| K = 1 singletons below 2^26 | window sum 39 | 39 | yes |
| K = 0 and K = 1 windows to 2^26 (n, mean255, c0..8) | explore/*_windows.txt | identical, 26 + 26 rows | yes |
| odd M_2 (N >= 3) to 2^20; above 2^15 | 2,024; 1,752 | 2,024; 1,752 | yes |
| odd M_2 (N >= 3) to 2^24; above 2^15 | 10,055; 9,783 | 10,055; 9,783 (1,835,066,963 prefixes visited) | yes |
| word-singletons above 1 with 2 <= K <= 8 (HYP-9494 (C)) | none claimed | M_7(1897) = 1; no others (N <= 2048 directly; to 10^9 by the depth argument in M4) | **no (M4)** |
| F_<=0/1/2 misses to 10^8 (largest) | 25,095,920 / 224,380 / 140 (99,999,995 / 99,998,133 / 3,780,967) | identical | yes |
| multiplicity windows below 2^26 (26 rows: n, min, min_mod3, mean255, sat255, c0..8) | rc_1e9.txt | identical | yes |
| F_<=2 misses above 10^5; N > 10^5 with M_2 <= 2 | 15 values; 39 entries | identical | yes |
| odd N in (2^26, 10^9] with M_2 = 0 (P1) | 0 | 0, so exactly 140 misses below 10^9 | yes |
| window minima of M_2, 2^26..10^9 (P2) | 10, 8, 8, 14 | 10, 8, 8, 14 | yes |
| minima by N mod 3 (P4) | (34,10,10) (42,8,16) (40,12,8) (48,14,20) | (>=15,10,10) (>=15,8,>=15) (>=15,12,8) (>=15,14,>=15); 4-bit cap | yes (P4 outcome reproduced) |
| M_2 = 1 above 1 up to 10^9 (P3) | none | none | yes |
| M_2 >= 4 for odd N in (2^22, 10^9) | yes | yes | yes |
| M_2 <= 12 in (2^22, 2^28): count, mod 3, minima | 130; 104/26/0; 4 at 5,661,961, 8,203,243, 11,340,277, 20,195,929, 21,375,313 | identical. Beyond the packet, in [2^28, 10^9]: only 277,438,727 (8) and 301,189,927 (12) | yes |
| F_<=k misses at 10^6, k = 3..7 | 44 / 33 / 31 / 31 / 29 (largest 27,781 / 3,523 / 2,647 / 2,647 / 2,647) | identical, with lists | yes |
| F_<=3 at 4*10^6; F_<=2 list; depth-3 count | 44, largest 27,781; 140; 96 | identical; 140; 96 | yes |
| no word at all (any number of J), N <= 30,001 | Session197's 29 | the same 29 | yes |
| exact depth >= 4 (no cap on J), N <= 30,001 | 11 at depth 4, 667 and 3523 at 5, none at 6, 1877 and 1897 at 7 | identical; no word with <= 6 J has value 1877 or 1897 | yes |
| literal witnesses (18 words): depth, #J, F, vertices, (H, C) | LITERAL_WITNESSES.txt | identical. 1877 = ZJJZJJJJJOO and 1897 = ZZJJJJJJJOO give 14 vertices and (H, C) = (N, 1); my own minimal witnesses (for example 1877 = ZZJJJJJOJJO) also give (N, 1) | yes |
| builder against state machine, all 29,524 words of length <= 9 | (packet: witnesses only) | H(w) = D, H(w - root) = d, H(wJ) = F, C(wJ) = 1; 0 failures | yes |
| mean M_2 by v3(N) on [2^23, 2^24) | 205.8, 332.0, 457.3, 510.7, 556.8, 614.2, 618.7 | identical (exact 32-bit; max M_2 below 2^24 is 3,792, so the packet's 16-bit counters did not wrap) | yes |
| mean M_2 by N mod 3 / mod 9 | 382.16 / 199.96 / 211.73, ... | identical, all windows from 2^20 | yes |
| rc_census.c against the registered file; outcome evaluation | registered | sha256 cb2dc660... equal; evaluate_prereg.py reproduces PREREG_OUTCOMES.txt | yes |
| totient statements (2*7^k, the 2p rule, five of nine doubles, about 30%, A005277 start, A_sigma(3) = 1) | stated | all true (carmichael_checks.out) | yes |

## Task notes

**1. Proofs.**

- **Theorem 1.** It is tautological once one notes that Z and J both set D' = F(s), and it is correct.
- **Census word set.** The census counts the multiset {empty word} together with {wZ : #J(w) <= K} and
  {wJ : #J(w) <= K - 1}.
  - `rc_census.c`'s `tail` marks h = D(J(s)) and the final blocks that end in Z.
  - `block` enumerates every intermediate O/Z continuation, including a trailing O before the next J.
  - Words ending in O are never counted.
  - My enumerator counts that multiset directly. It matches every window statistic of rc_census, including saturated
    means and sat255, at 2^20, 10^8 and 10^9. So the closed-form tail implements M_K exactly.
- **Theorem 2.** It is correct. The dagger acts on the prefix w, whatever its last letter, and the counted word is wZ.
  (adj) concerns prefixes, so the case where w ends in Z needs no special treatment. The dagger of the counted word,
  O w^dagger, is never used.
- **Theorem 3.** It is correct, including the shape B J C C^dagger J B^dagger with B nonempty.
- **Lemma 4 and Theorem 5.** Both are correct. The dependence O_K(sqrt X (log X)^floor(K/2)) is right, with
  C_K = (K + 1) log2(1/c0)/2, about 0.59(K + 1), or C_K = 0 with the W-bound. "Sharp for K = 0" rests on data only
  (M1).

**2-3. Computations and the corollary.**

- **The data ranges combine without a gap.**
  - The K = 2 census to 10^9 covers (3,780,967, 10^9].
  - The K = 3 census to 4*10^6 has no miss above 27,781.
  - The K = 7 census to 10^6 has no miss above 2,647.
  - Hence every odd N in (2647, 10^9] lies in F_<=7.
- **The 29 are not C = 1 counts.** They are at most 2,647, below 5,001, so Session197's complete classification
  excludes them. My uncapped search also finds no word at all for them.
- **Realisation.** Each value is realised by the J-closed word, which my 29,524-word builder check supports.
- **Citation.** Session199 Theorem 1 is an acceptable citation, but the readout and the 50,001 list are Session198's
  (m3).

**4. Pre-registration.** It was handled correctly (m12 only asks to record the evidence).

- 6d0e88f0e8 is an ancestor of 58a5a8a3d1.
- The reflog shows it pushed at 11:00:57 and the outcome commit pushed at 11:42:13.
- `git diff 6d0e88f0e8 HEAD -- rc_census.c` is empty.
- The raw scratch output (`scratchpad/s14/rc_1e9.txt`, mtime 11:26:28) is identical to the packet's rc_1e9.txt
  modulo CR.
- The census prints nothing to stdout before the enumeration ends.
- The outcomes follow the registered definitions, including the P4 tie at (34, 10, 10).

**5. Literature.** See m10. The elementary totient claims are all true.

**6. Priority.** See M2 and m11. I searched session199-206, THM-4884, HYP-9477, HYP-9480, THM-4888 and the S13 packet.

- **No parity law for word multiplicities was found anywhere.** Session200 section 4 only warns that the dagger "may
  fix a word".
- **New for K >= 1:** the parity law, the dagger-fixed chamber with J, the sqrt(X) (log X)^floor(K/2) density bound
  and the depth table of the full C1 readout through 10^9.

**7. Exploratory numbers.** All reproduce exactly, from an exact 32-bit census to 2^24 and the near-miss lists. The
wording should stay EXPLORATORY (m9).

## Appendix A: a self-contained proof of (phys)

The set of states (D,E,U,d) satisfying the following five inequalities contains t and is closed under O, Z and J:

    E >= 0,   E <= d,   D - E >= d >= 1,   D - E + d <= U,   U >= 2d.

D <= U follows from the fourth and second inequalities.

- **t = (1,0,2,1):** 0 <= 1, 1 >= 1, 2 <= 2, 2 >= 2.
- **O, giving (D, 0, 2U, d):**
  - D' - E' = D >= d;
  - D + d <= U + U/2 <= 2U;
  - 2U >= 2d.
- **Z, giving (D-E+U, d, 2U-d, d):**
  - E' = d' = d;
  - D' - E' = (D - E) + (U - d) >= d;
  - D' - E' + d' = D - E + U <= 2U - d, by the fourth inequality;
  - 2U - d >= 2d because U >= 2d.
- **J, giving (D-E+U, D-E, D+U, D):**
  - E' = D - E <= D = d';
  - D' - E' = U >= D = d';
  - D' - E' + d' = D + U = U';
  - D + U >= 2D.

The parities D, d odd and E = U (mod 2) are preserved similarly.

## Appendix B: the height bound and Theorem 5 with C_K = 0

This is Session200 root-weighted-height section 1. Put W = U - d. The three moves act on W by

    O: W' = 2W + d,   Z: W' = 2W,   J: W' = W + d,

and W(t) = 1. Hence W(t.w) >= 2^(n(w)), and

    F = D - E + U >= d + U = W + 2d >= 2^(n(w)) + 2.

So a chamber word A c A^dagger of value <= X has n(A) < (1/2) log2 X. The number of such A with exactly j = floor(K/2)
letters J is at most

    2 sqrt(X) C((1/2) log2 X + j, j).

## Files created (all in this folder)

| File | Contents |
| --- | --- |
| `chamber_audit.c` | adjoint, Q identities, (phys), growth, Theorem 1, parity law (two Fix methods), Theorem 3, singletons and windows |
| `census_audit.c` | segmented OpenMP letter-by-letter census: per-level coverage and saturating M_K counters |
| `depth_literal_audit.c` | depth with no cap on J; literal builder; Hamiltonian path and cycle counts |
| `carmichael_checks.py` | the elementary totient statements |
| `compare_with_packet.py` | the automated comparison; output in `COMPARISON.txt` |
| `k_singletons.py` | the word-singletons for K = 2..8 (M4); output in `k_singletons_2048.out` |
| raw outputs | `chamber_audit_2p20.txt`, `chamber_audit_K01_2p26.txt`, `chamber_audit_K2_2p24.txt`, `census_K2_1e8.txt` (+ `.err.txt`), `census_K2_2p26_to_1e9.txt` (+ `.err.txt`), `census_K{2..7}_1e6.txt`, `census_K2_4e6.txt`, `census_K3_4e6.txt`, `census_K2_2p24_exact32.txt`, `depth_literal_30001.txt`, `carmichael_checks.out` |

Reproduce, with binaries outside the repository:

```
gcc -O2 -o ca chamber_audit.c -lm && ./ca 1048576 2 16 > chamber_audit_2p20.txt
./ca 67108864 1 0 > chamber_audit_K01_2p26.txt
./ca 16777216 2 0 > chamber_audit_K2_2p24.txt
gcc -O2 -fopenmp -o cen census_audit.c
./cen 100000000 2 24 8 1 100001 > census_K2_1e8.txt
./cen 1000000000 2 26 4 0 100000001 512 67108864 > census_K2_2p26_to_1e9.txt
for k in 2 3 4 5 6 7; do ./cen 1000000 $k 20 8 1 1 > census_K${k}_1e6.txt; done
./cen 4000000 3 22 8 1 1 > census_K3_4e6.txt
./cen 4000000 2 22 8 1 1 > census_K2_4e6.txt
./cen 16777215 2 24 32 0 > census_K2_2p24_exact32.txt
gcc -O2 -o dl depth_literal_audit.c && ./dl 30001 9 1877:ZJJZJJJJJOO 1897:ZZJJJJJJJOO ... > depth_literal_30001.txt
python3 -I carmichael_checks.py > carmichael_checks.out
python3 -I k_singletons.py 2048 > k_singletons_2048.out
python3 -I compare_with_packet.py > COMPARISON.txt
```

The full argument list for the last `dl` run is the 18 witnesses of LITERAL_WITNESSES.txt, given as N:word.
