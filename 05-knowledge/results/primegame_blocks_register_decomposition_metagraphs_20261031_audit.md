# Independent adversarial audit of opus-2026-10-31-S71 (note, THM-4854, index rows)

Auditor: independent (Claude), 2026-10-31. Worktree audited: `/c/tmp/wt-s71` at commit `5277b4505b`
(parent `5431278adc` = session129 checkpoint 3). Nothing in the worktree was modified. Every number below was
reproduced with the auditor's own code (`s71_primegame_check.py`, `s71_metagraph_check.py`,
`s71_bijection_scan.cpp`, `s71_linkcheck.py` in this directory; outputs `primegame_check.out`,
`metagraph_check.out`, `bijection_scan.out`, `linkcheck.out`). The author's script
`primegame_blocks_metagraph_20261031.py` was read but not run.

Files audited: `05-knowledge/results/primegame_blocks_register_decomposition_metagraphs_20261031.md` (the note),
`01-canon/theorems/THM-4854-primegame-four-blocks-...-quotient-of-complement-graph.md` (the theorem),
`04-computation/experiments/primegame_blocks_metagraph_20261031.py` / `.out`, the S71 segment at the end of line 2 of
`05-knowledge/hypotheses/INDEX.md`, the top row of `05-knowledge/results/INDEX.md`, the S71 bullet (line 803) of
`00-navigation/PROBLEM-LEDGER.md`.

## Verdict: PASS WITH CORRECTIONS

Every finite computation in the package reproduces exactly (PRIMEGAME table, flag sets and intersections,
bipartition, the unique reciprocal pair, the product census 2 / 78 / 6, the 19 entries and their triangular and
Fibonacci members, the semiprime lemma to 10^5, the merged metagraphs at orders 3, 4, 5 line by line against
THM-4729, 18 adjacencies, 15 of 20 by 40 bijections, the random calibration). The derivation of the three-move
restatement from the composition law is correct. Ming 1989 is cited correctly. But the package has three false
sentences (one of them the headline "fibres" claim of THM-4854 (A)), a lemma mis-numbered in six places, and
substantial missing credit: the same-day session129 packet (committed 13:19-13:34, i.e. before S71's 13:41 commit
and in S71's own commit ancestry) already states the per-fraction control/data table, the bipartite structure of
the flag graph, the block products, the order-3 and order-4 verdicts and the triangular identities; Conway 1987
itself names the data primes "storage registers" and gives the macro pairs; THM-4815 (C) already drew the
THM-4810 / THM-4732 parity link S71 presents as its own correction; the semiprime lemma is an OEIS comment. The
"PROVED (inspection)" typing of the register decomposition should be CITED, and a MISTAKES.md entry is warranted
(the repeated same-day-concurrency pattern: MISTAKE-644 / 677 / 684 / 694 / 697 / 698).

## MAJOR findings

### M1. THM-4854 (A): the blocks are NOT the fibres of "which data registers change" (false as stated)

Quote (THM-4854 (A)): "the owner's blocks are the fibres of \"which data registers move\": none (B1), {5, 7} (B2),
{3, 7} (B3), {2, 3, 5} (B4)." Also note section 0: "(1) The owner's four blocks are the register decomposition of
Conway's program (section 1)"; hypotheses INDEX: "PROVED: the owner's four blocks are Conway's register
decomposition"; results INDEX: "The owner's blocks are Conway's registers"; ledger: "The owner's decomposition is
the register decomposition of Conway's program".

Evidence (`primegame_check.out`): the per-fraction map "set of data registers changed" has SIX fibres, not four:
`{}` = {1/17, 11/13, 13/11}; `{2}` = {23/38}; `{2,3,5}` = {78/85, 15/2}; `{3}` = {19/51, 29/33};
`{5}` = {95/23, 55/1}; `{7}` = {17/91, 77/29, 77/19, 1/7}. So 1/17 (no data change) sits in B2 with the data
movers, register 7 is split across B2 and B3, register 5 across B2 and B4. What the author's code printed (and what
is true) is the per-BLOCK UNION of registers touched: B1 none, B2 {5, 7}, B3 {3, 7}, B4 {2, 3, 5} - four distinct
sets, but overlapping (7, 5, 3 shared pairwise), so they do not decompose the registers either. "Register
decomposition" and "fibres" are both wrong words.

Correction: "each block touches a distinct set of data registers (none, {5, 7}, {3, 7}, {2, 3, 5}); the blocks are
not the fibres of a per-instruction register map (1/17 moves no register; registers 7 and 5 are each split between
two blocks)". The only intrinsic characterisation of the owner's split in the repository is session129's height
grading (partition/README.md: Phi = 2e_7 + 2U + 3V - 4T, every X1 instruction lowers Phi by 2, every X2 raises it
by 1, X0/X3 are horizontal, X0 = the horizontal moves changing no data counter) - cite it (M2).

### M2. Missing credit: session129 (same owner directive, same four blocks, committed before S71's note)

`05-knowledge/results/session129_primegame_tournament_partition_20261031/` (checkpoints cb0f03e551 13:19:31,
fbc383f8e9 13:24:53, 5431278adc 13:34:24; S71's note commit 5277b4505b 13:41:19 has 5431278adc as parent). The
note says "Three lanes answered the first PRIMEGAME directive on 2026-10-31 and are cited for everything they
contain" and never mentions session129 (grep count 0 in note and theorem). Session129 already contains:

- the per-fraction control-arc / data-gain table (bridge_audit/README.md section 2, letters A..N: "A | 17/91 |
  13 -> 17 | 1/7", ..., "N | 55 | 0 -> 11 | 5") = the note's section-1 table;
- the bipartite structure of the flag graph (same file): "Its underlying simple graph, after deleting loops, is
  `K_(2,3)` on bipartition {11,17} | {0,13,19}, with a pendant `29` at `11` and a pendant `23` at `19`" - exactly
  the note's parts {11, 17, 23} / {13, 19, 29, none};
- the block products (same file, section 1): "a = (11/13)(13/11) = 1, b = 55(1/7)(1/17)(17/91) = (5/49)(11/13), c =
  ... = (49/9)(11/17), d = ... = (45/2)(13/17)" and the transitive triangle 17 -> 13 -> 11, 17 -> 11 they define;
- B1 as the no-data moves (partition/README.md: "J,K are exactly the horizontal instructions leaving all four data
  counters unchanged") and the contraction of the reciprocal pair "but only with a guard" (Theorem 4);
- the order-3 verdict (tournaments/REPORT.md, item 1 and section 2): "Order 3 has **one complement line**, not one
  node. This repairs the user's first proposed identification without discarding the reciprocal pair" and "The
  fraction pair `{11/13,13/11}` is likewise a two-element set with a fixed-point-free reciprocal involution ...
  identifies its reciprocal pair with the single element of `L_3`. It does **not** identify it with a single merged
  node" = the note's section-4 DICTIONARY reading of B1;
- the order-4 verdict (same report, section 3): "four literal complement lines: one blue `A--C`, one blue loop at
  C, and two black `B--C` lines ... the three node fibres have sizes **1,2,5**, not 4,4,4. Thus the three
  four-instruction blocks can be assigned as *names of nodes*, but cannot be bijected instruction-by-instruction";
- the triangular identities (same report, section 6): "`55=5*11=T_10=F_10`, `91=7*13=T_13`, `78=2*3*13=T_12`,
  and `15=3*5=T_5`", with the earlier repository file `05-knowledge/results/triangular_fibonacci_heron_boat_kps.out`
  (line 13: "5·11=55=T_10 Fib=True") credited there;
- the encoding-dependence typing (encoding/README.md, summarised in README section 5): "prime names can be changed
  ... but 55 then need not remain triangular or Fibonacci. Numerical identities are retained as facts about
  Conway's chosen encoding".

Correction: add session129 to section 0 and to the CITED line of the status; retype the section-1 table, the
bipartiteness, the B1 verdict, the order-4 verdict and the triangular entries as CITED (session129, with the
S71 additions stated as such: the fibre check of M1 done correctly, the pairwise flag intersections, the product
census to length four, the p(2p +- 1) reading, the calibrated order-5 test, the THM-4810 restatement). Record a
MISTAKES.md entry.

### M3. The register decomposition is Conway's own and session126's / mac-mini's / THM-4852's: CITED, not PROVED

Quote (note, status line): "PROVED (elementary: the register decomposition and the hub flag; ...)"; section 1
title "The register decomposition (PROVED by inspection ...)"; typing table "register decomposition, hub flag 17,
bipartite flag machine, unique reciprocal pair | PROVED (inspection)".

Evidence. Conway 1987 (gwern scan, OCR text extracted by the auditor, `conway_dl/conway1987.txt`), section 7:
"The different primes that arise in the numerators and denominators of the various fractions may be regarded as
storage registers, and in a state in which the current working integer is N = 2^a 3^b 5^c 7^d ..., we say that
register 2 holds a, ... register 7 holds d". Section 5: "To help in Figure 2, we have labeled the fractions:
A B C D E F G H I J K L M N ... and we note that AB = 2x3/(5x7), EF = 7/3, DG = 5/2" - Conway's own macro pairs,
with Figure 2 run on 5^n 7^d 13 (the control primes 11, 13, 17, 19 as the phase markers). Session126's FRACTRAN
report, section 3 "Four counters and seven control states": "(a,b,c,d)=(nu_2,nu_3,nu_5,nu_7) are unbounded data
counters and the remaining primes select one of seven control states {1,11,13,17,19,23,29}" with the full
per-state transition table (three rows for control 17). Mac-mini's note section 4: "data registers 2, 3, 5, 7;
control states 11, 13, 17, 19, 23, 29 and 'none'; state digraph arcs (1,1), (1,11), (11,13), ..., (29,11)".
THM-4852's setting: "data registers 2, 3, 5, 7 and state registers 11, 13, 17, 19, 23, 29".

Correction: type the data/flag split and the per-fraction table CITED (Conway 1987 sections 5 and 7; Guy 1983;
session126 section 3; mac-mini section 4; THM-4852; session129 bridge audit section 2); keep PROVED (inspection)
only for what is S71's: the pairwise intersections {11,17}, {13,17}, {17,19}, the unique reciprocal pair, the
length-<= 4 product census. The "hub 17" is implicit in session126's table (three rows at 17) and in session129's
K_(2,3) (17 is one of its two hubs); say so.

### M4. THM-4732 "Lemma 2" is Lemma 1; the parity link to THM-4810 is already in THM-4815 (C); the two-move
criterion is a different statement

Quote (note section 3): "THM-4732 Lemma 2 (\"an even number L of black moves needs n <= 2L + 2 and C(n-1, 2)
even; every black distance is odd for n = 0, 3 mod 4\") is the detour form of THM-4810's parity invariant (one move
toggles C(n-1, 2) pairs), and THM-4810's two-move criterion is its L = 2 case ... (a self-citation S63 missed)".
Repeated in the note's section 0 item (4) and section 5, THM-4854's title and (D), the hypotheses INDEX ("THM-4732
Lemma 2 credited as the precursor"), the results INDEX and the ledger.

Evidence. (i) THM-4732 section 2: "**Lemma 1.** K_Q flips the pairs outside E(Q) ... Even L needs C(n-1, 2) <= L
(n - 1), i.e. n <= 2L + 2, and C(n-1, 2) even ..., i.e. n = 1, 2 mod 4; so black distances are odd for n = 0, 3
mod 4 and at least 3 for n >= 7. **Lemma 2.** All non-pure-blue classes lie in one black component iff (a) ... and
(b) ...". The quoted content is Lemma 1; Lemma 2 is the black-core connectivity criterion. (ii) THM-4815 (C)
(mac-mini-2026-10-27-dipole, four days before S71): "For n ≡ 0, 3 (mod 4) it alternates (THM-4810's toggle count is
odd there; THM-4732 Lemma 1 uses this for detours), so the labelled move graph is bipartite" and, in its
concurrency paragraph, "the parity law (C; implicit in THM-4810 Theorem P's odd toggle count and used in THM-4732
Lemma 1)". The cross-reference S71 presents as a missed self-citation is already canon, with the right lemma
number; moreover THM-4810's Theorem P is the EVEN case (n = 1, 2 mod 4: invariance; "For n = 0, 3 (mod 4) ... no
such invariant exists"), while THM-4732 Lemma 1 is the ODD case (distance parity forced), which THM-4810's own
addendum credits to THM-4815 (C) (MISTAKE-682). (iii) THM-4810 (C): "**Theorem C2 (two moves).** d(T, transitive
class) <= 2 iff some total order r satisfies D_r(T) = E(P) xor E(Q)" is about the transitive class; the L = 2 case
of THM-4732 Lemma 1 is about the line partner (E_1 xor E_2 = A \ E(P), impossible for n >= 7 by size). Both are
corollaries of the composition law (THM-4810 Lemma 1 = the first sentence of THM-4732 Lemma 1) but they are not
the same statement.

Correction: "Lemma 1" everywhere (six places); replace "a self-citation S63 missed" by a citation of THM-4815 (C);
replace "THM-4810's two-move criterion is its L = 2 case" by "the L = 2 case of THM-4732 Lemma 1 (no two-move
detour for n >= 7) and THM-4810's Theorem C2 are both instances of the composition law".

### M5. False sentence: "L is constant along a blue line"

Quote (note section 3): "Score monotonicity (THM-4810, Theorem L) adds only that L is constant along a blue line
and changes by at most one per black move, which does not constrain a three-move detour."

Evidence (`metagraph_check.out`, L = (1/2) sum_i |s_(i) - i| as in THM-4810; L(TT) = 0, L(regular 5) = 3 checked):
order 5: all 8 blue lines have |Delta L| = 1 (e.g. T5 (L = 3) to Q5 (L = 2)); order 7: all 256 blue lines have
|Delta L| = 1; order 6: 20 of 32 blue lines have |Delta L| = 1, 12 have 0; order 4: T4 - S4 has |Delta L| = 1.
Black lines have |Delta L| in {0, 1} (order 5: 14 and 10; order 7: 7812 and 8316). Theorem L says |Delta L| <= 1
along every move; nothing makes a blue line score-preserving (kappa t is tau_P(t), not the converse of t).

Correction: delete "L is constant along a blue line"; "L changes by at most one along every line, blue or black".

### M6. False sentence: the Fibonacci tiling counts sit "on the five nodes with a nontrivial automorphism or a
merged pair"

Quote (note section 4): "The ten merged nodes' tiling counts 1, 1, 2, 3, 5, 9, 9, 10, 11, 13 contain the Fibonacci
numbers 1, 1, 2, 3, 5 on the five nodes with a nontrivial automorphism or a merged pair, and 13 on the node with
H = 13".

Evidence (`metagraph_check.out`, THM-4729 numbering): tiling counts 1, 1, 2, 3, 5 belong to nodes 0 (TT5: |Aut| =
1, self-converse - neither condition), 3 (1 -> C3 -> 1, |Aut| = 3), 1 (the (1,1,C3) merged pair), 9 (T5, |Aut| =
5), 4 (Q5, |Aut| = 3); the nodes with |Aut| > 1 or a merged pair are 1, 2, 3, 4, 9 with counts 2, 10, 1, 5, 3 - the
merged (1,S4) pair (node 2) has 10 tilings. (Typed NUMEROLOGY, so no consequence, but the sentence is false.)
Also: the count list "1, 1, 2, 3, 5, 9, 9, 10, 11, 13" is already S70 section 4's count coincidence; cite it.

### M7. The semiprime lemma is a recorded OEIS comment (and session129 lists the identities): CITED, not a new lemma

Quote (note section 2 / status): "**Lemma (semiprimes).** For primes p < q, pq is triangular iff q = 2p + 1 or
q = 2p - 1" typed "PROVED (elementary ...)". Evidence: OEIS A068443 "Triangular numbers which are the product of
two primes" (6, 10, 15, 21, 55, 91, 253, 703, 1081, 1711, 1891, 2701, ...) carries the comment "These triangular
numbers are equal to p * (2p +- 1)"; A156592 "Product p*q of two primes with q = 2*p + 1" is a subsequence
(Zumkeller 2009). The lemma is correct (my check: 0 violations below 10^5, 27 triangular semiprimes, the first
being 6, 10, 15, 21, 55, 91, 253, 703) but should be typed "elementary, OEIS A068443 (CITED)". The identities
15 = T_5, 55 = T_10, 78 = T_12, 91 = T_13 are session129's (M2).

## MINOR findings

- m1. "17, the only flag with three outgoing instructions" (note section 1; THM-4854: "the only flag with three
  instructions"): the idle state "none", which the note itself uses as a vertex of the flag digraph, also has three
  instructions (15/2, 1/7, 55/1; my count: {13: 2, 17: 3, 19: 2, 11: 2, 29: 1, 23: 1, none: 3}). Say "the only
  control prime with three instructions".
- m2. "Hence every closed walk of flag moves has even length" (section 1) and "PRIMEGAME's flag machine is
  bipartite and has no odd closed walk at all" (section 4); THM-4854 (A) "so every closed walk of flag moves is
  even": false with the two state-free loops included (15/2 and 1/7 are closed walks of length 1 at "none", and
  Conway's phase word L^n M^(d-1) traverses them). Qualify: "every closed walk avoiding the two state-free loops".
- m3. "re-marking along another Hamiltonian path is exactly a black K move of THM-4732" (section 3): a K move
  K_Q(T) = T^op xor E(Q) is kappa taken AFTER re-marking along Q (K_Q(T) is isomorphic to kappa of the tiling
  re-marked along Q); re-marking alone changes nothing in the quotient. And it is "black" only when Q is a black
  word (endpoint_three README: "a word is blue when its positional reversal is an anti-automorphism, and black
  otherwise").
- m4. The three-move restatement (section 3 and THM-4854 (D): "the base path as the symmetric difference of three
  successive legal paths") drops the blackness of Q_1, Q_2, Q_3. The derivation is right: tau_Q3 tau_Q2 tau_Q1 T =
  T xor A xor E_1 xor E_2 xor E_3 and tau_P T = T xor A xor E(P), so equality iff E_1 xor E_2 xor E_3 = E(P). But
  HYP-9347 (1) and THM-4732 (3) quantify over BLACK words ("three successive black K moves"), and the black distance
  is defined over black words only; without "black" the restatement defines a different (not larger) distance. Add
  "black" (positional reversal not an anti-automorphism of the current tournament).
- m5. "(D) PROVED identification" / "the merged tiling line metagraph ... is the quotient ... (PROVED
  identification)": this is definitional. THM-4809 (S62, same author) already defines X_n with "tau complements
  every free pair" and studies the reset cells ([x], [tau x]) (the unmerged metagraph lines; E = 7, 50, 732,
  24654), and THM-4729 Proposition K already calls the merged metagraph "the convention-free quotient". Retype as
  a remark (CITED THM-4809 / THM-4729 K) rather than a new PROVED item in the theorem title.
- m6. Broken link (note section 0): `../hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-chirality-character-on-black-twin-pairs.md`;
  the file is `HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md`.
  All other 10 links in the note and the single link in the theorem resolve (`linkcheck.out`).
- m7. Semiprime proof (section 2): "2pq = k(k + 1) with k, k + 1 coprime forces {k, k + 1} = {p, 2q} or {2p, q}"
  omits the trivial splits k in {1, 2, pq, 2pq} (k(k+1) = 2, 6, pq(pq+1), 2pq(2pq+1), none equal to 2pq); add
  half a sentence. The case p = 2 is covered by {2p, q} = {4, q} (6 = 2 x 3, 10 = 2 x 5).
- m8. "77, 91 are the two composites with the factor 7 and both are strong-only Hamiltonian-path values at order 7
  (THM-4810 (A))" (section 2): true (order-7 sweep: 77 has 10 classes, all strong; 91 has 15, all strong; 91 is the
  most populous value, as the summary output says), but every value >= 49 at order 7 is strong-only, so the remark
  carries no information; say so, or drop it.
- m9. "the program has no macro law beyond the control toggle" (section 1; THM-4854 (B)): true under the note's
  narrow definition (product equal to 1 or to a program fraction), but the natural macro laws exist and are Conway's
  own (section 5: AB = 6/35, EF = 7/3, DG = 5/2 - test-and-act pairs that restore the flag) and session126's word
  W_d = (AB)^d J (EF)^d K; the sentence should state the definition and cite Conway's pairs.
- m10. Index wording: results INDEX "The owner's blocks are Conway's registers" and hypotheses INDEX "the owner's
  four blocks are Conway's register decomposition" (and the ledger) inherit M1 / M3; fix with them, and fix
  "Lemma 2" (M4) in the hypotheses INDEX, the results INDEX and the ledger.
- m11. HYP-9449: released by deleting the RESERVED stub; the only reference left is the hypotheses INDEX line 2
  ("no hypothesis, HYP-9449 released"). Earlier releases (HYP-9402 by S57, HYP-9412 by S58) kept a stub with
  `status: RELEASED` and "The ID may be reused". No index change is required, but for consistency restore a
  RELEASED stub (otherwise the number is a silent gap in the sequence).
- m12. "The line-degree of a node equals its tiling count" (section 3): correct only with loops counted twice (two
  line-ends); verified at orders 3, 4, 5 with that convention (e.g. node H = 13: 1 + 2 + 1 + 2 + 2x2 + 2 + 1 = 13).
  Say "loops counted twice".
- m13. Calibration (section 4, THM-4854 (E)): 20 random graphs is a small sample; with 600 (three seeds x 200) the
  maxima are 14-17 with means 15.155, 15.105, 15.170 and 15 is the mode (389 of 600). The conclusion (inside the
  random range) stands; quote the larger sample or say "20 trials".
- m14. Section 4, "B2, B3, B4 against the three nodes at n = 4": the REFUTED verdict is correct (K_3 on three
  blocks against two non-loop adjacencies plus a loop), but session129's tournaments report already gave the
  equivalent verdict (fibres 1, 2, 5 against 4, 4, 4) - cite (M2).

## Numbers reproduced with the auditor's own code (exact outputs)

PRIMEGAME (`primegame_check.out`):
- per-fraction table: identical to the note's (flag moves 13->17, 17->13, 17->19, 19->23, 11->29, 29->11, 23->19,
  19->11, 17->none, 13->11, 11->13, none->none, none->none, none->11; data changes as tabulated); one-hot flag
  structure holds for every fraction; all 28 entries squarefree.
- per-fraction register fibres: 6 (M1); per-block unions: B1 [], B2 [5, 7], B3 [3, 7], B4 [2, 3, 5]; overlaps
  B2&B3 [7], B2&B4 [5], B3&B4 [3].
- flag sets: B1 [11, 13], B2 [11, 13, 17], B3 [11, 17, 19, 29], B4 [13, 17, 19, 23]; intersections B2&B3 [11, 17],
  B2&B4 [13, 17], B3&B4 [17, 19]; instructions per source: 13: 2, 17: 3, 19: 2, 11: 2, 29: 1, 23: 1, none: 3.
- bipartite without the two loops: True; parts [[11, 17, 23], [13, 19, 29, 'none']]; 12 non-loop arcs, all crossing.
- reciprocal pairs: [('11/13', '13/11')] only.
- ordered products equal to 1 or a fraction: length 2: 2, length 3: 78, length 4: 6, all of the toggle-inserted
  form; B2 x B3 x B4: none.
- distinct entries: [1, 2, 7, 11, 13, 15, 17, 19, 23, 29, 33, 38, 51, 55, 77, 78, 85, 91, 95] (19); triangular
  {1: T_1, 15: T_5, 55: T_10, 78: T_12, 91: T_13}; Fibonacci [1, 2, 13, 55]; composites [15, 33, 38, 51, 55, 77,
  78, 85, 91, 95] (10).
- semiprime lemma below 10^5: 0 violations; 27 triangular semiprimes, first (2,3,6), (2,5,10), (3,5,15), (3,7,21),
  (5,11,55), (7,13,91), (11,23,253), (19,37,703); for 2pq = k(k+1) the admissible k is 2p or q as the lemma says.
- composite-sharing graph: 20 edges, degree sequence [2, 2, 3, 4, 4, 4, 5, 5, 5, 6].

Metagraphs (`metagraph_check.out`; own conventions, canonical forms by all permutations, H and Aut by brute
force; sanity checks passed: sigma an involution commuting with kappa, kappa fixed-point-free, Camion closure
law c > 0 iff strong, THM-4729 L2 closure numbers, THM-643 fibre law tc |Aut| = H):
- order 3: 2 tilings, 2 classes, 2 nodes (H = 1 and H = 3, |Aut| 1 and 3, one tiling each), 1 line, blue; both
  nodes pure-blue.
- order 4: 8 tilings, 4 classes, 3 nodes: (H, |Aut|, t) = (1, 1, 1), (3, 3, 2) merged, (5, 1, 5) with 3 symmetric
  and 4 closed tilings; lines T4 - S4 blue, pair - S4 black x2, S4 - S4 blue loop; pure-blue: T4 only.
- order 5: 64 tilings, 12 classes, 10 nodes, 32 lines (8 blue, 24 black); node data (H, |Aut|, t, s, c) =
  (1,1,1,1,0), (3,3,2,0,0) merged, (5,1,10,0,0) merged, (3,3,1,1,0), (9,1,9,1,5), (9,1,9,3,5), (11,1,11,3,5),
  (15,3,5,1,5), (13,1,13,3,10), (15,5,3,3,2); pure-blue TT5, 1->C3->1, T5; under the renumbering to THM-4729's
  table (forced by the invariants) the 8 blue lines equal THM-4729's list exactly and the 12 doubled black lines
  equal THM-4729's list exactly (also S70 section 4's lists); 18 distinct non-loop adjacencies, loops at THM nodes
  5 (H = 11) and 8 (H = 13); two black lines between the H = 11 and H = 13 nodes; line-degree (loops twice) equals
  the tiling count at every node of orders 3-5.
- |Delta L| along lines: order 5 blue (1): 8; black (0): 14, (1): 10. Order 6 blue (0): 12, (1): 20; black (0):
  264, (1): 216. Order 7 blue (1): 256; black (0): 7812, (1): 8316. (M5.)

Bijection scan (`bijection_scan.out`, C++, all 10! bijections): composites vs the THM-4729 adjacency set (my
recomputed set equals it): max preserved 15 of 20, attained by 40 bijections; degree sequence of the metagraph
support [1, 1, 1, 3, 3, 5, 5, 5, 6, 6]; random 20-edge graphs, 200 per seed: seed 12345 histogram 14:20 15:133
16:43 17:4 mean 15.155; seed 777: 14:26 15:130 16:41 17:3 mean 15.105; seed 20261101: 14:21 15:126 16:51 17:2 mean
15.170 (the author's 20 trials: 14-17, mean 15.25: consistent).

Order-4 verdict: K_3 (three shared-flag edges) against a star with two non-loop adjacencies plus a loop: no
bijection; REFUTED is right (and is session129's fibres 1, 2, 5 verdict).

Canon data checked against the repository's own outputs (not recomputed): `h_spectrum_sweep_20261027_n7.out`:
H = 77: 10 classes, 10 strong; H = 91: 15 classes, 15 strong; `h_spectrum_summary_20261027.out`: "order 7: most
populous values (classes): [(91, 15), ...]"; f(6) = 15 (THM-4810 (A)). The endpoint_three README: "PROVED
structural theorem for every n>=11 ... d_T(a)=n-4, d_T(z)=3" - the note's "(n-4, 3) layer theorem proved for
n >= 11" is quoted correctly; HYP-9347 (1) "every closed blue line between two non-PB classes has black distance
exactly 3 for n >= 8" quoted correctly; THM-4732 (3)'s exhaustive orders 6, 8, 9, 10 quoted correctly.

## Prior art and credit (with sources)

- J. H. Conway, FRACTRAN: a simple universal programming language for arithmetic, in Open Problems in
  Communication and Computation (1987), pp. 4-26, DOI 10.1007/978-1-4612-4808-8_2; scan
  https://gwern.net/doc/cs/computable/1987-conway.pdf (OCR text extracted: `conway_dl/conway1987.txt`). Theorem 1
  = the owner's 14 fractions (verified). Section 5 labels the fractions A..N and gives the macro pairs AB, EF, DG
  and Figure 2 on 5^n 7^d 13; section 7 calls the primes "storage registers" (register 2 holds a, ...). Source
  of the data/flag reading: CITED.
- Luo Ming, On triangular Fibonacci numbers, Fibonacci Quart. 27 (2) (1989) 98-108,
  https://www.fq.math.ca/Scanned/27-2/ming.pdf (index https://www.fq.math.ca/27-2.html). Verified: "Theorem 2: The
  Fibonacci number u_n is triangular if and only if n = +-1, 2, 4, 8, 10" (values 1, 3, 21, 55), "verifies the
  conjecture of Hoggatt". The note's citation is correct.
- OEIS A068443, Triangular numbers which are the product of two primes, https://oeis.org/A068443 (comment "These
  triangular numbers are equal to p * (2p +- 1)"); A156592 (p q with q = 2p + 1). Prior art for the semiprime lemma.
- Repository, same day: session129 packet (paths and commits in M2); session126
  `session126_primegame_symmetry_20261031/primegame/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md` section 3 (four
  counters, seven control states, full transition table; the first README line "four data counters and six control
  primes, with first-enabled priority"); mac-mini `groups_of_five_primegame_merged_metagraph_20261031.md` section 4
  (data registers, control states, the 13 state-digraph arcs) and section 5 (the order-5 table and lines in the
  owner's gauge); S70 `primegame_variant_selector_automaton_metagraph_20261031.md` section 4 (the order-5 table,
  lines, tiling counts 1, 1, 2, 3, 5, 9, 9, 10, 11, 13 as a count coincidence) and THM-4852 (setting line).
- Repository, earlier: THM-4815 (C) (THM-4810 Theorem P / THM-4732 Lemma 1 link, 2026-10-27); THM-4809 (X_n, tau,
  reset cells ([x], [tau x])); THM-4729 Proposition K (convention-free quotient);
  `05-knowledge/results/triangular_fibonacci_heron_boat_kps.out` (55 = T_10, Fibonacci).
- Not prior art for the register split, but relevant: R. K. Guy, Conway's prime producing machine, Math. Mag. 56
  (1983) 26-33 (the 15/14 list with its flowchart; cited by THM-4852) - not re-read here.

## Answers to the seven audit questions

1. PRIMEGAME facts: all reproduced (table, blocks' register unions, flag sets, intersections, bipartition, 12
   crossing arcs, unique reciprocal pair, 2 / 78 / 6 all toggle-inserted, none from B2 x B3 x B4). Corrections:
   the blocks are not register fibres (M1); "17 the only flag with three instructions" needs "control prime" (m1);
   "every closed walk is even" needs the loops excluded (m2); the register split is Conway's (M3).
2. Triangular / Fibonacci: 19 entries, T_1, T_5, T_10, T_12, T_13, Fibonacci 1, 2, 13, 55 reproduced; lemma proof
   correct modulo the trivial splits (m7), 0 violations below 10^5; Ming 1989 exists and says exactly what is
   claimed. The lemma is an OEIS comment (M7); the identities are session129's (M2).
3. Metagraphs: orders 3, 4, 5 reproduced and identical to THM-4729 / S70 line by line; line-degree = tiling count
   with loops counted twice (m12). The blue-line L claim is false (M5); the Fibonacci-node sentence is false (M6).
4. Identification: correct but definitional (m5); THM-4732's parity lemma is Lemma 1, not Lemma 2, and the link to
   THM-4810 Theorem P is THM-4815 (C)'s (M4); "two-move criterion = L = 2 case" conflates two statements (M4);
   HYP-9347 (1) and the (n-4, 3) n >= 11 theorem are represented correctly; the three-move restatement is derived
   correctly but omits "black" (m4); "re-marking is a black K move" is imprecise (m3).
5. Calibrated test: 20 edges, 18 adjacencies, 15 of 20 by 40 bijections, random maxima 14-17 (mean about 15.1-15.2
   over 600 trials) all reproduced; order-4 verdict correct (and session129's).
6. Prior art: session126 (control table), mac-mini (register/state split and arcs), THM-4852 (setting) state the
   register decomposition; session129 states the bipartite K_(2,3) structure, the block products, the order-3/4
   verdicts and the triangular identities; Conway's own text is the source of the register reading (CITED, not
   PROVED). The hub 17 and the reciprocal-pair uniqueness are not stated as such elsewhere (session129 has the
   product a = 1 and the contraction theorem; session126 has the 11 <-> 13 balanced two-cycles). Typings:
   DICTIONARY / NUMEROLOGY / REFUTED are honest; "PROVED (inspection)" for the register decomposition and "PROVED
   identification" for (D) overclaim novelty; "no macro law" is true only under the note's narrow definition (m9).
   HYP-9449: no other index references; a RELEASED stub would match the S57 / S58 practice (m11).
7. Links: one broken (HYP-9347 filename, m6); the other eleven resolve.

## What could not be checked

- Conway's Figure 2 beyond its labels (the OCR of the figure is garbled); Guy 1983's wording was not re-read.
- THM-4732's and the endpoint lanes' censuses and proofs (cited by the note, not re-run here).
- The author's script was not executed (by instruction); its `.out` agrees with the auditor's numbers everywhere
  the two overlap, including the 40 maximising bijections and the 20-trial calibration list.
- Whether the owner's four blocks have a further intrinsic characterisation beyond session129's height grading.
