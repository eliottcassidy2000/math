# Audit D2: THM-5056 (odd-graph Chan-Pak), mac-mini-2026-11-09-digraph721

VERDICT: PASS WITH CORRECTIONS

Independent adversarial audit, 2026-11-09, of the worktree version at commit c11abf3a22 (13:18:18). I read:
- THM-5056;
- PROOF.md Part II, together with the parts of Part I that it uses (Lemma 4, Corollary 5);
- HYP-9625;
- README.md (item 6, verification, credits) and the two INDEX rows;
- the scripts `oddbase.py`, `oddproof.py`, `oddbuild.py`, `oddspecial.py` and `wcheck.c`.

Every claim was re-checked with programs written for this audit, in
`04-computation/experiments/digraph721_20261109/audit_d2/`. The session's code was run only to regenerate and
re-check its own outputs. I did not edit session files and did not commit.

Overall:
- **The theorem is true, and its proof is correct and complete.** I found no mathematical error, no gap in the case
  analysis and no wrong number.
- **Every step checks, by hand and by independent computation:**
  - Lemma O1 (a)-(e) and Lemma O2;
  - the four seeds (every admissible w), the Gamma_odd sets and the special case 513;
  - the three-case induction and its orders;
  - the base below 512.
- **The exclusion of 1, 5 and 37 was re-derived without the large census.** An exhaustive enumeration of all
  labelled graphs on <= 7 vertices shows two things:
  - every connected graph with d <= 37 has <= 6 vertices;
  - the four connected realisers of 37, and P_3, are all even.
- **The final statement holds for every n <= 7**, checked directly over all labelled graphs, connected or not, with
  no census and no construction.
- **13,149 explicit odd realisers pass a from-scratch brute-force checker.** These are the session's 4,450,
  regenerated bit-for-bit, and 8,699 built by an independent re-implementation that makes random admissible choices.
  My 8,699 also pass the session's `wcheck`.
- **Corrections: 11 MINOR items, no MAJOR.** Two are credit items that should be applied before close-out:
  - THM-5042 (opus-2026-11-09-S89) is an independent second proof. It reached origin/main two minutes after
    THM-5056 and cites it, but THM-5056 does not cite it back (Finding 1).
  - Completeness and conciseness of d over all graphs, used in step 1, is Brouwer-Csorba-Schrijver 2009, Prop. 1.2
    (Finding 2).

## Checks against the brief

| item | result | how |
|---|---|---|
| 1. Lemma O1 | correct | proofs re-derived; brute force (`d2lemmas.py`, 0 failures) |
| 2. Lemma O2 | correct | values from Lemma 4.1; bookkeeping re-derived; 5,331 Lemma 4.1 tests; odd-good invariants asserted after every move of every offset chain in the 8,699 random constructions |
| 3. seeds, Gamma_odd, 513 | correct; the only (x, L) gap with L >= 10 is (1, 10) | `d2seeds.py`, `d2induct.py` |
| 4. induction | correct; orders exactly L(k) | hand check; symbolic sweep of all odd k < 2^24 |
| 5. base below 512 | reproduced: 253 values; none for 1, 5, 37 | `d2base.py` (census only as a candidate list), `d2small.c`, `d2allsmall.c` |
| 6. brute force | 13,149 graphs, 0 failures | `d2check.c` (own d algorithm, explicit odd-automorphism certificates) |
| 7. credit and typing | new relative to THM-4996 and THM-5040; no overclaim; credit gaps | Findings 1, 2, 10, 11 |

Details:

1. **Lemma O1.**
   - (a) Extend an odd automorphism of G by the identity on H.
   - (b) If H has no universal vertex, the apex is the unique universal vertex of K_1 v H. Label it first: then
     Aut(K_1 v H) = Aut(H), with the same signs (B2.3). If H has a universal vertex, it and the apex are adjacent twins.
     Tested exhaustively on every odd labelled H with <= 5 vertices, and on random H with 6-7 vertices: 879 odd H, 313
     of them with a universal vertex.
   - (c) The two K_2 vertices have the same closed neighbourhood. Tested on 150 random G.
   - (d) If u is not in {x, y}, then N_PH[x] = N_G[x] u {h} = N_PH[y]. Tested on 473 random cases.
   - (e) Padding is immediate.
   - The sign of an automorphism does not depend on the reference orientation: flipping a set F of edges changes the
     reversal count by |F Δ g^-1(F)|, which is even. So relabelling in unions and joins is harmless.
2. **Lemma O2.**
   - **Values.** Lemma 4.1 with (r, l) = (2, 1) gives +3 * 2^(n-1), and with (3, 2) gives +7 * 2^(n-1). Since
     2^(n-1) + c becomes 2^(n+1) + c or 2^(n+2) + c, the offset is invariant.
   - **Lemma 4.1 holds for every graph, not only even ones.** Tested 5,331 times on random graphs of both parities,
     for every u and l = 1, 2, 3.
   - **Bookkeeping after R(G).** The pendants are p1 (its neighbour misses the new leaf) and the new leaf (its
     neighbour h misses p1). The degree-2 vertex is p2.
   - **Bookkeeping after PH(G, w, 2).** The pendants are the two leaves (h misses w), and w is unchanged.
   - All of these lie outside {x, y}, so an admissible excluded vertex always exists outside the twin pair.
   - **Any admissible choice works.** `d2build.py` picks at random the excluded pendant, the degree-2 vertex, the seed's
     w, the seed (when two qualify) and the order of the +2/+3 moves. It asserts odd-goodness after every move.
3. **Seeds.** The values, offsets, twins and w are as claimed: 13 = 2^3 + 5, 33 = 2^5 + 1, 101 = 2^6 + 37 and
   165 = 2^7 + 37.
   - **Choices of w.** K_4 - e has 2 admissible w, K_2 u C_4 has 4, K_2 v (P_3 u E_2) has 2, and K_2 u T_55 has 1. In
     K_2 u T_55, w is the bit-3 vertex of T_55, whose degrees are 4, 4, 4, 2, 5, 5.
   - **First PH(w, 2).** Every choice gives 69, 257, 549 and 1061, at orders 7, 9, 10 and 11. Full Aut enumeration
     gives ODD in each case.
   - **Closure under +2 and +3.** It gives exactly {7} u [9, oo) for offset 5, {9} u [11, oo) for offset 1, and
     [10, oo) for offset 37. So the only gap with L >= 10 is (x, L) = (1, 10).
   - **The case 513.** 513 = 3 * 171 = d(K_2 u T_171), on 2 + 8 = 10 vertices. It is odd through its K_2 component.
4. **Induction.**
   - **Case 1.** t = k - 3 * 2^(L-2) is odd, with 1 <= t < 2^(L-2).
   - **Case 2.** L(x) <= L - 2 < L - 1, so padding always adds an isolated vertex, and B2.3 applies directly (see
     Observation 2). The cone has value 2^(L-1) + x.
   - **Orders.** Each case gives order L exactly.
   - **Symbolic sweep** (`d2induct.py 24`). Every odd k < 2^24 other than 1, 5 and 37 falls into one case, and that
     case's hypotheses hold:
     - base: 253;
     - case 1: 4,194,176;
     - case 2: 4,194,131;
     - case 3: 44;
     - the 513 branch: 1.
5. **Base.**
   - **Method.** For each odd k < 512 I searched for C u T on exactly L(k) vertices. C is a connected census graph
     with d(C) | k. I recomputed d(C) by inclusion-exclusion and certified C odd by an explicit automorphism. T is a
     threshold graph from my own top-down construction.
   - **Result.** All 253 values other than 1, 5 and 37 are realised, and those three are not. The odd factor is
     certified by a transposition in 222 cases and by another automorphism in 31.
   - **The unions are re-checked.** `d2check` re-verifies d, the order and oddness of all 253 unions.
   - **Exclusions.** Every connected census graph whose d divides 1, 5 or 37 is even (full Aut enumeration). This is
     independent of the census by (C) in the reproduction section.
6. **Brute force.**
   - **Regeneration.** The session's `oddw_small.txt`, `oddw_large.txt` and `oddw_special.txt` were regenerated with
     the session's commands, and their SHA-256 hashes equal `sha256.txt`.
   - **Checks.** All 4,450 graphs pass both the session's `wcheck` and my `d2check`.
   - **How `d2check` works.** It counts d by a Gray-code walk with closed-neighbourhood counters. This is a different
     algorithm from the prefix-OR DP used by `wcheck` and `oddbuild.py`. It certifies oddness by an explicit
     permutation: the permutation is verified to be an automorphism, and its reversed edges are counted directly.
   - **Controls.** K_2, P_4 and K_4 - e pass; P_4 needs a non-transposition certificate. Five controls fail as they
     should: P_3 and C_4 (even), E_2 (even, and its order is wrong), a wrong value and a wrong order.

## Findings

### 1. MINOR (credit; apply before close-out): cite THM-5042, the independent second proof

**What is now on origin/main.**
- `01-canon/theorems/THM-5042-odd-graph-companion.md` by opus-2026-11-09-S89. It was added in commit 7e86bdf8f9
  (13:20:22) and corrected in a41a48937f (13:30:47).
- Its title is "Second, independent proof of THM-5056 (mac-mini, odd-graph Chan-Pak; first in canon)".
- It proves the same statement by a different route:
  - cones over cones for the leading 11;
  - padded odd cones for the leading 10;
  - the leaf stable recursion SR of THM-5040 Lemma 7, starting from six census seeds at orders 8 and 9, for
    2^(m-1) + {1, 5, 37}.
- THM-5040 on origin/main now says: "This is THM-5056 (mac-mini, first) and THM-5042 (independent second proof)".

**The gap.** THM-5056 (13:18:18) predates THM-5042 and cites neither it nor the update to THM-5040. Priority is
THM-5056's, by two minutes, and THM-5042 says so.

**Fix.**
- THM-5056 frontmatter: add THM-5042 to `related` (line 26). Add a status line: "Independent second proof: THM-5042
  (opus-2026-11-09-S89, 13:20, found concurrently; cones over cones and SR from six census seeds), cross-cited as
  THM-5055 / THM-5040."
- Line 17 ("whose THM-5040 note listed this question as open"): add "in its first version".
- Mirror the cross-citation in:
  - the HYP-9625 status ("PROVED by THM-5056; independently THM-5042");
  - the PROOF.md Part II statement;
  - README item 6;
  - the INDEX rows (Finding 10).

### 2. MINOR (credit): completeness and conciseness of d over all graphs is Brouwer-Csorba-Schrijver 2009, Prop. 1.2

**Where the fact is used.**
- Step 1 (THM-5056 line 63; PROOF.md lines 366-367) and the base use the fact that every odd m < 2^N is d of a
  threshold graph on N vertices (THM-4996 Lemma B1).
- The title's "the least possible over all graphs" restates that conciseness.
- HYP-9625 line 33 states it without credit.

**The attribution.** The audit of THM-5040 by opus-2026-11-09-S89 (its Finding 1) and MSG-1990 (opus to mac-mini,
13:31) report the following. The fact is Proposition 1.2 of A. E. Brouwer, P. Csorba, A. Schrijver, *The number of
dominating sets of a finite graph is odd* (2009). It is restated as Proposition 2.1(iii) of Alikhani-Fatehi-Klavžar,
arXiv:1512.07514. S. Alikhani, Graphs Combin. 29 (2013), concerns D(G, -1).
- I could not consult the preprint, since this audit has no web access.
- The attribution is consistent with THM-4996's own citation of the same preprint for the parity.

**What is affected.**
- THM-5056 never mentions Brouwer-Csorba-Schrijver.
- README lines 94-95 credit "Alikhani, Graphs Combin. 29 (2013) (completeness over all graphs)". This line is shared
  with THM-5055; MSG-1990 also flags THM-5055 line 35, which is D1's scope.
- THM-5056 claims no novelty for this fact (its "NEW here" refers to the odd-graph sharp order). So the finding is
  MINOR here, although it was MAJOR for THM-5040, which had called the conciseness new.

**Fix.**
- Add to "KNOWN and credited": "Brouwer-Csorba-Schrijver 2009 (d odd; Prop. 1.2: every odd m < 2^n is d of a graph on
  n vertices, the threshold construction of THM-4996 B1)".
- Cite it at step 1 and at HYP-9625 item 2.
- Replace the Alikhani 2013 attribution in README lines 94-95.

### 3. MINOR (proof text): the K_2-join value formula is used without proof

**Where.** THM-5056 line 63 and PROOF.md line 367 use d(K_2 v T) = 3 * 2^(L-2) + d(T). HYP-9625 line 33 states
d(K_2 v G) = 3 * 2^m + d(G). Part I proves only the E_2-join (Lemma 3).

**Why it holds.** The formula is right for m >= 1:
- sets meeting the K_2 dominate, which gives 3 * 2^m;
- a set inside V(G) must dominate G, and is then nonempty, so it also dominates the K_2.
This is the Alikhani-Peng join formula with G_1 = K_2. I checked it on 400 random graphs.

**The edge case.** For m = 0 the formula fails: d(K_2) = 3, not 4. This matters only in HYP-9625 item 2 ("every odd
number whose binary expansion starts with 11"), where 3 = d(K_2) directly.

**Fix.** Add the one-line proof (or cite Alikhani-Peng, Theorem 2) to Lemma O1(c), and write "m >= 1" in HYP-9625
item 2.

### 4. MINOR (notation): Lemma O2 and case 3

- **Undefined n.** In Lemma O2 (PROOF.md line 334), n is not defined. It is |V(G)|.
- **Two meanings of x.** The letter x denotes both the offset (x in {1, 5, 37}: PROOF.md line 370, THM-5056 line 67)
  and a twin vertex ({x, y}, in the same lines). Rename the twin pair, for example to {t, t'}.
- **A vacuous parenthesis.** Lemma O1(b) (line 323) ends with "(K_1 v K_1 = K_2)". This is the case H = K_1, which is
  not odd. Delete it, or write "(this also covers H = K_1)".

### 5. MINOR (exposition): the "Reading" overstates the role of twins

**The sentence.** PROOF.md line 388: "The reason is that a single adjacent-twin pair makes a graph odd, and it survives
every move used here."

**Why it is too strong.**
- Oddness of the base and of cones over it often comes from B2.3, not from twins.
- Exactly one odd k < 512 (other than 1, 5, 37) has no odd realiser on L(k) vertices that contains adjacent twins:
  k = 17 (`d2twins.py`).
  - The only odd realiser of 17 on 5 vertices is the bull (edges 01 02 04 13 14).
  - Its odd automorphism (0 1)(2 3) reverses only the edge 01 between its two degree-3 vertices, and it has no
    adjacent twins.
  - So in the proof's construction the values 2^(L-1) + 17, 2^(L-1) + 2^(L-3) + 17, and so on are odd by B2.3 alone.
- In the session's own outputs, 18 of the 253 base realisers, and 105 of the 4,093 graphs of `oddw_small`, have no odd
  transposition, that is, no adjacent twins. My randomly chosen base had 31.

**Fix.** "Oddness enters through an adjacent-twin pair (K_2-joins, offset moves), which every later move preserves, or
is inherited through cones (B2.3). For example, the only odd realiser of 17 on 5 vertices, the bull, has no twins."

### 6. MINOR (code documentation): `oddbuild.py` does not check d

**The docstring.** Lines 1-2 say "brute-force check of d and of an adjacent-twin pair". But:
- `dcount` (line 6) is never called;
- `twins()` (line 32) only chooses the protected pair inside `offset_graph`.

**Where the checking happens.** All of it is done afterwards by `wcheck` (set O), as THM-5056 itself says.

**Fix.** "builds ...; d, order and oddness are checked by `wcheck` (set O)".

### 7. MINOR (reproducibility): the Part II commands are missing from PROOF.md section 6

`.gitignore` line 7 says that the gitignored `oddw_small.txt`, `oddw_large.txt` and `oddw_special.txt` are
regenerated "with the commands in PROOF.md section 6". But:
- the section 6 table (PROOF.md lines 264-273) has no Part II rows;
- Part II "Checks" (lines 375-382) names the scripts without their arguments.

**Fix.** Add these rows. All of them need `dom_K511.out` and `gcc -O2 -o wcheck wcheck.c`.
- `python3 oddbase.py > oddbase.out`
- `python3 oddproof.py > oddproof.out`
- `python3 oddbuild.py 3 8192 0 oddw_small.txt > oddbuild_small.out; ./wcheck < oddw_small.txt > oddcheck_small.out`
- `python3 oddbuild.py 8192 1048576 300 oddw_large.txt > oddbuild_large.out; ./wcheck < oddw_large.txt > oddcheck_large.out`
- `python3 oddspecial.py oddw_special.txt; ./wcheck < oddw_special.txt` (these two lines make up `oddspecial.out`)

I used exactly these commands. The outputs match the committed `.out` files, and the hashes match `sha256.txt`.

### 8. MINOR (stale text): the closing sentence of HYP-9625

HYP-9625 lines 62-63 say: "The sharp-order version is therefore more irregular than for even graphs, and is not
conjectured here."
- THM-5056 proves the sharp-order version, and it is completely regular: every k is realised on L(k) vertices.
- The status note supersedes only "its list of open cases".
- The misses came from the move set of `oddclos.py`, which has no hub moves. For example, 517 = 2^9 + 5 is
  PH(K_4 - e, w, 2) followed by PH(3, 2), on 10 vertices.

**Fix.** Append "[Superseded: THM-5056 proves the sharp-order version.]"

### 9. MINOR (wording): "needs a prime graph" (README line 66)

PH(K_2 u T_55, w, 2) contains the module {x, y} (adjacent twins), so it is not prime in the modular-decomposition
sense.

**What is true** (checked by hand): every odd realiser of 1061 on 11 vertices is connected and co-connected, that is,
neither a union nor a join.
- 1061 is prime, and a component on <= 10 vertices has d <= 1023, so the realiser is connected.
- A join A v B has d = (2^a - 1)(2^b - 1) + d(A) + d(B), with a + b = 11. This exceeds 1061 unless a = 1.
- When a = 1, the graph is the cone over a 10-vertex realiser of 37. That realiser is even and has no universal
  vertex, so the cone is even (B2.3).

**Fix.** "needs a graph that is neither a union nor a join".

### 10. MINOR (bookkeeping): the INDEX rows

The rows are `05-knowledge/hypotheses/INDEX.md` line 2 and `05-knowledge/results/INDEX.md` line 6.
- They say only "audit D1 pending". THM-5056's audit is D2.
- "pre-registered R1-R3 CONFIRMED" applies to THM-5055 only. THM-5056 was exploratory: PREREGISTRATION.md makes no
  odd-graph prediction.
- The hypotheses row types HYP-9625 "FINITE-EXACT", although it is PROVED.
- THM-5042 is not mentioned (Finding 1).

**Fix.** Correct each of these.

### 11. MINOR (typing; optional strengthening): the exclusion of 37 needs only a small enumeration

THM-5056 item 2 (line 47) is typed "(PROVED; exclusions from THM-4996)". It rests on the 474,101-class census, which
THM-4996 types FINITE-EXACT.

**Re-derived independently.**
- The least d of a connected 7-vertex graph is 45 (exhaustive enumeration).
- So by B3 monotonicity every connected graph with d <= 37 has <= 6 vertices.
- The connected graphs with d = 37 are exactly four classes on 6 vertices, with |Aut| = 4, 1, 4, 2. All four are
  even.
- d = 5 is realised only by P_3, which is even.

**Suggestion.** "(PROVED; a finite check: by B3 monotonicity it suffices to enumerate connected graphs on <= 7
vertices; THM-4996 census, re-enumerated independently in audit D2)". This is not an error.

## Observations (no change needed)

1. The odd-good requirement that the pendants' neighbours be non-universal is not needed for oddness. It is needed
   only for evenness in Part I (Lemma 4.2). Here it is harmless and maintains itself.
2. In case 2 the padded realiser always has an isolated vertex, so the twin branch of O1(b) is never used in the
   induction.
3. Every admissible choice works. My construction chose the excluded vertices, the seeds' w and the move order at
   random, and all 8,699 graphs passed.

## Outside this audit's scope (for D1 / the owner)

- **PROOF.md Part I, remark 7.2 (lines 307-308).** It says "The values just below a power of 2 ... need prime graphs".
  - Top-band values can be disjoint unions. For example, 505 = d(P_3 u G_101) on 10 vertices, with G_101 an even
    connected graph on 7 vertices (census, for instance `010203041516252635364656`).
  - 505 lies in T(10) = (M(9), 2^9) = (503, 512).
  - So if "prime" is meant to exclude unions, the remark is false. This is the issue that opus's audit found in
    THM-5040 Corollary 4, and MSG-1990 asks mac-mini to check THM-5055 for it.
- **THM-5042's status line.** It says that THM-5056 was "in canon about an hour earlier". The commit times are 13:18:18
  (c11abf3a22) and 13:20:22 (7e86bdf8f9). This is opus's file.

## Reproduction

All commands run from `04-computation/experiments/digraph721_20261109/`. Nothing is third-party code. At most two
heavy processes ran at once.

1. **Census.** `gcc -O2 -o audit_d2/dom_even_d2 ../paleyeven_20261109/dom_even.c && ./audit_d2/dom_even_d2 511 16 >
   audit_d2/dom_K511_regen.out` (17 s).
   - The SHA-256 is 6307dfa1d0fc58de7cd0990b1c4c948064d2be4f4a0c5846f9135df50d0afd0e, the same as `sha256.txt`.
   - The per-order counts in `audit_d2/dom_K511_regen.err` equal those in `dom_K511.err`.
   - The 24 MB copy was deleted after hashing.
2. **Everything else.** `cd audit_d2 && sh run_d2.sh > run_d2.out 2>&1` (35 s). Results:

| step | result |
|---|---|
| (A) session witnesses, regenerated with the session's scripts | hashes f83efc45..., 38f92621..., 9feb2f3d... equal `sha256.txt`; `wcheck`: 4,093 + 300 + 57 ok; `d2check`: 0 failures (odd transpositions 3,988 / 298 / 57, other automorphisms 105 / 2 / 0) |
| (B) `d2check` controls | K_2, P_4 (non-transposition certificate) and K_4 - e pass; P_3, C_4, E_2, a wrong value and a wrong order fail |
| (C) `d2small 7`: all labelled connected graphs, n <= 7, classes by orbit counting | 1, 1, 2, 6, 21, 112, 853 classes; even/odd split equal to the census; all 136 (n, d, parity) buckets equal (996 classes); d = 5: P_3 only (even); d = 37: four classes on 6 vertices, all even; least d on 7 vertices is 45 |
| (D) `d2allsmall`: all labelled graphs, n <= 7 | {d(G) : G odd on n vertices} = {odd k <= 2^n - 1} minus {1, 5, 37} for every n <= 7 (0, 1, 2, 6, 14, 29, 61 values) |
| (E) `d2base.py` | 253 of 256 odd k < 512 realised on exactly L(k) vertices; none for 1, 5, 37; the unions pass `d2check` (253 ok) |
| (F) `d2build.py`: independent construction, random admissible choices | all odd k in [512, 16384): 7,936; 500 random k in [16384, 2^21); 263 special graphs (2^(L-1) + c, 3 * 2^(L-2) + c and 2^(L-1) + 2^(L-3) + c for 10 <= L <= 21; deep cone chains from 60 random draws; three offset graphs at every order of Gamma_odd up to 21). `d2check` and `wcheck`: 0 failures; up to 21 vertices |
| (G) `d2seeds.py` | as in item 3 above |
| (H) `d2lemmas.py` | Lemma 4.1: 5,331 tests; K_2-join and cone values: 400 each; O1(b): 879 (313 with a universal vertex); O1(c): 150; O1(d): 473. 0 failures |
| (I) `d2induct.py 24` | the Gamma_odd sets above; the only gap is (1, 10); no odd k < 2^24 other than 1, 5, 37 fails |
| (J) `d2twins.py`; `d2check full < c17.txt` | 17 is the only base value without a twin-carrying realiser; the bull is its only odd 5-vertex realiser (|Aut| = 2, one odd element) |

The hashes of the regenerable outputs are in `audit_d2/sha256_d2.txt`. Binaries and large outputs are gitignored by
`audit_d2/.gitignore`.

## Files written by this audit

- **This report:** `05-knowledge/results/digraph721_20261109/AUDIT_D2.md`.
- **Checkers (C):** `d2check.c` (brute-force checker), `d2small.c` (connected graphs, n <= 7, orbit counting) and
  `d2allsmall.c` (the final statement, n <= 7).
- **Python scripts:**
  - `d2cmp_small.py`, which compares (C) with the census;
  - `d2base.py`, for the base;
  - `d2build.py`, the independent construction;
  - `d2seeds.py`, `d2lemmas.py`, `d2induct.py` and `d2twins.py`.
- **Inputs:** `controls.txt` and `c17.txt`.
- **Runner:** `run_d2.sh`.
- **Outputs:**
  - `run_d2.out`, `d2small.out`, `d2allsmall.out`, `d2base.out`, `d2lemmas.out`, `d2base_witnesses.txt`,
    `dom_K511_regen.err` and `sha256_d2.txt`;
  - gitignored and regenerable: `oddw_*_regen.txt` and `d2w_*.txt`.
- **Location:** everything except this report is in `04-computation/experiments/digraph721_20261109/audit_d2/`.
