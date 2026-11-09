# Audit 2: THM-5042, THM-5041 and HYP-9610 (opus-2026-11-09-S89)

VERDICT: PASS WITH CORRECTIONS

Independent adversarial audit, 2026-11-09, of the worktree `C:/tmp/wt-s89` at HEAD `a2db8851ab`. Nothing outside this
directory was edited.

**What I read in full.**
- The three items: THM-5042, THM-5041, HYP-9610.
- The canon they rest on: THM-5056, THM-5040, THM-4996, THM-5071, THM-830 section 1, and Session234 `PROOFS.md`
  (DG1-DG8).
- Session89 `encoding/ENCODING_AND_ACTIONS.md`.
- HYP-9625, MISTAKE-747, and the messages MSG-2473 and MSG-2474 (mac-mini to opus) and MSG-1989 (opus to mac-mini).
- The scripts `odd_companion.py`, `tile_additivity.py`, `cover.py`, `e0search.c` and `egcp_lib.py`.

All code in this directory is new and shares nothing with `egcp_lib.py` or `e0search.c`. It has:
- two dominating-set counters (a subset DP in C, and a direct subset-by-subset check in Python);
- two automorphism engines (C backtracking, and networkx VF2);
- a Hamiltonian-path counter that enumerates every vertex ordering.

## Overall

**The mathematics is correct in all three items.**
- **THM-5042.** The theorem and its proof are right.
  - Confirmed exhaustively for n <= 8, over all 2^28 labelled graphs at n = 8.
  - I rebuilt from the proof text every value that needs exactly n vertices for 8 <= n <= 13 (4,032 graphs) and
    verified each one.
- **THM-5041.** Every identity checks out by brute force, including the reversal formula against THM-830 (1.3).
- **HYP-9610 (A).**
  - The "proved for all n" parts are right.
  - FINITE-EXACT holds for n <= 8, confirmed independently.
  - Containment holds at n = 10, 11, 12 and 13, verified directly.
  - The exclusion listed as open (M0(n) - 2 for odd n >= 11) can be proved now (Note 1).

**Corrections.**
- **Three MAJOR:**
  - two are credit errors (Findings 1 and 2);
  - one is a false verification claim caused by a buffer overflow in `e0search.c` (Finding 3). Its conclusion is
    nevertheless true.
- **Seven MINOR** and **ten NOTES**.

No correction changes a mathematical statement of THM-5042 or THM-5041, or the conjecture of HYP-9610.

## Findings

### 1. MAJOR (credit, THM-5042): the priority paragraph misdescribes both proofs and the timeline

**Location.**
- THM-5042, lines 5-9 (status) and line 23 (banner).
- Repeated in `05-knowledge/results/egcp_20261109/README.md` ("Companion results"), the INDEX row ("by cones over
  cones and SR seeds"), MISTAKE-747 (3) and MSG-1989.

**(a) "The two proofs share the cone step and differ in the leading-11 step": false.**
- K_1 v (K_1 v G) and K_2 v G are the same graph.
- Lemma A is THM-5056 step 1 verbatim: the same graph, the same value 3 * 2^(n-2) + d(T), and the same
  adjacent-twin argument.
- Lemma B is THM-5056 step 2: a cone over a padded odd realiser, odd by B2.3.
- Both lemmas are already "proved reductions" 2 and 3 of HYP-9625 (mac-mini, commit `153409f791`, 13:08:19).

**(b) SR is not a difference either.**
- THM-5056's offset moves are R = PH(pendant, 1) and PH(w, 2) (`digraph721_20261109/PROOF.md`, Part II, Lemma O2).
- THM-5040 section 6 already identifies R with SR.

**What genuinely differs is Lemma C:**
- **Oddness invariant.** THM-5042 carries an odd automorphism that fixes the SR leaf. THM-5056 carries a persistent
  pair of adjacent twins.
- **Seeds.** THM-5042 uses six census graphs at orders 8 and 9 and SR only. THM-5056 uses four small explicit seeds,
  shifts of +2 and +3 (R and PH(w, 2)), and the special case 513 = d(K_2 u T_171).

**(c) "in canon about an hour earlier": false.**
- THM-5056 landed in `c11abf3a22` at 13:18:18.
- THM-5042 landed in `7e86bdf8f9`, authored 13:19:19 and committed 13:20:22.
- The gap is about one minute. MISTAKE-747 (3) itself dates mac-mini's message to 13:17.

**(d) "found without knowledge of it" cannot cover Lemmas A and B.**
- MSG-2473 (mac-mini to opus, sent 13:09, commit `9015b7995a`) states: "K_2 v G covers every leading 11; cones over
  odd graphs are odd; ... The first open cases are 2^m + {1, 5, 37}".
- MISTAKE-747 says S89 last fetched at about 13:10. MSG-1989 replies to MSG-2473 by number.
- Whatever the reading order, Lemmas A and B and the reduction to the offsets {1, 5, 37} predate THM-5042 by about ten
  minutes, so they must be credited.

**Fix.**
- Rewrite lines 5-9 along these lines:
  > Second proof of THM-5056. Lemmas A and B are THM-5056 steps 1-2, already stated as HYP-9625 reductions 2-3
  > (mac-mini; MSG-2473, 13:09); note K_1 v (K_1 v G) = K_2 v G. The independent part is Lemma C: SR carries
  > 2^(m-1) + {1, 5, 37} upward while preserving an odd automorphism that fixes the SR leaf, from six census seeds at
  > orders 8 and 9. THM-5056 instead preserves a pair of adjacent twins through R (= SR) and PH(w, 2), from four small
  > seeds, plus the special case 513.
- Change "about an hour earlier" to "about one minute earlier (13:18 vs 13:19)".
- Add the 13:08/13:09 items to MISTAKE-747 (3).
- Update the README paragraph and the INDEX row to match.

### 2. MAJOR (credit, THM-5041): both items called "New" have repository prior art

**Location.**
- THM-5041 lines 5-9. The status says: "New: the pair-of-tilings coordinates for arbitrary digraphs, and the tile-wise
  additivity with its metagraph-flip corollary".
- Lines 27-28, where Session89 is mentioned only as a different split of n^2.

**(a) The coordinates are Session89's.** The source is
`05-knowledge/results/session89_consecutive_tournament_pairs_20261028/encoding/ENCODING_AND_ACTIONS.md` (2026-10-28).
- **Section 1** sets x_ij = M_ij, y_ij = M_ji (i < j) and d_i = M_ii. These are exactly (U, L, diag).
  - Each triangle orients a labelled tournament, which is a tiling of order n + 1 by THM-830 section 1.
  - "Imposing T = S retains only ... the symmetric matrices", i.e. undirected digraphs are the equal pairs.
- **Section 6:** "Matrix transposition acts by tau(x,y,d) = (y,x,d)". This is the transpose law of (1d).
- **Section 4:** the discrepancy graph r = x + y. Its zero set is THM-5041's defect set Delta.

The only difference is that Session89 attaches the diagonal to the second tournament as root arcs.

**(b) Tile-wise additivity follows at once from Session234 DG7.**
- Adding an absent arc a -> b raises H by sum over S of E_D(S, a) B_D(V\S, b).
- That sum counts paths inside S (which contains a) and inside V\S (which contains b), so it never involves the arc
  b -> a.
- Hence H(1,1) - H(0,1) = H(1,0) - H(0,0) at every pair of every digraph. That is (2b), and (2c) is (2b) on a
  tournament.
- THM-5041 credits the displayed example 17 + 19 - 15 = 21, but not DG7.

**Fix.**
- Credit Session89 sections 1, 4 and 6 for the coordinates, the transpose law, symmetric = equal pairs and the
  discrepancy graph.
- Credit Session234 DG7 for (2b).
- Restrict "New" to what remains:
  - the tiling reading through THM-830, with the convention fixed (Note 2);
  - the class table;
  - the reversal formula;
  - the explicit flip corollary and the signed expansion.
- Consider typing THM-5041 as a dictionary, not as a new theorem. Every identity in it is verified (checks below).

### 3. MAJOR (false verification claim, HYP-9610 (A)): `e0search.c` overflows `found[16]`, so 6 of the 38 "realised" values were never found

**Location.**
- `04-computation/experiments/egcp_20261109/e0search.c` line 18: `int found[16]={0};`. Up to 64 targets are read,
  and the runs used 18 and 20.
- `e0search_n10.out` and `e0search_n11.out`.
- HYP-9610 lines 8-9 and 40-42.
- `04-computation/experiments/egcp_20261109/README.md` line 16.
- `05-knowledge/results/egcp_20261109/README.md` line 50.
- MSG-1989 ("all 38 closure misses at n = 10, 11 realised").

**Evidence** (`a2_e0verify.out`).
- `e0search_n10.out` prints a FOUND line, with a graph, for only 16 of its 18 targets. 773 and 779 (target indices 16
  and 17) are marked "realised" with no graph.
- `e0search_n11.out` does the same for 16 of 20. 1533, 1537, 1541 and 1547 (indices 16-19) have no graph.
- The unchanged source reproduces this.
  - With 1000 iterations, at n = 10 and at n = 11, every target at index >= 16 prints "realised" with no FOUND line.
  - With 0 iterations (n = 10), 773 and 779 print "not found".
  - So found[t] for t >= 16 reads past the array: this is undefined behaviour, and in this build it aliases the
    neighbouring stack array adj[].
- The 32 graphs that were printed are all correct: d, evenness and no universal vertex all check, 0 bad.

**The conclusion is nevertheless true.**
- I found even graphs with no universal vertex for all six values (`a2_six_witnesses.out`, edge lists there).
  - Five have a trivial automorphism group and one has |Aut| = 2.
  - Each was checked by both d-counters and both automorphism engines.
- More strongly, a direct search (`a2_e0hunt.c` with `a2_e0full.py`) realised every value of the conjectured E0(n)
  at four orders. Each witness was re-verified by `a2_check`.

  | n | values realised |
  | --- | --- |
  | 10 | 503 of 503 |
  | 11 | 1,013 of 1,013 |
  | 12 | 2,038 of 2,038 |
  | 13 | 4,084 of 4,084 |

**Fix.**
- Make the array `int found[64]` and rerun, or cite `a2_six_witnesses.out` and `a2_e0full_n10.out` to
  `a2_e0full_n13.out`.
- Correct "18 at n = 10 and 20 at n = 11" and "all 38" in the three files above.
- HYP-9610 can now say: "containment verified at n = 10-13 by direct search (audit2)".

### 4. MINOR (HYP-9610 lines 21-22): "These are the co-connected or padded realisers" is false

A graph with no universal vertex need not be co-connected or padded.
- The cocktail-party graph comp(perfect matching) attains M0(n) for even n, as the HYP itself says. Its complement is
  disconnected, and it has no isolated vertex.
- C_4 = comp(2K_2) is another example.

**Fix:** "These are the realisers with no universal vertex: exactly the inputs of an even cone (B2.3)."

### 5. MINOR (HYP-9610 table, line 53): A334335 counts connected even graphs

- OEIS A334335 is "Inverse Euler transform of A000568": 1, 0, 1, 2, 8, 43, 398, 6413, 184596, ... These are the
  connected even graphs.
- The number of even graphs on n vertices is A000568, the tournament count, by RPGFD.
- My census reproduces both sequences for n <= 8.
- Also, RPGFD's equinumerosity is for isomorphism classes, while the tournament row counts labelled tournaments
  (2^C(n,2)). The labelled counts differ: at n = 2 there are 2 labelled tournaments and 1 labelled even graph.

**Fix:** write "A000568 (unlabelled; connected: A334335)".

### 6. MINOR (HYP-9610 table, lines 49, 51, 52): sources and status of three exceptional sets

- **Tournaments (line 49).** The exclusion of 7 and 21 is THM-1370 (h-spectrum-omits-7-21-all-n). The row cites
  THM-830 and Session234, and Session234 itself cites THM-1370.
- **One reciprocal pair (line 51).**
  - THM-5071 (4b) proves that 1, 4 and 7 are never attained.
  - That nothing else is missing is only FINITE-EXACT, below 198 and for orders <= 7 (4c).
  - Add "(complete below 198)".
- **Base 3 (line 52).**
  - 10 and 31 are excluded by THM-4940 (3) and (4).
  - Completeness of the rest is HYP-9530 (A), which is OPEN; it is verified to 10^6 by THM-5070 (4).
  - Write "conjecturally complete", as the tournament row does ("complete iff Chan-Pak").

### 7. MINOR (THM-5041 line 81): "sum over A in R" should be "sum over A ⊆ R"

The signed expansion runs over subsets A of R. I verified it in that form on 120 random semicomplete digraphs.

### 8. MINOR (THM-5041 line 100): "misses 7 and 21, and conjecturally nothing else"

Tournaments also miss every even number. **Fix:** "misses 7 and 21 among the odd numbers".

### 9. MINOR (THM-5041 line 102): credit for {1, 5, 37}

This line credits THM-5042 alone for the odd-graph exceptional set. **Fix:** cite "THM-5056 (first), THM-5042", as
HYP-9610 and THM-5040 do.

### 10. MINOR (THM-5041 line 42): "U and L = 0 at every tile" reads as "both are 0"

Read literally, that is the empty digraph. **Fix:** "U_p L_p = 0 (never both 1) at every tile".

## Notes

**Note 1 (positive, HYP-9610 line 43).** The exclusion of M0(n) - 2 holds for every n, so the "Open" item for odd
n >= 11 can be closed. The checks are `a2_excess.py` and `a2_excess_atlas.py`.

*Setting.* Let G have no universal vertex, and put H = comp(G). Then H has no isolated vertex, and
D(H) = 1 + n + e(H), where e(H) counts the sets of size >= 2 inside some open neighbourhood. G is even iff H is
co-even (B2.5).

*Facts about e.*
- e is additive over components (Brown-Nowakowski), and e(K_2) = 0.
- Adding a vertex to a connected graph raises e by at least 1, and deleting a non-cut vertex does not raise it. So
  every other component C has e(C) >= |C| - 2.

*Co-even components.*
- A disjoint union is co-even only if every component is: an automorphism acting on one component has the same
  edge-sign and permutation sign there.
- The connected graphs with e <= 5 are 12 classes, all on <= 7 vertices. I listed them twice: by growth with
  isomorphism rejection, and from the networkx atlas.
- Of these, only three are co-even: K_3 (e = 3), P_5 (e = 3) and P_6 (e = 4).

*Conclusion.*
- No sum of these e-values equals 1, 2 or 5. So e(H) is never 1, 2 or 5 for co-even H, and
  d = 2^n - n - 1 - e is never 2^n - n - 2, 2^n - n - 3 or 2^n - n - 6.
- For odd n, d is odd only when e is odd. So M0(n) comes from e = 3 (P_5, and also K_3), and M0(n) - 2 would need
  e = 5, which is impossible.
- For even n, e must be even, and e = 2 is impossible.
- Hence (A)'s exclusion side, values above M0(n) and M0(n) - 2, is PROVED for every n. Only containment remains open.

**Note 2 (THM-5041 (1d), line 52).**
- **Define s.** Write (sX)_ij = X_(n+1-j, n+1-i). This is THM-830's staircase reflection sigma, (a, b) ->
  (n+1-b, n+1-a), transported to gap pairs.
- **State the convention.** Under THM-830 (1.1), (G_t)_(i,j) = t_(i+1,j), so the tiling of the tournament A is L(A),
  not U(A); U(A) is the tiling of A^op.
- **Checks** (`a2_thm830.py`), over all tilings for n = 4, 5, 6 with 0 failures:
  - (1.3) holds with this sigma;
  - (U, L) -> (sL, sU) under i -> n+1-i holds.
- **The blue-line statement** is correct. Blue means sU = U, which says the reversal is an anti-automorphism.

**Note 3 (THM-5041 line 63).** Define "defect cycles":
- cycles of D through a (1,1) tile;
- cycles of the loopless complement through a (0,0) tile (THM-5071 (2)).

The parity claim was verified on 300 random general digraphs, by explicit packing counts.

**Note 4 (THM-5041 line 73).** "T a tournament (a line node)" mixes two sorts of the metagraph. In THM-830, nodes
are fibres of Phi and lines are complement pairs. Say "an ordered line".

**Note 5 (THM-5042 line 101; THM-5040 line 320).** "No map carries H to d (THM-4996 (4), THM-4531)" needs two fixes.
- THM-4996 (4) is the I(G, 2) census, not a no-map statement.
- For value-preserving bijections at a fixed order, the claim is already settled at n = 3: tournaments have
  H in {1, 3}, while even graphs have d in {1, 5}.
- For "no natural map", cite THM-4996 "NOT CLAIMED" and THM-4531.

**Note 6 (HYP-9610 lines 61-64).** The twin explanation in "Reading" does not cover 37.
- Its four connected realisers all have 6 vertices. Their automorphism groups have orders 4, 1, 4 and 2
  (`a2_real37.py`).
- One realiser is asymmetric. Another has an even involution that is not a twin swap.
- THM-4996's own table attributes {7, 21} to strong-component floors, not to twins.

**Note 7 (HYP-9610 lines 29-33).** Two simplifications for odd n:
- D(H) = 2^n - d is odd, so D = n + 1 and D = n + 3 are excluded by parity alone. Only D = n + 2 needs the co-even
  argument.
- comp(K_3 u sK_2) also attains M0(n), since K_3 is co-even.

**Note 8 (THM-5042 exclusions and base).** Both can be shortened.
- **Exclusions without the d <= 511 census.**
  - A connected realiser of 37 has at most 7 vertices. Count the dominating sets of size n - 3 and above:
    d >= 1 + n + (C(n,2) - L) + (C(n,3) - (n - 2)L - L_2), where L and L_2 count the vertices of degree 1 and 2. This
    gives d >= 43 at n = 8, d >= 65 at n = 9, and d >= 2 + C(n,2) >= 47 for n >= 10.
  - The exhaustive census here then gives exactly four realisers, all even.
- **Smaller base.** The induction can start at n = 8 from the base n <= 7 plus the six seeds. That is how
  `a2_thm5042.py` rebuilt every value.

**Note 9 (`odd_companion.py`).** The script's main check reads `dom511.out`, which is not committed.
- Without that file it re-checks only the seeds, and still prints "ALL: PASS".
- Suggest printing "census check skipped" when the file is missing.

**Note 10 (THM-5041 section 3).** Session234 sharpens DG2. A semicomplete digraph with H = 7 has exactly two
reciprocal pairs, and one with H = 21 has one, two or three. I confirmed this exhaustively for semicomplete n <= 5:
reciprocal-pair counts {2} at H = 7, and {2, 3} at H = 21 for n = 5.

## Independent checks (commands and key outputs)

All commands run from `05-knowledge/results/egcp_20261109/audit2/`. The executables are not committed; rebuild them
with `gcc -O3 -fopenmp -o a2_census.exe a2_census.c`, `gcc -O3 -o a2_check.exe a2_check.c` and
`gcc -O3 -o a2_e0hunt.exe a2_e0hunt.c -lm`.

**1. Exhaustive census of all labelled graphs, n <= 8.**
- Commands: `./a2_census.exe n > a2_census_n<n>.out` for n = 1..8 (n = 8: 268,435,456 graphs, 53 s on 12 threads),
  then `python a2_census_analyse.py 8`.
- The engine computes the full automorphism group of every graph, with the sign of each automorphism.
- **Validation.**
  - Classes: 1, 2, 4, 11, 34, 156, 1044, 12346 (A000088).
  - Even classes: 1, 1, 2, 4, 12, 56, 456, 6880 (A000568, RPGFD).
  - Connected classes: 1, 1, 2, 6, 21, 112, 853, 11117.
  - Connected even classes: 1, 0, 1, 2, 8, 43, 398, 6413 (A334335, THM-4996).
- **Results** for every n <= 8 (output: `a2_census_analyse.out`, all True):
  - E_odd(n) = {odd k <= 2^n - 1} minus {1, 5, 37};
  - E(n) = O(n) (THM-5040);
  - E(n) and E0(n) equal `base_sets.json`;
  - for 5 <= n <= 8, E0(n) equals the HYP-9610 formula.
- **Exclusions.**
  - No odd graph at 1, 5 or 37.
  - The connected realisers of 37 are four classes, all even, all on 6 vertices.
  - The only connected realiser of 5 is P_3.
- `base_sets.json` E0 at n = 9 also equals the formula. That is a consistency check, not an independent one.

**2. THM-5042, rebuilt from scratch.** Command: `python a2_thm5042.py 13`.
- **Seeds.** All six check out:
  - d by both counters;
  - a leaf q;
  - no universal vertex;
  - an odd automorphism fixing q, from C and from networkx.
- **Lemma C chain.** For m = 8..16 and c in {1, 5, 37} (27 graphs), each has d = 2^(m-1) + c, is odd, has the leaf
  fixed by an odd automorphism, and has no universal vertex. 0 failures.
- **Every new value at orders 8..13.** Built by Lemmas A, B and C from my census base (n <= 7) and the seeds, with
  64, 128, 256, 512, 1024 and 2048 values. 0 failures. A random dozen per order for n <= 10 were re-checked with
  networkx and the direct counter.
- **Re-runs of the packet's scripts.** `python odd_companion.py` (seed part) gives ALL: PASS.
  `python tile_additivity.py` gives 300 bad 0.

**3. THM-5041.** Command: `python a2_thm5041.py`; H is computed over all n! orderings. 0 failures throughout.
- Class table and counts, exhaustive at n = 3.
- Transpose, complement J - A and reversal (sL, sU): 400 random digraphs each.
- Loops irrelevant, and H(A^T) = H(A).
- (2b) at every pair: 2,942 tile tests on 400 random digraphs with n <= 6.
- H = 2 hp on 200 symmetric digraphs with n <= 7.
- THM-5071 (3) by explicit defect-cycle packings: 300 random general digraphs.
- (2c) at every tile of 150 random tournaments with n <= 7: 1,508 tests each of the identity, the parity statement and
  "difference = HPs through p".
- The signed expansion: 120 random semicomplete digraphs.
- The Session234 examples in (U, L) coordinates:
  - H = 7 with two (1,1) tiles;
  - H = 21 with two (1,1) tiles;
  - H = 21 with one (1,1) tile, whose resolutions and missing state give 17, 19 and 15.
- DG2, exhaustive over semicomplete digraphs with n <= 5.
- Also `python a2_thm830.py`: THM-830 (1.1)/(1.3) conventions and (1d), all tilings for n = 4, 5, 6.

**4. HYP-9610.**
- `python a2_e0verify.py`: the 32 printed e0search graphs are all correct. The values claimed realised with no graph
  are [773, 779] and [1533, 1537, 1541, 1547], at target indices 16-19.
- **Overflow reproduction.**
  - Commands: `gcc -O3 -o e0s ../../../../04-computation/experiments/egcp_20261109/e0search.c`, then
    `./e0s 1000 11 7 <the 20 targets>`, `./e0s 1000 10 3 <the 18 targets>` and `./e0s 0 10 1 <the 18 targets>`.
  - With 1000 iterations, 1533, 1537, 1541, 1547 (n = 11) and 773, 779 (n = 10) are "realised" with no FOUND line.
  - With 0 iterations, 773 and 779 print "not found".
- `python a2_six_witnesses.py`: the six witnesses, PASS.
- `python a2_e0full.py n` for n = 10..13: every conjectured value realised and re-verified, with 0 failures (see the
  table in Finding 3).
- `python a2_excess.py` and `python a2_excess_atlas.py` give the 12 connected classes with e <= 5. The only co-even
  ones are K_3, P_5 and P_6. For E = 1, 2 and 5 there are 0 co-even multisets (Note 1).
- `python a2_real37.py`: the four realisers of 37 and their automorphism groups (Note 6).
- `python ../../../../04-computation/experiments/egcp_20261109/cover.py 11 --show` reproduces
  `e0_closure_misses.txt`: 18 values at n = 10 and 20 at n = 11.
- OEIS A334335 was read on oeis.org: "Inverse Euler transform of A000568", the connected even graphs.

## Files in this directory

| File | Content |
| --- | --- |
| `a2_census.c`, `a2_census_n1..8.out`, `a2_census_analyse.py`, `a2_census_analyse.out` | exhaustive labelled census (n <= 8) with full automorphism groups |
| `a2_check.c`, `a2_lib.py` | independent graph checker and helpers |
| `a2_thm5042.py`, `a2_thm5042.out` | THM-5042 rebuilt and verified, orders 8-13; seeds; Lemma C chain to order 16 |
| `a2_thm5041.py`, `a2_thm5041.out`, `a2_thm830.py`, `a2_thm830.out` | THM-5041 identities; THM-830 conventions |
| `a2_e0verify.py`, `a2_e0verify.out` | re-check of the e0search outputs (Finding 3) |
| `a2_e0hunt.c`, `a2_e0full.py`, `a2_e0full_n10..13.out`, `a2_e0full_n10..13_witnesses.txt` | direct search of E0(n), n = 10-13; each witness line is `d` followed by the n neighbour bitmasks of G |
| `a2_six_witnesses.py`, `a2_six_witnesses.out` | the six values without a printed graph |
| `a2_excess.py`, `a2_excess.out`, `a2_excess_atlas.py`, `a2_excess_atlas.out` | the finite reduction of Note 1 |
| `a2_real37.py`, `a2_real37.out` | the realisers of 37 |
