# Follow-up audit of commit a5081c2031 (opus-2026-11-09-S87)

Independent auditor, 2026-11-09. Scope:
- the new THM-4951 (3e) and (3f), with the consequences for HYP-9540 (A) and (C);
- the corrections for M1, M2 and m1-m17 of [AUDIT.md](AUDIT.md), in THM-4950, THM-4951, HYP-9540, HYP-9530, the
  README, both INDEX rows, PREREGISTRATION.md, PREREGISTRATION-2.md and nu_lines.c;
- MISTAKE-742.

All new checks use code written for this follow-up, in this directory: `runwindow_audit.c`,
`runwindow_n3k_classify.py`, `runwindow_n3k_runs.py`.

## Verdict: PASS WITH CORRECTIONS

- The new clauses (3e) and (3f) are correct, and so are their proofs. They prove HYP-9540 (A) and (C) at every
  n != 0 (mod 3), as claimed.
- All corrections for M1, M2 and m1-m17 are applied correctly and completely.
- MISTAKE-742 describes the findings accurately.
- One new MAJOR finding (F1): a false finite claim about n = 3k, repeated in four places. It is a one-line fix in
  each place and touches no PROVED clause. There are also four minors (F2-F5).

## 1. Findings

### F1 (MAJOR). The run-window tiling at n = 3k does not fail "exactly when an isolated 0 separates two 1-runs"

**Where it appears.**
- THM-4951 (3f), third bullet group (line 165): "It fails exactly when an isolated 0 separates two 1-runs
  (`run_window_check.out`)."
- THM-4951 Remarks (line 378): "The words where an isolated 0 joins two 1-runs need other constructions."
- HYP-9540 (A), "The open case n = 3k" (line 52): "It fails exactly when an isolated 0 separates two 1-runs, so the
  windows overlap."
- The docstring of `run_window_check.py` (line 9).

**Why it is false.** `run_window_check.out` lists only the failures, so it cannot support "exactly".
- The THM's own count "5 of 5 at n = 6" includes the word 1010. That word has an isolated 0 between two 1-runs, and
  its windows [0, 2] and [2, 4] overlap. The run-window tiling still realises it: nu = (1, 1), and the word is right.
- At n = 9 the words 1010000, 1111010, 1110101 and 1010111 also have an isolated 0 between two 1-runs, and all are
  realised. The only failure at n = 9 is 1110111.

**The data** (`runwindow_audit.c n3k`, `runwindow_n3k.out`), tight colourings with x_0 = 1:

| n | tight colourings | realised by the run-window tiling | with an isolated 0 between 1-runs (realised) | failures |
|---|---|---|---|---|
| 6 | 5 | 5 | 1 (1) | 0 |
| 9 | 17 | 16 | 5 (4) | 1 |
| 12 | 50 | 47 | 15 (12) | 3 |
| 15 | 136 | 128 | 40 (32) | 8 |
| 18 | 352 | 332 | 100 (80) | 20 |

The 6, 9 and 12 columns agree with the session's 5/5, 16/17, 47/50.

**The correct pattern** (`runwindow_n3k_runs.out`), for every n = 3k <= 18:
- The tiling fails exactly when some isolated 0 lies between two 1-runs that both have length >= 2.
- No 1-run of length exactly 2 is ever adjacent to an isolated 0 in these tight colourings, so ">= 2" here means
  ">= 3".
- Every failing word contains 1110111.
- Every word whose isolated 0s each touch a 1-run of length 1 is realised.

**Fix.**
- Replace the four sentences with:

  > It fails exactly on the tight colourings in which an isolated 0 lies between two 1-runs of length >= 2 (1, 3, 8,
  > 20 words at n = 9, 12, 15, 18; checked for n = 3k <= 18).
- In the Remarks, write "the words with an isolated 0 between two long 1-runs (such as 1110111 at n = 9) need other
  constructions".
- Label the n = 3k counts in (3f) FINITE-EXACT. They are a census, not part of the proved theorem (see F3).

### F2 (minor). THM-4950 status: the HYP-9541 sentence contradicts (5)'s convention note

- The status (line 26) says "HYP-9541 (4) conjectures that 64 is the only value 1 (mod 3) attained but never by a
  strong tournament of order >= 3".
- With order >= 3 the value 1 qualifies too, as (5)'s own note says: "the values attained but never by a strong
  tournament of order >= 3 are exactly 1 and 64".
- HYP-9541 (4)'s sentence itself has no order qualifier.
- Write instead "... that 64 is the only value v > 1, v = 1 (mod 3), attained but never by a strong tournament
  (HYP-9541's title, which leaves order 1 out, lists 1 and 64)".

### F3 (minor). THM-4951 status: "(3a)-(3f)" is listed as PROVED

- (3f) also contains the n = 3k census (5/5, 16/17, 47/50, and the failure pattern). That census is FINITE-EXACT,
  not part of the theorem.
- Either move it to a remark or list it under FINITE-EXACT. The FINITE-EXACT list already holds the n = 3k count of
  (3c).

### F4 (minor). THM-4950 title, (5): add "v > 1"

- The body states (5) for v > 1. The title still reads "the values v = 1 (mod 3) that are H_3 of no strong
  tournament are exactly 10, 31 and 64". Under HYP-9541's convention, which leaves order 1 out, that would also
  include 1.
- Add "v > 1", as in the body.

### F5 (minor). Pending markers

- THM-4951 status "Their follow-up audit is PENDING (audit/FOLLOWUP.md)" and the README row "(3f) ... PROVED
  (follow-up audit)" should now cite this report: PASS WITH CORRECTIONS, with F1 applied.
- Record F1 in MISTAKE-742, or as a one-line addendum: a characterisation was typed "exactly" from a list of
  failures only.

## 2. The new clauses (3e) and (3f), checked line by line

**(3e), the word law at n = 3k + 1.**
- **Block words are tight.** Each block has diameter <= 2, so P(A) <= j and P(B) <= k - j. Subadditivity gives
  P(A) + P(B) >= P([0, L - 1]) = k. Correct.
- **Tight colourings are block words.** Induction on k.
  - k = 1: {0, 1} is monochromatic. Correct.
  - x_0 != x_1: P(A) + P(B) >= 2 + P([4, 3k - 2]) = k + 1. Correct; [4, 3k - 2] has 3k - 5 positions, so
    P = k - 1.
  - x_0 = x_1 = x_2: P(A) = 1 + P(A cap [3, ..]) and P(B) = P(B cap [3, ..]). The rest has 3(k - 1) - 1 positions.
    Correct.
  - x_0 = x_1 != x_2: P(A) = 1 + P(A cap [3, ..]) = 1 + P(A cap [2, ..]), because 2 is in B. So the colouring is
    tight iff its restriction to [2, 3k - 2] (3(k - 1) positions) is tight. Then (3a) shifted by 2 applies. Correct.
- **Count.**
  - A word determines its run lengths (one = 2 (mod 3), the rest = 0 (mod 3)), hence beta in {0, 1}^k and the run
    holding the short block. The pair (beta, run) determines the word.
  - So the count is sum over beta of runs(beta) = 2^k + (k - 1) 2^(k-1) = (k + 1) 2^(k-1), which is (k + 1) 2^(k-2)
    up to complement. Correct.
  - It agrees with W(k) = 2^(k-1) + 2 W(k-1) (AUDIT.md m12) and with brute force for k <= 7 (`hyp9540_a_blocks.out`).
- Every run of such a word has length 3m or 3m + 2 >= 2, so (3f) applies.

**(3f), the run-window theorem.**
- **Runs and P.** Runs of one colour are separated by at least two positions, so they lie at distance >= 3. Greedy
  from the left end of each run gives P >= sum ceil(l/3), and subadditivity gives <=. Correct.
- **The t side.**
  - A tile (a, a + 2) lies inside [p, q + 2] iff a is in [p, q]. The windows of distinct 1-runs are disjoint, since
    p' >= q + 3. So the triple word is w.
  - t_w is an order join of staircases of sizes l + 2 and single vertices, so nu(t_w) = sum floor((l + 2)/3) = P(A).
    Correct.
- **The u side, lower bound.** For x in B, the tile (x, x + 2) lies in no window, so {x, x + 1, x + 2} is a 3-cycle of
  u. Correct.
- **The u side, upper bound.**
  - **The points of F.** Each 0-run [p0, q0] contributes ceil(l/3) distinct points, all in [p0 + 2, q0 + 1].
    - Every B-window {x, x + 1, x + 2} with x in [p0, q0] contains one. The last point is q0 + 1 only when
      p0 + 3i_max >= q0 - 1, and then x >= q0 - 1.
    - The vertices p0 and p0 + 1 are the last two vertices of the previous 1-window. The vertices [p0 + 2, q0] lie in
      no 1-window. q0 + 1 is the first vertex of the next 1-window, or free if the 0-run is last.
    - So no point of F lies inside a 1-window except as its first vertex. Correct.
  - **Intervals of V - F.**
    - An interval of length >= 3 contains no B-window, so its positions form one stretch of 1s and the interval lies in
      that run's window. There u is transitive.
    - A window contains no free point of F. If it contains a point of F, that point is its first vertex, so the window
      minus F lies in one interval.
    - Forward arcs of u are path arcs or tiles inside windows, so every arc between different intervals points
      backwards.
    - Hence u - F is acyclic, and nu(u) <= tau(u) <= |F| = P(B). Correct.
- **Colour.** v -> n - 1 - v maps the window [p, q + 2] to [L - 1 - q, L + 1 - p], which is the window of the mirror
  run [L - 1 - q, L - 1 - p]. So t_w is blue iff the 1-runs are mirror-symmetric, i.e. iff w is a palindrome.
  Correct.
- **Computer check** (`runwindow_audit.c theorem`, `runwindow_theorem.out`).
  - Every colouring with all runs >= 2 was checked exhaustively for n = 4..26 (150,048 colourings), and 20,000 random
    ones for each n = 27..32. There were 0 failures.
  - At n = 5..13 there are 284 such colourings, the same number as `run_window_theorem_check.out`.
  - The checks cover the word, both packing numbers against the run sums and against P, |F| = P(B), the position of
    every point of F, that every 3-cycle of u meets F, and blue iff palindrome. The session's script checked neither
    the colour nor the position of F.
- **Consequences.**
  - (A) at n != 0 (mod 3): by (1), every tight line's word is a tight colouring. By (3a) and (3e), every tight
    colouring is a block word with runs >= 2. By (3f), its run-window tiling is a tight line with that word. Correct.
  - (C) at n != 0 (mod 3): a palindromic tight colouring gives a blue t_w, by (3f). Correct.
  - So HYP-9540 (A) and (C) are PROVED at n != 0 (mod 3). For n = 13 and 14 this supersedes the first-hit search.

## 3. The corrections from AUDIT.md

| item | where | status |
|---|---|---|
| M1 | HYP-9540 title, status, (B) statement, notes and falsifier; README; INDEX row | applied. The restated blue clause matches the audit's complete blue census at n = 4..12 exactly, including the tight pairs at n = 10, 11 (1,2) and n = 12 (2,2). It drops the redundant (1,1) exception and records the old clause as REFUTED. |
| M2 | HYP-9540 title, status, (D); THM-4951 (5)(e); README; INDEX rows | applied. The sign law (-1)^(c+1) matches the data: + at n = 11 and 15, - at n = 13 and 17. Its history is stated, including that U4 was pre-registered after U1-U3. |
| m1 | THM-4950 proof of (1) | applied |
| m2 | THM-4950 status and (3), HYP-9530 | applied: FINITE-EXACT; orders <= 10 suffice |
| m3 | THM-4950 (5), status, README | applied: method credited to THM-4952 (E), scope, convention, THM-4744, exact below 337. Its proof uses H_3 - 1 >= (3/2)(H - 1) and f(11) = 225, and the census shows only 10, 31, 64 below 337. (F2 and F4 are wording.) |
| m4 | THM-4950 Remarks, README | applied |
| m5 | PREREGISTRATION.md | applied. Outcome appended, predictions unedited; 4,031 irreducibles below 29,038 as in coverage10.out. |
| m6, m7 | THM-4951 status, Setting, README credits | applied |
| m8, m9 | THM-4951 proofs (2b), (3b) | applied with the suggested text |
| m10 | THM-4951 Setting (L, e, m), credit line | applied |
| m11 | THM-4951 status labels, (5) and (5d) | applied, including (P2) cited in the proof of (5d) |
| m12 | THM-4951 (3c), (3e); HYP-9540 (A) | applied and extended by (3e) and (3f) |
| m13 | HYP-9540 (C), README (2b) row | applied |
| m14 | THM-4951 (2a) | applied |
| m15 | `related` fields | applied |
| m16 | nu_lines.c | applied (comment removed) |
| m17 | both INDEX rows | applied |
| PREREGISTRATION-2 credit note | PREREGISTRATION-2.md | appended after the outcome, nothing above edited; credit correct (THM-4922 (4) audit, steps 3 and 5) |

The new THM-4951 (4) blue table at n = 12 and the (5)(e) table copy the audit's and the session's numbers correctly.
So do HYP-9540 (D)'s data: s_n and the D-densities 0.0898, 0.0274, 0.0238, 0.0092.

## 4. MISTAKE-742

- The entry is at the end of 01-canon/MISTAKES.md. MISTAKE-740 and MISTAKE-741 are mac-mini's and Session231's, as
  stated.
- It records M1 and M2 and their dates. Both had already been contradicted by the session's own results (3d), T4
  and U1.
- It records:
  - the swapped proof labels;
  - the uncredited method of THM-4950 (5) and its over-wide confirmation claim;
  - the Session26 and staircase credits;
  - the binary lane list;
  - the unsourced R1 proof in PREREGISTRATION-2.
- The repair list matches what was applied, and the lessons are apt.
- It is accurate. F1 could be added as one more item (F5).

## 5. Files added for the follow-up (audit/)

| file | content |
|---|---|
| runwindow_audit.c | `theorem` mode: (3f) exhaustive n <= 26, sampled n <= 32; `n3k` mode: the run-window tiling on all tight colourings at n = 3k |
| runwindow_theorem.out | the (3f) check |
| runwindow_n3k.out | n = 6, 9, 12, 15, 18: every failure and every word with an isolated 0 between 1-runs |
| runwindow_n3k_classify.py / .out, runwindow_n3k_runs.py / .out | the failure pattern (F1) |
