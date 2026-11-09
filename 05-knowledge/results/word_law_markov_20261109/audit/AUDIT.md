# Independent audit of opus-2026-11-09-S88 (THM-5020, THM-5021, HYP-9590, HYP-9540, MISTAKE-745)

Independent auditor, 2026-11-09. Pre-audit checkpoint 2fe6b181c5 on branch opus-2026-11-09-S88 (worktree
C:/tmp/wt-s88). Scope:
- THM-5020 (word law for tight lines at every order) and THM-5021 (tight colourings as a Markov chain);
- HYP-9590, and this session's changes to HYP-9540 (diff against a3267e1c21);
- the packet `word_law_markov_20261109/` (README, PREREGISTRATION.md, all scripts and outputs);
- MISTAKE-745.

All code in this directory was written for this audit. No session program was reused. The nu solver was first
validated against the complete census of THM-4951 (4).

## Verdict: PASS WITH CORRECTIONS

- **THM-5020 is correct.** Every proof was checked line by line, and (1)-(5) hold for all n.
  - The certificate case analysis in (2) has two cosmetic slips (m1).
  - Independent exact computation confirms (2), (3), the junction bound and the proof lemmas themselves on every
    word with n <= 22 (2,097,150 words).
  - It confirms the word law (4) on every tight colouring with n <= 36 (251,904 colourings).
- **THM-5021's counts, generating functions and junction counts are correct and proved.** There is one MAJOR finding
  (M1): the spectral headline "x^14 (x - 2)^4, a single Jordan block of size 4 at 2" is an artefact of the 18-state
  encoding. Six states never occur, and one state that does occur can never accept. The trim automaton has
  x^2 (x - 2)^3 and a Jordan block of size 3, which is what the generating functions see.
- **HYP-9590 is understated, and HYP-9540 (B) can be upgraded.**
  - HYP-9590 (C), and so the exclusion of all three families, has a short proof (m5).
  - Explicit lines prove HYP-9590 (B) for two of the three families and for part of the third. A structured family
    realises every remaining pair for n <= 26 (m6).
  - With the session's census, HYP-9540 (B) then holds for every n <= 26.
- **Prior art.**
  - No earlier statement was found of the junction lemma, the certificate argument, the complement realisation or
    the word law at n = 3k.
  - Two credits are missing (m4).
  - Session233's THM-5000 and THM-5001 are unrelated. THM-5002 is related but does not overlap.
- **Process.** The pre-registered predictions were not edited after the run, and MISTAKE-745 is accurate.

## 1. Findings

### M1 (MAJOR). THM-5021 (1): the Jordan block of size 4 and x^14 (x - 2)^4 are artefacts of the encoding

**Where.**
- THM-5021: title, statement (1) (lines 59-61), the status line on the Jordan form, Remarks line 93.
- README section 2 ("18 states", "a single Jordan block of size 4") and section 3 ("nilpotent Jordan tail").
- Both INDEX rows (hypotheses/INDEX.md line 2, results/INDEX.md line 3).

**What is true** (`aud_markov.py`, `aud_markov_states.py`; outputs `aud_markov.out`, `aud_markov_states.out`).
- **The stated numbers reproduce.** The 18 x 18 matrix built exactly as THM-5021 defines it has characteristic
  polynomial x^14 (x - 2)^4, with Jordan blocks 4 (at 2) and 3, 1^11 (at 0).
- **Six states never occur.** Only 12 of the 18 states (s_A, s_B, E) are used: (1,1,*), (1,2,0), (2,1,0) and
  (2,2,*) never occur at any position.
- **One of the six block-boundary states is dead.** Only six states occur at block boundaries (m = 0 mod 3), and one
  of them, (3, 3, 1), can never reach E = 0.
  - From (3, 3, 1) the first letter of every block is chosen. A second chosen element in the same block would give
    E = 2, so every block is monochromatic.
  - It returns to (3, 3, 1) with weight 2 (111 and 000), and E never returns to 0. An example of a word reaching it
    is 100111.
- **The dead state makes the size-4 block.** It lengthens the Jordan chain at 2 from 3 to 4:
  dim ker (T - 2)^3 = 3 and dim ker (T - 2)^4 = 4, and the top vector of the chain lives on the E = 1 states.
- **The trim automaton has a block of size 3.** Its states are the reachable and co-reachable block-boundary states
  (3,3,0), (3,1,1), (3,2,1), (1,3,1), (2,3,1), and

      T_trim = [[2,1,2,1,2], [0,1,0,1,1], [0,0,1,0,1], [0,1,1,1,0], [0,0,1,0,1]].

  It has characteristic polynomial x^2 (x - 2)^3 and a single Jordan block of size 3 at 2. It gives the same three
  generating functions (checked).
- **Size 3 is the invariant.** The maximal pole order is 3 = slack 2 + 1, which is exactly what the Remarks' reading
  needs ("each unit of slack adds one power of k ... That is the Jordan block"). A live block of size 4 would allow a
  k^3 2^k term, and none occurs.

**Fix.**
- State (1) for the trim automaton: 5 states, x^2 (x - 2)^3, one Jordan block of size 3 at 2.
- Or keep the 18-state numbers with the caveat that 6 states are unreachable, that (3, 3, 1) is dead, and that the
  size-4 block is not an invariant of the tight colourings.
- Update the title, the README and the INDEX rows. "18 states suffice" is fine as an upper bound; add that 12 occur.

### MINOR

**m1. THM-5020 (2) proof, the Claim (lines 149-161).**
- The case x(D) = x(D') is not mentioned. It is immediate, because x + 1 lies in both D and D'.
- "D contains x + 3 or (for T4 with c = x + 2) x + 2" (line 157): the alternative is vacuous.
  - In T4, (b, c) is a window tile, so c >= b + 2 = x + 3.
  - The window is [x + 1, x + 3], so c = x + 3 for T2 and T4 alike, and D always contains x + 3 = x' + 1, which lies
    in D'.
- Delete the alternative and its sub-bullet ("In the second, ...").
- The argument is otherwise complete. Every case was re-derived (section 2). The classification and the separation
  were also checked exhaustively for n <= 22:
  - 71,303,470 3-cycles of R(w) + 1 classified, and 742,385,889 disjoint pairs checked, with 0 violations;
  - a deliberately wrong T2 certificate is caught (`AUDIT_MUTATE=1`), so the check is not vacuous.

**m2. THM-5020 (4), line 83: "At n != 0 (mod 3) it reproves THM-4951 (3a)-(3f)."**
- (2)-(4) reprove the realisation half there, (3b), (3f) and the word laws.
- They do not reprove:
  - the block characterisations (3a) and (3e);
  - the counts (3c), which THM-5021 (2) reproves;
  - (3d)'s list of blue pairs.
- Write "at n != 0 (mod 3) it gives a second proof of the word laws (THM-4951 (3f))".

**m3. Symmetry remarks (THM-5020 Remarks lines 202-204; README section 3 lines 57-61).**
- **Complementation.** "The construction commutes with sigma and with complementation" is true for sigma and false
  for complementation.
  - sigma(R(w)) = R(reverse w) holds for every word with n <= 16.
  - R(complement w) = R(w) + 1 holds only for the two constant words (`aud_symmetry.out`).
  - Failing to commute with complementation is exactly the one asymmetric choice.
- **Blue lines and the choice.** "So blue lines never need an asymmetric choice" is false.
  - The palindrome 1110111 at n = 9 has a central 0-junction, and only R(0001000) realises it.
  - Every palindromic tight colouring with a central 0-junction needs the complement construction: 1, 2, 4, 8, 16 of
    them at n = 9, 15, 21, 27, 33. The line is still sigma-symmetric.
- **Defects.** README: "Blue lines are therefore the lines whose single possible defect sits on the mirror" and
  "Black lines are the same construction with the defect off-centre". Restrict both sentences to tight lines and to
  junctions.
  - With defects read as short blocks (THM-5021's own reading), a palindromic tight colouring can carry two
    symmetric defects off the mirror: 1100000011 at n = 12, with blocks 2 + 3 + 3 + 2 and pair (2, 2).
  - A black tight line need not have any defect.

**m4. Credits.**
- **Central-pair lemma (THM-5020 (5), NEW list).**
  - Its n = 3k + 2 case is THM-4951 (3d): a palindromic block word of even length has an even number of 1s.
  - The pattern "the blue tight pairs are the {j, k - j} with j even or k odd" was FINITE-EXACT for n <= 26 in the
    S87 audit (`ternary_metagraph_20261109/audit/palindromic_pairs.out`; HYP-9540 (B) Notes).
  - S87's audit also proved the n = 12 case, by listing the three palindromic tight colourings. README line 36
    ("which S87 had only observed") understates this.
  - The general lemma is new. Credit both.
- **THM-5021's automaton.**
  - It is the "greedy-state DP" of S87's `staircase_check.py`, part (d): the gaps since the last greedy pick in A
    and in B, plus P_A + P_B. That DP produced THM-4951 (3c)'s fitted count.
  - New in THM-5021: the prefix bound E <= 1 that makes it finite, the exact generating functions, and the proof of
    the count. Credit the S87 DP in KNOWN.
- **THM-5002 (Session233).** It is a related transfer: it covers the complementary run-window endpoint B(m) for runs
  >= 2. THM-5020 lists it under `related` without saying how it relates. Add one sentence that there is no overlap:
  - THM-5002 counts odd-cycle packings of B(m) when all runs have length >= 2;
  - THM-5020's new content is the case of runs of length 1;
  - THM-5021 counts words.

**m5. All three families are excluded by proof; HYP-9590 (C) is PROVED.**

Where: HYP-9590 status lines 6-13 and section (C); THM-5020 (6) lines 99-100; MISTAKE-745 ("Both gaps are
structural").

Proof that no palindrome has P(A) = P(B) = k at n = 3k + 1 (L = 3k - 1):
- A k-element 3-separated subset of [0, 3k - 2] has total slack 1. So it is
  S_j = {3(i - 1) + [i > j] : 1 <= i <= k} for some 0 <= j <= k.
- Suppose P(A) = P(B) = k, with S_j in A and S_j' in B.
  - Their least elements must differ, so exactly one of j, j' is 0, say j' = 0. Then S_0 = {1, 4, ..., 3k - 2}.
  - Disjointness forces j = k, so S_k = {0, 3, ..., 3k - 3}.
- Then w_0 != w_(3k-2), which is impossible for a palindrome.

At n = 3k the session's argument (the unique maximum packing) already excludes every word. So no palindrome has
(k, k) at n = 0, 1 (mod 3). This was checked for all palindromes of length <= 32 (`central.out`).

Fix:
- Upgrade HYP-9590 (C) to PROVED, and the "minus" part of (A) to PROVED for all three families.
- Add the argument to THM-5020 (6) and to MISTAKE-745.

Observation: for n <= 28 the missing pairs are missing for all words, admissible or not (`pairsx_*.out`). So in this
range admissibility costs nothing: the three families are exactly the pairs that are the P-pair of no word (black)
or of no palindrome (blue).

**m6. HYP-9590 (B) is largely provable; HYP-9540 (B) then holds for every n <= 26.**
- **Nested mirror tiles (families 2 and part of 3; `aud_nested.c`).** Let t_m reverse exactly the sigma-fixed tiles
  (a, n - 1 - a) for 0 <= a < m <= k, with n >= 6.
  - In t_m the backward arcs are the n - 1 - a -> a. They are pairwise vertex-disjoint, and every cycle uses one, so
    nu(t_m) <= m.
  - The triangles {a, m + a, n - 1 - a} are disjoint 3-cycles (since 3m <= n), so nu(t_m) = m.
  - In t_m + 1 the only forward tiles are the m mirror tiles, and none of them has length 2 when n >= 6 and m <= k.
    So every consecutive window is a 3-cycle, and nu(t_m + 1) = k.
  - So the blue pair (m, k) occurs for every 0 <= m <= k:
    - m = k gives family 2, blue (k, k), at every n >= 6;
    - odd m gives the members (i, k) of family 3 (even n, k odd).
  - m = 1 is a single sigma-fixed tile, (0, n - 1), giving blue (1, k) at every n >= 4 (`aud_onek.out`, n = 4..40).
    Until now HYP-9540 had (1, k) only in black, through the tile (0, 3).
- **Black (k, k) at n = 3k >= 6 (family 1).** Take t_k together with the tile (0, 3).
  - The packings {a, x_a, n - 1 - a} (x_a in [k, 2k - 1], x_0 != 3) and {3b, 3b + 1, 3b + 2} both avoid the arc
    between 0 and 3. So both survive, and nu = (k, k).
  - (0, 3) != sigma(0, 3), so the line is black.
- **Checks.** Exact nu confirms both constructions for n = 6..33 (`aud_nested.out`; the run was stopped during
  n = 34, where the exact confirmation becomes slow).
- **The rest of family 3 (both-odd (i, j) with j < k at even n).** All of these members are realised for every
  n <= 26 by a structured family (`aud_mirror.c`; `mirror_10_20.out`, `mirror_21_26.out`):
  - take R(w) for a palindrome w and flip at most three mirror tiles;
  - typical witnesses are w = 1^r 0^(L - 2r) 1^r with nested mirror tiles just inside the end windows;
  - every witness was re-verified with exact nu and colour.
- **A random search (`aud_residual.c`, stopped at n = 22) independently found witnesses for most members at
  n <= 19.** They include black (4, 4) at n = 12, for which no census existed.
- **Consequences.**
  - HYP-9590 (B) is PROVED for families 1 and 2 and for the j = k members of family 3, and FINITE-EXACT for
    n <= 26.
  - With the session's admissible-word census (black pairs at n = 10-12 included, `pairs_check_insample.out`) and
    S87's complete censuses, HYP-9540 (B) holds for every n <= 26.
  - THM-5020 (6)'s "There is no complete black census at n >= 10" can now point to these explicit lines.

**m7. THM-5021 (3), lines 71-74.**
- "317 reachable states" is the count when any phase is allowed at any step (markov_junctions.py's search). With
  the phase tracked, 103 states occur. Correct the number, or say how it was counted.
- "which is THM-5020 (1)'s 'J <= 1' seen by the automaton": the flags record only whether 11011 and 00100 occur.
  The automaton therefore sees "no tight colouring has both kinds", not J <= 1: two 0-junctions are invisible to
  it.

**m8. THM-5021 proof, line 79.**
- The first generating function written, v0 (I - xT)^(-1) T fin, equals 2/(1 - 2x) = sum N(3k) x^(k-1).
- Statement (2) uses x^k, that is 2x/(1 - 2x), and `markov_tight.out` prints -2/(2x - 1) under "generating
  function, L = 3k".
- Add the factor x, as in the other two, or note the shift.

**m9. "Markov chain" and the Ulam/Tao comparison (THM-5021 title, Remarks lines 98-103; README lines 92-93,
"just as").**
- T is a counting transfer matrix with spectral radius 2, not a stochastic matrix. It becomes a Markov chain only
  after normalisation: uniform letters, T/8 per block, killed at E = 2. Say so once.
- "Only nonzero eigenvalue 2, just as the Ulam and Tao chains are rank one plus nilpotent": here the part at 2 is a
  Jordan block, 3-dimensional in the trim automaton (4 in the 18-state matrix), not a rank-one projection.
  T = 2 Pi + N with Pi of rank 3. Keep this as ANALOGY, with the difference stated.

**m10. README section 3, Euler-product bullet (line 81).**
- "The product over cycles of (1 + lambda^|C|) equals I(Omega, lambda) only if cycles never interact" has the
  wrong exponent.
- For pairwise disjoint cycles I(Omega, lambda) = Prod_C (1 + lambda), because I counts packings by their number of
  cycles. With lambda^|C| the identity fails even for a single 3-cycle.
- The slip is inherited from 00-navigation/CONCEPT-MAP.md (line 992), whose next line uses Prod_C (1 + 2). Fix the
  README; the CONCEPT-MAP line deserves the same fix.

**m11. "No feedback set".**
- README line 29 ("No feedback set can work here"), THM-5020 Remarks line 205 ("a feedback set may not exist") and
  the hypotheses INDEX row ("no feedback set exists in general") should all say "no feedback set of size P(B)".
- For w = 010 the claims nu = 1 < tau = 2 and "three pairwise-intersecting 3-cycles with no common vertex" are
  right. The 3-cycles are {0,1,2}, {0,1,3}, {1,3,4} and {2,3,4}.

**m12. Stale cross-references.**
- **THM-4951** has no pointer to THM-5020 or THM-5021. It still calls the n = 3k count fitted (status list
  lines 16-17, (3c) lines 145-146) and the word law at n = 3k NOT CLAIMED (line 47). Add a forward pointer.
- **HYP-9540 status lines 15-16.** "Through (C), at n != 0 (mod 3) exactly those occur" now holds at every n.
- **HYP-9540 (B) Notes.** The exclusion of both-odd blue tight pairs is now PROVED at even n (THM-5020 (5)).
  - At odd n it is automatic, except at n = 3k + 1 with k even.
  - There the central run of a palindromic tight colouring has length 5 (mod 6), so both packing numbers are even.
- **HYP-9540 (A) body.**
  - "There are (k - 2) 2^(k-3) such words" counts both first letters. S87's failure counts 1, 3, 8, 20 are those
    with x_0 = 1, which number (k - 1) 2^(k-4) (confirmed here to k = 12). State the convention.
  - "That every 0-junction word fails is FINITE-EXACT only, for n <= 18": `tight_3_36.out` extends this to n <= 36.
    R(w) fails on all 12, 32, 80, 192, 448, 1024, 2304 and 5120 such colourings at n = 15, ..., 36.

**m13. The pasted text is not archived.** The pasted text on self-avoiding walks, the honeycomb lattice and Collatz
is not in the packet, so "The pasted text is wrong on this point" cannot be checked against it. Quote the specific
pasted claims in the README, or archive the paste. The corrections themselves are right (section 5).

**m14. Labels.** PREREGISTRATION.md (Consequences) and the HYP-9540 status say (B) is "PROVED for 13 <= n <= 40" on
the realised pairs. THM-5020's status labels the same census FINITE-EXACT. Use one label: a FINITE-EXACT census plus
the proved lemmas (2) and (3).

**m15. HYP-9590 ranges.**
- The gap law (exactly the three families) holds for every n = 6..40, the in-sample range included
  (`compare_pairs.out`). So "n >= 13" can read "n >= 6".
- (B)'s evidence list omits members found in complete censuses: blue (2, 2) at n = 6 and 7, blue (3, 3) at n = 9,
  and blue (1, 3) and (3, 3) at n = 10. "At every n >= 9" could read "n >= 4". m6 supersedes most of this.

## 2. Proof review, line by line

**THM-5020 (1), the junction lemma.**
- **Deletion step.** The shifted images of optimal packings of w' are 3-separated, because distances only grow, and
  they keep their letters. No B-position lies within distance 2 of z, since z +- 1 and z +- 2 are 1-positions, so z
  joins the B side. Correct.
- **Survival of the other junctions.**
  - A 1-junction's isolated 1 and its two 0-runs (each of length >= 2) are disjoint from R1, {z}, R2.
  - A 0-junction that shares R1 or R2 as a neighbouring run sees that run grow to |R1| + |R2| >= 4.
  - So every other junction survives, and none is created. Confirmed: J(w') = J - 1 for every junction of every
    word with L <= 20.
- **Induction.** The induction and the arithmetic ceil((L - J)/3) + J = ceil((L + 2J)/3) are correct.
- **Consequences.** 2J <= 3k - L; J = 0 at n != 0 (mod 3) and J <= 1 at n = 3k; a palindromic tight colouring has
  its junction at the centre. All correct.
- **Exhaustive check.** The bound holds for all 268,435,454 words with L <= 27 (`junction.out`), and it is attained
  by many words, so it is sharp in general.

**THM-5020 (2).**
- **Arcs.** The arc lists of t and u are correct. Two windows share at most one vertex, so each tile lies in at most
  one window, and the triple word of R(w) is w.
- **No reversed 3-cycles.** The orientation a -> c -> b -> a is impossible in u: a -> c forward makes (a, c) a
  window tile, and then b -> c.
- **The types are right.** T1-T4 are exhaustive and mutually exclusive (b = a + 1 or not, c = b + 1 or not), and
  their descriptions are correct:
  - T2 needs p >= 1 and c in [p + 2, q + 2];
  - T3 needs a in [p, q] and q + 3 <= n - 1, so q + 1 <= L - 1;
  - T4: W'' ends at b and W' starts at b, with a <= b - 2 and c >= b + 2.
- **Certificates.** They lie in B, and x(D) + 1 lies in D.
- **The Claim.** Every case of the Claim was re-derived:
  - x' = x + 1 excludes T2 and T4 for D, and T1 and T3 both contain x + 2 = x' + 1;
  - for x' = x + 2 with D of type T1 or T3, D contains x + 1 and x + 2. A D' of type T1 or T2 contains x + 2. A D' of
    type T3 or T4 has the run {x + 1} (since x is in B) and contains its member x + 1;
  - for x' = x + 2 with D of type T2 or T4, D contains x + 3 = x' + 1.
  - The only gaps are those of m1.
- **Lower bound.** Correct.
- **Edge cases.** Runs of length 1 (chains of shared vertices, as in 1010101), runs at the ends (p = 0 has no T2;
  q = L - 1 has no T3) and n = 3, 4, 5 are all covered by the exhaustive run.

**THM-5020 (3).**
- **Classification.** The 3-cycles of t are exactly the window triangles and the shared-vertex triangles. Every
  shared-vertex triple with a in W'' - {b} and c in W' - {b} is a 3-cycle.
- **The feedback set meets every 3-cycle.** Each a in A lies in [g_i, g_i + 2], by the greedy rule. With no
  0-junction, one run at the isolated 0 has length 1.
  - If the run before is {b - 2}, then g in {b - 4, b - 2}, so g + 2 is in {a, b}.
  - If the run after is {b}, then g in {b - 2, b}, so g + 2 is in {b, c}.
- **Conclusion.** t - F is acyclic. Correct, and confirmed exhaustively for n <= 22. On words with a 0-junction,
  R(w) does fail sometimes (1110111 is the first, at n = 9; nothing fails at n <= 8), so the hypothesis cannot be
  dropped.

**THM-5020 (4).**
- **Realisation.** Correct, including the complement realisation through R(w-bar) + 1.
- **Colour.** The mirror of the window of a run is the window of the mirrored run, so R(w) is sigma-fixed iff w is a
  palindrome.
- **Converse.** It follows from THM-4951 (1).
- **Convention.** "Blue iff t = sigma(t)" is unambiguous, because sigma(t) = t + 1 is impossible: the tile (0, n - 1)
  is sigma-fixed.

**THM-5020 (5).** Both bounds are correct: 2 P(X ∩ [0, h - 2]) <= P(X) <= 2 P(X ∩ [0, h - 1]). Parity then forces
h - 1, and by symmetry h, into X. This needs L >= 2, which even L supplies.

**THM-5020 (6).**
- **Proved.** The lines of admissible words come from (2) and (3). (k, k) is impossible for any word at n = 3k (the
  unique maximum packing; the argument is correct). The both-odd family follows from (5).
- **Census.** FINITE-EXACT, reproduced exactly (section 3).
- **Second family.** Left open in THM-5020, but provable (m5).

**THM-5021.**
- **Automaton.** The automaton is correct: the greedy rule is optimal per colour, and E changes by
  chosen - [m = 0 (mod 3)].
- **Prefix bound.** The proof of the prefix bound is correct. E >= 0 on prefixes is used implicitly; it is Session26
  applied to the prefix.
- **Upper bound.** "18 states suffice" is a valid upper bound (but see M1).
- **Generating functions.** The three generating functions and the closed forms are correct, including the
  partial-fraction identity 2[C(k+1,2) 2^(k-1) - C(k,2) 2^(k-2) - C(k-1,2) 2^(k-3)] = (k^2 + 9k - 2) 2^(k-3).
- **Junctions.** The junction generating function x^3/(1 - 2x)^2 is correct.

## 3. What was independently verified, and how

| check | method | range | result | file |
|---|---|---|---|---|
| nu solver | all lines {t, t+1} by unordered nu-pair and colour | n = 5..8 | reproduces THM-4951 (4) exactly (e.g. n = 8: 12 + 0, 5024 + 17, 1041492 + 2030) | `aud_words.c validate` |
| THM-5020 (2), (3), junction bound, deletion step, T1-T4, certificates, F meets every 3-cycle of t | exact nu (own branch and bound), packing by the DP f(i+1) = max(f(i), [i in X](1 + f(i-2))), not greedy | every word, n = 3..22 (2,097,150 words) | 0 failures; mutation test bites | `words_all.out` |
| junction bound | exhaustive | every word, L <= 27 | 0 violations | `junction.out` |
| word law and blue word law (4) | pruned DFS over tight colourings; realisation by R(w) or R(w-bar); exact nu of both ends, endpoint word, colour | every tight colouring, n = 3..36 (251,904) | 0 failures; J <= 1 always; R(w) alone fails on every tight 0-junction colouring | `tight_3_36.out` |
| central-pair lemma | every palindrome of even length | L <= 32 | 0 violations | `central.out` |
| tight-colouring counts | brute force over all words (L <= 27); DFS (n <= 36); own automaton on DP difference states (k <= 40) | | 2^k, (k + 1) 2^(k-1), (k^2 + 9k - 2) 2^(k-3) | `junction.out`, `tight_3_36.out`, `aud_markov.out` |
| generating functions | sympy, own automaton; 18-state and trim | exact | the three claimed generating functions | `aud_markov.out`, `aud_markov_states.out` |
| THM-5021 (1) | the 18-state matrix as defined; reachability and co-reachability; trim matrix | exact | x^14 (x - 2)^4 reproduced; trim x^2 (x - 2)^3, block of size 3 (M1) | `aud_markov_states.out` |
| junction counts (k - 2) 2^(k-3) | own junction automaton (k <= 30); brute force (k <= 9); DFS (k <= 12) | | confirmed; both kinds never together; none at n != 0 (mod 3) (L <= 60) | `aud_markov.out`, `junction.out`, `tight_3_36.out` |
| pairs census (THM-5020 (6), Q1-Q2) | own DP on (difference states, counts, last four letters, flags), palindromes enumerated | n = 6..40 | identical to pairs_check.out and pairs_check_insample.out (missing lists and admissible counts) | `pairsdp_6_40.out`, `compare_pairs.out` |
| the same, exhaustive, and for all words | full enumeration | n = 6..28 (spot checks at 13, 16, 18, 21 included) | identical; missing pairs are missing for all words | `pairsx_6_21.out`, `pairsx_22_28.out` |
| gap law (exactly three families) | from the census | n = 6..40 | holds at every n | `compare_pairs.out` |
| HYP-9590 (C) | every palindrome | L <= 32 | no (k, k) at n = 0, 1 (mod 3); proof in m5 | `central.out` |
| residual pairs (HYP-9590 (B)) | explicit constructions and structured / random search, exact nu | n <= 33 (constructions; m = 1 to n = 40), n <= 26 (all members) | all realised | `aud_nested.out`, `aud_onek.out`, `mirror_*.out`, `residual_9_24.out` |
| commutation with sigma / complementation | every word | n <= 16 | sigma yes; complementation only for constant words | `aud_symmetry.out` |

## 4. Prior art and concurrent work

**Searched.**
- Scope: git grep over 00-navigation, 01-canon, 04-computation, 05-knowledge and 07-reflections.
- Terms: 11011, 00100, 1110111; junction; isolated letter or isolated 0; central pair or central positions;
  palindrome together with both-odd parity; 3-separated; tight colouring; run-window; triple word; transfer matrix,
  automaton or generating function in this context; 9k - 2; Narayana.

**Found.**
- **The S87 packet (THM-4951, its audit and FOLLOWUP.md).**
  - The n = 3k failure pattern is FINITE-EXACT for n <= 18. There is no proof, and no complement construction. The
    pre-session HYP-9540 (lines 56-57 at a3267e1c21) says "a proof at n = 3k needs a construction for the 1110111
    pattern", and the THM-4951 Remarks say such words "need other constructions".
  - The palindromic parity pattern is FINITE-EXACT for n <= 26.
  - The greedy-state DP is there.
  - Credits for the last two: m4.
- **Session26** (session26-blue/RESULT.md) has the window counts nu_c only as lower certificates. THM-5020 credits
  it.
- No earlier junction lemma, certificate argument for nu(R(w) + 1), word law at n = 3k, central-pair lemma in
  general, or finite-state reduction with generating functions for tight colourings was found.
- **07-reflections/the-slope-ladder-of-pascal-and-the-growth-constants-kps.md** has the Narayana / hard-core
  transfer for 3-separated subsets, as THM-5021 cites.

**Concurrent work.** `git fetch`, then `origin/main` since 2026-11-09 10:00. Besides this session's commits there
are only:
- mac-mini-2026-11-09-paleyeven: Paley doubles and dominating sets, unrelated;
- Session233's final checkpoint (10:07).

**Session233.**
- THM-5000 (coprime interface composition) and THM-5001 (endpoint spectrum) are unrelated.
- THM-5002 is related but does not overlap, for these reasons:
  - it gives a three-state transfer for the odd-cycle packing polynomial of bridge-linked transitive blocks;
  - it contains the complementary endpoint of run-window lines with all runs >= 2;
  - it obstructs H = 11 for that family;
  - it says nothing about packing numbers of words with runs of length 1, about junctions, or about counting
    tight colourings.

## 5. Wording, status and process

**Labels.** PROVED, FINITE-EXACT and ANALOGY are used accurately, except for M1, m2, m3 and m14. The README's lens
section is typed ANALOGY. Its factual claims were checked:
- **Duminil-Copin and Smirnov**, Ann. of Math. 175(3) (2012) 1653-1665. The parafermionic observable satisfies half
  of the discrete Cauchy-Riemann relations, and convergence to SLE(8/3) is conjectural. The README is correct.
- **Tao**, arXiv:1909.03562, Forum Math. Pi 10 (2022) e12. The theorem is "almost all orbits attain almost bounded
  values" in logarithmic density. The method uses Syracuse random variables and a two-dimensional renewal process,
  and it does not prove the conjecture. The README is correct.
  - The repository's earlier refutation of a pasted "dyadic martingale renewal process / entropy decrement"
    description (`collatz_mod6_20260922_w6_tao_minus_sheet_ladder.md`, S2) is correctly cited.
- **V. R. R. Uppuluri**, "Prime Numbers Treated as a Markov Chain", Oak Ridge National Laboratory (1968; NTRL
  PB180746, 12 pp.). It embeds an ergodic Markov chain in the last digits of successive primes. It exists, and the
  README describes it correctly.
- **Lemke Oliver and Soundararajan**, PNAS 113(31) (2016) E4446-E4454. Correct.
- **"Marty's Primes as a Markov Chain".** My searches also found nothing by anyone named Marty. One lead for the
  owner, not an identification: arXiv:2605.00301 (Alexeev, Barreto, Li, Lichtman, Price, Shah, Tang, Tao, 2026),
  "Primitive sets and von Mangoldt chains", builds Markov chains with von Mangoldt weights on the divisibility
  structure.
- **Session140** (`session140_collatz_hexagonal_geometry_20261102`): a faithful directed honeycomb subdivision, with
  unbounded routing. Correctly summarised.
- **The other references exist and are cited correctly:** the quasirandom note
  (`quasirandom_arithmetic_structure_20261029/PROOFS.md`), THM-4940 (1) (multiplicative over order joins, that is,
  over strong components), THM-4473, THM-4821/4822, and the CONCEPT-MAP line (its formula is the slip of m10).
- **Other remarks.** "Literally such a rule" for the greedy packing is analogy wording. The hard-core
  identification is correct.

**Pre-registration.**
- PREREGISTRATION.md was created in 654969dcc5 (10:22:05). The outcome was appended in 4cd63082f1 (10:24:52).
- The only later change is MISTAKE-744 -> MISTAKE-745 in the outcome paragraph (ad3ed2bb81).
- The predictions (Q1-Q3) were not edited. The in-sample run postdates the push, as disclosed.
- The outcome's statements were reproduced exactly (section 3).

**MISTAKE-745.** Accurate: Q1 and Q2 were bound to fail for structural reasons, and the repair is as described. Two
additions are suggested:
- the third family, blue (k, k) at n = 1 (mod 3), is structural too (m5);
- all three families are realised by other lines (m6).

**Working-tree note.** Two untracked files appeared in the packet at 10:35, after the audited checkpoint (10:34) and
before this audit began: `palindrome_kk_check.py` and `palindrome_kk_check.out`.
- They are not part of the audited commit and were not modified here.
- They contain the same argument as m5 (the k-packings S_j at n = 3k + 1 and the mirror swap of S_0 and S_k), so the
  session reached it independently. Its output agrees with `central.out`.
- The docstring calls the script "added after the audit", although it predates this audit. State the provenance
  accurately when it is committed.

## 6. Files (audit/)

| file | content |
|---|---|
| `aud_words.c` | nu solver; `validate` (all lines, n <= 8); `words` (THM-5020 (1)-(3) and the proof lemmas, every word); `tight` (word law, every tight colouring); `AUDIT_MUTATE=1` mutation test |
| `words_all.out`, `tight_3_36.out` | outputs (every word, n = 3..22; every tight colouring, n = 3..36) |
| `aud_pairs.c` | junction bound and counts (all words), central-pair lemma and HYP-9590 (C) (all palindromes), pair census by full enumeration (`pairsx`) and by DP (`pairsdp`) |
| `junction.out`, `central.out`, `pairsx_6_21.out`, `pairsx_22_28.out`, `pairsdp_6_40.out` | outputs |
| `compare_pairs.py`, `compare_pairs.out` | auditor census against pairs_check.out; the gap law |
| `aud_markov.py`, `aud_markov.out` | own automaton (DP difference states), brute-force agreement, 18-state spectrum, generating functions, closed forms, junction automaton (reachable states with and without the phase) |
| `aud_markov_states.py`, `aud_markov_states.out` | used, unreachable and dead states; trim automaton and its spectrum and generating functions |
| `aud_symmetry.py`, `aud_symmetry.out` | R versus sigma and complementation |
| `aud_onek.c`, `aud_onek.out` | blue (1, k) by the single tile (0, n - 1), n = 4..40 |
| `aud_nested.c`, `aud_nested.out` | nested mirror lines (m, k) and black (k, k) at n = 3k, exact nu |
| `aud_mirror.c`, `mirror_10_20.out`, `mirror_21_26.out` | structured family R(w) XOR mirror tiles: every residual blue pair, n <= 26 |
| `aud_residual.c`, `residual_9_24.out` | random / greedy search for residual pairs (stopped during n = 22) |

Build: `gcc -O2 -march=native -o X.exe X.c`. Python 3 with sympy 1.14.
