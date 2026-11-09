# Audit of THM-4950, THM-4951, HYP-9540 and the HYP-9530 update (opus-2026-11-09-S87)

Independent auditor, 2026-11-09. Scope: the checkpoint adfb2402b7 plus a5a3b2793e (THM-4950 (5)), i.e.
THM-4950, THM-4951, HYP-9540, the HYP-9530 diff against 9c48b9ac1e, and the packet
`05-knowledge/results/ternary_metagraph_20261109/` (README, PREREGISTRATION, -2, -3, all sources and outputs).
The later commits de27f055b2, c2ac2a65c5, 8d36acba2e (PREREGISTRATION-4, U1-U4) were read only where they bear on
HYP-9540 (D). Every number below was recomputed with code written for this audit (this directory); no packet program
was reused.

## Verdict: PASS WITH CORRECTIONS

- Every number in THM-4950 and THM-4951 was reproduced independently, most of them by a method different from the
  session's (section 3). No numerical error was found.
- The PROVED clauses are correct. Three proofs have small gaps that need one or two sentences each (m1, m8, m9).
- Two MAJOR findings, both in HYP-9540:
  - M1: the blue clause of (B) is false. It is refuted at n = 12 by THM-4951's own results, and the audit's complete
    n = 12 blue census confirms it.
  - M2: (D) was refuted by the session's own later pre-registered test U1 and needs a status change. The audit's
    sampler confirms the refutation independently.
- Prior art and same-day work need more credit in several places (m3, m4, m6, m7, m11), but nothing NEW turns into
  KNOWN. THM-4950 (5)'s method is THM-4952 (E)'s. THM-4951 (1)'s independence corollary is Session26's.

## 1. Major findings

### M1. HYP-9540 (B), blue clause: false at n = 12, 13, 14 and at every n = 3k + 2 with k even >= 4

**Location.**
- HYP-9540 title: "(B) ... blue lines join the same pairs except (1, 1) from n = 6 on, plus the transitive line".
- HYP-9540 status line "(B): FINITE-EXACT for n <= 9 ... The boundary i + j = k is realised at every n
  (THM-4951 (2))". This holds for black lines only.
- HYP-9540 section (B), "Statement", second bullet.
- README section "Conjectures (HYP-9540)", (B) bullet.

**Why it is false.**
- (B) asserts a blue line joining nu-classes 1 and 3 at n = 12: i = 1 >= 1, i + j = 4 >= k = 4, and
  (1, 3) != (1, 1).
- Any line with nu-pair (1, 3) at n = 12 is tight, since 1 + 3 = k. By THM-4951 (1), its triple word is a tight
  colouring with {P(A), P(B)} = {1, 3}. A blue line's word is a palindrome, because sigma reverses the triple word.
- The palindromic tight colourings of [0, 9] (x_0 = 1) are exactly three:
  - 1100000011, with (P(A), P(B)) = (2, 2);
  - 1110000111, with (2, 2);
  - 1111111111, with (4, 0).

  THM-4951 (4) itself counts these 3 words. So no blue line has nu-pair (1, 3) at n = 12. This is a proof, not a
  census.
- The packet already contains the counterexample:
  - THM-4951 (4) and PREREGISTRATION-3 (T4) state that there are no blue (1, 3) tight lines at n = 12.
  - THM-4951 (3d) proves the same at n = 14 ("Example: at n = 14 there is no blue tight line with pair (1, 3)").

  HYP-9540 (B)'s own falsifier, "An allowed pair with no line at some n >= 10", is therefore met.
- **Complete census (audit).** `audit_lines.c` in mode `bluenu` enumerates all 2^29 = 536,870,912 blue lines at
  n = 12.
  Their nu-pairs are (0,4) 1, (1,4) 67, (2,2) 296, (2,3) 5,513, (2,4) 129,241, (3,3) 22,586, (3,4) 2,313,819 and
  (4,4) 534,399,389. The pair (1, 3) does not occur.
- **The pattern (`palindromic_pairs.py`).** For every n <= 26, the nu-pairs of palindromic tight colourings, which
  contain all blue tight pairs, are exactly the {j, k - j} with j even or k odd. This is the rule THM-4951 (3d) proves
  at n = 3k + 2. So (B)'s blue clause fails on the boundary at n = 12-14, 18-20 and 24-26, and by (3d) at every
  n = 3k + 2 with k even >= 4.

**Fix.**
- Restate the blue clause:

  > some BLUE line joins i <= j iff (i, j) = (0, k); or i >= 1 and i + j > k; or i + j = k and (i, j) is the pair
  > (P(A), P(B)) of a palindromic tight colouring (by (C); for n <= 26 these are the pairs with i even or k odd).
- The separate "(1, 1) from n = 6 on" exception is then redundant:
  - at n = 6, 7, 8, (1, 1) is a tight pair with k = 2 and i = 1 odd, so the parity rule already excludes it;
  - for n >= 9, i + j = 2 < k.

  THM-4951 (1)'s blue (1, 1) exclusion is the k = 2 case of the parity rule.
- Record the original clause as **REFUTED at n = 12 (proved)**.
- The audit's census supports the restated clause at n = 10, 11 and 12. It finds every non-tight allowed pair
  i >= 1, i + j > k in blue (`audit_lines.out`, `audit_blue12.out`).
- Change the status sentence to "the boundary is realised in black at every n (THM-4951 (2a)); in blue exactly the
  palindromic pairs are possible".
- Make the same change in the README.

### M2. HYP-9540 (D) is refuted, and the canon text still asserts it

**Location.**
- HYP-9540 title, "(D) the blue-line (1,3) excess of THM-4951 (5) persists at every odd n >= 11".
- HYP-9540 status "(D): one data point ... OPEN".
- HYP-9540 section (D).
- README "Conjectures (HYP-9540)", (D) bullet.
- 05-knowledge/hypotheses/INDEX.md line 2 ("(1,3) excess at odd n >= 11").

**The refutation.**
- The session's own pre-registered U1 (PREREGISTRATION-4, committed c2ac2a65c5 after the checkpoint) is REFUTED: at
  n = 13 the (1,3)-share is 1/2 - 0.000906 (24,000,000 samples, z = -53.6).
- The audit's sampler (`blue_sample.c`) agrees: 600,000 samples give a share of 1/2 - 0.00109 (se 0.00011,
  z = -10.1) and |D| / 2^30 ~ 0.0275. The session found 0.02744 and 1/2 - 0.000906; the two estimates are 1.7 se
  apart. The sampler uses only odd-cycle parities and reads e(y) directly, so it does not use (P2) or the mirror-path
  code.

**Fix.**
- Mark (D) REFUTED at n = 13.
- If the session wants to keep a conjecture here, the sign law (-1)^(c+1) (n = 2c + 1) may be stated as a NEW
  conjecture. Its data are n = 11 exact and n = 13, 15, 17 sampled; U4 at n = 17 was pre-registered and CONFIRMED.
- Remove "persists" from the title, README and INDEX.

## 2. Minor findings

**m1. THM-4950, proof of (1): the directions are swapped.**
- "Necessity is immediate, since every irreducible of M3 lies in M3" does not show that the irreducible is H_3 of a
  STRONG tournament. The strong-component argument placed under "sufficiency" is what proves necessity.
- Sufficiency needs only that M3 is generated by its irreducibles and that Spec_3 is multiplicative.
- Rewrite:

  > If Spec_3 = M3, an irreducible v equals H_3(T); v is the product of H_3 over the strong components (THM-4940 (1)),
  > each factor lies in Spec_3, which is contained in M3, and irreducibility leaves one factor != 1: v is H_3 of a
  > strong component. Conversely, every element of M3 is a product of irreducibles (induction on size), and Spec_3 is
  > closed under products.

  The mathematics is all present.

**m2. THM-4950 status: (3) is labelled "PROVED constructively".**
- THM-4940 (9), the same kind of statement, is labelled "FINITE-EXACT, by explicit witnesses".
- Use FINITE-EXACT (by explicit witnesses; pre-registered Q1-Q3 CONFIRMED).
- Add the audit's observation (`closure_check.out`, step 1): the strong classes of order <= 10 **alone** already cover
  every v = 1 (mod 3) up to 29,035 except 10 and 31. S86's pieces (ii) and (iii) contribute nothing below 29,038.
  They do not change m* (9679 with or without them). The witness list in (3) and in HYP-9530 can say so.

**m3. THM-4950 (5): credit, scope and convention.**
- (a) **Method.** The argument "H = I(Omega, 2) <= I(Omega, 3) = H_3, so Busch's floor bounds the order" is
  THM-4952 (E)'s (mac-mini, commit 844706ab1c at 23:00, before a5a3b2793e at 23:14). There it is applied at
  order <= 8 to exclude 64. Credit the method, not only the 64.
- (b) **Scope.** (5) confirms only the strong-spectrum sentence of HYP-9541 (4) ("64 is the only value 1 (mod 3)
  attained by a tournament but by no strong tournament"). It does not address the ribbon-completeness sentence.
- (c) **Convention.** HYP-9541's title writes the strong spectrum as "Spec_3 minus {1, 64}", so order 1 is not
  counted as strong. (5) counts order 1 as strong. State the difference when claiming confirmation.
- (d) **Citation.** Cite THM-4744 for Busch's formula f(3q) = 3 * 5^(q-1), f(3q+1) = 5^q, f(3q+2) = 9 * 5^(q-1)
  (EJC 13 (2006) N3, Cor. 1), next to THM-4870.
- (e) **Optional strengthening.**
  - H_3 - 1 = sum alpha_k 3^k >= (3/2) sum alpha_k 2^k = (3/2)(H - 1) for every tournament.
  - So a strong tournament of order >= 11 has H_3 >= 1 + (3/2) * 224 = 337, and (5) is exact (not only
    constructive) below 337.
  - With the weaker Busch bound 5^((n-1)/3), which needs no monotonicity, it is still exact below 322.
  - The audit's complete order <= 10 data show no further gap below 337 (`compare_census.out`).

**m4. THM-4950, Remark "Comparison with base 2": the binary lane list is incomplete.**
- "THM-4094's binary lanes are p (odd prime != 7) and 7p" omits the condition p != 3 in the 7p lane (7 * 3 = 21).
- It also omits the binary carry atoms 63 and 343 (THM-4094 status; S86 README section 5 table: the irreducibles of
  the odd numbers without 7 and 21 are p != 7, 7p with p not 3 or 7, and 49, 63, 343).
- 63 = 9 * 7 and 343 = 7^3 are the binary counterparts of the sporadic 40 = 4 * 10 and 29791 = 31^3. So base 2 has
  sporadic irreducibles too. Correct the remark and the sentence explaining why base 3 has more lanes.

**m5. PREREGISTRATION.md (Q1-Q3) has no appended outcome.**
- Its last line is still "Outcome: to be appended below without editing the predictions". The outcome exists only in
  coverage10.out, THM-4950 and the README.
- Append it, as was done in PREREGISTRATION-2 and -3:
  - m* = 9679, value 29,038;
  - Q1, Q2 and Q3 CONFIRMED;
  - Q3 holds automatically once the closure is right: an irreducible covered by the closure must be a single piece
    value.
- The file history shows the predictions were not edited. The scratch timestamps are consistent with
  "pushed before the run": ternary_order10.exe was built at 22:06:38, after commit 30c08c495b at 22:06:23, and the
  t10 outputs date from 22:28.

**m6. THM-4951 (1), third bullet, is Session26's corollary.**
- The bullet is "Nodes with nu <= d form an independent set, loops included, when n >= 6d + 3".
- Session26 states the same with tau, with the same proof (session26-blue/RESULT.md section 1: "The nodes tau<=d form
  an independent set whenever n>=6d+3"; session26-tournament/WIRING_REVIEW.md: "A loop needs
  tau>=ceil(floor(n/3)/2)").
- The nu form follows because the consecutive-triples proof exhibits disjoint 3-cycles.
- Credit it in the KNOWN list, and exclude it from "NEW: the tight-line analysis (1)-(3)".

**m7. Credit for the staircase lemma (Setting, "Staircases").**
- The staircase on [0, n-1] is the repository's full tiling tournament T*_n, the transitive node's line partner. It
  appears in:
  - THM-4731 (6), whose proof has the structure lemma "the out-neighbours of i are i+1 and all j <= i-2";
  - THM-337, as the "base-path staircase tournament", relabelled.
- The concurrent THM-4952 (C) records that its odd-cycle packings are the sets of disjoint odd intervals, i.e. that
  its cycles are the intervals.
- The lemma is elementary and the proof given is complete, but these should be cited. THM-4731 should also go into
  `related`.

**m8. THM-4951, proof of (2b).**
- Case c >= 2: the upper bound "since the other cycles lie in M - {a, a + c - 1}" covers only packings that have a
  cycle meeting L1 and a cycle meeting L2. Add:

  > If a packing has no cycle meeting L1 (resp. L2), its remaining cycles lie in M - {a + c - 1} (resp. M - {a}) or in
  > M, giving at most 1 + floor((c - 1)/3) or floor(c/3) cycles, both <= 2 + floor((c - 2)/3). One cycle meeting both
  > L1 and L2 contains both a and a + c - 1, which is covered by the same count.
- Case c in {0, 1}: "Every cycle meets both L1 and L2 or contains the single middle vertex" should read:

  > For c = 0 every cycle meets L1 (L2 is transitive), hence contains a; for c = 1 every cycle meeting L1 or L2
  > contains the middle vertex a, and no cycle lies in M = {a}. So all cycles share a.

  The conclusion nu(u) = 1 is right.

**m9. THM-4951, proof of (3b), "Colour".**
- "sigma reverses the triple word, so t_beta is sigma-fixed iff beta is a palindrome" proves only "blue implies
  palindrome".
- For the converse, add:

  > v -> n - 1 - v maps the window [3s, 3s + 3r + 1] of the run [s, s + r - 1] to the window of the mirrored run
  > [k - s - r, k - s - 1], so for a palindromic beta the set of tiles inside windows is sigma-invariant.
- The audit verified both directions for all 126 block words with k <= 6.

**m10. Notation.**
- In THM-4951, L is used in the proof of (3a) ("[5, L - 1]") and e in (4) ("2^(e-1) blue lines") without
  definition. Define L = n - 2 = 3k, and e = (m + f)/2 = N + c, the number of sigma-orbits of tiles.
- The credit line "The tiling frame, sigma, tau and blue/black lines: THM-643, ..." uses tau for THM-643's complement
  involution. THM-4951 uses tau for the feedback vertex number everywhere else, as in (2a). Write "the complement
  involution" there.

**m11. THM-4951 (5): labels and credit.**
- **Labels.** (5)(a) is THM-791 (a). (5)(b) is the THM-4922 (4) audit proof restricted to blue tilings, with the
  two-line line-count corollary. Label them KNOWN (the corollary PROVED), not PROVED. THM-4951's own body text already
  credits them as known.
- **Proof of (5d).** The constant 1 in f(t) + f(t + 1) = A(y) + A(y + 1) + 1 + x.(B(y) + B(y + 1)) is
  sum_i B_i(y + 1) = 1, i.e. (P2), which is step 5 of the THM-4922 (4) audit proof (|O| odd). Cite it at that line.
- **Credit.** THM-791's "odd-n contrast" is the first record of the odd-order tiling balance: at n = 5 the 16 blue
  tilings split 8/8 by H mod 4, and blue Delta H = 2 occurs. The THM-4922 audit's prior-art note says the same.

**m12. Upgrade available for THM-4951 (3c): the n = 3k + 1 count can be proved with (3a)'s induction.**
- **Claim.** The tight colourings of [0, 3k - 2] are exactly the words that are constant on the blocks of a
  decomposition of [0, 3k - 2] into one block of 2 and k - 1 blocks of 3, in any order.
- **Proof.**
  - If x_0 != x_1, then P(A) + P(B) >= 2 + P([4, 3k - 2]) = k + 1, so the word is not tight.
  - If x_0 = x_1 != x_2, then P(A) + P(B) = 1 + (P + P on [2, 3k - 2]). The rest is an n = 3(k - 1) + 2 instance,
    covered by (3a).
  - If x_0 = x_1 = x_2, the rest [3, 3k - 2] is an n = 3(k - 1) + 1 instance, and induction applies.
  - Base case: k = 1.
- **Count.** W(k) = 2^(k-1) + 2 W(k-1), with W(1) = 2. So W(k) = (k + 1) 2^(k-1), which is (k + 1) 2^(k-2) up to
  complement.
- **Check.** Brute force agrees for k <= 7 (`hyp9540_a_blocks.out`).
- **Text changes.** Upgrade the n = 3k + 1 count from "fitted" to PROVED. In HYP-9540 (A), "appear to be the block
  words" becomes "are the block words" (only realizability stays open).
- The n = 3k count (k^2 + 9k - 2) 2^(k-4) stays fitted. It was brute-forced here for k <= 7.

**m13. HYP-9540 (C), wording.**
- "Even j is realised at every n >= 6j + 4 by the double staircase" mixes the two j's. The double staircase with
  parameter j realises the pair {2j, k - 2j} for every n >= 6j + 4.
- The README row (2b), "blue tight lines with pair {2j, k - 2j}, for every n >= 10", needs the same qualifier: some
  blue tight line for every n >= 10, and the pair {2j, k - 2j} for n >= 6j + 4.

**m14. THM-4951 (2a), wording.**
- "Together with the transitive line (j = 0)" is redundant: j = 0 is in the range, and t_0 is the transitive tiling.
  Write "(j = 0 gives the transitive line)".

**m15. `related` fields.**
- THM-4950 cites THM-4952, HYP-9541 and Busch (via THM-4870) in its status, but its `related` field lists none of
  them. Add THM-4952, HYP-9541, THM-4744 and THM-4870.
- THM-4951: add THM-4731 (m7) and the Session26 packet path.

**m16. Code comment.**
- The header of `nu_lines.c` says it "also records which endpoint is strong (closed tilings ...; THM-4729)". It does
  not. Delete the sentence.

**m17. INDEX rows.**
- 05-knowledge/hypotheses/INDEX.md line 2 and 05-knowledge/results/INDEX.md line 5 say "audit PENDING". The first
  also repeats the refuted (D) ("(1,3) excess at odd n >= 11", M2).
- Update them when the corrections are applied.

**HYP-9530 diff.** It is accurate. The coverage to 29,035, the threshold 29,038 = 2 * 14,519, the lane list and the
falsifier range m <= 9678 were all reproduced. Optionally add m2's remark (orders <= 10 suffice).

## 3. What was verified independently, and how

All code in this directory was written by the auditor. Large intermediate files (class lists, order-10 witnesses,
the order-10 bitmap) are in the auditor's scratch folder
`C:/Users/Eliott/AppData/Local/Temp/claude/800365894db7d24c/70a6d0f6-6709-458c-a124-7d36217695ea/scratchpad/aud87/`
(`acls_1.txt` ... `acls_9.txt`, `my_t10.bin`, `my_t10_wit.txt`). They are regenerated by `gen_classes` and `ext10`.

### 3.1 THM-4950 (ternary Chan-Pak)

**Class generator (`gen_classes.c`, `gen_classes.out`).**
- One-vertex extension with the auditor's own individualisation-refinement canonical labelling. A key of the form
  "min code over relabellings" never merges non-isomorphic classes, so the lists are complete whatever the canonical
  form does.
- Class counts are 1, 1, 2, 4, 12, 56, 456, 6880, 191,536 for orders 1-9 (A000568). The strong counts are 1, 0, 1, 1,
  6, 35, 353, 6008, 178,133 (A051337).
- For every class, H_3 was computed by an independence-polynomial subset recursion, and OCF was checked against a
  Hamiltonian-path DP: 0 failures.

**Orders <= 9 (`compare_census.py`, `compare_census.out`).**
- The multiset of strong packing vectors is identical to S86's alpha_n.txt at every order 3-9.
- There are 1,909 distinct strong H_3 values (as in strong_spectrum10.out).

**Order 10 (`ext10.c`, `ext10.out`).** This is a different method from `ternary_order10.c`: no canonical forms and no
order-10 subset DP.
- For each of the 191,536 order-9 classes and each out-set O of the new vertex, the cycles through the new vertex on
  S + {v} are counted as sum over u in S cap O, w in S \ O of hp(S; u -> w).
- Then I(Omega, x) = Z_x(R) + x * sum_S cyc(S) Z_x(R - S).
- 94,224,222 strong extensions were processed, with repeats.
- **Values.** There are 10,070 distinct H_3 values; the minimum is 199 and the maximum 31,549. The set is identical to
  the union of the session's 12 bitmaps and to witnesses_order10.txt.
- **Checks.** 46,002 sampled extensions were recomputed directly at order 10 (H_3, OCF and Hamiltonian paths):
  0 failures.

**Witnesses (`check_witnesses.c`, `check_witnesses.out`).** All 10,070 codes were decoded and checked: strong, H_3 equal
to the stated value, and I(Omega, 2) equal to the Hamiltonian-path count. There were 0 failures.

**(5) Strong spectrum.** Below 225 the missing values v = 1 (mod 3) are exactly 10, 31 and 64. The first value >= 225
with no strong realiser of order <= 10 is 28,981. No other value is missing below 337.

**(3), (4) Coverage (`closure_check.py`, `closure_check.out`; `lanes_check.py`, `lanes_check.out`).**
- The closure code is the auditor's own.
- Strong orders <= 10 alone cover all 1 + 3m with m <= 9678, except m = 3 and 10.
- With S86's pieces (inventory; one-apex words recomputed from the stated formula with the auditor's own routine;
  rays), m* = 9679, value 29,038 = 2 * 14,519, prime. The next misses are 9699, 9730, 9750, 9751, ...
- 4,063 order-10 values are not given by S86's pieces.
- Calibration: m* = 370 (orders <= 8), 1800 (orders <= 9), 9679 (orders <= 10).
- 29,038 / 31,549 = 0.9204.
- The first missing irreducible in each lane is Q 29,038, P 29,191, 31P 30,721, 31Q 31,093 and 10P 31,210.
- The sporadics 40, 250, 310, 9610 and 29791 are all strong values; 29,791 is attained at order 10.

**(2) Lanes below 200,000 (`lanes_check.out`).**
- Irreducibility was tested by direct divisor search in M3, with the auditor's own lane classifier.
- There are 23,925 irreducibles: P 8987, Q 12882, 10P 1123, 31P 413, 31Q 515, sporadic 5.
- No irreducible lies outside the lanes, and every lane element is irreducible.

**Proofs.**
- (1): correct, with the labels fixed as in m1.
- (2): re-derived case by case (l(v) >= 4, 3, 2, 1, including the pairing argument for 40 and 250 and the uniqueness
  argument for 9610 and 29791). Correct.
- (5): correct (see m3 (e) for the stronger range).
- The bound H_3 <= (3/2)^K H, with K = 3 and max H(10) = 15,745 (A038375), gives order-10 values < 53,140.

**Logs.** The session's logs give 70,960,990 processed order-10 classes (with repeats) and 277,185 OCF checks with
0 bad, as stated.

### 3.2 THM-4951 (lines between packing-number classes)

**nu-pair census, all lines, n = 3-9 (`audit_lines.c`, `audit_lines.out`).**
- Own tile order and own nu routine (branch on the highest free vertex).
- Line totals are 2^(m-1), and the blue totals 2^(e-1).
- Every row of the THM-4951 (4) table and of nu_pairs.out is reproduced exactly, including n = 9 (134,217,728 lines).

**Blue lines, n = 4-12.**
- n = 4-11: nu-pairs plus exact H mod 4 and H mod 8 by Hamiltonian-path DP.
- n = 12: nu-pairs only, all 2^29 lines (`audit_blue12.out`).
- New data at n = 10 and 11: no blue (1, 1). Every allowed pair occurs (n = 10: (1,2) 9, (1,3) 25, (2,2) 409, ...;
  n = 11: (1,2) 18, (1,3) 535, (2,2) 184, ...).
- n = 12: see M1.

**Tight lines, n = 3-12 (`tight_audit.c`, `tight_audit_12.out`).** This search differs from `tight_lines.c`: it runs
word by word over the tight colourings, with the lookahead bound nu(T[0..v]) + P(A cap [v + 1, n - 3]) <= P(A); n = 12
takes 2 min 48 s on one core.
- **Tight lines.** 1, 1, 1, 57, 31, 13, 3673, 1802, 541, 563,609 for n = 3-12. Blue: 1, 1, 1, 1, 1, 1, 39, 10, 19, 297.
- **Pairs.** n <= 11: only (0, k) and (1, k - 1). n = 12: (1, 3) 96,254 black; (2, 2) 467,058 black + 296 blue.
- **Words.** Every tight colouring is realised (1, 1, 1, 5, 3, 2, 17, 8, 4, 50 with x_0 = 1). Every palindromic tight
  colouring (1, 1, 1, 1, 1, 1, 3, 2, 2, 3) is realised by a blue line.
- The n = 12 per-word histogram is identical to tight_lines.out (50 of 50 words).
- The n <= 9 values agree with the full census above.

**T5 witnesses.** All 28 first-hit lines at n = 13 and 14 were re-checked: word, nu-pair and colour
(`staircase_audit.out`).

**Tight-colouring counts.** Brute force for n = 3-21 (k <= 7) agrees with 2^(k-1), (k + 1) 2^(k-2) and
(k^2 + 9k - 2) 2^(k-4). The session's DP lists to k = 13 also match the closed forms term by term.

**Constructions (`staircase_audit.py`).**
- (2a): n = 3-21, all j. nu = (j, k - j). The proof's explicit feedback sets have sizes j and k - j and leave acyclic
  tournaments, so tau = nu at both ends. The line is black iff j >= 1.
- (2b): j = 1-3, c = 0-9 (n <= 26). Blue, with nu = (2j, k - 2j).
- (3b)/(3d): all 126 block words with k = 1-6. Triple word = beta; nu = (|beta|_1, |beta|_0); u - F is acyclic; blue
  iff beta is a palindrome.
- Single reversed tile (0, 3), n = 4-20: nu-pair (1, k), and the line is black for n >= 5 (HYP-9540 (B),
  `single_tile_check.out`).

**(5) Blue lines by H mod 4.**
- Even n = 4, 6, 8, 10: every blue line is (1, 1), by exact Hamiltonian-path counts.
- Odd n:

  | n | (1,1) : (1,3) : (3,3) | sources |
  |---|---|---|
  | 5 | 2 : 4 : 2 | exact Hamiltonian-path counts and the odd-cycle parity method |
  | 7 | 64 : 128 : 64 | same |
  | 9 | 8192 : 16384 : 8192 | same |
  | 11 | 4,014,080 : 8,749,056 : 4,014,080 | Hamiltonian paths on all 2^25 blue tilings; `blue_parity.c` (odd-cycle parity) |

- `blue_parity.c` computes H mod 4 = 1 + 2 (alpha_1 mod 2) from the parity of the number of directed odd cycles (OCF),
  with no Hamiltonian paths and no mirror paths.
- The H mod 8 pair counts at n <= 9 equal blue_lines.out.

**(5b), (5d) at n = 11 (`blue_parity.out`).**
- H = 3 (mod 4) on exactly 2^24 of the 2^25 blue tilings.
- The affine law f = A(y) + x.B(y) holds on all 2^25 (y, x).
- (P2) holds for all 2^20 values of y.
- |D| = 94,208, |D_1| = 69,632 and |D_0| = 24,576.
- The formula 2^(N+c-2) + 2^(c-2)(|D_1| - |D_0|) gives 8,749,056.
- D is empty at n = 5, 7, 9. (At n = 3, N = 0 and D is trivially {empty word}; the formula still holds.)

**HYP-9540 (D) at n = 13 (`blue_sample.c`, `blue_sample.out`).**
- The sampler was validated at n = 11: 400,000 samples give 1/2 + 0.02111 (se 0.00023), against the exact 0.02148.
- n = 13: 600,000 samples give 1/2 - 0.00109 (se 0.00011, z = -10.1), with |D| / 2^30 ~ 0.0275. U1 is refuted
  independently.

**Proofs checked line by line.** THM-4951 (1), (2a) (including tau, j = 0, j = k - 1 and b = 1), (2b) (gap m8), (3a),
(3b) (gap m9), (3d), (5b) (the involution flips only mirror tiles, which are never path arcs) and (5d) (needs (P2),
m11). All are correct.

### 3.3 Pre-registrations and prior art

**Pre-registrations.**
- The predictions in PREREGISTRATION.md, -2 and -3 are unchanged since their first commits (`git diff` 30c08c495b,
  d89a19995b, 32a7f2e85a against HEAD). Only outcome sections were appended.
- The outcomes match the predictions as written: Q1-Q3 CONFIRMED (not appended, m5); R1 CONFIRMED; R2 and R3 REFUTED;
  T1-T5 CONFIRMED.

**Prior art searched** (git grep over 01-canon, 00-navigation, 07-reflections, 05-knowledge, 04-computation).
- Terms: staircase, packing number, disjoint 3-cycles or triangles, tight line, tight colouring, block word,
  multi- / double / prefix staircase, nu-pair, floor(n/3), tau(T) + tau(T'), 3673, 563,609, 4,014,080, 8,749,056,
  1:2:1, 29,035, 29,038, 9610, 29791, M3, 10p / 31p / 31pq, irreducible with mod 3 or monoid, ternary lane,
  not factorial.
- No earlier statement was found of:
  - the tight-line analysis, the staircase families, the block law or the counts;
  - the nu-pair census;
  - the line-level mod-4 splits;
  - the base-3 lane reduction.

  S86's README section 8 proposed the lane reduction as "the natural route".
- Credits to add: m3, m4, m6, m7, m11.
- **Concurrent mac-mini lane** (THM-4952, THM-4953, HYP-9541, ternary_ribbon_metagraph_lines_20261108.md):
  - no overlap with THM-4951, apart from the staircase lemma (m7);
  - THM-4953 (A) is a complementary line bound for the number of strong components;
  - for THM-4950 (5), see m3.

## 4. Files (audit/)

| file | content |
|---|---|
| gen_classes.c / .out | own class generator, orders 1-9, strong flags, H_3, OCF |
| ext10.c / .out | order-10 strong H_3 set by incremental extension (no canonical forms) |
| my_t10_values.txt, strong_values_le10.txt | the auditor's order-10 values; all strong values of orders <= 10 |
| check_witnesses.c / .out | re-check of all 10,070 witnesses |
| compare_census.py / .out | comparison with S86 and the session's data; strong spectrum |
| closure_check.py / .out | multiplicative closure, m* with and without S86's pieces |
| lanes_check.py / .out | lane counts, first missing irreducibles, 4,063, calibration |
| audit_lines.c / .out, audit_blue12.out | all-line nu census n <= 9; blue lines n <= 11 (nu, H mod 4/8); blue n = 12 |
| tight_audit.c, tight_audit_12.out | tight lines n <= 12 by word-by-word lookahead search |
| staircase_audit.py / .out | constructions (2a), (2b), (3b), (3d); tight-colouring brute force; T5 witnesses |
| single_tile_check.py / .out | the (1, k) construction of HYP-9540 (B) |
| palindromic_pairs.py / .out | possible blue tight pairs n <= 26 (M1) |
| hyp9540_a_blocks.py / .out | n = 3k + 1 block characterisation (m12) |
| blue_parity.c / .out | odd-n blue lines by odd-cycle parity: splits, affinity, (P2), D, D_0, D_1 |
| blue_sample.c / .out | sampled sign of the (1,3) excess (n = 11 validation, n = 13) |

**Build and run.**
- Build each C file with `gcc -O3 -march=native -fopenmp`.
- Run `gen_classes OUTDIR 9`, then `ext10 OUTDIR/acls_9.txt OUTDIR/my_t10`.
- The Python scripts take the paths shown in their docstrings.
