# Independent audit of kind-pasteur-2026-11-13-S26 (THM-5151, HYP-9678)

**Verdict: PASS WITH CORRECTIONS.**

Every PROVED item is correct: Propositions 1, 2, 6, 6', 7 and 8. Every FINITE-EXACT number was reproduced with
independent code. That covers the order <= 9 census, the capacity formula, Proposition 3, the uniquely Hamiltonian
counts to order 12, the order-10 collisions, the 317 pair, the clock values and primality. I also replicated the
pre-registered search S3 and its H = 333 refutation of HYP-9678 (A).

The corrections are about wording, scope and credit:

- a one-apex direction convention that makes several sentences false under the usual reading;
- an unranged "2647 is the largest" claim;
- missing or misattributed credits;
- one false side remark about Session 310's A.

Chan–Pak is not claimed anywhere; every file says it remains open.

Audited state: commit c4a76394db (checkpoint 3), including the coordinator's addenda: Proposition 6', Proposition 8,
the outcomes of S1, S2 and S3, ysearch_bic.c and yhash.c. All audit scripts are in the auditor's scratch directory,
not in the repository (listed at the end).

## MAJOR findings

**M1. The one-apex direction is unstated in the canon and misdescribed in prose.** The mathematics is right under
PROOF.md §4's convention: body i -> j for i < j, so body vertex 1 is the body's source. The prose reads the other way.

- In O(0011111) the apex beats body vertices 1 and 2. Their scores inside TT7 are 6 and 5, so they are the two
  *top* vertices. The natural reading of "apex beats the two lowest", O(1111100), is **not strong**: vertex 1 is then a
  source. I checked this by computation.
- Sentences affected:
  - README §"What was found": "B is the transitive tournament on seven vertices plus one apex that beats the two
    lowest."
  - PROOF.md §4: "Its apex beats the two lowest vertices of TT7".
  - PROOF.md §4 Credit: "For O(m,k) (apex beats the bottom m, loses to the top k)".
  - PROOF.md Proposition 6': "Omitting the lowest is admissible iff w_2 = 0". This is correct only if "lowest"
    means label 1.
  - verify_s26.py lines 174 and 177: "apex beats the two lowest / bottom two of TT7".
- verify_s26.py line 6, the docstring, contradicts its own check. It says "Session 310's B is isomorphic to the
  one-apex tournament O(5,2) (apex beats the bottom five of TT7)", but line 175 asserts that B is *not* isomorphic to
  O(0000011) = O(5,2).
- THM-5151 (6) defines O(w) with "transitive body 1..m, with w_i = 1 iff i beats x" and gives no direction. Then
  "strong iff w_1 = 0 and w_m = 1" and "B ≅ O(0011111)" are underdetermined: under the opposite direction they
  become w_1 = 1, w_m = 0 and O(1111100). HYP-9678 (B) inherits the same ambiguity.
- Session 191's one_apex/PROOF.md uses the complementary convention: x -> i iff b_i = 1, strong iff b_1 = 1 and
  b_n = 0. The packet's O(w) is Session 191's T(complement of w). The same word therefore names different
  tournaments in the two packets.

**Fix.**

- In THM-5151 (6) and HYP-9678, write: "body 1 -> 2 -> ... -> m (i -> j for i < j; vertex 1 beats the rest of the
  body), w_i = 1 iff body vertex i beats x".
- Replace "two lowest / bottom m / top k" by "the first (top-ranked) two / first m / last k".
  - B: "the apex beats the two top-ranked vertices of TT7 and loses to the other five".
- Fix verify_s26.py line 6 to read "O(2,5) = O(0011111) (apex beats the first two body vertices)".
- Add one line: "O(w) = Session 191's T(w̄)".

**M2. "2647 is the largest C = 1 miss" is stated without its range.** The canonical status is FINITE-EXACT through
10^9 only. HYP-9494 (D) is FINITE-EXACT through 10^9, and THM-5148 says "the largest 2647: KNOWN (HYP-9494 (D), to
10^9)". The all-order statement is open.

- Sentences affected:
  - README §Tangents: "2647 is also the largest value no uniquely Hamiltonian tournament attains."
  - HISTORY_TANGENTS §A: "…= 2647 is the largest odd number that no uniquely Hamiltonian (C = 1) tournament
    attains." The next sentence mentions 10^9, but the claim itself is unqualified.
  - 05-knowledge/results/INDEX.md row: "Tangents: 2647 = largest C=1 miss".
- **Fix.** Write "the largest odd value below 10^9 that no C = 1 tournament attains (HYP-9494 (D), FINITE-EXACT);
  whether it is the largest overall is open."

**M3. Missing and misattributed credit.**

(a) **Proposition 7 / THM-5151 (8).**

- The mechanism already appears in Session 309's
  [MIGRATING_CLOCK.md](../../session309_changing_module_returns_20261113/MIGRATING_CLOCK.md) §2: "Modulo 7 the ratio
  1611/19 is 3 of order 6 … Thus 7|H_n iff n≡0 mod6". That is the same Π/Γ clock family and the same day.
- "3/2 is a primitive root mod 7" is also recorded in THM-4800.
- The packet says only "noticed by this session's history sweep".
- **Fix.** Credit Session 309 MIGRATING_CLOCK §2 in PROOF §5, THM-5151 status (8) and HISTORY_TANGENTS §A. Drop
  "Rabinowitsch-type" as a name, or mark it as an analogy only: there is no class-number content.

(b) **One-apex credits.**

- THM-5151 status (6) and PROOF §4 cite "Session 191 one_apex/PROOF.md (1)" for "odd cycles of a one-apex
  tournament are increasing chains through the apex". Formula (1) of that file is the H append-recurrence, and the
  file never mentions cycles. The citation appears to be copied from THM-5070's recurrence credit.
- The chain/parity description of cycles through feedback vertices predates THM-5070. It is in
  feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md §4 (Oct 20), and THM-5134 (G) states "every odd cycle
  of a one-apex tournament passes through its apex".
- PROOF §4's "H = 2 + (2^m − 1)(2^k − 1)" for O(m,k) is Session 191's U(o,z), formula (5), and it is uncredited.
- **Fix.** Replace the Session 191 (1) citation with THM-5070 (2), THM-5134 (G) and intercyclic_rays.md §4. Credit
  Session 191 (5) for the O(m,k) value.

(c) **Profile-collision censuses.** Session 261
[constructors/MODULE_REBASE.md](../../session261_creative_chanpak_20261111/constructors/MODULE_REBASE.md) §3 already
ran full cover-profile collision censuses at orders 5–6: 56 classes and 24 profiles at order 6, which matches my
count. It did not separate strong from non-strong, so THM-5151 (4) remains new. Cite it as the precursor in PROOF §3.

**M4. A false statement in README §S1.** The README says: "Session 310's partner A has this shape ({5,6} × {0,2})."

- That tournament is TT7 with 5 -> 0, 5 -> 2, 6 -> 0, 6 -> 2 reversed. It is **not** isomorphic to Session 310's A
  = (2,28,25,81,97,79,7).
- Session 310's A has no biclique presentation of any sizes, although its minimum number of backward arcs is 4.
- Two *other* classes in A's five-class profile group are bicliques: ({0,2},{5,6}) and ({0,1},{3,5,6}) for one, and
  ({0,1},{4,6}) and ({0,1,3},{5,6}) for the other.
- The remark about the order-9 H = 317 partners is correct: neither class is a biclique of any sizes.
- **Fix.** Write "one of the five order-7 partner classes in A's profile group (not Session 310's literal A) has
  this shape".

## MINOR findings

**m1. README "What was found": a stale heading.** The heading "The complete stock is tiny and uniquely Hamiltonian"
is now contradicted two lines below by the S3 result. Retitle it "Through order 9 the stock is tiny and uniquely
Hamiltonian".

**m2. Proposition 7 notation.** PROOF §5 says "Session 309's clock (I = 499, c = 83,404,800)".

- In the stated form H_k = I + c 2^(r−k) 3^k, Session 309's clock has c = 1,303,200 and r = 6. The number
  83,404,800 is c·2^r.
- The next sentence, "Primes l with −I/c outside the subgroup", should read −I/(c 2^r).
- Fix both. The listed primes 13, 19, 29, 43, 47, 53 are correct, and they are *all* such primes below 60. The next
  ones are 61, 67, 71, 97, 101.

**m3. The divisibility-by-7 condition is incomplete.**

- README: "any run of six consecutive transfers meets a multiple of 7 unless 7 divides the intercept". The THM-5151
  title (8) says "unless 7 | I".
- The correct condition is "unless 7 | c·I": when 7 | c and 7 ∤ I, no term is divisible by 7. For the 317 clock this
  happens whenever 7 | P_M(3).
- "Six consecutive transfers" spans seven values. Say "any six consecutive values H_k".
- The body bullets of (8) and Proposition 7 are already correct.

**m4. S1 scope sentence.** The README says: "Detachable partners there, if any, need at least 6 backward arcs and
are not bicliques."

- S1 only targeted Q_X/t for **uniquely Hamiltonian** order-12 X.
- The biclique run only covered |L|, |U| <= 5.
- **Fix.** Write "Order-11 partners of uniquely Hamiltonian order-12 members, if any, have no vertex order with at
  most 5 backward arcs and are not bicliques with both sides of size at most 5."

**m5. Pre-registration timing.** The README says "[PREREGISTRATION] was pushed at b1b011e4d2 before these ran."

- S1's output file was created at 11:06:02 −0700, and commit b1b011e4d2 is at 11:10:10. S1 was therefore already
  running at registration.
- Its k <= 4 progress line, all null, was probably already printed, because k <= 4 is about 10% of the work.
- S2 started at 11:10:37, after the commit. S3 started after 5d04e33349 (11:30:06): yhash.exe was built at 11:31:07.
- PREREGISTRATION.md's own wording, "before the two decisive searches below finish", is accurate.
- **Fix.** Use that wording in the README and state S1's start time.

**m6. THM-5151 status block.** It says "PRE-REGISTERED (b1b011e4d2), both restricted to order-11 tournaments…" and
then lists S3, which is order 10 and was pre-registered at 5d04e33349. Fix the commit reference and drop "both".

**m7. Proposition 5(d) can be sharpened.** "Vertex-deletion fingerprints give at least 2 isomorphism classes" can
become "exactly 2 classes, a converse pair, among the 120 labelled H = 333 matches". My canonical forms give 2 classes
each at H = 153, 219, 317 and 333.

**m8. HYP-9678 needs updating after the refutation.**

- (B)'s list of block-type members omits H = 333 (order 10, C = 5).
- The status bullet "no (3,9) or (5,7) block partner exists within C = 1" belongs under block pairs, not under (B),
  which is about detachable pairs.
- The H1 heading still reads "exact pair members are uniquely Hamiltonian". Update it to the (A') form. The filename
  can stay for link stability.

**m9. Conditional "must" sentences.**

- HISTORY_TANGENTS §D.2: "…must already hold an even component with a partner, because spent odd cores have no
  hidden capacity (Propositions 2-3)".
- README §"What this means": the bullet pair that says a return "must" alternate a fresh detachable even component
  with a split.
- Both depend on HYP-9678 (C) beyond order 9. Proposition 3 covers only two-core merges with |Y1| + |Y2| <= 10.
  Proposition 8 alone gives only "an unbounded return needs non-neutral moves".
- Prefix each with "Under HYP-9678 (C)".
- The README sentence "spent cores cannot be refilled (Propositions 3 and 3′)" also needs "through order 9".

**m10. "Odd merge" should include merges into an even component.** In PROOF §2 "Reading" and §6, "Merge odd strong
components into fewer ones" should explicitly allow merging three odd cores into one even component: R_X = R_Y1 R_Y2
R_Y3 gives +1 site. Proposition 3 does not cover that case, though Proposition 3' does through order 9. Without this,
§6's "Proposition 3 excludes them in small orders" reads as broader than it is.

**m11. Selective sweep line.** HISTORY_TANGENTS §A says "Ramanujan's ternary form x^2 + y^2 + 10 z^2 misses 7, 21,
43, 67 among small values". Its odd misses below 400 are 3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253,
307, 391. List the first few in full, or say "among its misses are …".

**m12. Reproducibility.**

- scripts/README has no step for S3 (yhash.c).
- No packet script writes the S1, S2 and S3 target files. The construction exists only in the author's scratch
  directory.
- Add the target generator: wouldY(Q_X) for S1, Q(Y5 => Y7)/t for S2, and the 48,908 two-block profiles for S3.
- scripts/README step 6 points to "PROOF section 7 notes" for the construction, but §7 has no construction.

**m13. Small wording nits.**

- PROOF Proposition 6: "Camion's theorem gives strongness". The Hamiltonian cycle gives strongness trivially;
  Camion is the converse.
- THM-5151 title (10): ">= 3" needs n >= 6. At n = 4 the value is 2, and the "never uniquely Hamiltonian" conclusion
  still holds.

## What I recomputed independently

All of the following used my own code: C plus Python with sympy and networkx.

**Engine.**

- Class lists come from my own one-vertex extension with an invariant-refined backtracking canonical form.
- Profiles come from the **backward-adjacency distribution**: F_j = Σ_b N_b binom(n−1−b, j−1−b) over all vertex
  orders. This is not the packet's subset convolution. For the literal graphs I also used a third route, the
  odd-cycle packing expansion (I2).
- Hamiltonian cycles were counted by direct DP, and isomorphism was tested with networkx VF2.

**Class counts.** Orders 1–9 give 1, 1, 2, 4, 12, 56, 456, 6880, 191,536 classes (A000568). Strong classes are 1, 0,
1, 1, 6, 35, 353, 6008, 178,133.

**Proposition 4 / THM-5151 (4).**

- There are 52,432 profile groups at orders 3–9; the per-order counts are 2, 3, 8, 24, 156, 1728, 50,511.
- Groups mixing strong and non-strong classes:

| order | H | strong classes | C | non-strong classes |
|---:|---:|---:|---:|---|
| 6 | 15 | 2 | 1 | 2, shapes (1,5) and (5,1) |
| 8 | 45 | 3 | 1 | 2, shapes (3,5) and (5,3) |
| 8 | 95 | 2 | 1 | 10, A → TT1 and TT1 → A over 5 A-classes |

- There are none at orders 4, 5, 7 or 9.
- Proposition 1's shape held in every mixed group.

**Proposition 3' / HYP-9678 (C).** The formula cap = s + d holds in all 52,432 groups through order 9. The detachable
even profiles are exactly H = 15 at order 6 and H = 95 at order 8.

**Proposition 3.** I checked against **all** classes Z, not only strong ones: 389 strong odd block pairs at order 10
against all 191,536 order-9 classes, plus the (3,3) and (3,5) cases. There are 0 merges.

**Uniquely Hamiltonian slice.**

- Extending **all** classes of order m, so without using Moon or Session 197's Lemma 1, gives 1, 3, 8, 21, 55, 144,
  377 C = 1 classes at orders 4–10.
- C = 1 extension then gives 987 and 2584.
- Order 10 was checked against **every** non-strong order-10 tournament: the last strong component of any order with
  any first part, so this does not rely on Proposition 1. Exactly 6 classes collide:
  - H = 153 (2) and H = 219 (2) with partners C3 => Y7;
  - H = 317 (2) with partners Y9 -> TT1.
- The partner groups have 4 strong order-7 classes (H = 153), 6 (H = 219) and 2 order-9 classes (H = 317).
- Order 12: 0 collisions with (3,9) or (5,7) block sums. There are 1,291 distinct S1 targets; the S2 targets are 655
  distinct among 2,118 pairs. Both counts match the packet.

**The six literal pairs, plus the H = 333 pair.** Each has F(X) = F(Z) and Q_X = Q_Z symbolically, X strong, and the
(I2) route agreeing. C(X) = 1 for the six pairs and 5 for H = 333. The Q polynomials:

| pair | Q |
|---|---|
| order 6, H = 15 | t^6+8t^4+6t^2 |
| order 8, B | t^8+20t^6+60t^4+14t^2 |
| order 8, H = 45 | t^8+12t^6+24t^4+8t^2 |
| order 10, H = 153 | t^10+16t^8+58t^6+66t^4+12t^2 |
| order 10, H = 219 | t^10+20t^8+82t^6+100t^4+16t^2 |
| order 10, H = 317 | t^10+26t^8+134t^6+140t^4+16t^2 |
| order 10, H = 333 | t^10+20t^8+106t^6+162t^4+44t^2 |

**The 317 pair (c).** F(X) is exactly as claimed. Q_X = t·Q_Y holds symbolically. X is strong with C = 1; Y is strong
with C = 8; Q_Y(2) = 9280. Neither is one-apex.

**B (b).**

- Among all 128 words, only O(0011111) is isomorphic to B, under PROOF's convention.
- B is not isomorphic to O(0000011).
- A is strong with C(A) = 7 and cycle counts {3: 10, 5: 18, 7: 7}. It is not one-apex.

**Proposition 6.**

- The closed form matches on all 1,022 apex words of length 1–9.
- "Strong iff w_1 = 0 and w_m = 1, then C = 1" matches on all of them.

**Proposition 6'.**

- c_(n−1)(O(w)) = (n−3) + [w_2 = 0] + [w_(n−2) = 1] by direct (n−1)-cycle counts on all 168 strong words at
  n = 6, 8, 10, and via the closed form on 10,920 words for even n up to 16.
- C(Y5, code 152) = 3, C(Y5', code 739) = 2, C(A) = 7.
- C(Y7) is 3 at H = 153, 4 at H = 219 and 11 at H = 333; C(Y9) = 8.
- c_9 + 2d = 6, 8, 8 and 22 match 2·C(C3)·C(Y7) or C(Y9) in every case.

**Proposition 8.** I checked the proof by hand, including the reordering step. It is correct and unconditional.
Proposition 2's realizability is automatic, since s = β − Σ ord(S_j) is forced by degree.

**The H = 333 refutation (addendum 3).** X = (1002,…,16) is strong with C(X) = 5, H = 333, and Q_X = Q(C3 ⇒ Y7). Y7
(code 43104) is strong with H = 111 and C = 11, and it is unique in its profile.

**My own replication of S3 (s3.c, pull-DP profiles).**

- 48,908 target profiles; 1,385,979 candidates, of which 507,458 are strong; 862 matches. All of these equal the
  packet's numbers.
- Matches split 440 / 272 / 30 / 120 at H = 153 / 219 / 317 / 333.
- The detachable (1,9) tag occurs only at H = 317.
- The H = 333 matches form exactly 2 classes, a converse pair.

**Clocks (d).**

| run | values | status |
|---|---|---|
| M = TT3, r = 2 | 558,047,332,489; 837,070,948,489; 1,255,606,372,489 | prime; 79·1321·2053·3907; 7 × 179,372,338,927 with the cofactor prime |
| M = TT1, r = 2, k = 2 | 139,511,908,489 | prime |
| M = TT1, r = 3, k = 3 | 3,884,008,766,575,013 | prime |

- The recurrence holds.
- The B^r clock is divisible by 5 in all tested cases.
- 3/2 ≡ 5 has order 6 mod 7.
- Session 309's values reproduce, and 187,661,299 = 7 × 26,808,757.
- sympy agrees with the packet's Miller–Rabin.

**HISTORY_TANGENTS *checked* items.**

- The C = 1 frame values Q_X(2) + H(X) through order 9 have 85 entries at most 2647, with 45 distinct values. The
  misses are exactly 231 (the order-6 member) and 2647 (B). The expected number of hits for a random subset is 0.99.
- The Heegner tournament is a tournament with H = 71.
- 17, 73 and 107 are the 7th, 21st and 28th primes.
- (740, 760) is the unique pair of primes of the form n² + 1 that differ by 30000.
- 2·3^k + 1 for k = 0..4 gives 3, 7, 19, 55, 163, and 55 is a miss.
- The circulant on Z_7 with connection set {1,2,3} has H = 175 = 2^8 − 3^4.

**Packet verifier.** `python verify_s26.py` and `python -O verify_s26.py` both pass, as a cross-check after my own
runs.

**Novelty sweep.**

- Not found anywhere in the canon: a prior strong/non-strong exact-pair census, B identified as one-apex, the
  R-factor / cap(Q) bookkeeping, the 317 pair, 9280, 558,047,332,489, or 179,372,338,927.
- Session 311's WREATH_ATLAS.md uses B and A -> TT1 as leaves of C3 wreath towers but never mentions an apex. The
  "new identification" stands.
- The prior items found are those listed in M3.

Auditor scripts, in C:\Users\Eliott\AppData\Local\Temp\claude\2f1de76a73e81b8a\a7c1f371-6162-4db2-99bd-68121f2beaaf\scratchpad\audit26:

- engines: aud.c, s3.c;
- censuses: census_audit.py, c1_collisions.py, prop3_audit.py;
- literal graphs: pairs_audit.py (with pairs_audit_lib.py), h333.py, partner_counts.py;
- one-apex: biclique_check.py, oneapex_audit.py, prop6p_audit.py;
- arithmetic and tangents: clock_audit.py, tangents_audit.py.
