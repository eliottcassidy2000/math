# Audit of opus-2026-11-08-S86: THM-4940, HYP-9530 and `packing_bases_repunits_20261108`

Independent hostile-but-fair referee, 2026-11-08.

- **Audited revision:** 725d3bb4cf (the checkpoint named in the brief). The worktree HEAD moved to 1d57e6b23c during
  the audit. The three later commits (ebe3096762, 1a1606cc66, 1d57e6b23c) touch only `PREREGISTRATION-2.md`,
  `coverage_products.*` and `ternary_attainment.*`. They are outside this audit, and the audited files did not change.
- **Method:**
  - Every PROVED claim was read line by line against its cited sources.
  - Every computation was redone with the audit's own code (`aud_*`, written from scratch). The only reused tool is
    kind-pasteur's `classgen.c`, used for the class lists of orders <= 9.
  - No session script was imported or copied. The session's *output files* were read only as data to compare against.

## Verdict: PASS WITH CORRECTIONS

- **The mathematics is sound.**
  - Every proof of THM-4940 (1)-(8) is correct.
  - Every number in THM-4940, HYP-9530 and the README tables that the audit recomputed agrees exactly, including the
    exact gap tables, rescue sets, the {4, 8} set, the repunit lists, the 26 singletons, the fibre-move counts and the
    corrected S/C ledger.
  - The audit's enumerator differs from the session's in three ways: its canonical form, its cycle counter, and its
    stopping rule (it runs until a level is EMPTY). It reproduces the 5,324-class inventory and the 1,461-class
    intercyclic inventory class for class.
- **The pre-registration ordering is confirmed** (section 7).
- **Corrections needed:**
  - **Two MAJOR items, both fixable by text:**
    - M1: the base-2 half of statement (8) goes beyond its proof. The claim is true, and the audit's witnesses close
      the gap.
    - M2: a dictionary row typed PROVED repeats the old false claim "alpha_1 = 3 is impossible".
  - **Fifteen MINOR items:**
    - several MISTAKE-559-type sentences that state open completeness (Chan-Pak, HYP-9477 (A)) as fact;
    - incomplete credits for earlier evaluations of I(Omega, x) at x != 2;
    - a stale merge of the concurrent kind-pasteur S17 lane;
    - three false or overstated README sentences (even repunits, "prove that no map exists", Sylvester's sequence).

## 1. Proofs, line by line (task 1)

**(1) Basic properties.** All four are correct.

- **Multiplicativity.** Every directed cycle lies inside one strong component, so Omega(A => B) is the disjoint union
  of Omega(A) and Omega(B), and independence polynomials multiply over disjoint unions.
- **Converse.** Reversing every arc reverses every cycle, so Omega(T^op) is isomorphic to Omega(T).
- **Mod b.** alpha_0 = 1, and every other term is divisible by b.
- **Deletion.** Omega(T - v) is an induced subgraph of Omega(T), so alpha_k(T - v) <= alpha_k(T) for every k.
  Hence H_b(T - v) <= H_b(T) for every b >= 0.
- **Credit.** Multiplicativity at every x is already THM-467 (d)1, and the congruence is HYP-886 / HYP-892 (MINOR m6).

**(2) Block lemma.** Correct.

- If sum_{k>=1} alpha_k b^(k-1) = a < b, then every alpha_k with k >= 2 is 0 and alpha_1 = a. The converse is
  immediate.
- The phrase "Spec_b and A agree on [1, b^2)" is loose (m4). The README's "below b^2, Spec_b equals 1 + bA" is the
  right statement.

**(3) The base-free 7-hole.** Correct, including the real-x extension.

- **Reduction.** Divide by x > 0. Then I(Omega, x) = 1 + 3x means alpha_1 + alpha_2 x + alpha_3 x^2 + ... = 3.
  - Any alpha_k >= 1 with k >= 2 forces alpha_2 >= 1, hence alpha_1 >= 2.
  - If alpha_1 = 2, then alpha_2 <= 1 and alpha_3 = 0. So 2 + alpha_2 x = 3, which forces alpha_2 = 1 and x = 1.
  - If alpha_1 >= 4, the sum exceeds 3.
  - If alpha_1 <= 1, then all higher alpha_k vanish and the sum is at most 1.
  - So for every real x > 0 with x != 1, only alpha = (3, 0) remains.
- **(3, 0) is impossible.** Three pairwise intersecting cycles lie in one strong component C, and alpha_1(C) = 3.
  - A strong tournament on 3 vertices has 1 odd cycle, and on 4 vertices it has 2.
  - On n >= 5 vertices it has at least m(5) = 4.
  - This is THM-343's S5 proof (Steps 2-4), which really proves that three pairwise intersecting odd cycles force a
    fourth. THM-029 records the same repaired argument. THM-343's older "Key Lemma ... open for n >= 8" paragraph is
    stale and is superseded by its own S5 proof, so cite the proof, Steps 2-4.
- **x = 1 and x = 0.**
  - At x = 1, C3 => C3 has Omega = 2K_1, so I(Omega, 1) = 4 = 1 + 3. The gap fails, as stated.
  - At x = 0, I(Omega, 0) = 1 = 1 + 3*0 for every T. The gap fails trivially.
  - The annotation of HYP-1147 / HYP-1109 ("ALL x >= 0" is false at x = 0 and x = 1) is accurate.

**(4) The 21-hole set {4, 8}.** Correct.

- **b >= 11.** An alpha-sum of 10 < b forces alpha = (10, 0), i.e. an intercyclic tournament with H = 21. That is
  excluded at all orders by THM-1370, which proves "h = 21 is never attained" for every n (checked: proof steps 1-5).
  Equivalently, 10 is not in A (Session 16). The same argument already covers b = 9, 10, since (1, 1) and (0, 1) fail
  the chain constraint.
- **2 <= b <= 10.**
  - Every realisation has alpha_1 <= 10.
  - Each strong component therefore has alpha_1 <= 10. Since m(8) = 6 + 4 + 2 = 12 > 10, it has order <= 7.
  - The packing vectors of all tournaments with alpha_1 <= 10 are the products of strong vectors with alpha_1 <= 10.
    The audit recomputes these as 108 vectors, the same count as the session.
  - Only (6, 1) at b = 4 and (2, 1) at b = 8 are realised (`aud_spectra.out`).
- **Completeness argument.** Valid.
- **The decomposition list in the statement is incomplete for b = 2 (m3).** The seven chain-compatible decompositions
  are (10, 0), (8, 1), (6, 2), (4, 3), (2, 4), (6, 0, 1) and (4, 1, 1); the statement lists four. All seven are
  unrealisable, so the conclusion stands.
- **Witnesses.**
  - The unique strong 6-class with alpha = (6, 1) has out-masks (2, 4, 9, 35, 15, 23), scores (1, 1, 2, 3, 4, 4),
    H = 17 and H_4 = 41, and is self-converse.
  - C3 => C3 has alpha = (2, 1) and H_8 = 81.

**(5) Repunit criterion.** Correct.

- 111_b = 1 + b(b + 1).
- alpha_3 = 0 because b^2 > b + 1, and alpha_2 <= 1 because 2b > b + 1.
- alpha_2 = 1 would force alpha_1 = 1, which the chain constraint forbids (THM-338 for k = 2; HYP-1015 is the (1, 1)
  case).
- So alpha = (b + 1, 0), which is realisable iff b + 1 is in A.
- **Recomputed lists.**
  - Failing bases: 2, 9, 13, 16, 19, 20, 26, 27, 28, 32, 34, 36, 37, 40, 43, 45, 49, ...
  - Failing repunit primes with b <= 100: 7 (b = 2), 421 (20), 757 (27), 2971 (54), 4423 (66), 6481 (80),
    8191 (90). The list is complete through b = 100.
  - The bases with 111_b prime and b <= 100 are 2, 3, 5, 6, 8, 12, 14, 15, 17, 20, 21, 24, 27, 33, 38, 41, 50, 54,
    57, 59, 62, 66, 69, 71, 75, 77, 78, 80, 89, 90, 99. Here 81 and 91 (b = 80, 90) are holes of A, while 79, 90 and
    100 (b = 78, 89, 99) lie in A.
- **Density.** Session 16's bound (4) gives at least (1/8 - o(1)) 2^(L-1) holes of A below 2^(L-1). So the failing
  bases have positive upper density, as stated.

**(6) Gap tables: completeness.** Correct.

- **The three ingredients:**
  - a component of a realisation has alpha_1 <= m <= 80;
  - Moon's sharp count bounds its order: m(19) = 17 + 15 + ... + 1 = 81 > 80;
  - extension completeness follows from Moon's deletion lemma together with alpha_1(T - v) <= alpha_1(T).
- **Citation checked.** J. W. Moon, "On subtournaments of a tournament", Canad. Math. Bull. 9 (1966) 297-301,
  Corollary 2.1: the minimum number of k-cycles in a strong n-tournament is n - k + 1. THM-343 cites the same
  corollary.
- **The deletion lemma** follows from vertex-pancyclicity (Moon 1966): an (n - 1)-cycle spans a strong
  (n - 1)-subtournament.
- **Floor values.** m(8) = 12, m(9) = 16, m(10) = 20 and m(19) = 81 are all correct.
- **The audit's checks:**
  - On the census, the minimum number of l-cycles over strong n-classes equals n - l + 1 for every odd l <= n <= 9.
  - Over the inventory, the least alpha_1 at each order 3-18 equals m(n), so the floor is attained.
  - Most importantly, the audit's own extension step finds no strong 19-class with alpha_1 <= 80. Completeness
    therefore also follows from the deletion lemma alone (UPGRADE U1).

**(7) Only bases 2 and 3 can have gap set {3, 10}.** Correct.

- For b >= 15, the gap 14 follows from (2).
- For 4 <= b <= 14, it is read off the complete tables of (6), which the audit reproduced. Every base in that range
  contains 14, 17 or 27, or (b = 4, 8) misses 10.
- It is PROVED in the sense of a complete finite computation (m5).

**(8) Carry-free singletons.**

- **Completeness is correct.** For a strong intercyclic n-tournament, alpha_1 >= 2^(n-4):
  - One-apex reduced words start with 1 and end with 0. The pair (first, last) alone contributes c(n - 2) = 2^(n-4).
  - The two-feedback templates of Session 16 (1^a20, 2^a10, 2^a1 2^b, 2^a3^b0, their duals, 310 and 320) have
    alpha_1 = 3 * 2^(n-4) or 9 * 2^(n-5) - 2 at order n. The audit checked each template's order-to-length relation.
  - So alpha_1 < 1024 forces n <= 13.
  - Independently, the audit's intercyclic enumerator capped at alpha_1 <= 1023 finds order 14 EMPTY.
- **The 26 values, their self-converse status and A below 1024 all check out.**
- **The binary half has a gap: MAJOR M1.**

## 2. MAJOR findings

### M1. Statement (8) claims more than its proof gives (title and status of THM-4940)

**Problem.**

- The title says: "in base 2 carries leave only those [carry-free singletons] with a <= 7, which are H = 3, 5, 11, 13".
  The status block types all of (8) "FINITE-EXACT ... below 1024".
- The proof establishes the base-2 half only for a < 312, by citing THM-4884's census below H = 625.
- For the twelve singletons 323 <= a <= 933 (H = 647, 739, 829, 859, 1089, 1209, 1361, 1399, 1493, 1677, 1719, 1867),
  no second strong realiser is exhibited.
- The body sentence is correctly scoped ("in the binary range below 625"). The title and status are not.

**The claim is true.**

- Every one of the 22 carry-free singletons with 8 <= a < 1024 has a strong tournament with H = 1 + 2a and two
  vertex-disjoint cyclic triangles. Witnesses for the last twelve have order 9; see `aud_compare.out`, section 4,
  which lists the counts by order and the least-order alpha vectors.
- The fix therefore costs nothing.

**Replacement text.**

- **Title (8):** "(8) Carry-free singletons: below 1024 exactly 26 values of A have a unique strong intercyclic
  realiser; in base 2 each of the 22 with a >= 8 acquires a second, non-intercyclic strong realiser of H = 1 + 2a
  (THM-4884 below 625; explicit order-9 witnesses for 323 <= a <= 933, audit), so among them only a <= 7, i.e.
  H = 3, 5, 11, 13, stay binary singletons".
- **Status block:** replace "(8) below 1024 (by Session 16's classification, ...)" with "(8) below 1024 (by Session
  16's classification, intercyclic strong tournaments with alpha_1 < 1024 have order <= 13; the base-2 twins are
  THM-4884's census below H = 625 and explicit order-9 witnesses for 625 < H < 2048, audit/aud_compare.out)".
- **Body, last bullet of (8):** "So carries give every carry-free singleton with 8 <= a < 1024 a second strong
  realiser in base 2: for a < 312 by THM-4884's census, and for the twelve values 323 <= a <= 933 by explicit strong
  order-9 tournaments with two disjoint cyclic triangles (audit)."
- **Alternative if the audit data is not adopted:** title "... ; in base 2, below H = 625, carries leave only those
  with a <= 7, which are H = 3, 5, 11, 13".

### M2. README section 5 types a false statement PROVED

**Problem.**

- The row "three pairwise intersecting odd cycles force a fourth" gives its role as "alpha_1 = 3 is impossible" and
  types it PROVED (THM-343; Session 16).
- That is false. alpha_1 = 3 occurs: C3 => (the strong 4-tournament) has alpha = (3, 2) and H = 15.
- THM-029 was corrected for exactly this claim in March 2026 ("alpha_1 = 3 IS achievable").
- What THM-343 excludes is alpha = (3, 0).

**Replacement text** (role column): "alpha = (3, 0) is impossible (three odd cycles, pairwise intersecting); alpha_1 = 3
itself occurs, e.g. C3 => strong 4-tournament has alpha = (3, 2) (THM-029 correction)".

## 3. MINOR findings (with replacement text)

### m1. README section 2, "Binary statement" is false as written

Even repunits, such as 1111_3 = 40 and 1111_5 = 156, are never H-values (Redei).

**Replace with:** "Within the verified coverage (2^32 by THM-4875 (D); 2^35 by the HYP-9478 addendum), the only *odd*
repunits with at least three digits in some base that are not H-values are 7 = 111_2 and 21 = 111_4 (even repunits
such as 1111_3 = 40 are never H-values, by Redei). This is a restatement of Chan-Pak's two exceptions, not a new fact."

### m2. Open completeness stated as fact (the MISTAKE-559 rule: cite canon at the strength of the file)

Chan-Pak (binary completeness) and HYP-9477 (A) (binary singletons only 3, 5, 11, 13) are OPEN. The following
sentences state them as facts. THM-4940's own "NOT CLAIMED" block is correct.

**(a) README section 0, item 1, last bullet.**

- Current: "So 7 and 21 are exactly the two holes of the carry-free layer that base-2 carries cannot fill, and base 3
  (pre-registered, confirmed) has the same two."
- Replace with: "So 7 and 21 are the two holes of the carry-free layer that base-2 carries provably cannot fill. Every
  other hole a < 2^34 is filled in base 2 (all-odd coverage to 2^35), and every hole is filled if Chan-Pak holds. Base
  3 has the same two unfillable holes and no other gap up to alpha-sum 80 (pre-registered P1, confirmed)."

**(b) README section 0, item 3.**

- First change: "The carry-free layer has singletons up to 933 and beyond (26 below 1024)." becomes "The carry-free
  layer has 26 singletons below 1024, the largest 933 (64 below 8192 by the audit's extension)."
- Second change: "In base 2, carries both fill holes (keeping only the gaps 7, 21) and create twins (keeping only the
  singletons 3, 5, 11, 13). Both exceptional sets live in the carry-free range a <= 10."
- Replace that with: "In base 2, carries both fill holes (only the gaps 7, 21 are known, and Chan-Pak says they are
  the only ones) and create twins (only the singletons 3, 5, 11, 13 are known below 625, and HYP-9477 (A) says they
  are the only ones). The known exceptional values of both kinds lie in the carry-free range a <= 10."

**(c) README section 1, rescue sets, last bullet.** Replace with: "So, up to a = 80, the carries rescue every hole of
the carry-free layer except 3 and 10 in both bases 2 and 3."

**(d) README section 3, "Binary carries".**

- Replace "That is why the binary list stops at 13." with "That is why the binary list below 625 stops at 13."
- Add: "(The audit extends the twin statement to every carry-free singleton with 8 <= a < 8192.)"

**(e) README section 3.** Replace "The binary exceptions of both problems live in the carry-free range a <= 10" with
"The known binary exceptions of both problems (all of them, if Chan-Pak and HYP-9477 (A) hold) live in the carry-free
range a <= 10".

**(f) HYP-9530 (A), "Comparison with base 2".**

- Replace "Base 2 (Chan-Pak) and base 3 have the same gap alpha-sums {3, 10}." with "Base 2 and base 3 have the same
  known gap alpha-sums {3, 10}: base 2 to alpha-sum 2^34 (all-odd coverage), base 3 to alpha-sum 80. Chan-Pak and (A)
  assert this at every alpha-sum."

**(g) HYP-9530 (C), "Reading".**

- Replace "carries both fill holes (Chan-Pak keeps only 3, 10) and destroy uniqueness (Carmichael-type singletons keep
  only a <= 7)" with "carries both fill holes (Chan-Pak, OPEN: only 3 and 10 survive) and destroy uniqueness
  (HYP-9477 (A), OPEN: only a <= 7 survive; verified below 625, and for the carry-free singletons below 8192 by the
  audit)".

**(h) README section 4.2.**

- Replace "first creators of new values are strong (HYP-9481 (b)), so they are produced by extensions, not joins" with
  "first creators of new values are strong through order 11 (THM-4922, THM-4942; in general HYP-9481 (b), OPEN), so
  there they are produced by extensions, not joins".

### m3. THM-4940 (4) and README section 1: the decomposition lists for alpha-sum 10 are incomplete

**THM-4940 (4).**

- Current: "Every other base has every decomposition of the alpha-sum 10 unrealisable".
- Replace with: "In every other base, every decomposition of the alpha-sum 10 that satisfies the chain constraint is
  unrealisable (the list omits those already excluded by alpha_2 <= C(alpha_1, 2) or by alpha_3 >= 1 =>
  alpha_2 >= 3, namely (2, 4), (6, 0, 1), (4, 1, 1) at b = 2 and (2, 2) at b = 4)".

**README section 1.**

- Current: "The vectors are (alpha_1, alpha_2) and respect the chain constraint."
- The listed b = 2 row omits (2, 4), which also respects that constraint.
- Replace with: "The vectors are (alpha_1, alpha_2). Listed are those that also satisfy alpha_2 <= C(alpha_1, 2) and
  alpha_3 >= 1 => alpha_2 >= 3; `base_spectra.out` lists all seven chain-compatible vectors at b = 2."

### m4. THM-4940 (2), wording

Replace "So Spec_b and A agree on [1, b^2)." with "So Spec_b ∩ [1, b^2) = {1 + a b : a ∈ A, a < b}."

### m5. THM-4940 status, (7)

Replace "PROVED: (1), (2), (3) (all orders), (4) ..., (5), (7)." with "... (5), (7) (for b >= 15 by (2); for
4 <= b <= 14 by the complete finite tables of (6))."

### m6. Credits and prior art: missing or imprecise

No NEW item is invalidated. The all-order proofs, the base set {4, 8}, the tables, the repunit criterion and the
singleton list are new in the repository and, as far as a brief literature search shows, outside it. But the following
earlier work should be credited.

**Add to "PRIOR OBSERVATIONS credited":**

- "`04-computation/ip_hierarchy.py` (opus-S71e, March 2026; the script behind HYP-1147) also found I(Omega, 3) = 31
  absent at n <= 6 and asked whether 10 and 31 are 'gaps analogous to H = 7 and H = 21'. That is the ternary pair of
  HYP-9530 (A) and of (4) at b = 3, observed at n <= 6."
- "`07-reflections/cubic-family-x6-coloring.md` (oracle-2026-05-16-S2): I(Omega, 6) = 19 = 1 + 6*3 is never attained,
  by the same alpha = (3, 0) obstruction. That is the b = 6 case of (3) and the block-lemma mechanism. For n <= 5 the
  note also maps the gaps 7, 21 to 19, 61 by H -> 3H - 2."
- "`07-reflections/fugacity-axis-five-points.md` (oracle-2026-05-17-S1): the degree-1 alpha_1 gaps
  {3, 10, 14, 17, 20, 21, ...} at n = 7 (5,000 samples), an early sampled view of the holes of A."
- "HYP-886 (kind-pasteur-S62, x = 3) and HYP-892 (opus-S72, every x): I(Omega, x) = 1 (mod x)."
- "THM-467 (d)1: I(Omega(T), x) is the product over the strong components for every x. These are the first and third
  items of (1)."
- "HYP-1144 = HYP-1106 (opus-S71e): the base-2 row of (4) ('all six decompositions of T = 10 are blocked', sampled at
  n = 7; proved by THM-1370)."

**Fix two existing entries:**

- "THM-467: generating function of I(Omega, t)" becomes "THM-467: I(Omega(T), t) = eps([det(I + XA)/det(I - XA)]^(t/2))
  for even t, and (d)1 the strong-component factorisation at every x".
- "07-reflections/the-tournament-alphabet.md: H_q as ANALOGY" becomes "07-reflections/the-tournament-alphabet.md
  defines the same object H_q(T) = I(Omega(T), q) (the 'q-tournament') and notes that it is odd iff q is even. Its
  readings are ANALOGY."

**README section 5:**

- Credit kind-pasteur S17's concurrent answer to the same directive item, THM-4942's note section 6 "the 7/21 proof
  objects are the generation operations (typed analogy)", pushed at 19:49, before this session's pre-registration.
- In the projective-carrier row, the join law a + a' + 2aa' on (H - 1)/2 predates Session 210. It appears as "a star b"
  in `paired_prime_collatz_20261003.md` and `tiling_modular_atoms_20261003.md`, and Session142's jacobian bridge uses
  q = (h - 1)/2 with forbidden values 3, 10.

**Search scope.**

- **Repository.** git grep with untruncated listings for I(Omega/Ω/CG, 3); "moving gap"; 1 + 10x / K10; 111_b and q^2 +
  q + 1; Sylvester; Klee; the 7/21 dictionary terms.
- **Literature.** Searched for evaluations of the tournament odd-cycle independence polynomial at x != 2 and for
  Chan-Pak. Nothing earlier was found for H_b with b != 2, for the set {4, 8}, or for the repunit criterion.
- **Unverified.** "Ford proved Sierpinski's conjecture for all k" (README) versus "Ford proved it for odd k"
  (CARMICHAEL_MAP) is a difference of wording only.

### m7. README section 3 overclaims Session218 and Session219

Replace "Session218 and Session219 prove that no gap-preserving map from totient fibers to tournament fibers exists."
with "Session218 and Session219 record that the analogy supplies no map: halving the nontotients 14 and 26 gives 7 and
13, but 13 is an H-value, and exact totient multiplicity is not a map to tournament-class multiplicity."

### m8. README section 3, the Sylvester sentence is false

Sylvester's sequence appears in at least 20 places outside the LRC/Heegner lane. Examples:

- THM-2005 (eq. 38);
- S33's `number_tournament_bridges_20261014.py`;
- `bernoulli_chain_s116n.py` (von Staudt-Clausen);
- `sigma_phi_uniqueness_s116n.py`, which builds the chain 2, 3, 7, 43 in a sigma/phi context;
- `all_telescopes_s90co.py` and others.

Replace "Sylvester's sequence appears in the repository only in the LRC/Heegner lane." with "Sylvester's sequence
already appears in several repository lanes: the covering-min and Heegner hypotheses HYP-3724, HYP-3730 and HYP-4079;
THM-2005; S33's number-tournament bridges; and the s90/s116 scripts, including `sigma_phi_uniqueness_s116n.py` with the
chain 2, 3, 7, 43. None attaches a tournament mechanism."

### m9. README section 4 is a stale merge of the concurrent kind-pasteur S17 lane

**What happened.** S17 pushed THM-4942 / HYP-9532 checkpoints at 19:49 (14487ef208) and 19:57 (f870c2abf1). S86 then
pushed its pre-registration (20:08) and its checkpoint (20:27). The README still describes S17 as "live and
exploratory" with only "R1-R11 pre-registered". This is the MISTAKE-677 pattern: re-read concurrent READMEs.

**Replace the section 4 introduction sentence and the "Extensions" bullets with:** "kind-pasteur-2026-11-08-S17
(THM-4942 / HYP-9532, pre-audit; checkpoints 14487ef208, f870c2abf1): order 11 pre-registered R1-R5 all CONFIRMED. Each
of the 36,248 values new at order 11 is a one-vertex extension of an order-10 first creator, so 63 stays the only
orphan value. Through order 11 every tournament is a first creator, a first creator plus a vertex, an ordinal sum, a
substitution, or an ordinal sum re-connected by at most two arcs."

### m10. The fibre-move counts depend on the truncation at order 8

**Facts.**

- The counts 3/56, 34/456 and 632/6,880 are components of the move graph restricted to orders <= 8. They can only fall
  when higher orders are added.
- The README's qualitative "Consequence" is nevertheless true at all orders, by the H = 9 argument (UPGRADE U3).

**README section 0, item 4.** Replace "632 of the 6,880 order-8 classes are cut off." with "In the move graph on
orders <= 8, 632 of the 6,880 order-8 classes lie in no component with a first creator of their value. The H = 9 fibre
has exactly two components at every order (audit), so the obstruction is not an artefact of the truncation."

**README section 4.2, "Result" and "Consequence".** Add: "(components of the move graph on orders <= 8; the H = 9
fibre is disconnected at every order: see the audit's proof)".

### m11. README section 5: two rows to fix

**The floor row.**

- Current: "Busch/Moon floors f(n) = min 3^a 5^b ... PROVED (THM-1370)".
- The formula as written is incomplete, and the general formula is Busch's theorem, not THM-1370. THM-1370 proves only
  f(5..8) and monotonicity.
- Replace with: "Busch's floor f(3q) = 3*5^(q-1), f(3q+1) = 5^q, f(3q+2) = 9*5^(q-1) (Busch, EJC 13 (2006) N3,
  Cor. 1, CITED in THM-4744; f(5..8) exhaustive in THM-1370) | ... | KNOWN (Busch); the gap picture is ANALOGY".

**The monoid row.**

- The irreducible list is correct and matches THM-4094's lanes; the audit re-derived it.
- But the value monoid is not factorial, so "unique factorisation in a free monoid" describes the tournament side
  only.
- Append: "Unique factorisation holds for tournaments (the free monoid of strong classes under ordinal sum, THM-1862).
  H maps it onto the odd numbers without 7 and 21, a monoid that is atomic but not factorial: 35 * 77 = 49 * 5 * 11
  and 63^2 = 3^4 * 49."

### m12. INDEX row wording (05-knowledge/hypotheses/INDEX.md, line 2)

Replace "(fails b = 2, 9, 13, 16, ...; 91 = 111_9 is 21's repunit, not 111_4)" with "(fails b = 2, 9, 13, 16, ...;
the repunit failure tied to 21's alpha-sum 10 is 91 = 111_9, while 21 = 111_4 is a base-4 value)".

### m13. "Independent" confirmation in README section 1

- Current: "an independent exhaustive confirmation of Session 16's spectrum below 1024. Two different methods were
  used".
- The two generators share one canonical-form routine, adapted from `classgen.c`.
- Replace with: "an exhaustive confirmation of Session 16's spectrum below 1024 by two generators (sharing one
  canonical-form routine); the audit's independent enumerator confirms it below 8192".

### m14. Moon citation in `base_spectra.py` docstring (cosmetic)

"MISTAKE-560 / THM-115" should read "Moon, Canad. Math. Bull. 9 (1966) 297-301, Cor. 2.1 (as cited in THM-343;
MISTAKE-560)". THM-115 uses the weaker pancyclicity count.

### m15. HYP-9530 (C) status line

After the audit's extension, it can read: "(C): FINITE-EXACT below 8192 (64 carry-free singletons, all self-converse;
every one with a >= 8 has a strong non-intercyclic binary twin; audit). The binary half beyond that is HYP-9477 (A)."

## 4. UPGRADES (proved or verified by the audit)

**U1. Completeness without Moon's cycle floor.**

- The audit's enumerator extends level by level until a level is EMPTY: order 19 at budget 80, order 14 at
  intercyclic cap 1023, order 17 at cap 8191.
- By the deletion lemma and monotonicity (or heredity) alone, no larger order contributes.
- Suggested sentence for (6): "Independently, the extension step finds no strong 19-class with alpha_1 <= 80, so
  completeness also follows from Moon's deletion lemma alone."

**U2. Carried twins for every carry-free singleton with 8 <= a < 8192.** FINITE-EXACT: the singleton list is complete
and every witness is explicit and re-verified.

- **The 22 singletons below 1024:** witnesses at orders 6-9 come from the census (`aud_compare.out`, section 4). This
  closes M1.
- **The 38 singletons in [1024, 8192):** found by randomised search at order 10, plus order 11 for a = 7845 (H = 15691
  is near max H(10) = 15745). Each witness was re-verified in C and independently in Python: H by subset DP, strong,
  two disjoint cyclic triangles (`aud_twins_hi.out`).
- **Consequence:** HYP-9530 (C), second half, holds for every carry-free singleton below 8192.

**U3. The H = 9 fibre has exactly two components at every order.**

- **Setup.** Every tournament with H = 9 is of one of two types:
  - (i) a padding of one of the two strong 5-classes with H = 9 (f(6) = 15);
  - (ii) a padding of two C3 components.
- **H-preserving deletions** remove only singleton components (Session133: H(T) - H(T - v) >= 2 otherwise), so they
  preserve the type.
- **A flip from type (ii) to type (i) is impossible.** Let X and Y be the two C3 components of T, and let T' be T with
  one arc reversed.
  - **The flipped arc lies inside neither X nor Y.** Then X and Y both remain cycles of T'. They cannot share a strong
    component: such a component would have at least 6 vertices, hence H >= f(6) = 15 > 9. So they lie in two
    nontrivial components, and the result is type (ii).
  - **The flipped arc lies inside X** (the case of Y is symmetric). Y stays a cycle, so it lies in the 5-component D of
    a type-(i) result.
  - Every vertex outside Y beats all of Y or loses to all of Y: this holds in T because Y is a component, and the
    arcs at Y are unchanged.
  - So D strong needs w, w' in D \ Y with w -> Y -> w', and w' -> w in T'.
  - In T the arc was w -> w', since otherwise Y ∪ {w, w'} would be strong in T. So the flipped arc is {w, w'}, which
    lies inside X.
  - But X is a component of T lying entirely before or entirely after Y. That contradicts w -> Y -> w'.
- **Conclusion.** The flip direction (i) -> (ii) is the reverse move, so it is impossible too. Within each type,
  deletions reach the creators.
- This proves the README's "Consequence" (same-H elementary moves cannot generate fibres from first creators) at all
  orders.

**U4. The carry-free layer below 8192.**

- A below 8192 has 4,207 values, equal to Session 16's formula (one-apex values by an independent chain-path DP, plus
  the rays).
- There are 64 singletons, all self-converse: 26 below 1024 and 38 in [1024, 8192).
- They sit at orders 3, 4, 5 (two), 6, 7, 8, 9 (three), 11 (eight), 13 (17) and 15 (29). Above order 9 they occur only
  at odd orders.
- This supports HYP-9530 (C), first half.

**U5. Moon's sharp count is attained.**

- On the census (n <= 9), the minimum number of l-cycles over strong n-classes is n - l + 1 for every odd l.
- In the inventory, the minimum alpha_1 is m(n) at every order 3-18.

**U6. The corrected S/C ledger (README section 4.1).**

- Recomputed from the census for all 77 values born at orders <= 7. It is identical to Session224's
  `birth_fibers.csv`.
- Every entry quoted for orders 4-6 is correct.
- The claims about 35, 39, 61 (born at order 7), 41 (one class), 25 and 27 (one converse pair each) and 11 (two
  non-strong classes at order 6, a converse pair) are all correct.
- "124 of 7,663" matches kind-pasteur S16's `aud_fibers_le10.txt`.

**Out of scope at first, noted.** The post-checkpoint `ternary_attainment.out` claims constructive attainment of every
value congruent to 1 (mod 3) up to 5398 except 10 and 31. Combined with THM-4940 (3) and (4), that would make Gap_3
exact to alpha-sum 1799. This is now audited in the section "Follow-up: post-audit commits" at the end.

## 5. Independent computations

| # | universe | method (audit code) | result | agrees with |
|---|---|---|---|---|
| 1 | all tournament classes of orders 1-9 (198,948; classgen lists) | own packing vectors (explicit odd-cycle list), own canonical form, own self-converse test | OCF H = sum alpha_k 2^k equals classgen's DP H on every class; class counts = A000568; strong counts = A051337; own codes injective; strong intercyclic counts at orders 3-9: 1, 1, 6, 11, 20, 37, 70 | README section 6.1 and 6.4 |
| 2 | strong classes, n <= 9 | census | min #l-cycles = n - l + 1 (all odd l); min alpha_1 = m(n) = 1, 2, 4, 6, 9, 12, 16 | Moon (1966) Cor. 2.1 |
| 3 | strong tournaments with alpha_1 <= 80, all orders | own one-vertex extension enumerator; order 19 EMPTY | 5,324 classes; per order 1, 1, 6, 35, 353, 1480, 1115, 838, 625, 401, 245, 131, 57, 25, 8, 3 (orders 3-18); one class with alpha_6 = 1 (order 18, H = 25281) | `few_cycle_strong_80.txt` (identical multiset of (order, alpha, H)); `few_cycle_strong_44.txt` (901, identical); census code sets at orders <= 9 identical |
| 4 | gap tables, m <= 80, bases 2-90 and 100, 127, 128, 500, 1000, 10^6 | own multiplicative knapsack over strong values | THM-4940 (6) exactly; bases >= 81 give the 31 holes of A; bases 70, 73, 76, 77, 79, 80 also have no rescue below 81 | `base_spectra_80.out` |
| 5 | rescue sets of the holes a <= 80 | same | 3: none; 10: {4, 8}; 14: {2, 3, 5, 8, 12}; 17: {2, 3, 4, 5, 6, 7, 8, 11, 15}; every other hole rescued in 2 and 3 | THM-4940 (6), README section 1 |
| 6 | decompositions of the alpha-sums 3 and 10 | all products of strong vectors with alpha_1 <= 10 (108 vectors) | only (6, 1) at b = 4 and (2, 1) at b = 8 realised; 7 chain-compatible decompositions at b = 2 | THM-4940 (4) (see m3) |
| 7 | strong intercyclic, alpha_1 <= 1023 | own enumerator; order 14 EMPTY | 1,461 classes; A below 1024 = Session 16 formula (565 values, 459 holes); the 26 singletons as listed, all self-converse | `intercyclic_strong_le13.txt` (identical multiset of (order, alpha_1, self-converse)); `carry_free_layer.out` |
| 8 | strong intercyclic, alpha_1 <= 8191 | same; order 17 EMPTY | 11,433 classes; A below 8192 = Session 16 formula (4,207 values); 64 singletons, all self-converse | UPGRADE U4 |
| 9 | carried binary twins | census n <= 9 (a < 1024); randomised search n = 10, 11 with re-verification in C and Python (1024 <= a < 8192) | all 60 singletons with 8 <= a < 8192 have a strong non-intercyclic tournament with H = 1 + 2a | closes M1; UPGRADE U2 |
| 10 | one-apex words | brute-force cycle enumeration (length <= 9) vs chain-path DP (length <= 16) | c(d) = 1, 2^(d-2) confirmed; min alpha_1 per reduced length 1, 2, 4, 7, 12, 22, 40, 76, 144, 280, 544, 1072, 2112, 4192, 8320; 3,992 values below 8192 not attained | `single_apex_spectrum.out` |
| 11 | repunits | tables (b <= 79), A below 1024 | criterion agrees for 2 <= b <= 79; failing bases and failing repunit primes as in section 1 (5) | THM-4940 (5) |
| 12 | witnesses | inventory | (6, 1): out-masks (2, 4, 9, 35, 15, 23), H = 17, H_4 = 41, self-converse, unique; C3 => C3: H_8 = 81, I(Omega, 1) = 4 | THM-4940 (3), (4) |
| 13 | same-H move graph, orders <= 8 | own canonical form on classgen lists | 1,272 deletion edges and 6,791 flip edges; classes cut off 3/56, 34/456, 632/6,880; H = 9 fibre: 2 components (sizes per order 2, 4, 6, 8 and 1, 3, 6) | `fibre_moves_le8.out` |
| 14 | birth fibers with n(h) <= 7 | census | 77 values identical to Session224's `birth_fibers.csv` | README section 4.1 |
| 15 | canonical-form invariance | 20 random relabellings of each of the 5,324 + 1,461 classes | 0 failures | none (internal check) |
| 16 | irreducibles of the odd numbers without 7 and 21 | hand derivation (cases 7^e m) | odd primes other than 7; 7p with p not 3 or 7; 49, 63, 343; not factorial | README section 5, THM-4094 |

## 6. Scripts (all in this folder) and how to run them

- `aud_fewcycle.c`: the audit's enumerator and checker. Its modes:
  - `aud_fewcycle M 0 HCHK`: every strong class with alpha_1 <= M;
  - `aud_fewcycle M 1 HCHK`: every strong intercyclic class with alpha_1 <= M;
  - `aud_fewcycle 0 2 n cls_n.txt`: census mode, recomputing a classgen list;
  - `aud_fewcycle 0 3 R file`: canonical-form invariance under R random relabellings;
  - `aud_fewcycle 0 4 N dir`: the same-H move graph on orders <= N.
- `aud_spectra.py FC IC`: gap tables, rescue sets, decompositions, A against Session 16's formula, singletons,
  repunits, witnesses. Writes `aud_spectra.out`.
- `aud_compare.py WORK SESSIONDIR S224CSV`: comparison with the census and the session's output files, carried twins
  below 1024, the extended intercyclic run, birth fibers. Writes `aud_compare.out`.
- `aud_single_apex.py`: one-apex brute force against the DP, minima, non-attained values. Writes
  `aud_single_apex.out`.
- `aud_twins_search.c` and `aud_singles_hi.py`: carried-twin search for the singletons in [1024, 8192). Writes
  `aud_twins_hi.out`.
- `aud_followup.py WORK SESSIONDIR`: the follow-up audit of the post-checkpoint commits (see the last section). Writes
  `aud_followup.out`.
- **Logs:** `aud_census.log`, `aud_fewcycle_80.log`, `aud_intercyclic_1023.log`, `aud_intercyclic_8191.log`,
  `aud_invariance.out`, `aud_fibre_le8.out`.
- **Run everything:** `bash run_audit.sh WORKDIR` takes about 2 minutes (gcc and Python 3, standard library only). It
  rebuilds kind-pasteur's `classgen.c`, writes the large inventories to WORKDIR and the summaries here. A rerun
  reproduced every `.out` file byte for byte.

## 7. Pre-registration ordering (confirmed)

- **Push order.** The remote-tracking reflog shows 18463fd3cf (`PREREGISTRATION.md` with P1-P5, the only file in the
  commit) pushed at 20:08:06.
- **Later file times.** The budget-80 inventory `few_cycle_strong_80.txt` was last written at 20:12:29 and
  `base_spectra_80.out` at 20:12:51, both after the push.
- **The fix commit.** `few_cycle_strong.c` was modified at 20:08:53, which is the 192-bit fix and rerun the outcome
  describes.
- **Append-only outcome.** The outcome section was added in 725d3bb4cf (pushed 20:27:40). Its diff against 18463fd3cf
  only appends, so the predictions are unedited.
- **One disclosed post-run change.** The alpha_6 loader fix in `base_spectra.py` (20:12:51) came after the run. It is
  disclosed in README section 6 and does not touch the predictions.
- **The in-sample table** (gaps to 44) matches the audit's tables. P1-P5 are CONFIRMED as stated. P2, P4 and P5 are
  consistency checks, as the file says.

## Follow-up: post-audit commits

Scope: the three commits the first audit did not review, at the coordinator's request.

- ebe3096762 adds `PREREGISTRATION-2.md`, with P6-P8.
- 1a1606cc66 adds `coverage_products.py` and `.out` and appends the outcome to `PREREGISTRATION-2.md`.
- 1d57e6b23c adds `ternary_attainment.py` and `.out`.

The audit's own script is `aud_followup.py`, with output `aud_followup.out`. Its inputs are the audit's census, its
alpha_1 <= 80 inventory and its exhaustive intercyclic inventory below 8192. It reads the session's `.out` files only as
data.

### Verdict: PASS WITH CORRECTIONS (minor)

**The proposed statement (9) is correct, and the audit reproduced it independently.**

- **Base 3 is even stronger than claimed.** No products are needed: every value congruent to 1 (mod 3) up to 5398,
  except 10 and 31, is already H_3 of a tournament of order <= 9.
- **The ordering is clean.**
- **No gap is claimed anywhere in these commits.**
- **The printed witnesses are genuine.**

**Two things need fixing:**

- A truncation bug in `ternary_attainment.py`. It does not touch (9), but it falsifies one reported count (1698 should
  be 1701).
- One stale sentence in the appended outcome of `PREREGISTRATION-2.md`.

### (i) Ordering

- **The pre-registration came first.** ebe3096762 contains only `PREREGISTRATION-2.md`. It was committed at 20:29:03
  and pushed at 20:29:09, according to the remote-tracking reflog.
- **The coverage script did not exist at that commit.** `coverage_products.py` is absent at ebe3096762 (checked with
  `git cat-file`). It first appears in 1a1606cc66, committed at 20:30:30 and pushed at 20:30:36. Its file time, 20:29:35,
  and that of its output, 20:29:41, are both after the pre-registration push.
- **The outcome only appends.** Between ebe3096762 and 1a1606cc66 the diff of `PREREGISTRATION-2.md` is pure appending,
  so the predictions are unedited.
- **The attainment script came last.** `ternary_attainment.py` and `.out` first appear in 1d57e6b23c, committed at
  20:31:26 and pushed at 20:31:33. That is after the outcome. Its docstring correctly calls the run "an exploratory
  calibration ... not a pre-registered test".
- **One deviation, which makes no difference.**
  - The pre-registered universe uses one-apex reduced words of length <= 23; the script uses length <= 21.
  - The one-apex values up to 100,000 are the same for both bounds (audit check).
  - So P6-P8 tested exactly the registered universe.

### (ii) Independent verification of the base-3 range

**Products of the audit's own pieces.**

- **The pieces:**
  - all 184,537 strong classes of orders 3-9 (census), with 64,196 distinct packing vectors;
  - all 5,324 strong classes with alpha_1 <= 80 (2,072 distinct vectors);
  - all 11,433 strong intercyclic classes with alpha_1 <= 8191 (4,206 distinct values of alpha_1), a complete set,
    since order 17 is empty.
- **The closure** is a divisor-pair closure, a different algorithm from the session's forward multiplication.
- **Base 3.** Every value 1 + 3m with m < 1800, other than m = 3 and 10, is attained. The first uncovered alpha-sum
  other than 3 and 10 is m = 1800, the value 5401.
- **Base 2.** Every odd value up to 2925 except 7 and 21 is attained. The first uncovered alpha-sum is m = 1463, the
  value 2927.
- Both thresholds are identical to `ternary_attainment.out`.

**No products at all.**

- H_3 evaluated directly on all 198,948 classes of orders <= 9 attains every value congruent to 1 (mod 3) up to 5398
  except 10 and 31. The next misses are 5401, 5425, 5437, ...
- H_2 evaluated directly attains every odd value up to 2881 except 7 and 21. The next miss is 2883.
- Products of strong classes of order <= 9 alone reach alpha-sum 1799 in base 3, but only 1442 in base 2. The base-2
  range 2925 needs the pieces (b) and (c).

**The witness pieces were checked independently.** The pieces with H_3 = 3001 and H_3 = 5398 were re-derived in pure
Python (cycle enumeration) and agree with the C census.

- H_3 = 3001: alpha = (454, 164, 6), out-masks (136, 13, 25, 448, 43, 143, 55, 342, 119).
- H_3 = 5398: alpha = (755, 309, 13), out-masks (88, 385, 99, 102, 46, 323, 402, 61, 157).
- Their Hamiltonian-path counts, 1613 and 2851, were confirmed by a separate DP.

**Why 5401 is where the construction stops.**

- 5401 = 11 * 491, and both factors are congruent to 2 (mod 3). So 5401 has no factorisation into factors congruent to
  1 (mod 3).
- Any tournament with H_3 = 5401 would therefore have a single nontrivial strong component, with base-3 alpha-sum 1800.
- By the three complete inventories, that component would have order >= 10, alpha_1 > 80, and two disjoint odd cycles.
- Whether 5401 is in Spec_3 is OPEN.

### (iii) The witnesses, and the absence of any gap claim

**Every witness is a product of genuine pieces.**

- `ternary_attainment.out` prints eight witness factorisations: four in base 2 and four in base 3.
- Each one multiplies out correctly.
- Each factor is H_b of an explicit strong class in the audit's own census. Two examples: 43 = H_3 of an order-6 class
  with alpha = (11, 1), and 5398 = H_3 of an order-9 class with alpha = (734, 313, 14).
- No factor needs a truncated vector (see the bug below).

**Nothing claims a gap.**

- The `.out` reports "uncovered" alpha-sums only.
- The docstring says "Gaps are NOT proved here".
- No line asserts non-attainment.

**Bug: `ternary_attainment.py` truncates packing vectors.** It does not affect (9).

- **Where.** It builds its pieces with `{v[:4] for v in cp.load_inventory(inv_path)}`.
- **What is lost.** That drops alpha_5 and alpha_6 for the 39 inventory classes of orders 15-18 that have
  alpha_5 > 0 (23 distinct vectors). The truncated vectors are not packing vectors of those classes.
- **The fake values are too large to matter for (9).** The smallest is 3385 in base 2 and 10537 in base 3, above 2925
  and 5398.
- **The base-3 count is wrong.** The reported count to m = 10,000 in base 3 is 1698. The correct count for this
  universe is 1701.
  - The alpha-sums 4016, 7885, 8030 and 8876 are "covered" only through fake values.
  - The alpha-sum 7966 is covered only by a genuine alpha_5 > 0 piece and is lost.
- **The base-2 count is right only by cancellation.** It reports 1183, which is correct, but only because 3335, 3350
  and 5966 are gained through fake values while 3174, 6300 and 9523 are lost.
- **Fix.**
  - Use full vectors: replace `{v[:4] for v in cp.load_inventory(inv_path)}` with
    `{tuple(v) for v in cp.load_inventory(inv_path)}`, and pad the census vectors to the same length, e.g.
    `tuple(...) + (0, 0)`.
  - Rerun. Expect "base 3: uncovered alpha-sums m <= 10000: 1701".
  - The thresholds 1463 and 1800, and all printed witnesses, are unchanged.

### (iv) Wording of the `PREREGISTRATION-2.md` outcome

**Every number checks out (audit recomputation in the same universe):**

- P7: 18,038 uncovered alpha-sums in (80, 100000], the first being m = 179.
- P6: 24,059, the first being m = 340.
- P8, bases 4, 5, 6, 7, 8: 5,182, 5,613, 5,827, 6,561, 6,470 uncovered in (80, 20000], the first being 273, 318, 333,
  296, 351.
- Within m <= 80, the closure reproduces the exact gap tables in every base 2-8.
- The unrun in-sample analogue:
  - pieces with alpha_1 <= 20, plus the intercyclic families, leave 14 alpha-sums in (20, 80] uncovered in base 2
    (29, 33, ..., 80) and 14 in base 3 (35, 38, ..., 80);
  - pieces with alpha_1 <= 40 leave exactly m = 80, in each of the two bases.

**The "Why" paragraph is borne out.** In 1d57e6b23c, adding the strong classes of order <= 9 lifts coverage to
alpha-sum 1799 in base 3. Those classes have up to 954 odd cycles, and most have disjoint cycles.

**Suggested edits.** The outcome section is append-only, so the stale bullet is fixed by a dated addendum, not an edit.

- **Clarify one sentence.** Replace "Pieces with alpha_1 <= 40 leave one (m = 80)." with "Pieces with alpha_1 <= 40
  leave one, m = 80, in each of bases 2 and 3."
- **Supersede the stale bullet.** "Status of the ternary law. HYP-9530 (A) rests on the exact tables to alpha-sum 80
  (P1) and nothing more." was true when written and was superseded one minute later by 1d57e6b23c. Append:
  "**Addendum (after 1d57e6b23c; not part of the pre-registered outcome).** An exploratory, non-pre-registered
  constructive check (`ternary_attainment.py`, reproduced independently by the audit) attains every value
  = 1 (mod 3) up to 5398 except 10 and 31, already by tournaments of order <= 9. With THM-4940 (3) and (4) this makes
  Spec_3 ∩ [1, 5398] exact (THM-4940 (9)). Beyond 5398, HYP-9530 (A) remains OPEN, and P6 stands REFUTED as
  registered."

### Replacement text for the proposed THM-4940 (9)

**Statement.** "(9) **Constructive attainment** (FINITE-EXACT, explicit witnesses).

- Every value v = 1 (mod 3) with v <= 5398, other than 10 and 31, is H_3(T) of some tournament T of order <= 9 (the
  complete census of orders <= 9).
- Every odd v <= 2925 other than 7 and 21 is H(T) of an order join of explicit strong pieces. The pieces are the strong
  classes of order <= 9 or with alpha_1 <= 80, the strong one-apex tournaments, and Session 16's rays. (Every odd
  v <= 2881 except 7 and 21 is already attained at order <= 9; this base-2 half is a control, inside the known all-odd
  coverage to 2^35.)
- With (3) and (4): Spec_3 ∩ [1, 5398] = {v <= 5398 : v = 1 (mod 3)} minus {10, 31}, i.e.
  Gap_3 ∩ [0, 1799] = {3, 10}.
- The construction first fails at 5401 = 11 * 491. A tournament with H_3 = 5401 would have a single nontrivial strong
  component, of order >= 10, with alpha_1 > 80 and two disjoint odd cycles. Whether one exists is OPEN.
- No gap is claimed beyond (3), (4) and (6)."

**Status line.** "(9) FINITE-EXACT (explicit witnesses; exploratory, not pre-registered; reproduced independently by
the audit, which finds the base-3 range already attained at order <= 9)."

**Follow-on edits.**

- **HYP-9530 (A) status** becomes: "FINITE-EXACT and complete for alpha-sums <= 1799 (values <= 5398): the gaps there
  are exactly 10 and 31 (THM-4940 (3), (4), (9)); OPEN beyond."
- **HYP-9530 (A) falsifier.** "Every m <= 80 is already settled" becomes "Every m <= 1799 is already settled".

### Follow-up computations (added to the table of section 5)

| # | universe | method (audit code) | result | agrees with |
|---|---|---|---|---|
| F1 | products of all strong classes of order <= 9, all with alpha_1 <= 80, all intercyclic with alpha_1 <= 8191 | own pieces, divisor-pair closure | base 3: all alpha-sums m < 1800 except 3, 10 (values <= 5398); base 2: all m < 1463 except 3, 10 (values <= 2925) | `ternary_attainment.out` thresholds |
| F2 | all 198,948 classes of orders <= 9, no products | direct H_3, H_2 | every value = 1 (mod 3) up to 5398 except 10, 31; every odd value up to 2881 except 7, 21 | stronger form of (9), base 3 |
| F3 | the session's universe (a)-(c), correct vectors versus the `v[:4]` truncation | same closure | to m = 10000: base 3 has 1701 uncovered (truncated: 1698); base 2 has 1183 both ways, by cancellation | `ternary_attainment.out` (count wrong in base 3) |
| F4 | PREREGISTRATION-2 universe | same closure | P6, P7, P8 counts and first misses exactly as stated; in-sample analogue 14 / 14 and 1 / 1 | `coverage_products.out`, outcome text |
| F5 | witness pieces H_3 = 3001, 5398 | pure-Python cycle enumeration and path DP | alpha = (454, 164, 6) and (755, 309, 13); H = 1613 and 2851 | C census |
