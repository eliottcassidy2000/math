# Independent audit: kind-pasteur-2026-11-09-S19 (THM-5085, HYP-9655)

**Target.** Commit c322f6b577 (checkpoint, pre-audit):
- `01-canon/theorems/THM-5085-fixed-port-tails-closed-form-signed-offset-transport-and-coverage-to-a-million.md`
- `05-knowledge/hypotheses/HYP-9655-chan-pak-as-cofiniteness-of-finitely-many-fixed-port-tail-forms.md`
- this packet (README.md, PROOF.md, PREREGISTRATION.md, port_stats.c, fp_cover.c, fp_greedy.c, fp_check.py,
  fixedport.py, out/).

The later commit c5b81ee9eb (pre-registered addendum Q5-Q6) was also read. THM-5085, HYP-9655, README, PROOF and the
programs are unchanged in c5b81ee9eb and in the working tree at the time of writing. PREREGISTRATION gained only the
Q5-Q6 addendum, plus its outcomes, which were still uncommitted while this audit ran.

**Method.** Every number below was recomputed with new code in `audit/`. No S19 program was run. The class
generator, the counters, the two coverage engines, the Q4 programs and the brute-force checkers were all written from
scratch. Binaries were built outside the repository.

## Verdict: PASS WITH CORRECTIONS

**What is correct.**
- The mathematics of THM-5085 (1)-(3) is correct. I re-derived the Lemma 1 proof. The transfer, the closed form, the
  final port state and the invariance of h_u all hold exactly. Brute force on actual tournaments in bases 2, 3 and 4
  found no mismatch: 19,811 exhaustive instances and 2,100 random ones.
- Every coverage number in THM-5085 (4) reproduces exactly, including the search-state counts. Two structurally
  different engines agree.
- Q1 and Q2 are reproduced. The addendum predictions Q5 and Q6 also hold: the misses below 10^7 are exactly {7, 21} and
  {10, 31}.

**What must change before promotion.**
- **(M1) Novelty.** At x = 2, THM-5085 (1) and (2) are KNOWN from Sessions 195, 197, 198 and 205. THM-4888, a canon
  theorem from the same machine (kind-pasteur S13), states the same tail formula with the same Q. The new content is the
  base-x extension and the computations.
- **(M2) HYP-9655.**
  - (A) contradicts the heuristic of the same machine's open, audited HYP-9480 (C), which HYP-9655 does not cite.
  - (C) is HYP-9480 (A)'s root case, which already has a pre-registered conjectured limit of 1/2.
- **(M3) Wrong example.** The displayed example of (3) is wrong.
- **(M4) Over-claimed gloss.** The pure-dyadic gloss of (3) is over-claimed. A bridge-appending pure-dyadic ray with a
  negative offset exists.
- **(U1) Q4.** Q4 is decidable and is REFUTED. An LP dual certificate, verified exactly, shows that every core cover of
  [1001, 2*10^5] needs at least 44 core states.

No proof is wrong. All items are typing, credit, wording or number corrections.

## Confirmed (no correction needed)

| claim | audit result | evidence |
|---|---|---|
| Lemma 1 (S and B transfers, h_u invariant), all x | exact for x = 2, 3, 4 on tournaments of order <= 14: 19,811 exhaustive instances (all labelled cores of order 1-4, every port, every word of length <= 6 or 5) and 2,100 random ones (cores of order 5-7) | `out/check_transfer.out`, `out/check_transfer_seed99.out` |
| Theorem 2 closed form, final (R, Delta) incl. Delta0 at L = 0 | exact, same runs; Q pair sum = Q recursion Q(2m) = Q(m), Q(2m+1) = Q(m) + m + (m mod 2) | same; `out/dictionary_check.out` |
| (3) transport along S^s(pattern I), s >= 1; pure dyadic => abs(I) = 1, offset h0 | correct as formally stated, for every x (for odd x the hypothesis never holds) | proof below; sink rays covered by the exhaustive runs |
| PROOF Cor. (c): sink-ray offsets can be 7 or 21 | true, e.g. a one-apex core (h0, hu) = (9, 1), I = {1, 2}: offset 9 - 2 = 7 (380 instances) | `out/dyadic_rays.out` (c) |
| I_1 = SB, I_2 = SSB add nothing | they are words; the audit closure uses S and B only and reaches the same states (identical state counts) | `out/cover_*` |
| class counts, order <= 7 | 1, 1, 2, 4, 12, 56, 456 (A000568); 3,613 rooted classes | `out/class_counts.out` |
| distinct core states (h0, R0, Delta0, hu) | base 2: 46 / 273 / 2,831 (order <= 5 / 6 / 7); base 3: 46 / 278 / 2,901; identical to the sets from all 2,097,152 labelled order-7 tournaments (no isomorphism code) | `out/core_states.out`, `out/labelled_order7.out`, `out/cores_b*_le*.txt` |
| (4) base 2, order <= 7, < 10^6 | misses exactly 7, 21; 57,330,333 states (= S19) | `out/cover_b2_le7_1e6.out`, `out/cover2_b2_le7_1e6.out` |
| (4) base 3, order <= 7, < 10^6 | misses exactly 10, 31; 34,006,533 states (= S19) | `out/cover_b3_le7_1e6.out`, `out/cover2_b3_le7_1e6.out` |
| (4) base 2, order <= 6, < 10^6 | the same 20 misses (7, 21, 34869, ..., 949777); 20,372,148 states (= S19) | `out/cover_b2_le6_1e6.out` |
| README rows: order <= 6 base 3 < 2*10^4 (10, 31, 5401); order <= 5 base 2 "hundreds" < 2*10^4 | 10, 31, 5401; 320 (base 3: 193) | `out/cover_b3_le6_2e4.out`, `out/cover_b2_le5_2e4.out`, `out/cover_b3_le5_2e4.out` |
| Q4 greedy covers (81 base 2, 86 base 3) | both lists consist of genuine core states and cover every target; audit greedy: 82 / 85 | `out/q4_b2.out`, `out/q4_b3.out` |
| Q4 lower bound 17 / 20 (pairwise disjoint covering sets) | the argument is valid (each such target needs its own core); the numbers are reproduced and the disjointness re-verified | same |
| HYP-9655 (C) table, overall 49,728 / 99,500 | every entry reproduced; the table is exact (proof under U4) | `out/oneapex_share_B22.out` |
| HYP-9655 (A) => binary and ternary Chan-Pak | valid, given the exclusions THM-1370 (7, 21) and THM-4940 (3), (4) (10, 31); see m4 | - |
| Q5, Q6 (addendum, c5b81ee9eb) | misses < 10^7 exactly {7, 21} (573,476,037 states) and {10, 31} (340,144,258 states); the closed-form engine gives identical hit sets on [0, 10^7] | `out/cover_b*_le7_1e7.out`, `out/cover2_b*_le7_1e7.out` |

**Lemma 1, re-derived.**
- In T' = T + w (B letter), w's only out-neighbour is u. So every new cycle is u -> P -> v -> w -> u, with P a path from
  u in T. It is odd iff abs(V(P)) is even, and P = (u) is impossible because u does not beat w.
- All new cycles contain u, so a packing holds at most one. The remainder is T - V(P), which gives h' = h + xE.
- A path from u in T' is an old path, possibly followed by the new vertex, which can only be last. The single exception
  is (u, w). An unused new vertex is a sink of the remainder.
- Hence f'(z) = (1 + z) f(z) - [B] z^2 I(Omega(T - u), x). Removing u leaves the new vertex beaten by all of T - u, so
  h_u is invariant.
- Evaluating at z = 1 and z = -1 gives the stated laws.

**(3), pure-dyadic clause.**
- If 2^s (x/2) R0 N(I) is a power of 2, then for odd x it is impossible.
- For even x every positive integer factor must be a power of 2, so N(I) is one, abs(I) = 1, Q(I) = 0 and o(I) = h0.

**x = 2 dictionary.** E_2 - O_2 = -e_u, where e_u is the number of Hamiltonian paths ending at the port. Checked on
2,000 random tournaments (`out/naive_check.out`).

## MAJOR findings

### M1. THM-5085 (1) and (2) are KNOWN at x = 2; only the base-x extension is new

**The same construction exists earlier in the repository, for arbitrary rooted tournaments, with the same formulas.**
The letters are O = S (append a universal sink, keep the root) and Z = K = B (append a vertex beating only the
retained root).
- `session195/combined-growth/PROOF.md` §2 (2026-11-05, audited in `combined-growth-audit`) gives the laws
  K(h,e,R,d) = (h+R-e, d, 2R-d, d) and D(h,e,R,d) = (h, 0, 2R, d). They are proved for every rooted tournament.
- `session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md` (audited in one-core-review):
  - (1) gives the O/Z laws;
  - (6) states: "For an arbitrary rooted input (h,e,R,d), the scalar output of any O/Z suffix is
    h + t(R-d) + d(F(t)-1) + (t mod 2)(d-e)".
- `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md` §4-5 (audited in 4b1a06ec1e):
  - (7) is the same suffix law;
  - (8) is the two-sided prefix/core/suffix formula for an arbitrary rooted core;
  - §5 gives a complete finite target inverse.
- `session198/source-digit-blocks/DIGIT_BLOCK.md` (1):
  - D' = D + UA - E eps_0 - d C(A), where C(A) = sum_{h<i} eps_h eps_i 2^(i-h-1) + sum_{i>=1} eps_(i-1) eps_i is
    THM-5085's Q(I) verbatim;
  - U' = BU - dR and E' = dt, which are THM-5085's final port state.
- `session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md` (independently audited) gives
  H_s(m) = h0 + R0 m - d Q(m) - e0 (m mod 2), with Q(2m) = Q(m) and Q(2m+1) = Q(m) + m + (m mod 2). Its proof unrolls
  R_i = 2R_(i-1) - d b_i, exactly as in THM-5085 Theorem 2.
- **Canon.** The Setting of THM-4888 (kind-pasteur-2026-11-06-S13, this machine) states the same tail formula and the
  O/Z laws, credited to Session 205.

**The dictionary.**
- At x = 2: (h, e, R, d) = (h, -Delta, R, h_u), with e the number of Hamiltonian paths ending at the port.
- The four formulas coincide identically (`out/dictionary_check.out`: 0 mismatches for N < 2^17).
- E_2 - O_2 = -e_u is verified in `out/naive_check.out`.

**What the packet says instead.**
- THM-5085 types (1)-(3) as "PROVED".
- Its KNOWN and credited list names only the one-apex recurrence, Session236, THM-5040/5055 and THM-5070.
  README "Concurrent work", the
  PREREGISTRATION ("What the tails add ... they satisfy the closed form") and the INDEX rows present the closed form as
  this lane's contribution.
- None of Sessions 195, 197, 198, 205 or 206, THM-4888 or HYP-9480 is cited anywhere in the packet. MISTAKE-722
  recorded the same kind of omission about the same O/Z calculus.

**What is new.**
- The base-x version. Remainders are counted by I(Omega(T - V(P)), x), and Delta = E - O replaces -e, which has no
  base-x meaning. The proof is the same, and the audit verified it at x = 3 and x = 4.
- The base-3 family it produces.
- The coverage computations.

(3) is an immediate corollary, since Q(2m) = Q(m).

**Related finite coverage that should be cited next to (4).**
- Session197's one-core bank (cores of order <= 5 with binary prefixes and suffixes) misses exactly 7, 21 and 462661
  below 10^6.
- Session195 combined-growth support (order-5 gadgets with root switches) misses 7, 21, 67 and 173 below 2*10^4.

**Required.**
- Retype (1) and (2) as "KNOWN at x = 2 (the citations above); PROVED for general x (the base-x extension is new)".
- Mark (3) as a corollary.
- Add the citations to THM-5085, README, PROOF, PREREGISTRATION (as an appended note) and the INDEX rows.
- Add a MISTAKES entry.

### M2. HYP-9655 conflicts with, and partly duplicates, HYP-9480 (same machine, audited) without citing it

**(A) versus HYP-9480 (C).**
- HYP-9480 (kind-pasteur-2026-11-06-S13; audit PASS WITH CORRECTIONS, MISTAKE-727) is OPEN. Its (C) conjectures:
  "every finite physical bank misses a positive proportion of the odd integers".
- It adds that "under a Poisson heuristic, a bank of effective intensity Lambda misses a proportion of order
  exp(-Lambda)". Citing a one-apex x order-<=7 bank that covers every odd N <= 2^35 except 7 and 21, it concludes that
  this "shows only that a positive miss proportion can lie beyond any feasible census".
- HYP-9655 (A) asserts the opposite for the finite bank C_7 (2,831 tail forms). Its evidence is coverage below 10^6,
  which by HYP-9480's own analysis cannot discriminate between the two.
- The two are not formally contradictory. HYP-9480's banks consist of unit-start seeds, while C_7 has arbitrary ports.
  But the tail forms h0 + R0 m - d Q(m) - e0[m odd] are identical, so the heuristic is stated for exactly this kind of
  bank. HYP-9480 notes only that S82's family is far from Poisson, so its census does not test (C). C_7 is a one-sided
  O/Z tail bank, exactly the kind (C) addresses.
- HYP-9655 must cite HYP-9480 (C), state the scope difference, and either retype (A) as a competing conjecture with an
  explicit caveat or downgrade it to a question.
- Its falsifier, "a value v < 10^9", is weak at the Poisson scale.

**(A) versus Session 205 and Session 206.**
- At x = 2, Session205 already formulates the reduction "eventual odd coverage by a fixed finite bank's O/Z tails would
  suffice for the unrestricted spectrum".
- Session206 proves a necessary capacity condition and the non-cofiniteness of the one-apex (root) image.
- The audit evaluated the capacity bracket: 31.2 for C_7 and 8.74 for order <= 6 (`out/capacity_b2.out`). This gives
  no obstruction, but it belongs in HYP-9655.

**(C) is HYP-9480 (A)'s root case.**
- The one-apex image is the root image of Sessions 205/206 and HYP-9480. Two checks confirm it:
  - 249,077 odd values <= 999,999, the figure recorded by Session205 and Session197 (`out/oneapex_share_B22.out`);
  - share 0.4974 on [2^20, 2^21), HYP-9480's figure (`out/oneapex_share_B22.out`).
- HYP-9480 (A) already conjectures the limit 1/2. Its pre-registered P2 measures 0.498538 on [2^27, 2^28).
- Session206 proves that the image has at most (log 2) X/2 + O(sqrt X) elements below X.
- HYP-9655 (C)'s "open whether the share has a positive limit, tends to 0, or oscillates" should be replaced by these
  citations.

**Internal inconsistency.** The README says "So the core order needed may grow slowly with the target
(HYP-9655 (B))". That contradicts (A), in which order 7 suffices for every target.

**Supporting data for (B).** The order <= 6 bank misses 169 values below 10^7. Their density is stable at about
3.3*10^-5 per odd integer on [10^5, 10^7] (`out/b2_le6_density_1e7.out`). This is consistent with a positive miss
proportion, that is with HYP-9480 (C), and it strengthens (B).

### M3. The displayed example of (3) is wrong

**The error.**
- THM-5085 (3), README (3) and PROOF.md Cor. (c) all state: "T0 = {u} with I = {1, 2} gives 3 * 2^(s+1) - 1".
- The theorem's own formula gives o(I) = 1 - 2 = -1 and N(I) = 3, hence -1 + 3 * 2^s = **3 * 2^s - 1**.
- Brute force on the actual tournaments {u} + S^s BB gives H = 5, 11, 23, ..., 1535 for s = 1..9, which is
  3 * 2^s - 1 every time (`out/check_transfer.out`).

**The likely source.** The stated form is Session236's one-apex ray 6 * 2^m - 1 = 3 * 2^(m+1) - 1. That ray is a
different word family (appended bridges) with the same set of values.

**Required.** Correct all three places, and record the error in MISTAKES.

### M4. The pure-dyadic gloss of (3) is over-claimed; bridge-appending dyadic rays carry negative offsets

**What the formal statement covers.** The formal statement (3) concerns the sink-prepending rays S^s(pattern I) and is
correct.

**The glosses that go further.**
- THM-5085 title: "pure dyadic rays force one bridge and a positive H-valued offset".
- THM-5085 status: "Here (3) extends the barrier to multi-bridge pure-dyadic rays".
- README (3): "Pure dyadic leading parts force a single bridge. Then the offset is h0 > 0, an H-value (Session236's
  barrier, extended to all tails). So the values just below powers of 2 always need non-dyadic leading parts."
- PROOF.md §3: "The family therefore realises 2^K + c along a pure dyadic ray only for offsets c that are themselves
  values ... extends it to multi-bridge tails".
- INDEX: "pure dyadic rays keep Session236's barrier".

**The counterexample.**
- After a bridge, the state is (h, R, -h_u). Appending k more bridges gives
  h = [h - (x/2)(R - h_u)] + 2^k (x/2)(R - h_u).
- Take the order-7 core with rows [94, 84, 32, 6, 12, 27, 60] and port 0. Its base-2 state is
  (h0, R0, Delta0, h_u) = (71, 287, -5, 31).
- After one B the state is (353, 543, -31), and R - h_u = 512 = 2^9. Brute force gives
  H(core + B^k) = 2^(k+8) - 159 for k = 2..8, at orders 9-15 (`out/dyadic_ray_bridges_bruteforce.out`).
- This is a pure-dyadic ray with offset -159, realised inside the fixed-port family by multi-bridge words. It produces
  values 2^K - 159 for every K >= 10.
- Among the cores of order <= 7 and words of at most 6 letters ending in B, there are 144 such (core, word) pairs with
  offset <= 0, or with offset 7 or 21 (`out/dyadic_rays.out` (b)).

**Session236's barrier.** It is scoped to its single-final-bridge mechanism and remains correct. Its (3) also gives
offset 1 for pure-dyadic one-apex rays.

**Required.**
- Restrict every such sentence to sink-prepending rays.
- Add the bridge-appending formula and the example.
- Drop "extended to all tails" and "extends it to multi-bridge tails".

## Outcome retyping and audit upgrades

**U1. Q4 is REFUTED, not UNDECIDED.**
- Dropping targets only weakens the covering constraints. So the LP relaxation over the targets reached by at most 45
  core states bounds the full minimum from below.
- Its dual gives a certificate: 489 targets with weights y_t >= 0, sum 43.3707, and every core column sum <= 1.
- This was re-verified with per-core state closure and integer arithmetic, independently of the program that produced
  it (`out/q4_cert_verify_b2.out`, `out/q4_lp_b2_cert_K45.txt`).
- Hence every cover of the odd targets in [1001, 2*10^5] by order-<=7 core states has at least **44 > 40** members.
  The pre-registered Q4 ("at most 40") is therefore REFUTED.
- In base 3 every cover needs at least **50** (`out/q4_cert_verify_b3.out`).
- With the audit's greedy covers (82 and 85), the minimum covers lie in [44, 81] (base 2) and [50, 85] (base 3), not
  [17, 81] and [20, 86].
- The S19 bounds 17 and 20 are valid but weak.
- Required:
  - append "Q4: REFUTED by the audit's LP certificate" to the outcomes, without editing the original outcome line;
  - correct THM-5085 status, README, HYP-9655 Route and the INDEX rows.

**U2. Q5 and Q6 are confirmed independently.**
- The state-closure engine gives the same state counts as the S19 outcome that appeared in the working tree during the
  audit: 573,476,037 and 340,144,258.
- The closed-form engine gives identical hit sets on [0, 10^7].
- With this, ternary Chan-Pak holds to 10^7, and the binary result to 10^7 sits within the known 2^35.

**U3. HYP-9655 (B) evidence extended.**
- Order <= 6 in base 2: 169 misses below 10^7, at a stable density (M2).
- Order <= 6 in base 3: 15 misses below 10^6, namely 10, 31, 5401, 40597, 97282, ..., 904876
  (`out/cover2_b3_le6_1e6.out`).

**U4. The HYP-9655 (C) table is exact, not a lower bound.**
- H(2W+1) = H(W).
- For even W' = 2W >= 2, a(W') >= ceil(W/2), so H(W') >= 1 + W'/2.
- Every value v is therefore attained at W = 0 or at an even W' <= 2(v - 1).
- So W < 2^22 decides all values below 2^21, which covers all windows k <= 20.
- The core-{u} closed form gives the same value set (`out/oneapex_share_B22.out`).

## MINOR findings

- **m1. Q3 typing.**
  - "The intended statement is CONFIRMED" overstates the evidence. The misses below 2*10^5 (7, 21, 34869, ..., 153451)
    were listed under "Already seen" before the pre-registration, and any miss in [2*10^5, 10^6) is automatically
    above 30,000.
  - So the "above 30,000" clause carried no predictive content. The only genuine prediction was "at most 500 misses",
    and it holds (20 misses).
  - Suggested wording: "REFUTED as worded; its only predictive clause (at most 500 misses) holds".
- **m2. Missing MISTAKES entry.** PREREGISTRATION says "This wording slip is recorded in MISTAKES". No MISTAKES.md
  entry references S19, THM-5085 or HYP-9655, either at c322f6b577 or at origin/main d0c7de6204.
- **m3. HYP-9655 title.** "0.5006 at 2^16" should be 0.5001 (16386/32768 = 0.50006), as in its own table.
- **m4. HYP-9655 (A) citations.**
  - "Why it implies Chan-Pak" must name the exclusion theorems: THM-1370 for 7 and 21, THM-4940 (3), (4) for 10 and 31.
  - The `related` list should add THM-4940, HYP-9480 and THM-4888.
- **m5. Records in README.**
  - Binary is known to 2^35 (THM-4875 (D) addendum, one run), not only 2^32.
  - THM-5070's binary port-chain bound is 2*10^5. So the README's "a second explicit family, after THM-5070's port
    chains, that realises both Chan-Pak spectra to 10^6" is wrong for base 2.
  - "the ternary bound 10^6 equals THM-5070's" is also wrong. The S18 audit's union (order-5 chains to 10^6, plus the
    15 generators to 3*10^6 with the same ten small misses) reaches 3*10^6, as the S19 addendum itself states.
- **m6. Credit for negative offsets.**
  - THM-5085's "(3) ... shows that non-dyadic multi-bridge rays carry negative offsets" should credit Session236
    offsets/PROOF.md (3).
  - That section already exhibits negative one-apex offsets (6 * 2^m - 1, values 5, 11, 23, ... - the same values as
    S19's example) and offset 1 on pure-dyadic one-apex rays.
- **m7. Base-dependent wording.** "Positive H-valued offset" is vacuous at x = 3, where the leading coefficient is
  never a power of 2. Write "the value h0 of the core (x even)".
- **m8. Unproved claim in README.** "So there is no tournament cone" is unproved. Soften it to "no port move of this
  type is a cone".
- **m9. Credits to add in README and PROOF.** Session236's single-bridge ray with a singleton block is, at x = 2, also
  an instance of the earlier O/Z suffix law (M1). Keep the Session236 credit for the barrier reading, and add the
  earlier sources for the formula.
- **m10. INDEX rows.** Both INDEX rows need the M1, M3, M4 and U1 corrections: KNOWN typing, the example, the dyadic
  scope, and "Q4 REFUTED, 44-81".

## Suggested replacement text (for the author)

**THM-5085 status.**
- KNOWN at x = 2:
  - (1) is the K/D (= Z/O) law of `session195/combined-growth` §2 (also Session197 ONE_CORE (1));
  - (2) is the O/Z tail formula of Session205 `finite-oz-tail-bank` (also Session197 ONE_CORE (6),
    ROOTED_ORDERED_AMALGAMATION (7), Session198 DIGIT_BLOCK (1)), as restated in THM-4888;
  - dictionary e = -Delta, d = h_u.
- PROVED (new): the base-x versions of (1) and (2).
- (3) is an immediate corollary, restricted to sink-prepending rays.
- FINITE-EXACT: (4), to 10^6 (audit: 10^7).
- Pre-registration: Q1 and Q2 CONFIRMED; Q3 REFUTED as worded, its predictive clause holds; Q4 REFUTED by the audit's
  LP certificate (at least 44 core states are needed).

**(3), example.** "T0 = {u} with I = {1, 2} gives 3 * 2^s - 1 (s >= 1)."

**(3), dyadic clause.**
- "For sink-prepending rays, a power-of-2 leading coefficient forces abs(I) = 1 and offset h0."
- Then add: "Bridge-appending rays differ. From a state (h, R, -h_u) the word B^k gives
  [h - (x/2)(R - h_u)] + 2^k (x/2)(R - h_u). Pure-dyadic leading parts can then carry negative offsets: an order-7
  core gives 2^K - 159 for every K >= 10."

**HYP-9655.**
- (A) needs a paragraph on HYP-9480 (C), covering the conflict, the scope difference and the Poisson scale.
- The Route should read "between 44 and 81".
- (C) needs the HYP-9480 (A) root conjecture (limit 1/2, P2 = 0.498538) and Session206's bound.
- Fix the title number 0.5006 -> 0.5001.

**MISTAKES entry.** Record:
- (a) the x = 2 transfer and closed form typed as new, although they are KNOWN from Sessions 195/197/198/205 and
  THM-4888 (same machine); and HYP-9655 written without HYP-9480;
- (b) the exponent in the (3) example;
- (c) the dyadic gloss;
- (d) the Q3 wording slip promised in PREREGISTRATION.

Process rule: before typing a transfer law new, grep the session19x-20x trees and the machine's own canon for the
formula's signature (`2R - d`, `universal sink`, `O/Z`, `retain`).

## Independent checks run

All outputs are in `audit/out/`. Scripts and C sources are in `audit/`, and `run_all.sh` reproduces everything.

1. **tbrute.c and naive_check.py.**
   - The new subset-DP counter is compared with a naive counter (permutations, explicit cycle lists, explicit disjoint
     packings, DFS over paths) on 2,000 random tournaments of order 1-8: 0 mismatches.
   - OCF and E_2 - O_2 = -e_u also hold.
2. **check_transfer.py** (two seeds).
   - Lemma 1 and Theorem 2 (h, final R and Delta, h_u invariance) in bases 2, 3, 4.
   - Cases: 19,811 exhaustive cases plus 10 example cases plus 600 or 1,500 random cases; 0 mismatches.
   - Also the (3) example values.
3. **cores.c and dedup_cores.py.**
   - Classes by one-vertex extension plus a brute-force canonical form: counts 1, 1, 2, 4, 12, 56, 456.
   - Port states from the classes, and independently from all labelled tournaments of order 5, 6 and 7 with port 0.
   - The state sets are identical; the counts are as in the table.
4. **cover.c** (state closure, letters S and B only).
   - Every coverage row, at 2*10^4, 2*10^5, 10^6 and 10^7.
   - The state counts equal S19's.
5. **cover2.c** (closed-form parameter enumeration with the a-priori bound 0 <= Q(m) <= m - 1).
   - Same misses, and identical hit sets to cover.c on [0, 10^6] and [0, 10^7].
6. **q4.c, q4_lp.py and verify_cert.c.**
   - Validity of S19's greedy lists; independent greedy covers.
   - The disjoint-set bound.
   - LP dual certificates for 44 (base 2) and 50 (base 3), each re-verified by closure.
7. **oneapex.py.**
   - The HYP-9655 (C) table, the [1001, 2*10^5] share, the exactness proof, and the cross-check against the core-{u}
     closed form.
   - The k = 20 window, and the count 249,077 below 10^6.
8. **capacity.py.** Session206's capacity bracket for the banks of order <= 5, 6 and 7.
9. **dyadic_rays.py, find_core.c and dyadic_bridge_example.py.**
   - Single-bridge sink rays near powers of 2.
   - Bridge-appending dyadic rays with offsets <= 0, 7 or 21, plus the brute-forced order-7 example (M4).
   - Sink-ray offsets 7 and 21.
10. **dictionary_check.py.** THM-5085's Q coincides identically with Session205's Q recursion, Session197's
    F-form and Session198's C(A).
11. **Priority sweep** (grep):
    - the `2R - d` / `universal sink` / `O/Z` literature (session195-206, THM-4882/4883/4888, HYP-9480/9494);
    - Session236 offsets/PROOF.md and last_letter/PROOF.md;
    - SEAM.md, Session191 one_apex/PROOF.md, THM-4866 and THM-4875;
    - `git log origin/main` (latest d0c7de6204: no newer S19-relevant commits).

    No earlier source has the base-x (x != 2) fixed-port law or the order-<=7 coverage runs.

## Reproduction

`cd audit && BUILD=/some/dir/outside/repo sh run_all.sh`. It needs gcc and Python 3 with numpy and scipy (HiGHS) for
q4_lp.py. Everything else uses the standard library. On Windows, set BUILD to a C:/ path so Python can reach the
binaries.

The script was run once from scratch after the report was drafted. It exited 0, and every regenerated file in `out/`
was identical to the drafted outputs (`diff -r`, line endings ignored). The run takes several minutes; the two 10^7
state closures dominate, at about 900M hash inserts in total.
