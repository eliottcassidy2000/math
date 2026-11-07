# Audit of THM-4922 and HYP-9512 (first creators, converse fibers, genealogy)

Independent audit, 2026-11-07, of checkpoint 9d7c9f4735 (packet `first_creators_genealogy_20261107/`,
`THM-4922-...md`, `HYP-9512-...md`). The auditor wrote all code from scratch. The packet's code was read only to
compare definitions and was not run.

## Verdict: PASS WITH CORRECTIONS

Every count in the packet reproduces exactly, at orders 1-10, with two exceptions: the maximum fiber size quoted for
orders <= 9 (minor m1), and the named smallest counterexample (M1). The independent order-10 run took 119 s on 12
threads. Its per-value table equals order10.out
on all 6,143 values (fiber, #SC, #pairs, #fully critical), and its k distribution is identical.

The following must be corrected:
- the smallest counterexample, which is h = 9, not 23 (M1);
- "every first creator is strong", which is not proved in general (M2);
- credit to THM-4889 and to the same-prompt Session224, which pushed 10 s earlier (M3);
- the blindness of Q4 (M4);
- the scope of the orphan and no-hypo statements (M5).

One status goes up: HYP-9512 (D), the odd-order self-converse balance, is a theorem. A short proof is given in section 4.

## 1. Recomputed numbers

Method:
- `aud_enum.c` uses canonical augmentation (McKay), not a global hash set.
- Its canonical labelling is its own individualization-refinement:
  - initial colours (score, cyclic triangles through v);
  - refinement by in-neighbour counts per cell;
  - the smallest non-singleton cell is individualized;
  - the maximal leaf code is kept, with all leaves explored.
- The labelling gives |Aut| and the Aut-orbits. An extension is accepted when its new vertex lies in the orbit of the
  last canonical vertex.
- Validation:
  - the brute-force canonical form (all n! labellings) agrees at orders 2-7 (`aud_enum_brute7.out`);
  - the orbit-stabilizer identity sum n!/|Aut| = 2^C(n,2) holds at every order 1-10;
  - klein's census tables agree (below).
- `aud_analyze.py` derives everything below from the auditor's data.

| Quantity | Packet | Auditor | |
| --- | --- | --- | --- |
| classes, n = 1..10 (A000568) | 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056 | same; orbit-stabilizer exact at every n | match |
| SC classes (A002785) | (odd-order splits only) | 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784 | match |
| klein census, (h, c3, strong) and (h, #Ham. cycles) class counts | "reproduced, orders 3-9" | identical at orders 3-9; at order 10 the per-h totals, least c3 and least hc (with their counts) are identical | match |
| values with n(h) <= 9 | 1,520 "(h >= 3)" | 1,520, **including h = 1** (= \|S(9)\|) | label wrong (m1) |
| (s, p) = (1, 0) / (0, 1) | 44 / 25, lists in PROOF section 2 | 44 / 25, identical lists | match |
| failures, n(h) <= 9 | 1,451 | 1,451, of which 8 are two SC classes (9, 15, 17, 45, 145, 621, 3075, 3105) | match |
| smallest failure | h = 23 (0, 2), order 6 | **h = 9 (2, 0), order 5**; then 15 (2, 0) at order 5, 17 (2, 0), 23 (0, 2), 29 (2, 1), 33 (1, 2), 37 (1, 2) | **mismatch (M1)** |
| 29; 33 and 37 | (2, 1); (1, 2) | (2, 1); (1, 2) | match |
| largest fiber, n(h) <= 9 | 2,458 (PROOF section 2, PREREGISTRATION) | **257 (h = 1503)**; 2,458 is order 10 (h = 7799) | mismatch (m1, M4) |
| kinds, n(h) <= 9: SC-only / NS-only / mixed | 57 / 371 / 1,092 | 57 / 371 / 1,092 | match |
| blue/black, 8 cells of PROOF section 3 | 89/43, 0/124, 521/103, 538/101 | same | match |
| blue / black totals | 1,149 / 371 | 1,149 / 371 counting h = 1; the table's cells (h >= 3) sum to 1,148 / 371 | m2 |
| SC classes by H mod 4, odd n = 3, 5, 7, 9 | 1/1, 4/4, 44/44, 1376/1376 | same (n = 1: 1/0); even n: all = 1 mod 4, including 8,784 at n = 10 | match |
| median h/max H(n(h)): SC-only, NS-only, mixed | 0.918; about 0.55 | 0.918 (50 values); 0.545; 0.551 | match |
| first creators at orders 3-9 | 328, 5,457, 163,676 at 7-9; 169,501 at 4-9 | 1, 1, 6, 33, 328, 5,457, 163,676; 169,501 | match |
| k = n at orders 7, 8, 9, 10 | 78, 2,118, 82,553, 5,433,631 | same | match |
| orphans (k = 0) | only the two first creators of 63 | only those two at orders 4-10, **but the 3-cycle (h = 3) at order 3** | scope (M5) |
| hypo analogue | none | none at orders 3-10; **the order-2 tournament is one** (trivially) | scope (M5) |
| order 10: classes, first creators, hypo | 9,733,056; 8,647,232; 0 | same | match |
| order-10 k distribution | 0:0 1:0 2:41 ... 10:5433631 | identical | match |
| order 10: new values, dichotomy, SC-only / NS-only / mixed | 6,143; 55 in [14675, 15745]; 5 / 3,605 / 2,533 | same; SC-only = 14937, 14973, 15597, 15601, 15745 (the middle two are (2, 0)) | match |
| order 10, h = 3 mod 4 | all 3,067 black | 3,067 black, 0 blue | match |
| order 10, h = 1 mod 4 blue | 2,538 of 3,076 | 2,538 of 3,076 | match |
| max H(10) = 15745 | unique SC class | (s, p) = (1, 0) | match |
| values with > 2 classes at least order and least c3 | 5,691 of 7,663 | 5,691 of 7,663; by least number of Ham. cycles: 1,901 | match (m8) |
| genealogy thresholds n max H(n - 2) | 1,701 (n = 9), 6,610 (n = 10) | same; 671 of 1,200 and 4,360 of 6,143 new values lie above | match (m3) |
| non-strong first creators | "none" (claimed PROVED) | 0 at every order <= 10; order 11 (klein census_11): no non-strong class has a new value | finite only (M2) |
| the order-8 first creators of 63 | THM-344's two classes | a converse pair, strong, \|Aut\| = 1, c3 = 8, deletion counts {1, 11, 25, 25, 33, 33, 37, 37} | match |

Raw outputs:
- `aud_enum.out`, `aud_enum_brute7.out`, `aud_analyze.out`;
- `aud_fibers_le10.txt` (h, n(h), fiber, s, p for all 7,663 values);
- `aud_o10_byH.txt` (order-10 per-h table);
- `aud_sc_mod4.out`.

## 2. Major findings

**M1. The smallest counterexample is h = 9, not h = 23.**
- h = 9 is first attained at order 5 by two self-converse classes. They are not converses of each other, so (s, p) =
  (2, 0) and the owner's "one class or two converse classes" fails.
- h = 15 (order 5) and h = 17 (order 6) fail the same way.
- The packet's own dichotomy lists exclude 9, 15 and 17, so the packet contradicts itself.
- Session224 proved the h = 9 case (`session224_first_occurrence_20261107/proof/PROOF.md` section 1), pushed in
  2054513960 at 15:25:57, before 9e07033ca4.

Replacements:
- THM-4922 body (3), first bullet, "The smallest failure is h = 23: two converse pairs at order 6." becomes:
  "The smallest failure is h = 9: two self-converse classes at order 5, neither the converse of the other (Session224,
  PROVED). h = 15 (order 5) and h = 17 (order 6) fail the same way, and h = 23 is the smallest value with two converse
  pairs."
- README Results 1, "The smallest failure is h = 23: two converse pairs at order 6. Fibers reach 2,458 classes."
  becomes: "The smallest failure is h = 9: two self-converse classes at order 5 (proved first by Session224). h = 23 is
  the smallest value with two converse pairs. Fibers reach 257 classes through order 9 (h = 1503) and 2,458 at order
  10 (h = 7799)."
- PROOF section 2, "Smallest counterexamples", replace the three bullets by:
  - "h = 9 and h = 15: two SC classes at order 5 (Session224 proved h = 9 first). h = 17: two SC classes at order 6.
    Eight values with n(h) <= 9 fail this way: 9, 15, 17, 45, 145, 621, 3075, 3105."
  - "h = 23: two converse pairs at order 6."
  - "h = 29: two SC classes and one pair."
  - "h = 33 and h = 37: one SC class and two pairs (33 is Session224's recovered H = 33 fiber)."

**M2. "Every first creator is strong" is not proved.**
- The argument removes a one-vertex strong component. A non-strong first creator could instead be an ordinal sum
  A => B of strong components, all of order >= 3.
- Excluding that needs H(A)H(B) in S(|A| + |B| - 1). This is exactly klein's HYP-9481 (b) in its reading "products of
  smaller counts are never new", which is OPEN.
- The claim also implies HYP-9481 (b) as stated (N(n) ⊆ Ss(n)), since every new value has a first creator.
- What is proved: there is no one-vertex component, and each strong component is a first creator of its own count, by
  replacement. This is Session224 section 3, already KNOWN.
- What is finite-exact: the claim holds through order 11. Order 10 was recomputed here (0 non-strong first creators);
  at order 11, klein's census_11 has no non-strong class with a new value (THM-4889 (2)).

Replacements:
- PROOF section 1, item 1 becomes: "**1. No one-vertex strong component, and closed under converse.** A first
  creator of h > 1 has no one-vertex strong component. H is the product of the counts of the strong components, and
  deleting a one-vertex component would realise h at a smaller order. For the same reason, each strong component C is
  itself a first creator of H(C): replacing it by a smaller tournament with the same count keeps H (Session224 section
  3). So a first creator is either strong, or an ordinal sum of at least two strong first creators, each of order at
  least 3. The second case never occurs through order 11. Order 10 was recomputed by the audit. At order 11, klein's
  census_11 has no non-strong class with a new value, as the audit checked (cf. THM-4889 (2)). Excluding it in
  general is HYP-9481 (b), which is OPEN. The first creators are closed under converse, since H(T^op) = H(T), so they
  are s SC classes plus p converse pairs."
- THM-4922 body (1) becomes: "**(PROVED)** The first creators of h are closed under converse. For h > 1, none has a
  one-vertex strong component, and each strong component is a first creator of its own count (Session224).
  **(FINITE-EXACT, orders <= 11; OPEN in general, HYP-9481 (b))** Every first creator is strong."
- THM-4922 status, the line "PROVED (elementary): strong and converse-closed fibers; the endpoint injection; the
  genealogy above n(h) max H(n(h) - 2)." becomes: "PROVED (elementary): converse-closed fibers with no one-vertex
  strong component (Session224 section 3); the genealogy above n(h) max H(n(h) - 2), from THM-4889's counting identity
  and padding; the odd-order self-converse balance (audit proof, section 4 of the audit). FINITE-EXACT, orders <= 11:
  every first creator is strong (HYP-9481 (b) in general, OPEN)."
- In the THM-4922 title, replace "are strong and closed under converse" with "are closed under converse and have no
  one-vertex strong component (strong through order 11)".

**M3. Priority and credit.** The following are KNOWN and were not credited.
- **The endpoint inequality.** H(T) <= sum over v of H(T - v) is THM-4889, listed under PROVED (elementary) as "the
  counting identity h(T) <= sum_v h(T - v)". It is also in klein's packet, section 4 ("the trivial bound is n").
  - The genealogy corollary (section 1.3) is new, but it is a one-line consequence of that identity and padding.
- **The dichotomy question and its refutation.** Session224 (`session224_first_occurrence_20261107`, commit
  2054513960, same owner prompt) contains:
  - the smallest counterexample h = 9;
  - the complete (s, p) birth ledger through order 8 (320 values, 285 with several converse orbits; the audit's count
    of dichotomy values with n(h) <= 8 is 35 = 320 - 285);
  - the strong-component min-plus structure of minimal realizations.

  S16's commit abfe525263 (15:26:07) was rebased onto 2054513960 at 15:26:12 and pushed as 9e07033ca4. So S16 had
  Session224 locally 14 minutes before its checkpoint, and Session224 is first-pusher on these points.
- **The fiber sizes.** For every h with n(h) <= 10, the number of first creators is in klein's least_order_10.txt (the
  strong class count at n(h), e.g. "7799 10 10 0 10[2458]"). What is new in S16 is the (s, p) split at orders 9 and 10.
- **The even-order black cell at class level.** "a_n(h) is even whenever n is even and h = 3 (mod 4)" is Session29
  converse congruences, section 4. "A uniquely attained value belongs to a self-converse class, and is = 1 (mod 4) at
  even n" is in `h_spectrum_multiplicities_20261021.md`. Both follow from THM-791, which is credited.

Replacement:
- Append to the KNOWN list in THM-4922's status (and to "Credited prior work" in PROOF.md):
  - "the counting identity H(T) <= sum over v of H(T - v) (THM-4889, klein-2026-11-06, section 4)";
  - "the refutation of the one-orbit dichotomy at h = 9, the (s, p) birth ledger through order 8 and the
    strong-component structure of minimal realizations (Session224, 2054513960, same owner prompt, pushed first)";
  - "the first-creator fiber sizes for n(h) <= 10 (klein's least_order_10.txt)";
  - "the class-level parity a_n(h) even for even n and h = 3 (mod 4) (Session29 converse congruences)".
- Add to `related`: THM-4889, HYP-9481, Session224's README, THM-4748, THM-466.
- Say what is new:
  - the (s, p) fibers at orders 9-10;
  - the genealogy data k(T), with the orphan pair of 63;
  - the no-hypo analogue;
  - the threshold corollary;
  - the blue/black value table;
  - the odd-order self-converse balance (now proved).

**M4. Q4 was not a blind prediction.**
- PREREGISTRATION.md (9e07033ca4, unchanged since) says, under "Data before this registration (orders <= 9)": "Fiber
  sizes reach 2,458".
- 2,458 is the order-10 maximum (h = 7799). The order-<= 9 maximum is 257 (h = 1503). So klein's order-10 fiber sizes
  (census_10 or least_order_10, committed 2026-11-06) were consulted before the registration.
- Those sizes already decide Q4. Only 57 of the 6,143 new order-10 values have at most two classes (0.93% < 5%), and
  the final count 55 is 57 minus the two (2, 0) fibers.
- Q1, Q2, Q3 and Q5 need deletion or self-converse data that klein's census does not have (klein's deletion data stop
  at order 9, dels_6..9, and the census records no converse type), so they stand as genuine tests.

Replacements:
- In THM-4922's status, "Q1-Q4 CONFIRMED, Q5 FAILED, explained by THM-791" becomes "Q1-Q3 CONFIRMED (blind); Q4
  CONFIRMED but pre-determined by klein's order-10 fiber sizes, which the registration quoted (2,458); Q5 FAILED
  (forced by THM-791)".
- In the PROOF section 4.1 table, the Q4 outcome becomes "**CONFIRMED, not blind**: the registration quoted klein's
  order-10 maximum fiber (2,458), and klein's fiber sizes bound the count by 57 (0.93%). Final count 55 of 6,143, all
  in [14675, 15745]; max H(10) = 15745 has a unique SC class."
- README Results 5 changes the same way.

**M5. The orphan and no-hypo statements are false as written at h = 3 and at order 2. Scope them, and bound "only
orphan" to order 10.**
- The 3-cycle is the first creator of 3, and it is an orphan: both its deletions have H = 1, first created at order 1.
- The order-2 tournament is a hypo analogue, trivially: it is not a first creator of 1, but its one-vertex deletions
  are.
- The packet's scripts start at order 4, but the titles do not say so.

Replacements:
- THM-4922 title, "the only orphan first creators" becomes "the only orphan first creators at orders 4-10". Also
  "no non-first-creator has every deletion a first creator" becomes "no non-first-creator of order 3-10 has every
  deletion a first creator".
- THM-4922 (3), Genealogy bullet becomes: "At orders 4-10, every first creator except the converse pair of first
  creators of 63 (order 8) has a vertex whose deletion is a first creator. So among the values h >= 5 with n(h) <= 10,
  63 is the only one whose first creators are not grown from first creators. (The 3-cycle, at order 3, is also an
  orphan: its deletions have H = 1.)"
- THM-4922 (3), No-hypo bullet becomes: "No class of order 3-10 that is not a first creator has every vertex deletion
  a first creator. (At order 2 the unique class trivially is one.)"
- HYP-9512 title and (A): "every odd h other than 7, 21 and 63" becomes "every odd h >= 5 other than 7, 21 and 63".
  "and 63 is the only orphan value" becomes "and 63 is the only orphan value with h >= 5".
- HYP-9512 (B): "no non-first-creator ever has every vertex deletion a first creator" becomes "no non-first-creator of
  order n >= 3 has every vertex deletion a first creator".
- README Results 3: "So 63 is the only value not grown from smaller first creators." becomes "So, among values h >= 5
  with n(h) <= 10, 63 is the only one not grown from smaller first creators."

## 3. Minor findings

- **m1.** PROOF section 2 table:
  - "values with n(h) <= 9 (h >= 3) | 1,520" becomes "values with n(h) <= 9 (all h >= 1) | 1,520".
  - "everything else (several classes) | 1,451 values; fiber sizes up to 2,458" becomes "everything else (two SC
    classes, or three or more classes) | 1,451 values (8 with two SC classes); fiber sizes up to 257 (h = 1503; 2,458
    at order 10, h = 7799)".
- **m2.** PROOF section 3: "Through order 9 there are 1,149 blue and 371 black values." becomes "Through order 9
  there are 1,149 blue values (1,148 with h >= 3, as in the table) and 371 black values."
- **m3.** PROOF section 1.3:
  - The proof uses padding (S(m) ⊆ S(m + 1)) without saying so. "Hence H(T - v) is not attained at order n - 2, and
    its least order is n - 1. QED." becomes "Hence H(T - v) is not attained at order n - 2, nor, by padding with a sink
    (S(m) ⊆ S(m + 1)), at any smaller order. T - v attains it at order n - 1, so T - v is a first creator. QED (n >=
    3)."
  - Section 4, "Proved part": "1701 at n = 9, 6610 at n = 10." becomes "the threshold is 1,701 at n = 9 (671 of the
    1,200 new values lie above it) and 6,610 at n = 10 (4,360 of 6,143)."
- **m4.** The "exact statement" overclaim:
  - PROOF section 4: "This is the owner's "generate larger structures from smaller ones", now an exact statement."
    becomes "This is the owner's "generate larger structures from smaller ones": finite-exact through order 10, proved
    above n max H(n - 2) (section 1.3), conjectured in general (HYP-9512 (A))."
  - README Results 3 heading: "(the owner's generation idea, made exact)" becomes "(the owner's generation idea,
    tested through order 10)".
- **m5.** The terminology is inverted. In THM-4870 a vertex whose deletion stays in the class is *non-critical*.
  - HYP-9512 (B), (C): "A vertex is critical when its deletion stays in the class of first creators, in THM-4870's
    sense of non-critical reductions." becomes "A vertex is non-critical, in THM-4870's sense, when its deletion stays
    in the class of first creators; k(T) counts the non-critical vertices."
  - Replace "fully critical" with "fully non-critical" everywhere: README 4, PROOF section 4 table header, HYP (C) and
    title, and the Q3 wording.
- **m6.** The blue shares, (E) and "about 83%". Among values not forced black, the blue share is 6/7, 42/58 = 0.72,
  82/124 = 0.66, 1012/1200 = 0.84 and 2538/3076 = 0.825 at orders 6-10. At even orders <= 8 it is 89/132 = 0.67.
  - HYP-9512 (E): "stays near 0.83" becomes "is 0.72, 0.66, 0.84 and 0.825 at orders 7-10 (82% in aggregate); a limit
    near 0.83 is EXPLORATORY".
  - README: "Elsewhere about 83% are blue." becomes "Elsewhere 82% are blue in aggregate (0.72, 0.66, 0.84, 0.825 at
    orders 7-10)."
  - THM-4922: "Otherwise about 82-84% of values are blue." becomes "Otherwise 82% are blue in aggregate, 0.66-0.86 by
    order (orders 6-10)."
- **m7.** HYP-9512 (A) body: "the new values of each order are exactly the values of one-vertex extensions of the
  previous order's first creators, plus the single orphan 63" becomes "the new values of each order n are exactly the
  new values among the one-vertex extensions of the first creators of order n - 1, plus 63 at n = 8". Extensions also
  produce old values, for example by adding a source.
- **m8.** PROOF section 2, the refinement sentence becomes: "5,691 values still have more than two classes at their
  least order and least level c3, and 1,901 at their least order and least number of Hamiltonian cycles, out of 7,663
  values with n(h) <= 10 (klein's census; recomputed by the audit)."
- **m9.** PROOF status: "**PROVED**: sections 1 and 4.1." becomes "**PROVED**: section 1 (as corrected), the 'Proved
  part' of section 4, and the odd-order SC balance (audit proof). **FINITE-EXACT**: sections 2, 3, 4 and 4.1."
  Section 4.1 is a computation.
- **m10.** The statement that no tournament is hypohamiltonian or hypotraceable is not in THM-4870's statement (its
  (C) concerns square-sum graphs). It is in THM-4870's source packet,
  `thomassen_merge_hypohamiltonian_groups_20261104.md` section 2, and is classical (Rédei; Camion). Cite it that way.
- **m11.** "It holds only at small values and near each order's maximum" becomes "mostly at small values and near
  each order's maximum (at order 7 also at 61, 85, 87, 125, 127, which lie at 0.32-0.67 of max H(7))".
- **m12.** Blue and black here are value colours. They differ from the line and class colours of the tiling
  metagraph (THM-791, HYP-4997, THM-2444-2454). Add after the definition: "These value colours are distinct from the
  metagraph's line colours."
- **m13.** Add these credits for 63:
  - THM-4748: the two first creators of 63 are the Fibonacci fan F_6 and its converse. Each is a transitive
    7-tournament plus an apex, which is why one deletion has H = 1.
  - Session24 multiplicities: they are the only strong H = 63 classes at any order, and they form a converse pair.
  - HYP-9359 item 4 already asks why 63 is an order-7 hole.
- **m14.** PREREGISTRATION.md is frozen. Add an erratum note in PROOF section 4.1: "The registration's data item 1
  has two errors: the smallest counterexample is 9, not 23; and 2,458 is the order-10 maximum fiber (order <= 9: 257)."
- **m15.** Cosmetic: the header comment of order10.c is copied from classgen.c ("usage: classgen NMAX outprefix").

## 4. Status upgrade supplied by the audit: HYP-9512 (D) is a theorem

**Theorem.** For every odd n = 2m + 1 >= 3, exactly half of the self-converse isomorphism classes of order n have
H = 1 (mod 4). (n = 1 is the trivial exception.)

*Proof.*

1. **Anti-involutions.** Let T be self-converse with anti-automorphism a. Then a^2 is an automorphism of odd order r,
   so sigma = a^r is an anti-automorphism with sigma^2 = 1. Sigma fixes at most one vertex, since an arc between two
   fixed vertices would be reversed. So at odd order it fixes exactly one vertex v0 and pairs the rest as {a_i, b_i}.
   Any two anti-involutions sigma and tau are conjugate in Aut(T): put g = tau sigma, of odd order r. Then
   sigma g sigma = g^-1, and h = g^((r+1)/2) gives h sigma h^-1 = tau (Session29 converse congruences, section 1).
2. **Orbit model.** Fix sigma0 of cycle type 1·2^m on V. Let X be the set of labelled tournaments on V with sigma0 as
   an anti-automorphism (|X| = 2^(m(m+1))), and G the centralizer of sigma0 (order 2^m m!). Then T -> [T] is a
   bijection from X/G to the self-converse classes of order n.
   - Surjective: every self-converse T has an anti-involution of the same cycle type.
   - Injective: if pi T1 = T2 with T1, T2 in X, then pi^-1 sigma0 pi is an anti-involution of T1. By step 1 it equals
     alpha sigma0 alpha^-1 with alpha in Aut(T1), and pi alpha lies in G and maps T1 to T2.
3. **The second binary digit.** By OCF (THM-466), H = 1 + 2 c_odd (mod 4). The map C -> reverse(sigma0 C) is an
   involution on the directed odd cycles.
   - Its fixed cycles have sigma0-invariant odd support, so they pass through v0. They are exactly
     v0 -> x_1 -> ... -> x_k -> sigma0(x_k) -> ... -> sigma0(x_1) -> v0, with x_1, ..., x_k in distinct pairs and
     x_k -> sigma0(x_k); the backward half is forced by anti-invariance.
   - Hence H = 1 + 2F (mod 4), with F = sum over x of P(x)[x -> sigma0 x]. Here P(x) is the number of directed paths
     v0 -> ... -> x with at least one arc that meet each pair at most once. This is Session29 section 3 modulo 4.
4. **The involution.** P never uses an arc inside a pair. Let O = {i : P(a_i) + P(b_i) odd}, and let Phi reverse the
   internal arcs of the pairs in O.
   - Phi maps X to X and fixes P, so Phi^2 = 1.
   - Phi commutes with G.
   - F(Phi T) - F(T) = |O| (mod 2).
5. **|O| is odd.** |O| = sum over x of P(x) (mod 2).
   - Group the paths by their sequence of pairs i_1, ..., i_k. The number of paths with that sequence is
     e^T M_(i1 i2) ... M_(i(k-1) ik) 1. Here M_ij is the 2x2 arc matrix from pair i to pair j, and e marks the unique
     vertex of pair i_1 that v0 beats.
   - Anti-invariance gives M_ij = [[x, y], [1 - y, 1 - x]], whose row sums are both = x + y = eps_ij (mod 2). So the
     count is the product of the eps along the sequence (mod 2).
   - Also eps_ji = 1 + eps_ij (mod 2). So eps defines a tournament Q on the m pairs, and the sum is the number of
     directed paths of Q, which is sum over nonempty U of H(Q[U]) = 2^m - 1 = 1 (mod 2) by Rédei.
6. **Conclusion.** Phi is a G-equivariant involution on X that exchanges H = 1 and H = 3 (mod 4). By step 2 it gives
   a bijection between the two kinds of self-converse class. QED.

Checks:
- `aud_sc_mod4.py`:
  - sum P odd, H = 1 + 2F (mod 4), Phi flips H mod 4, and Phi^2 = 1 hold on all of X for n = 3, 5, 7 (4, 64 and 4,096
    tournaments);
  - they also hold on random samples at n = 9 and 11.
- `aud_sc_orbits.py`: the G-orbits equal the self-converse classes (8 and 88) and split 4/4 and 44/44 at n = 5, 7.
- The census gives 1/1, 4/4, 44/44 and 1376/1376 at n = 3-9. At n = 11 the theorem predicts an even split. The OEIS
  value A002785(11) = 279,968 would give 139,984 each; that value was not recomputed here.
- Prior art: not found. THM-791's "odd-n contrast" records only an 8/8 split of the 16 grid-symmetric tilings at
  n = 5. Session29 gives H = 1 + 2F + 4I (mod 8) without the balance.

Recommended text:
- HYP-9512 (D) status becomes "PROVED for every odd n >= 3 (audit proof, THM-4922 (4))".
- THM-4922 gains "(4) **(PROVED)** At every odd order n >= 3 the self-converse classes split exactly evenly between
  H = 1 and H = 3 (mod 4): an equivariant involution reverses the internal arcs of the pairs of an anti-involution
  selected by a path-count parity."
- Remove "not found stated in the repository" from PROOF section 3 and add the proof.

## 5. The proofs in PROOF.md section 1

- **1.1.** Closure under converse is correct. Strongness has a gap (M2): only the one-vertex-component case is
  argued.
- **1.2.** The endpoint injection is correct: (v_1, ..., v_n) maps to (v_n, (v_1, ..., v_(n-1))). It is KNOWN (M3).
- **1.3.** The corollary is correct, but it needs padding. "Not attained at order n - 2" gives "least order n - 1"
  only because S(m) ⊆ S(m + 1) (m3). With that sentence the corollary holds for all n >= 3.

## 6. Order-10 code and pre-registration

**Logic of order10.c (read, not run).**
- The first-creator test is `!in9[h]`. The k test is `in9[g] && !in8[g]`. Because S(m) ⊆ S(m + 1) (pad with a sink),
  these are exactly "least order 10" and "least order 9".
- HM = 20000 exceeds max H(10) = 15745.
- The hypo count is taken over non-first-creators with k = 10. "Lineage" means some first creator of h has k > 0.
- The full independent recomputation agrees on every output line.

**Ordering (`prereg_check.txt`).**
- The registration commit abfe525263 (15:26:07) was rebased onto 2054513960 and pushed as 9e07033ca4 (15:26:15), and
  is on origin/main.
- order10.c, classgen.c, analyze.py, hypo.py and PREREGISTRATION.md are in that commit and unchanged since.
- The checkpoint 9d7c9f4735 is at 15:40:35. A 750 s run fits between the two commits, but git cannot prove when it
  ran: order10.out and PROOF.md share the mtime 15:39:37.
- The outcomes follow from order10.out as registered:
  - Q1: 0 orphans, no ORPHAN lines;
  - Q2: 0;
  - Q3: 62.84% > 50.44%;
  - Q4: 0.90%;
  - Q5: mixed 41.2%, FAILED.
- Q4's blindness is M4.

## 7. Prior art: what is known

- **Endpoint injection.** KNOWN: THM-4889 and klein section 4.
- **Genealogy threshold corollary.** NEW, elementary.
- **Every first creator is strong.**
  - KNOWN as FINITE-EXACT to order 11 (THM-4889 (2); klein: "no such product was ever new").
  - OPEN in general (HYP-9481 (b)).
  - The partial structure is KNOWN (Session224 section 3).
- **The dichotomy question.** Session224, same prompt, pushed first:
  - h = 9 counterexample, PROVED;
  - (s, p) ledger to order 8;
  - H = 33 = (1, 2), recovered from THM-4810 and S50.

  Related all-order strong-fiber results, which concern every order rather than the least order:
  - THM-4810 (A): strong fibres are singletons exactly at 3, 5, 11, 13 for v <= 74;
  - THM-4884 (2);
  - HYP-9477 (B): every h > 63 is attained by two non-isomorphic, non-converse strong tournaments.

  Other related work:
  - Session24 multiplicities (N, S, M per (n, H));
  - `h_spectrum_multiplicities_20261021.md`.

  NEW in S16: the (s, p) fibers at orders 9-10, and the counts 69/1,520 and 55/6,143.
- **Fiber sizes at least order.** KNOWN (klein's least_order_10.txt).
- **Odd-order self-converse balance.** Not found (see section 4). Now PROVED.
- **Even-order black cell.** KNOWN (THM-791 Cor. 2; Session29 section 4).
- **No hypo analogue and orphans.** NEW (finite). The tournament fact behind the analogy is classical and is in
  THM-4870's packet, section 2. That 63 is an order-7 hole is klein's and HYP-9359's.
- **Metagraph blue/black** (THM-2444, THM-2450, THM-2453, THM-2454, HYP-4997). These colour lines and classes, not
  values, so there is no overlap (m12).
- **Session219.** Its w(h), the maximum number of strong components, controls the degree in n of fixed-H class counts.
  It is unrelated to least-order fibers, so there is no overlap.

## 8. Files (audit/)

| File | Role |
| --- | --- |
| `aud_enum.c` | independent enumerator to order 10: canonical augmentation, own canonical form, \|Aut\| and orbits; built-in orbit-stabilizer and brute-force self-tests |
| `aud_analyze.py` | every comparison in section 1 (auditor data against packet claims, order10.out and klein's census_3..11) |
| `aud_sc_mod4.py`, `aud_sc_orbits.py` | checks of the section 4 proof |
| `aud_enum.out`, `aud_enum_brute7.out`, `aud_analyze.out`, `aud_sc_mod4.out` | outputs |
| `aud_fibers_le10.txt`, `aud_o10_byH.txt` | per-value fiber table to order 10; order-10 per-h table |
| `prereg_check.txt` | git provenance for section 6 |

Reproduce:

    gcc -O2 -fopenmp -o aud_enum aud_enum.c && ./aud_enum brute 7 && ./aud_enum 10 RUNDIR    # 2 min on 12 threads
    python3 aud_analyze.py RUNDIR .. ../../path_spectrum_two_angles_metagraph_20261106/data
    python3 aud_sc_mod4.py 200 && python3 aud_sc_orbits.py
