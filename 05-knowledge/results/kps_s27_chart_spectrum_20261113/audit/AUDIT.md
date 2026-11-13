# Independent audit: kps-S27 (THM-5152, HYP-9679, SELECTION_MAP) and kps-S28 (THM-5153)

Independent adversarial audit, 2026-11-13. Every repository file was read-only except this report. Nothing was
stashed, committed or pushed. **Chan–Pak remains OPEN.** Nothing in either packet claims otherwise.

## Verdicts

| object | verdict |
|---|---|
| THM-5152 and kps_s27_chart_spectrum_20261113 (README, PROOF, PREREGISTRATION, scripts) | **PASS WITH CORRECTIONS.** Two MAJOR corrections block canon as written. |
| HYP-9679 | PASS (one minor wording fix) |
| SELECTION_MAP.md | PASS WITH MINOR CORRECTIONS |
| THM-5153 and kps_s28_exact_pairing_sc_census_20261113 | PASS (minor corrections) |

**THM-5152 in one paragraph.** The computation is right: the 7,048 grammar classes and the 125-value miss list
through order 11 were reproduced exactly with independent code. Every one of the 18 listed values really is omitted
at every order. Proposition 4, however, is proved in only one direction, and the other direction is false, because
the chart class is not closed under ordinal sums. The packet itself observes this (13 at order 5 but not at order 6)
and then wrongly dismisses it. As a result the "exactly" list below 375 is incomplete. **215 = 5·43, 267 = 3·89 and
305 = 5·61 are also omitted at every order.** So "all except 341 and 343 are primes or 7p" and "below 225 every
omission is prime or 7p" are false.

---

# Part A: THM-5152 / kps-S27

## MAJOR findings

### M1. Proposition 4 is only an inclusion: the chart class is not closed under ordinal sums

**Exact sentences.**
- PROOF.md §1, Proposition 4: "The set of odd values below B attained by C is the multiplicative closure, inside
  [1, B), of the values attained by strong members of C of order at most n_0." The proof is "Combine Lemmas 2 and 3."
- THM-5152 title, (1): "whether C attains an odd value below B is decided by its strong members of order <= n_0,
  where f(n_0 + 1) >= B."
- THM-5152 body (line 45): "Hence the odd values below B attained by C form the multiplicative closure, within
  [1, B), of the values attained by strong members of order <= n_0, where f(n_0 + 1) >= B."
- README line 19: "Three facts make this decisive below 375: Busch's minimum for strong tournaments, heredity of the
  class, and multiplicativity of H."
- PROOF.md lines 84-85: "Proposition 4 is unaffected, because it only deletes singleton components."

**The problem.** Lemmas 2 and 3 prove only that the attained set lies inside the closure.
- Lemma 2 deletes singleton components and factorizes H over strong components.
- Lemma 3 bounds the order of each component.

The reverse inclusion needs more. Each product h_1 ⋯ h_k of values of strong chart members must be attained, so some
ordinal sum S_1 ⇒ … ⇒ S_k of such members must lie in C. C is not closed under ordinal sums. The scripts implement
the false direction: closure.py treats every product of attained values as attained.

**Evidence.** All tests use the auditor's own membership test, which runs over every switch with an independent
modular-decomposition test for G.
- For every strong order-5 chart tournament with H = 13 (18 labelled representatives), both T ⇒ z and z ⇒ T lie
  outside C.
- The values 43, 61 and 89 occur in C only at orders 6, 8 and 7 respectively. A strong member with one of these
  counts never survives an added singleton on either side.

**Fix.**
1. State Proposition 4 as an inclusion: the odd values below B attained by C lie in the multiplicative closure of
   the values of strong members of order ≤ n_0. Hence every odd v < B outside the closure is omitted at every order.
2. Add an exactness step. A value in the closure that is missed through order n_0 is attained iff two things hold
   for some factorization v = h_1 ⋯ h_k:
   - there are strong chart members S_i with H(S_i) = h_i;
   - some ordinal arrangement of the S_i lies in C.

   This is a finite check, since the total order is at most Σ n(h_i). By heredity it fails as soon as some S_i
   cannot take a singleton on the needed side.
3. In THM-5152 (1), replace "decided by its strong members" with: "omission is decided by its strong members of
   order ≤ n_0; attainment of a composite value additionally needs ordinal-sum membership checks."
4. Delete PROOF.md's sentence "Proposition 4 is unaffected …".

### M2. The omission list below 375 is incomplete: 215, 267 and 305 are omitted at every order

**Exact sentences.**
- THM-5152 title, (3): "Consequently C omits at every order exactly the odd values 7, 21, 31, …, 341, 343 below
  375; all except 341 = 11*31 and 343 = 7^3 are primes or 7p."
- THM-5152 §(3), its list, and "other: 341 = 11 · 31 and 343 = 7^3".
- PROOF.md §2 (lines 58-62): "Below f(12) = 375 the omissions are decided … 2 other values: 341 = 11 · 31 and
  343 = 7^3."
- PROOF.md §3 (lines 74-75): "Below 225 every omission is prime or 7p. Below 375 two more appear, 341 = 11 · 31 and
  343 = 7^3. In both, every factorization uses an omitted value (31 or 7), as Proposition 4 requires."
- README lines 21-25.
- PREREGISTRATION, Outcome (line 45).

**What was missed.** Through order 11, the odd values below 375 that are missed are the 18 listed values plus 125,
215, 267 and 305. The auditor reproduced this independently. closure.py marks the last four as attained, because each
is a product of chart values: 125 = 5·25 = 5³, 215 = 5·43, 267 = 3·89 and 305 = 5·61.

**Correct status of the four values.** In each proof, first delete singleton components; heredity keeps the result
in C.

| value | status | reason |
|---:|---|---|
| 125 | **attained**, only at order 12 | S4 ⇒ S4 ⇒ S4 is a grammar member, since G is closed under ordinal joins, and has H = 125 (confirmed). The sums S4 ⇒ T25 and T25 ⇒ S4 are never in C (42 representatives tested). The packet's conclusion is right, but its reason is not. |
| 215 = 5·43 | **omitted at every order** | A strong member with H = 215 would have order ≤ 10 (f(11) = 225), and none exists through order 11. Otherwise the components have H = 5 (order 4) and H = 43. Strong chart members with H = 43 exist only at order 6, since f(8) = 45 > 43 and order 7 has none. So the member has order 10, inside the data, which has no 215. Directly: 48 sums S4 ⇒ T43 and T43 ⇒ S4 were tested and none is in C. |
| 267 = 3·89 | **omitted at every order** | No strong member has H = 267 (order ≤ 11). Strong chart members with H = 89 exist only at order 7, with none at orders 8 or 9 (f(10) = 125). So C3 plus T89 has order 10, inside the data, which has no 267. Directly: 384 sums were tested and none is in C. |
| 305 = 5·61 | **omitted at every order** | No strong member has H = 305 (order ≤ 11). Strong chart members with H = 61 exist only at order 8, so S4 plus T61 has order 12, beyond the data. But deleting three vertices of the S4 would leave T61 with a singleton: an order-9 chart member with H = 61. 61 occurs only at order 8, so this is impossible. Directly: all 80 order-12 sums S4 ⇒ T61 and T61 ⇒ S4 were tested, and none is in C. |

**Corrected statement.** Below f(12) = 375, C omits at every order exactly these 21 odd values:

    7, 21, 31, 35, 53, 113, 127, 149, 173, 181, 203, 215, 217, 251, 263, 267, 277, 287, 305, 341, 343

- primes (11): 7, 31, 53, 113, 127, 149, 173, 181, 251, 263, 277;
- 7p (5): 21, 35, 203, 217, 287;
- other (5): 215 = 5·43, 267 = 3·89, 305 = 5·61, 341 = 11·31, 343 = 7³.

Every other odd value below 375 is attained through order 11, except 125, which is attained at order 12.

**Consequences to edit.**
- For 215, 267 and 305 every factor is a chart value. These omissions come from failed ordinal sums, not from an
  omitted factor. This is a qualitatively new kind of composite omission and worth a sentence in §3.
- "Below 225 every omission is prime or 7p" is false, because of 215.
- With 19 counts (21 values minus 7 and 21), the title becomes "omits 31 and eighteen more counts".
- **Unaffected:**
  - the prime omissions, the 7p omissions, 341 and 343, since each has no factorization into chart values;
  - all of HYP-9679's evidence;
  - SELECTION_MAP §3's sentence listing ten omitted values, which is true as stated.
- **P1** is still CONFIRMED: its prediction of at least one more omission in (225, 375) holds even more strongly, with
  eight. Append to the Outcome that the reported list was computed with the invalid closure rule and is completed
  here.

## MINOR findings

- **m1. Wrong order bound for 8191.** PROOF.md §3, lines 77-78: "… which needs strong members of order up to n with
  f(n) <= 8191, that is n <= 18." It should be n ≤ 17, since f(17) = 5,625 ≤ 8,191 < f(18) = 9,375. HYP-9679 already
  says 17.
- **m2. Title count.** "31 and sixteen more counts" was already off by one for the packet's own list, which has 16
  counts once 7 and 21 are removed: that is 31 and fifteen more. After M2 it should be "eighteen more". The file name
  will then disagree with the title; rename the file or add a note.
- **m3. Provenance of the independent check.** PROOF.md lines 67-68 cite "(sweep27c/sw.c)", and the THM status says
  the order-10 spectrum was "reproduced independently by this session's read-only sweep". That code is not in the
  packet or the repository; it exists only in a session scratch directory. Either commit it under `scripts/` or label
  the check as uncommitted.
- **m4. Credit for 13.** THM-5152 lines 13-15 and PROOF.md lines 14-15 say Session 268's census "names 13 as outside
  it". Session 268 names 13 only as the hostile for its spectrum lane's atom bank (5), h(b,c) = 2 + (2^b − 1)(2^c − 1),
  and its product closure (spectrum/PROOF.md §6). The absence of 13 from the whole grammar is true and can be read
  from Session 268's modular H histograms. The auditor's unswitched spectrum misses 13 through order 12, and since 13
  is prime this decides it. Rephrase accordingly.
- **m5. The mechanism of Proposition 4 is KNOWN.** Strong floor, multiplicativity and singleton deletion are already
  used to turn finite data into all-order statements in:
  - THM-1370's own proof for 21 ("21 = 3·7 needs the impossible 7");
  - THM-4873 (iii);
  - THM-5129 (A) and (C);
  - chanpak_55_105_golden_20261110/targets (P1-P4);
  - decomposition16_odd_spectrum_20261110 (J1).

  Credit these in THM-5152's status. README line 49 ("The nearest existing mechanism is the Busch-bound finiteness
  used in THM-5152") should name THM-5129 (C) first, as SELECTION_MAP §5 already does. What is new here is the
  application to the chart class and, after M1, the one-sided form with membership checks.
- **m6. The tricolour rule is a candidate, not the owner's construction.** SELECTION_MAP §4, lines 49-56, presents
  the lowest-Zeckendorf-digit colouring ("This is the Wythoff AA/B/AB split") as "The owner's construction". The
  repository records this rule as a near-fitting candidate decoding of the owner's colour word: it matches the first
  34 colours and fails at colour 35. Sources: collatz_crossings_20260926_potential_and_seeds.md §1.1,
  crossroads_crossing_20260926_colour.md (C1)/(C2), and Session 313 PREFIX_SELECTION. Rephrase as "the repository's
  candidate decoding (Wythoff AA/B/AB), which agrees with the supplied colours for d ≤ 34 but not at 35". README line
  45 ("your September 26 diagonal colouring") needs the same caveat.
- **m7. Composite omissions need the membership step.** SELECTION_MAP §5, lines 89-91: "this packet's chart-class
  omissions are decided below Busch's bound." Change this to: "its prime and 7p omissions are decided by Busch's
  bound alone; composite values need the ordinal-sum membership step." The final paragraph's claim, that containing
  or omitting a prime p is decidable, is correct as stated.
- **m8. HYP-9679 title and body disagree.** The title states (A) as "relative density strictly less than 1 among
  primes, and omits infinitely many 7p values". The body states (A) as "infinitely many" and gives density < 1 as
  the "Stronger form". Align them, for example as (A) infinitely many and (A′) density < 1.
- **m9. Range source for HYP-9661.** The SELECTION_MAP HYP-9661 row says "r(N) >= 1 suffices; FE to 382,687". The
  range comes from THM-5129 (C) and THM-4865 (e): the two-step family misses only 7, 11, 13, 21, 47 and 67 below
  382,688, which HYP-9661 restates. Cite those. HYP-9661's own FINITE-EXACT item (A) is r(N) ≥ N/64 on [129, 16384).
- **m10. verify_s27.py uses the same closure rule.** Check 3 is still correct below 75: with order ≤ 8 data no value
  below 75 is closure-only, and the same holds below 125 with order ≤ 9 data. Add a comment, or better, assert that
  the closure-only set is empty, so the verifier does not appear to certify the general rule. Check 4's message "so
  31 is a count outside the chart class" relies on check 3. A direct membership test would make it self-standing;
  the auditor's test confirms this tournament is not in C.

## Checks that passed

1. **Heredity (Lemma 1).** Session 312's argument holds: a subset of a point configuration is a point configuration.
   Equivalently, G is hereditary (an induced piece of X ⇒ Y or C3[X,Y,Z] is again a join of induced pieces, or a
   single piece) and switching commutes with restriction.

   Sanity check of the characterization: 3,000 random configurations in Q_3² (orders 5-12) all gave members of C.
   Uniform random tournaments mostly did not (order 8: 12 of 124; order 9: 3 of 119; order 10: 0 of 135). Session
   312's explicit chart x = (12,1,3,9,5,0,2,27) gives a strong grammar member with H = 95, whose Q equals that of
   Session 310's B.
2. **Lemma 2 and multiplicativity** (THM-3166 §1) are correct.
3. **Busch bound.**
   - min{3^a 5^b : 2a + 3b = n − 1} equals Busch's f(3q) = 3·5^(q−1), f(3q+1) = 5^q, f(3q+2) = 9·5^(q−1) for n ≥ 3.
     The values are 3, 5, 9, 15, 25, 45, 75, 125, 225, 375 and 625 for n = 3..13, strictly increasing.
   - Edge cases: n = 3 is C3 with H = 3; n = 4 is the unique strong 4-tournament with H = 5. THM-4744 cites Busch
     2006 and Moon 1972.
   - Off-by-one: for B = 375 = f(12), n_0 = 11, since f(11) = 225 < 375 ≤ f(12). The run is through order 11, so
     this is correct. THM's "f(n_0+1) ≥ B" and PROOF's "n_0 largest with f(n_0) < B" agree because f is increasing.
   - The chart class's own strong minima at orders 9, 10 and 11 are 81, 153 and 289, above Busch. This is no
     contradiction.
4. **gtrees.py enumeration** is correct, both by argument and by computation.
   - Gallai's decomposition is unique.
   - A non-strong tournament's flattened L node has its strong components, in order, as children.
   - A C3 node's children are exactly its maximal proper modules. Any module meeting two parts of a C3 quotient is
     the whole vertex set.
   - Rotation-only canonicalisation over all ordered triples includes both orientations.

   Computationally, the 7,048 tournaments printed by `gtrees.py 11` canonicalise to 7,048 distinct classes. This is
   the same set as the auditor's own generator, with no misses and no duplicates.
5. **swspec.c.** Switch sets avoid vertex 0, which is valid because S and V∖S give the same tournament. The pull DP
   reads only entries already written for the current tournament. Index bounds are fine for n ≤ 12 and H < bound.
6. **The order-6 example.** H(36,17,2,7,13,30) = 31. The rows form a valid tournament, which is strong and not in C:
   none of its 32 switches is a grammar member.
7. **Scope.** Chan–Pak is stated OPEN in README, PROOF and THM. FINITE-EXACT claims carry ranges (order ≤ 11, values
   below 3,000 or 375). HYP-9679 (B) is labelled "curiosity, weakly supported" and "rests on three data points".
   Nothing above 375 is claimed.
8. **Novelty and credits.**
   - No prior computation of the switching-closure spectrum was found. Searches covered "switching closure",
     "chartable", "chart class", "3-adic chart" and switch mentions in Sessions 268 and 310-315 and in
     symmetry_port_return.
   - The class is Session 312's, from the same day. Session 313 does not compute its spectrum.
   - Credits to Session 268 (the spectrum of G) and Session 312 (the class, characterization, heredity remark and
     D-boundary) are accurate, apart from m4.
   - The miss-list arithmetic in HYP-9679 is correct: 73 odd primes below 375, 11 omitted. The omitted primes meet
     every unit class mod 3, 4, 8 and 9. 63, 255, 511, 1023 and 2047 are attained.
9. **SELECTION_MAP citations.** Thirteen were spot-checked against the cited files.

   | citation | result |
   |---|---|
   | THM-4094 Theorem 6.1 | (21) P_odd∖{7} ⊆ S_str and (22) 7(P_odd∖{3}) ⊆ S_str; matches the map |
   | THM-4875 (D) | ≤ 2^32 except 7, 21; 2^35 in one unreplicated run (addendum) |
   | HYP-9494 (D) | 29 misses (largest 2647); F_≤2 misses 140 below 10^9, largest 3,780,967 |
   | THM-5148 (E) | L = 5, 63 blocks, no miss in (10391, 2^19] |
   | THM-5148 blocks | the map's claims also re-derived: for n ≥ 1, det M_w = [P even][Q even]·2^n ∈ {0, 2^n}; the empty word gives charpoly t³ − t² − t − 1 |
   | HYP-9656 | N(4095) = 145,433 (audit) |
   | THM-5129 (C) | omits 7, 11, 13, 21, 47, 67 at every order; complete on [3, 225) |
   | inverse65 ZECKENDORF.md | §Z2 "Prefix-correct append-only compilation is impossible"; (Z1) = (X−1)(X² − L_q X + (−1)^q) |
   | reset_20260926_colours.md §2 | two or three representations |
   | ZERO-SLOT.md | three or four representations |
   | THM-5126 | c ≤ 129 |
   | THM-5130 (E) | seed and anti-seed for every odd c ≤ 4097 |
   | THM-5151 (9) | cap(Q_N0) transfers |
   | Session 310 POSITIVE_CLOCK | mod 96 |

   Exceptions are m6 and m9.

## Recomputed independently

All auditor code is new. It lives in
`C:\Users\Eliott\AppData\Local\Temp\claude\2f1de76a73e81b8a\a7c1f371-6162-4db2-99bd-68121f2beaaf\scratchpad\audit27\`
and was not committed.

| check | method (auditor code) | result |
|---|---|---|
| G class counts, orders 1-12 | `agen.c`: closure under X ⇒ Y and C3[X,Y,Z], deduplicated by a generic adjacency canonical form (partition refinement plus full individualisation, least relabelled adjacency); no decomposition trees | 1, 1, 2, 4, 9, 24, 65, 182, 528, 1,558, 4,674, 14,225; equals the packet (orders 1-11), Session 268 (14,225 at order 12) and the species equation A = S/(1−S), S = x + (A³ + 2A(x³))/3 |
| canonical-form soundness | `acanon.c`: 5 random relabellings of each of 7,048 classes | 35,240 tests, 0 failures |
| gtrees.py class set | packet's `gtrees.py 11` output canonicalised by `acanon` | 7,048 distinct classes, identical to the auditor's set |
| chart spectrum, orders 1-11, every switch | `aspec.c`: forward ("push") subset DP; switch sets avoid the **last** vertex (the packet avoids vertex 0); Landau strongness; 5,747,451 switched tournaments; order 11 in 46 s on 12 processes | odd values below 3,000 missed through order 11: **identical** to `scripts/miss11_order_le_11_below_3000.txt` (125 values) |
| brute force, no class enumeration | `memb.c brute`: every labelled tournament of orders 5, 6 and 7 (1,024; 32,768; 2,097,152), membership tested over all switches with a modular-decomposition test for G | 1,024; 24,320; 727,040 members. Attained and strong-attained H sets equal the enumeration's at each order |
| multiplicative-closure rule | `analyze.py` and the packet's `closure.py` | the packet's rule reproduces its 18 values. Closure-only values below 375 are {125, 215, 267, 305}. Only these four have a factorization into strong chart values; every listed omission has none |
| ordinal-sum membership | `targeted.py` plus `memb test` | none of these is in C: T13 ⇒ z and z ⇒ T13 (18 each); T43, T61 and T89 with a sink or source (24, 40 and 192 each); S4 ⇒ T43 and T43 ⇒ S4 (24 each; H = 215); C3 ⇒ T89 and T89 ⇒ C3 (192 each; H = 267); S4 ⇒ T61 and T61 ⇒ S4 (40 each, order 12; H = 305); S4 ⇒ T25 and T25 ⇒ S4 (42 each; H = 125). S4 ⇒ S4 ⇒ S4 **is** in C (it is in G) with H = 125 |
| example with H = 31 | `memb test` | H = 31, strong, not in C |
| packet's own replay | `verify_s27.py` normally and with `-O`; `swspec.c` on `gtrees.py` output to order 9 (bound 375); `gram.py 9` against the auditor's unswitched spectrum (`unsw.c`) | all pass and agree. Unswitched G misses 13, 19, 29, … through order 12 |

---

# Part B: THM-5153 / kps-S28 (addendum requested by the coordinator)

Verdict: **PASS** with minor corrections. Every number checked is right.

## Recomputed independently

The auditor's own census, `tclass.c` plus `census5153.py`, covers orders 1-9 and was run end to end through
order 9:
- **Classes** by one-vertex extension with the auditor's canonical form. Counts are 1, 1, 2, 4, 12, 56, 456, 6,880
  and 191,536 (A000568).
- **Self-converse flag** by testing canon(T^op) = canon(T), with the converse class found by lookup.
- **|Aut|** as the number of optimal leaves.
- **Profiles** q_r, the number of ordered r-tuples of disjoint nonempty paths covering V, by the auditor's own subset
  convolution. These are fprof's F_r: for TT5 both give (1, 30, 150, 240, 120).
- **Strongness** by Landau's criterion.
- **Built-in checks**, all of which pass: q(T) = q(T^op), q_1 = H and q_n = n! for every class.

**1. Census table, orders 3-9.** All ten columns agree exactly at every order:

| n | classes | SC | merged | Q-classes | ≥2 merged | SC-SC | SC-NSC | strong-twin | with SC member | SC alone |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| 4 | 4 | 2 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 2 |
| 5 | 12 | 8 | 10 | 8 | 2 | 1 | 1 | 1 | 1 | 5 |
| 6 | 56 | 12 | 34 | 24 | 9 | 3 | 2 | 5 | 3 | 5 |
| 7 | 456 | 88 | 272 | 156 | 69 | 15 | 28 | 49 | 27 | 37 |
| 8 | 6,880 | 176 | 3,528 | 1,728 | 827 | 34 | 69 | 707 | 73 | 53 |
| 9 | 191,536 | 2,752 | 97,144 | 50,511 | 21,217 | 158 | 1,175 | 19,528 | 1,161 | 1,348 |

The task asked for orders 3-8; order 9 was also recomputed in full.

**2. Self-converse gaps.** The odd values below 75 attained by no self-converse class of order ≤ 9 are 7, 19, 21,
23, 35, 39, 43, 47, 53, 55, 59, 63. This equals Session 315's list (SPECTRUM_AND_SIEVE.md; scsmall.py asserts the
same set).

**3. THM-5151 (4).** The Q-classes mixing strong and non-strong classes through order 9 are:

| order | H | strong classes | self-converse | non-strong classes |
|---:|---:|---:|---|---:|
| 6 | 15 | 2 | none; a converse pair, rigid | 2 |
| 8 | 95 | 2 | none; a converse pair | 10 |
| 8 | 45 | 3 | one, plus a converse pair | 2 |

The seven strong members thus include exactly one self-converse class.
- Decoding classgen code 167804930 (upper-triangle bits) gives a strong order-8 tournament with H = 45. It is
  self-converse and rigid (|Aut| = 1), and it is the auditor's self-converse H = 45 class.
- Its Q, (45, 1350, 13140, 60120, 147600, 200160, 141120, 40320), equals Q(C3 ⇒ Y5') with Y5' = code 739.
- The H = 95 Q-class contains Session 312's chart for B, which has the same profile.

**4. Credits.**
- The S26 audit (kps_s26 audit/AUDIT.md, "Proposition 4 / THM-5151 (4)") records "2, 3, 8, 24, 156, 1728, 50,511".
- Session 261 MODULE_REBASE §3 records "The order6 census has 56 unlabelled classes and 24 full cover profiles".
- Session 315 SPECTRUM_AND_SIEVE.md gives the global self-converse gap list.

All three citations are accurate. No earlier self-converse refinement of the profile census was found in Sessions
261, 310 or 314-315 or in kps-S26, so the NEW claim is plausible.

## MINOR findings (THM-5153)

- **t1. The SC-SC and SC-NSC counts overlap.** The THM-5153 title says "Of these, 0, 0, 1, 3, 15, 34, 158 hold two or
  more self-converse classes, and 0, 0, 1, 2, 28, 69, 1175 mix a self-converse class with a converse pair". The
  README table headers "of which SC-SC | of which SC-NSC" read as a partition. The two counts overlap: 0, 0, 0, 1, 9,
  22 and 101 Q-classes (orders 3-9) do both. Say "(overlapping)" or add a "both" column.
- **t2. Credit for "invisible to every uniform exterior".** README lines 11-12 credit Session 310 NEUTRAL_TAILS for
  "every pairing is invisible to every uniform exterior". The contextual congruence (equal q implies equal H in every
  uniform substitution context) and its converse are proved earlier, in Session 268 modular/PROOF.md §1 (Corollary)
  and §5, dated 2026-11-11. NEUTRAL_TAILS itself credits Session 306 for the run interpretation. Add the Session 268
  credit.
- **t3. Unpinned replay inputs.** scpairs.py reads `cls_n.txt` and `prof_n.txt`, which are produced outside the
  packet by first_creators classgen and kps_s26 fprof; the kps_s26 `scripts/README.md` pointer does exist. Pin the
  inputs by stating the exact build commands and the expected line counts (2, 4, 12, 56, 456, 6,880 and 191,536 for
  n = 3..9) or checksums.

