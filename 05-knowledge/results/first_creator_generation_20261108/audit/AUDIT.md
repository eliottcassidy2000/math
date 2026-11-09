# Independent audit of THM-4942 and HYP-9532 (first_creator_generation_20261108)

Audit of kind-pasteur-2026-11-08-S17, 2026-11-08. Two phases:

- **Phase 1** (orders <= 10): [PREREGISTRATION.md](../PREREGISTRATION.md) (af0112b596),
  [PREREGISTRATION-2.md](../PREREGISTRATION-2.md) (c201ec7b34), and PROOF.md sections 1, 3.1-3.3, 4.1, 4.2 and 5.
- **Phase 2** (order 11): `witnesses_order11.txt`, `order11_summary.txt`, `lambda3_witnesses_order11.txt`,
  `lam11_summary.txt`, `lambda3_classes_order11.txt` and `lineage_le10.txt`. Also the order-11 text of PROOF.md
  (sections 2, 3.1 and 3.3), README.md, THM-4942 and HYP-9532, in the versions on disk at 20:00 on 2026-11-08.

**Independence.** All verification code in this folder was written for this audit. The packet's programs were read
only to compare definitions and were never run: ext.c, ext_wit.c, gen.c, lam.c, lam11.c, lam11_all.c, canon11.c,
dump10.c, mindeg10.c and the Python scripts.

One block is reused, as the brief allows. It is the canonical labelling (`cctx`, `refine`, `search`, `canon`,
`decode`, `converse`) of the previous independent audit's
`first_creators_genealogy_20261107/audit/aud_enum.c`. It is copied verbatim between marked comments and shares
nothing with the packet. Everything that bears on a claim is new:
- the Hamiltonian-path DP, the attachment and deletion forms, arc-strong connectivity and modules;
- the categories, statistics and witness checks;
- the isomorphism certificates.

The auditor's order-10 classes are its own enumeration, not the packet's `classes10.bin`. Its code format differs
from the packet's: the pair (0,1) is the most significant bit, where the packet uses bit 0.

## Verdict: PASS WITH CORRECTIONS

Every count, floor, witness and prediction outcome in the packet reproduces exactly at orders 3-11 (sections 2 and 3),
with one exception: the generation threshold k*(n) is 0, not 1, at orders 3-5 (m3). The attachment form and the three
short proofs (sections 1, 3.2 and 4.1) are correct.

Four MAJOR corrections are needed (section 4):
- **M1.** The "born at the Busch bound" list is FINITE-EXACT only for h <= 80,451, or more generally for h in S(11). It
  is not exact "through 95,095": 3,635 values in [80,453, 95,093] are undecided.
- **M2.** "The first creators are the members of Strong(h) at the least order" is not proved, yet it stands in a
  PROVED section. It is HYP-9481 (b), exact only through order 11. This repeats the slip corrected as M2 of
  MISTAKE-735.
- **M3.** Lemma 3.2 is KNOWN. THM-4811 (B) (2026-10-27) proves lambda = delta^0 whenever n <= 4 delta^0 + 1, with the
  same cut count. It also already records the per-lambda class counts through order 10 that section 3.3 rebuilds.
- **M4.** The attachment form of section 1 is KNOWN: HYP-2879 / THM-4097, sections 2-3, and THM-4114, section 6, all
  from 2026-08-25.
  - THM-4097 (4) is the same insertion formula.
  - Its eq. (5) builds the slot matrix Q_ab by exactly the packet's subset DP, so Q_ab = W[a][b].
  - Its cut-field identity (9) and THM-4114's (32)-(34) give the quadratic, submodular form with pair differences
    -[Q(u,v) + Q(v,u)] = -M_uv.
  - The coordinator found this after the audit's first draft, through the concurrent Session230. The auditor read the
    three sources and confirms it.

Eighteen minor corrections follow in section 5.

**One status can go up.** Class-level genealogy at order 11 is listed as NOT CLAIMED, with the remark that it "needs
canonical forms at order 11". It does not, and the audit settled it exhaustively (section 3.6).
- **Every first creator of order 11 has a first-creator deletion.** There are 0 class-level orphans among the
  506,912,096 (non-first-creator parent, pattern) pairs that give an order-11 first creator.
- Every order-10 first creator lies on a chain of first creators down to C3, and first creators are strong through
  order 11. So every order-11 first creator lies on a Moon chain of first creators, one member at each order 3-11.
- HYP-9532 (B) and the class-level form of HYP-9512 (A) therefore become FINITE-EXACT through order 11.

## 1. Verdicts

| Item | Claim | Verdict |
| --- | --- | --- |
| 1 | PROOF section 1: the attachment form and its three consequences | **CONFIRMED** (correct proof, 0 numerical failures) but **KNOWN** (M4). Minor: m1 (the bound on M_ab needs one more line), m2 (neutral-pattern convention) |
| 2 | Generation census at orders 3-10, and R6-R9 | **CONFIRMED**, every number exact |
| 3 | Floors f_1-f_4 for n <= 10; lambda >= 3 implies FC; non-FC counts with lambda = 2 and min degree 2 | **CONFIRMED**, every number exact. Minor: m3 (k*(n) at n <= 5) |
| 4 | Lemma 3.2 and its sharpness | **CONFIRMED** (correct proof, sharp at n = 14) but **KNOWN** (M3) |
| 5 | Section 4.1: antichains, W(h) in [3, B(h)] | **CONFIRMED** for the theorem. **CORRECTED** last sentence (M2); minor citation m4 |
| 6 | Section 4.2: the eleven values through 95,095 | **CORRECTED** (M1): exact through 80,451 |
| 7 | Section 5 dictionary and credits | **CONFIRMED** in substance. Credit gaps M3 and M4 (THM-4811; HYP-2879, THM-4097, THM-4114); minors m6-m9, m18 |
| P2-a | the 36,248 R1 witnesses | **CONFIRMED**, every line checked by direct DP |
| P2-b | R1-R5 at order 11 | **CONFIRMED** exhaustively, on an independent engine and independent classes, plus a brute-force sample |
| P2-c | lambda >= 3 at order 11: f_3(11) = 11,475, R10, R11 | **CONFIRMED** exhaustively and independently |
| P2-d | the order-11 grammar: 11 classes, SUB/EXT, C3[R5, R5, 1] | **CONFIRMED**. The logic is airtight (section 3.5). Minor wording: m10, m11 |
| P2-e | lineage ledger `lineage_le10.txt` | **CONFIRMED**: rebuilt from the auditor's classes, all 7,661 values agree (section 3.7) |
| new | class-level genealogy at order 11 (NOT CLAIMED by the packet) | **SETTLED by the audit**: every order-11 first creator has a first-creator deletion. Exhaustive, with positive controls (section 3.6) |
| docs | README, THM-4942, HYP-9532 | **CORRECTED**: M1, M3 and M4 propagate; minors m3, m5, m6, m9-m15, m18 |

## 2. Phase 1 recomputation

### 2.1 Enumeration controls (`aud_gen.c`, `aud_gen.out`)

- Classes at orders 1-10: 1, 1, 2, 4, 12, 56, 456, 6,880, 191,536, 9,733,056 (A000568).
- The orbit-stabilizer identity sum n!/|Aut| = 2^C(n,2) holds at every order. At order 10 the sum is
  35,184,372,088,832.
- The auditor's S(n), n = 1-10, is identical to the values with n(h) <= n in klein's `least_order_11.txt`. The sizes are
  1, 1, 2, 3, 7, 19, 77, 320, 1,520 and 7,663, and the maxima are 1, 1, 3, 5, 15, 45, 189, 661, 3,357 and 15,745
  (`aud_busch.out`).
- The order-10 run took 113 s on 12 threads.

### 2.2 Item 1: the attachment form (`aud_attach.c`, `aud_attach.out`)

**The proof is correct.** The expansion of x_a (1 - x_b) is checked in full in the source comments. The position
classes (v first, v last, v between a and b) are disjoint and exhaustive, and sum_b g[V][b] = H(T).

**Numerical check.**
- Coverage: 34,937 tournaments and 2,397,258 pattern checks, every pattern of every tested tournament:
  - all labelled tournaments of orders 1-6;
  - random tournaments of orders 7-10 (600, 300, 120, 50), uniform and biased toward transitive, so that 1-10 strong
    components occur.
- Method: direct subset DP of the (n + 1)-vertex tournament, against the form. **0 failures.**
- The same run checks:
  - L_a = H(T + v_{e_a}) - H(T) >= 0;
  - the insertion lemma, H(T + v_x) >= H(T);
  - the H-neutral patterns: exactly the x whose support {a : a -> v} is a union of initial strong components, and
    #SCC + 1 of them;
  - W[a][b] = #Hamiltonian paths of T + v_x containing a, v, b consecutively, for x with x_a = 1 and x_b = 0. This
    uses a pruned-digraph DP independent of f and g, on every labelled tournament of order <= 4 and one random
    tournament in ten.
- Largest coefficients seen at n = 10: W 5,794, L 34,104, M 6,312, against max H(11) = 95,095.

**The bound "every coefficient is at most max H(n + 1)" is true**, but the reason given covers only W and L. A proof
for M_ab is in m1.

**Literature and canon.**
- **The form itself is KNOWN (M4).** THM-4097 eq. (4) (from HYP-2879) is
  H(T + x_S) = sum_(b in S) s_b + sum_(a notin S) e_a + sum_(a notin S, b in S) Q_ab, with S the set the new vertex
  beats. Its eq. (5), Q_ab = sum_(A: a in A, b notin A) E_A(a) B_(V\A)(b), is W[a][b].
- With x_a = 1 - [a in S], eq. (4) expands to the packet's (1). THM-4097 (6a) is the second difference -2 w_ij =
  -M_ij. THM-4114 (32)-(34) give the directed-cut form, the pair differences -[Q(u,v) + Q(v,u)] and the vanishing of
  all higher differences, that is, degree two.
- The split products E_a(S + a) B_b(R + b) are also the building blocks of the transfer matrix of THM-030 and of
  THM-050. When a -> b, W[a][b] is the number of Hamiltonian paths of T through the arc a -> b (THM-050, case 1,
  without signs).
- The packet's own contribution in section 1 is the Gray-code engine and its exhaustive use at order 11. THM-4097
  already used the same Q-compiler exhaustively over order-8 strong parents.

### 2.3 Item 2: the generation census (`aud_gen.out`)

| n | classes | FC | EXT | SUM | SUB | REST |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 2 | 1 | 0 | 1 | 0 | 0 |
| 4 | 4 | 1 | 2 | 1 | 0 | 0 |
| 5 | 12 | 6 | 2 | 4 | 0 | 0 |
| 6 | 56 | 33 | 14 | 9 | 0 | 0 |
| 7 | 456 | 328 | 90 | 37 | 1 | 0 |
| 8 | 6,880 | 5,457 | 1,162 | 216 | 33 | 12 |
| 9 | 191,536 | 163,676 | 23,855 | 2,489 | 1,010 | 506 |
| 10 | 9,733,056 | 8,647,232 | 962,686 | 49,755 | 42,514 | 30,869 |

All identical to the packet. The definitions are taken from the text: EXT uses S(n - 2) for the deletions, and modules
are found by brute force over all subsets, a different method from the packet's pair closure.

**REST.**
- Pendant classes: 12/12 at order 8, 505/506 at order 9, and 30,525/30,869 = 98.8856% at order 10.
- At order 10, the 344 non-pendant REST classes split into:
  - 4 with lambda = 1, the "single-back-arc grafts";
  - 340 with lambda = 2, of which 120 have a non-strong deletion and 220 do not.
- REST classes with a non-strong deletion: 12/12, 506/506 and 30,649/30,869. The other 220 all have lambda = 2.
- REST lambda histogram at order 10: lambda = 1: 30,529; lambda = 2: 340.

**The order-9 exception** has H = 385 and is self-converse, with |Aut| = 1, lambda = 2 and min degree 2. Its unique
score-4 vertex v is the middle vertex: T - v = P => Q, with P and Q strong 4-tournaments (H = 5 each, H(T - v) = 25),
and v beats exactly two vertices of each (`aud_h385.out`). It is the only class attaining f_2(9) = 385.

**Pre-registered outcomes.**

| ID | Outcome |
| --- | --- |
| R6 | 30,869 / 1,085,824 = 2.8429%: CONFIRMED |
| R7 | 98.8856%: FAILED |
| R8 | 344 > 100: FAILED |
| R9 | 962,686 / 1,085,824 = 88.6595%: CONFIRMED |

**The registration's "data seen" (orders <= 10) also reproduce.**
- Value genealogy: the only new value without a first-creator parent one order down is 63, at order 8 (orders 4-10).
- Class genealogy: k = 0 occurs only for C3 at order 3 and for the 63 pair at order 8. At order 9, 152 first creators
  have a creator of 63 as a deletion, and all 152 have another first-creator deletion (`aud_o11_fert_root63.out`).
  Every first creator through order 10 is strong.
- 0 hypo analogues at orders 3-10.
- The only-first-creator share reproduces by two independent methods: the auditor's extension engine, and the
  census's fully non-critical fibers. The values are 3/12, 18/58, 85/243 = 0.3498, 510/1,200 = 0.4250 and
  2,886/6,143 = 0.4698 at orders 6-10.
- Extending all classes of orders 3-9 reproduces S(4)-S(10) exactly.
- Fertility at orders 3-9:
  - no first creator is sterile;
  - the only first creators below half are the 63 pair, at 77/256;
  - the minima are 6/8, 14/16, 26/32, 51/64, 75/128, 77/256 and **299/512** (order 9, H = 615).

### 2.4 Item 3: the connectivity floors (`aud_gen.out`)

| n | f_1 | f_2 | f_3 | f_4 | max H(n - 1) | non-FC, lambda = 1 | non-FC, lambda = 2 | non-FC, min degree 2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 3 | | | | 1 | 0 | 0 | 0 |
| 4 | 5 | | | | 3 | 0 | 0 | 0 |
| 5 | 9 | 15 | | | 5 | 0 | 0 | 0 |
| 6 | 15 | 41 | | | 15 | 2 | 0 | 0 |
| 7 | 25 | 99 | 171 | | 45 | 25 | 0 | 0 |
| 8 | 45 | 195 | 621 | | 189 | 551 | 0 | 0 |
| 9 | 75 | 385 | 2,055 | 3,159 | 661 | 14,333 | 124 | 124 |
| 10 | 125 | 727 | 5,753 | 14,877 | 3,357 | 668,316 | 40,401 | 40,418 |

All equal the packet's values.
- No class with lambda >= 3, or with min degree >= 3, is a non-first-creator at orders 3-10.
- For every class at orders <= 10, (lambda >= 3) = (min degree >= 3). The lambda = 3 and 4 counts equal the min-degree
  3 and 4 counts: 10,237 and 15 at order 9, 1,400,468 and 13,333 at order 10.
- The 40,418 non-first-creators with min degree 2 at order 10 are:
  - the 40,401 with lambda = 2;
  - 16 with lambda = 1;
  - 1 with lambda = 0, which is R5 => R5.
- The per-lambda class counts at orders 5-10 equal those already listed in THM-4811's status (M3). For example,
  2,817,450, 5,124,698, 1,400,468 and 13,333 at order 10.

### 2.5 Item 4: Lemma 3.2 (`aud_lemma.py`, `aud_lemma.out`)

**The proof is correct.** For |X| = k, the arcs into X number sum indeg - C(k, 2) >= k(7 - k)/2 >= 3 for 1 <= k <= 6.
For n <= 13, one side of every cut has at most 6 vertices.

**The threshold is sharp at n = 14.** Take R7 => R7, with R7 the circulant on Z_7 with connection set {1, 2, 3}. It
has min degree 3 and lambda = 0. With one back arc lambda = 1, and with two back arcs lambda = 2, while min degree
stays 3.

The general form is THM-4811 (B): lambda = delta^0 whenever n <= 4 delta^0 + 1, and the first separation is two
R_(2 delta^0 + 1) blocks at n = 4 delta^0 + 2. 137 random tournaments with n <= 4 delta^0 + 1 at orders 5-13 all have
lambda = delta^0.

### 2.6 Item 5: section 4.1

**The proof is correct.** B strong with |B| >= 3 has no singleton strong component. THM-1900 (B), equality iff the
added vertex is a singleton component, gives H(B - v) < H(B). The insertion lemma (THM-1370 (2)), iterated, gives
H(A) <= H(B - v). The orders of Strong(h) are bounded by B(h) by the definition of B and f.

The citation and the last sentence need correction: m4 and M2.

### 2.7 Item 6: section 4.2 (`aud_busch.py`, `aud_busch.out`)

**The list reproduces.** Among the odd h with n(h) <= 11, those with n(h) = B(h) are exactly 3, 5, 9, 11, 13, 17,
19, 23, 35, 39 and 63. Equivalently, they are klein's length-1 strong windows below 375 that end at B(h): klein has
14 length-1 windows, and the other three (47, 55, 67) end at B(h) - 1.

**The range is wrong.**
- The odd values below 80,453 other than 7 and 21 all lie in S(11).
- The 3,635 odd values in [80,453, 95,093] that are absent from S(11) all have B(h) = 22 (f(22) = 78,125) and
  n(h) >= 12.
- For them, n(h) = B(h) = 22 is not excluded by any census. See M1.

## 3. Phase 2 recomputation

### 3.1 P2-a: the R1 witnesses (`aud_o11_check.c` mode `wit`, `aud_o11_witness.out`)

All 36,248 lines were checked in the packet's bit order:
- H(T) not in S(9), by direct DP of T;
- H(T + v_x) = h, by direct DP of the 11-vertex tournament;
- h is new at order 11 (klein's n(h) = 11, which agrees with S(11) minus the auditor's S(10)).

There are no duplicates, and no new value lacks a line. **0 failures.**

### 3.2 P2-b: R1-R5, exhaustive (`aud_o11.c`, `aud_o11.out`)

The auditor's engine is the attachment form plus deletion forms built from the same f/g tables restricted to V \ u.
It ran on the auditor's own 9,733,056 order-10 classes, all 1,024 patterns each: 9,966,649,344 extensions, 439 s on
12 threads. A first run gave identical results; two report-only counters had overflowed a 32-bit `long`, and the
committed run is after that fix.

| Claim | Auditor |
| --- | --- |
| S(11) equals klein's | yes: 43,911 values, 0 missing, 0 extra |
| R1: 36,248 new values, none without a first-creator parent | 36,248; 0 |
| R2: 0 hypo pairs | 0 |
| R3: 0 sterile first creators | 0 |
| R4: min fertile 564 at H = 2,901 | 564, attained by exactly two classes: a converse pair with \|Aut\| = 1, H = 2,901. These are the packet's codes 35572154624 and 18176436867200, canonicalised and brute-forced (564 each). The packet's third code 4956527526144 has 568, the next value of the distribution. 0 first creators below 512 |
| R5: 18,661 only-first-creator values | 18,661 / 36,248 = 0.5148 |
| multiplicities 1 to 445,304, median 251,107 | min 1, max 445,304. The two middle values are 251,103 and 251,107 (36,248 is even) |
| three values with a single first-creator pair | exactly 92,411, 93,467 and 95,095 (A = F = 1); next come 93,275 and 93,955 (3 each) |

**Brute-force sample.** Every 2,500th class was also brute-forced: 3,894 classes, 3,447 of them first creators.
- Every pattern was checked by direct DP of the 11-vertex tournament: H and all 11 deletions, 3,987,456 pattern checks.
- 0 mismatches. The fertile counts and the hypo flags are identical, and the sample contains 0 hypo pairs.

### 3.3 P2-c: lambda >= 3 at order 11 (`aud_o11.out`, `aud_o11_lambda3.txt`)

**The search is exhaustive, and the audit found the same 97 pairs.**
- It extends all 6,538,516 classes with min degree >= 2 by every pattern that gives min degree >= 3.
- This parent count equals the packet's "parents used" total. It is THM-4811's 6,538,499 2-arc-strong classes plus
  the 17 classes with min degree 2 > lambda.
- The degree-valid patterns number 2,834,126,741.
- The least H over **all** of them is 11,475, with no 15,745 cutoff. So f_3(11) = 11,475 exactly.
- With H <= 15,745 the search finds 97 pairs, by value:

  | H | 11,475 | 14,787 | 15,213 | 15,273 | 15,561 | 15,637 | 15,703 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | pairs | 3 | 7 | 18 | 18 | 18 | 11 | 22 |

- Every pair has lambda = 3 by brute force over all 2^11 - 2 cuts, and H confirmed by direct DP.
- 11,475, 15,213, 15,273 and 15,561 lie in S(10). R10 and R11 are CONFIRMED.

**Why exhaustive.** Any deletion of a min-degree->=3 11-tournament has min degree >= 2, so it is one of the
enumerated classes. The search does not even need Lemma 3.2: lambda <= min degree always holds, and lambda is tested
directly.

**Protocol deviation (m12).** PREREGISTRATION-2 registered parents with lambda >= 2. That scope would not have been
provably exhaustive, because of the 17 classes. The executed scope, min degree >= 2, is.

### 3.4 P2-d: the eleven classes (`aud_o11_check.c` mode `lam3`, `aud_o11_lambda3_classes.out`)

The packet's 97 pairs and the auditor's 97 pairs give the same 11 canonical classes, with identical pair counts per
class. Every merge was certified by an explicit isomorphism, the composition of the two canonical labellings checked
arc by arc: 0 failures.

| H | status | classes | \|Aut\| | SC | module | first-creator deletions |
| --- | --- | --- | --- | --- | --- | --- |
| 11,475 | non-FC | 1 | 25 | yes | size 5 | **0** (deletion values 2,985 x 10 and 225, all in S(9)) |
| 15,213 | non-FC | 2 | 3, 3 | no | size 3 | 8 each |
| 15,273 | non-FC | 2 | 3, 3 | no | size 3 | 8 each |
| 15,561 | non-FC | 2 | 3, 3 | no | size 3 | 8 each |
| 14,787 | FC | 1 | 9 | yes | size 3 | 8 |
| 15,637 | FC | 1 | 1 | yes | prime | 8 |
| 15,703 | FC | 2 | 1, 1 | no | prime | 8 each |

**C3[R5, R5, 1]**, built from its definition (A => B => v => A, with A and B the circulant R5), has H = 11,475,
lambda = 3, min degree 3, |Aut| = 25 and is self-converse. It is certified isomorphic to the H = 11,475 class.

### 3.5 The logic of the order-11 grammar (requested check)

The statement: every order-11 class is FC, EXT, SUM, SUB, or has lambda <= 2. It holds, and the argument is airtight:

1. A non-first-creator of order 11 has H in S(10), so H <= 15,745.
2. If moreover lambda >= 3, then min degree >= 3, since lambda <= min degree. So it lies in the exhaustive list of 3.3
   (H <= 15,745, min degree >= 3), which contains exactly the 7 non-FC classes above.
3. All 7 have a nontrivial module, so they are SUB; six are also EXT.
4. **lambda <= 2 means "an ordinal sum re-connected by at most two arcs", in both directions and for strong and
   non-strong classes alike.**
   - lambda = 0 is exactly "not strong", which is an ordinal sum X => V \ X with both parts nonempty.
   - For lambda in {1, 2}, take a cut X with exactly lambda arcs from V \ X into X. Reversing them gives X => V \ X.
   - Conversely, reversing at most two arcs of an ordinal sum leaves a cut crossed backward by at most two arcs.

The content of the grammar at order 11 is therefore precisely "every non-first-creator with lambda >= 3 is a
substitution", and that is verified. Two wording issues remain (m10, m11).

### 3.6 Class-level genealogy at order 11 (new; `aud_o11_cls.c`, `aud_o11_cls.out`)

Let T' be an order-11 first creator with no first-creator deletion. Then T = T' - v is a non-first-creator for every
v, so T' arises from one of the 1,085,824 non-first-creator classes of order 10. Extending all of them by all
patterns, and testing every deletion of each extension that is a first creator, is therefore exhaustive. No order-11
canonical form is needed.

**Positive controls**, with the same program compiled at smaller orders:
- parents of order 7 find exactly the 63 pair at order 8: 16 (parent, pattern) pairs, all with H = 63;
- parents of order 9 find nothing at order 10. Both agree with THM-4922.

**Result.**
- The run covered 1,085,824 non-first-creator parents, and 506,912,096 (parent, pattern) pairs give an order-11 first
  creator. **0** of them lack a first-creator deletion (414 s on 4 threads).
- Hence every first creator of order 11 has a first-creator deletion.
- Section 2.3 shows that every order-10 first creator has a chain of first creators down to C3. At order 9 the chain
  avoids the 63 pair by the 152-check.
- First creators are strong through order 11 (THM-4922 (1)). So every order-11 first creator lies on a Moon chain
  inside the minimal set, one member at each order 3-11.

This makes HYP-9532 (B) FINITE-EXACT through order 11, and the class-level form of HYP-9512 (A) holds at order 11. It
implies R1 at order 11.

The result rests on the auditor's engine. That engine was validated by the 3,894-class brute-force sample of 3.2 and by
the two controls. The result was not pre-registered.

### 3.7 P2-e: the lineage ledger (`aud_o11_check.c` mode `lineage`, `aud_lineage_le10.txt`, `aud_lineage_cmp.out`)

**Method.** The ledger was rebuilt from the auditor's own classes: orders 4-9 regenerated, order 10 from the auditor's
codes. For each value h with 4 <= n(h) <= 10 it records:
- the number of first-creator classes;
- the number of those with a first-creator deletion;
- the set of H values of those deletions.

**Result.**
- All 7,661 values agree with `lineage_le10.txt`.
- 196 lines are compared exactly. The other 7,465 are abbreviated in the packet (head, "...", tail, "(N values)") and
  are compared by head, tail and count. There are 0 disagreements.
- The only orphan is 63. It is also the only value whose first creators are not all grown from first creators.

## 4. MAJOR issues

### M1. Section 4.2 is exact through 80,451, not through 95,095

The list depends on n(h), and n(h) is known only for h in S(11). Odd h <= 95,095 fall into three cases:
- **h < 375 = f(12).** B(h) <= 11 is decided by the censuses.
- **h in S(11) with h >= 375.** Then n(h) <= 11 < 12 <= B(h).
- **The 3,635 values of [80,453, 95,093] absent from S(11).** Here B(h) = 22 and n(h) >= 12, and no census decides
  whether n(h) = 22.

HYP-9481 (a) would give n(h) <= m(h) + 1 = 12. So for these values the statement is conditional, exactly as it is
above 95,095.

**Replacement text** for PROOF section 4.2, second sentence:

> Among the odd h <= 80,451 (every one except 7 and 21 lies in S(11)), and more generally among all odd h in S(11),
> these are exactly 3, 5, 9, 11, 13, 17, 19, 23, 35, 39, 63. Equivalently, these are klein's length-1 strong windows
> below 375 that end at B(h); the other three length-1 windows, 47, 55 and 67, end at B(h) - 1. The 3,635 odd values
> in [80,453, 95,093] absent from S(11) have B(h) = 22 and n(h) >= 12, so for them n(h) != B(h) is conditional on
> HYP-9481 (a), which gives n(h) = 12.

The same scope applies wherever the claim appears. Replace "through 95,095" or "up to 95,095" by "through 80,451
(conditional on HYP-9481 (a) for the 3,635 values of [80,453, 95,093] outside S(11))" in:
- the THM-4942 title, status ("the eleven values born at their Busch bound through 95,095") and item 6;
- README item 5;
- the HYP-9532 (E) status and body ("Through 95,095, n(h) = B(h) holds exactly for ...").

The HYP-9532 (E) status then reads: "FINITE-EXACT through 80,451; conditional on HYP-9481 (a) for the 3,635 values of
[80,453, 95,093] outside S(11) and for all h > 95,095."

### M2. Section 4.1: first creators = Strong(h) at the least order is not proved

The section is typed PROVED, but its last sentence asserts that every first creator is strong. That is THM-4922 (1),
FINITE-EXACT through order 11, and HYP-9481 (b), OPEN in general. MISTAKE-735 already corrected the same slip (as M2
of that audit).

**Replacement text:**

> So every tournament with H = h is an ordinal sum of singletons and members of the finite sets Strong(d), d >= 3,
> over the ordered factorizations of h. Through order 11 the first creators of h are exactly the members of Strong(h)
> of order n(h) (THM-4922 (1), FINITE-EXACT); in general this is HYP-9481 (b), OPEN.

### M3. Lemma 3.2 is KNOWN (THM-4811 (B)); the floors table rebuilds THM-4811's counts

THM-4811 (B) (2026-10-27, independently audited) states and proves:

> lambda = min(delta^+, delta^-) whenever n <= 4 delta^0 + 1: a cut X of size k >= 2 has
> f(X) >= k delta^+ - C(k,2) and f(X) >= (n - k) delta^- - C(n - k, 2), and f(X) < delta^0 forces
> k >= 2 delta^0 + 1 and n - k >= 2 delta^0 + 1.

It also records:
- the first separation at n = 4 delta^0 + 2;
- exhaustive lambda = 2 < delta^0 searches at orders 14-17;
- the class counts by lambda at orders 5-10, which are the counts behind section 3.3;
- that lambda is exactly the distance to the non-strong locus.

The packet's lemma is the case delta^0 = 3, with the same proof, and is presented without credit:
- PROOF section 3.2 says PROVED;
- the THM-4942 status says "PROVED (elementary): for n <= 13, ...";
- README item 4 says "which is proved".

**Replacement text** for PROOF section 3.2:

> **Lemma (KNOWN: THM-4811 (B), the case delta^0 = 3).** For n <= 13, lambda(T) >= 3 if and only if every vertex has
> in- and out-degree at least 3. [proof as written] The bound is sharp: R7 => R7 (n = 14) has min degree 3 and
> lambda = 0, and with one or two back arcs lambda = 1 or 2. THM-4811 proves lambda = delta^0 whenever
> n <= 4 delta^0 + 1, and its census lists the class counts by lambda through order 10 reproduced in section 3.3.

Make the matching changes elsewhere:
- In the THM-4942 status, move the lemma from "PROVED" to "KNOWN and credited: THM-4811 (B)".
- Change README item 4 to "(THM-4811 (B): for n <= 13, lambda >= 3 iff min degree >= 3)".
- Add THM-4811 to the PROOF credit table and to THM-4942 `related`.

The new content of section 3 is unaffected: the floors f_2-f_4, "lambda >= 3 implies first creator" through order 10,
and its failure at order 11.

### M4. Section 1 (the attachment form) is KNOWN: HYP-2879, THM-4097 sections 2-3, THM-4114 section 6

The coordinator found this item after the audit's first draft, through the concurrent Session230
(`session230_creator_grammar_20261109`). The auditor read the sources and confirms it.

- **THM-4097 section 2** (2026-08-25) gives the insertion formula (4), from HYP-2879:

      H(T + x_S) = sum_(b in S) s_b + sum_(a notin S) e_a + sum_(a notin S, b in S) Q_ab,

  where x -> v iff v in S. It also gives the subset-DP compiler (5),
  Q_ab = sum_(A: a in A, b notin A) E_A(a) B_(V\A)(b), which is exactly W[a][b], computed in O(n^2 2^n) per parent.
  Substituting x_a = 1 - [a in S] turns (4) into the packet's (1).
- **THM-4097 section 3** gives the cut-field identity (9), H(T + x_S) = H(T) + cut_w(S) + sum_(i in S) h_i with
  w_ij = (Q_ij + Q_ji)/2 = M_ij/2. Its (6a) is the Boolean second difference -2 w_ij = -M_ij.
- **THM-4114 section 6** states that H(T + x_S) is the capacity of a directed s-t cut with capacities Start, End and
  Q (32). The pair differences are -[Q(u,v) + Q(v,u)] = -M_uv (33), and every difference of order >= 3 vanishes (34).
  So the response is a submodular quadratic: this is the "degree two" the packet calls folklore.
- **THM-4118** has the ear-response lattice.

None of these is cited in the audited checkpoint (HEAD f870c2abf1). Section 1 is typed PROVED, and THM-4942 lists
"the attachment form" first under "PROVED (elementary)".

**Replacement text.**
- Section 1 heading: "The attachment form (KNOWN: HYP-2879, THM-4097 sections 2-3, THM-4114 section 6; restated in
  pattern coordinates)".
- First sentence: "THM-4097 eq. (4)-(5) (from HYP-2879) give H(T + x_S) with the subset-DP slot matrix
  Q_ab = W[a][b]; in the coordinates x_a = [a -> v] = 1 - [a in S] this reads (1)."
- Replace "Degree two in x is folklore ... (Session189's response polynomials)" with "Degree two, with pair
  differences -M_ab, is THM-4114 (33)-(34); THM-4097 (9) is the cut-field form."
- State the packet's own contribution as: "new here: the Gray-code walk with O(n) updates and its exhaustive use on all
  9,733,056 order-10 classes, including the deletion forms." The neutral-pattern consequence stays credited to
  THM-1900 (B). The M_ab bound of m1 can stay as a remark.
- THM-4942: move "the attachment form" from PROVED to "KNOWN and credited: HYP-2879, THM-4097, THM-4114". Change the
  README item 1 wording "is an explicit quadratic form ..." to credit THM-4097.

## 5. Minor issues

**m1. Section 1, consequence 1.** W and L are bounded as stated; M_ab needs an argument. Add:

> For M_ab: if b -> a, take x with x_a = 1, x_b = 0. W[a][b] counts the paths of T + v_x containing a v b. Each pair
> counted by W[b][a] concatenates (b -> a) to a Hamiltonian path of T containing "b a". Inserting v at a feasible
> place never puts it between b and a, since v -> b. These paths are therefore distinct from the first family and
> from each other. So M_ab <= H(T + v_x) <= max H(n + 1). If a -> b, exchange the roles.

**m2. Section 1, consequence 3: the convention is the complement of THM-1900's.** In THM-1900 (B), P is the set the
new vertex beats, and a down-set is a suffix C_j u ... u C_k of terminal components. With x_a = [a -> v], P is the
zero set of x. Replace with:

> The pattern x is H-neutral exactly when the set {b : x_b = 0} that v beats is a condensation down-set (a union of
> terminal strong components), equivalently {a : x_a = 1} is a union of initial components; there are #SCC(T) + 1
> such patterns (THM-1900 (B)).

**m3. k*(n) at small orders.** At orders 3-5 every non-first-creator is non-strong (2.4), so k*(3) = k*(4) = k*(5) =
0. The statement appears in PROOF section 3.3 ("1 through order 8"), in the HYP-9532 status ("k*(n) = 1 for n <= 8")
and in the HYP-9532 (D) body. Replace with:

> k*(n) = 0 for n <= 5, 1 for n = 6, 7, 8, 2 for n = 9, 10, and at least 3 at n = 11.

The sentence "through order 8, within one [arc reversal]" stays correct.

**m4. Section 4.1 citation and edge case.**
- n <= B(h) follows from the definition of B. The finiteness of W(h) needs f(n) -> infinity, which comes from Busch
  2006 / THM-4744, or simply from Camion (a strong n-tournament has a Hamiltonian cycle, so H >= n). THM-1370 (3)
  gives only monotonicity, which does not give unboundedness.
- Say "for odd h >= 3". The order-1 tournament is trivially strong with H = 1, while klein's W(1) is empty.

**m5. Section 2 precision.**
- "at least 256/512" should be **299/512** (order 9, H = 615).
- "The minimum is at a first creator with H = 2,901" should read "at a converse pair of first creators (|Aut| = 1)
  with H = 2,901".
- The median 251,107 is the upper of the two middle values (251,103 and 251,107).

**m6. Section 5, Moon's-theorem row.**
- "HYP-9512 (A) says the first creators have Moon chains inside the minimal set ... (class level through order 10,
  value level at order 11)" misattributes. HYP-9512 (A) is the value-level conjecture. The class-level data are
  THM-4922 (3) plus this packet's 152-check.
- A value-level statement gives no chains of tournaments. Suggested text:

  > Every strong tournament is reached from C3 by strong one-vertex extensions (Moon chains). Through order 11 every
  > first creator other than the two creators of 63 has a Moon chain inside the minimal set: THM-4922 (3) and the
  > 152-check through order 10, and the audit's class-level test at order 11. HYP-9512 (A) is the value-level
  > conjecture.

**m7. Notation clashes.**
- f_2(n) is already HYP-9472's and klein's two-cycle floor: f_2(10) = 205, f_2(15) = f(15) + 8 x 125. The packet's
  f_2(10) is 727. Rename the connectivity floors, for example f^lambda_k(n), or add an explicit disambiguation.
- M_ab collides with the transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S + a) B_b(R + b) of THM-030 and THM-050. Those
  theorems are built from the same split products whose unsigned sum is W[a][b].

**m8. Credits to add to the PROOF table** (beyond M3):
- **THM-030 and THM-050:** the split products. When a -> b, W[a][b] counts the Hamiltonian paths through a -> b.
- **THM-4889:** klein's census is canon THM-4889; cite the ID.
- **THM-4744:** Moon's interval chains attain Busch. It is in THM-4942's `related`, but absent from the PROOF credits.
- **THM-4869, HYP-9459 and HYP-9472:** connectivity floors for Hamiltonian cycles. These include the minimum hc of
  2-strong tournaments, 2, 4, 8, 13, 19 and 30 at n = 5-10. They are the closest in-repo antecedent of the "graded
  floors".
- **Grinberg-Stanley** (arXiv:2412.10572), in the OCF row.
- **Thomassen 1980 and Busch 2006, Corollary 2:** exponential lower bounds on Hamiltonian cycles of 2-connected
  tournaments. They give literature context for f_k and are vacuous below about 32 vertices.

**m9. "Each operation has an explicit H-calculus"** (PROOF 3.1, README 3).
- The cited graft calculus of Sessions 225-229 is a rooted operation with the state (D, E, U, d), built on ordered
  two-path covers. It covers one back arc. For A => B with the single arc a -> b reversed:

      H = H(A) H(B) - e_A(a) s_B(b) + U_A(a) U'_B(b),

  where e and s count paths ending and starting at the vertex, U_A(a) counts ordered two-path covers of A whose
  second path starts at a, and U'_B(b) counts ordered two-path covers of B whose first path ends at b. The first path
  P_1 of A and the second path Q_2 of B may be empty. The audit checked this formula by brute force on 300 random
  cases (`aud_backarc.py`).

  **Erratum to this audit.** An early draft of m9 gave H = H(A) H(B) + U_A(a) U'_B(b). That omits the subtracted
  term for the paths destroyed at the reversed arc; the case A = {a}, B = {b} gives 2 instead of 1. The coordinator
  caught the omission independently (`graft_check.py`, 300 cases), and the identity above is the corrected one.
- With two back arcs, a Hamiltonian path can use both, which needs ordered three-part covers on each side. Either
  state that formula or say "one back arc: Sessions 225-229; two: the analogous three-cover transfer".

**m10. "So every tournament through order 11 is generated from the first creators"** (PROOF 3.1, README 3).
- This overreaches. SUM, SUB and the re-connection use no first creator, and K1 belongs to no category.
- Suggested text:

  > Every class of order 3-11 is a first creator, a one-vertex extension of a first creator, or is obtained from
  > strictly smaller tournaments by an ordinal sum, a substitution, or an ordinal sum with at most two arcs reversed;
  > the highly connected (lambda >= 3) classes are first creators or, at order 11 only, seven substitutions.

**m11. "C3[R5, R5, 1] ... built literally from first creators"** (PROOF 3.3; README 4: "all first creators"). The
singleton is not a first creator, and R5 is one of the two first creators of 15 (the other has lambda = 1). Say:

> built from the 3-cycle, the regular first creator R5 of 15 (twice), and a singleton.

**m12. PREREGISTRATION-2 protocol deviation.** The registration says "every order-10 class with lambda >= 2". The run
used min degree >= 2, which adds 17 classes. The executed scope is exhaustive; the registered one is not provably so.
Add one sentence to PROOF section 3.3.

**m13. THM-4942 title.**
- "this fails at order 11" now reads as if the five-operation description fails. Write instead: "the floor inequality
  fails at order 11 (f_3(11) = 11,475 in S(10)), but the description survives there: the seven non-first-creator
  classes with lambda >= 3 are substitutions".
- The title's Busch clause needs M1's scope.

**m14. Section 3.2 Corollary, second sentence.** Repeat the hypothesis: "Under the same hypothesis, every
non-first-creator of order n becomes an ordinal sum after at most two arc reversals." This half needs neither the
lemma nor n <= 13. Also, "Every non-strong class (lambda = 0) is a non-first-creator" (section 3.3) holds through
order 11 (THM-4922 (1)), not in general.

**m15. Section 2 "Scope" paragraph.** The class-level statement does not need order-11 canonical forms (section 3.6
of this audit). Replace the sentence with the audit's result:

> Every first creator of order 11 has a first-creator deletion: 0 class-level orphans among the 506,912,096
> first-creator extensions of the 1,085,824 non-first-creator classes of order 10 (independent audit, section 3.6).

Then move "class-level genealogy at order 11" in THM-4942 from NOT CLAIMED to FINITE-EXACT, with credit to the audit.
Change the HYP-9532 (B) status to "FINITE-EXACT through order 11".

**m16. Reporting.**
- "9.97 x 10^9 extensions" is correct (9,966,649,344).
- The packet's two 564-classes appear in different slices of `order11_summary.txt`. The summary line lists one. Say
  "two classes (a converse pair)".

**m17. Housekeeping.** The lineage section counts 7,661 values "with n(h) <= 10". That is all of S(10) except 1 and
3, so say "the 7,661 values h >= 5 with n(h) <= 10".

**m18. Related context for SUB and prime classes, and for the exhaustive order-11 run: THM-4169** (2026-08-26). It
imports Schmerl-Trotter (Discrete Math. 1993): at each odd order >= 5 the only critical prime tournaments are
T_(2k+1), U_(2k+1) and W_(2k+1), and there are none at even order. Hence every prime order-11 tournament outside
T_11, U_11, W_11 has a prime order-10 deletion.
- This is the natural tool for the REST / prime discussion of section 3.1, since every REST class is prime.
- It is also a precedent for an exhaustive one-vertex augmentation census at order 10 -> 11 (prime parents, 7,906,161
  prime classes at order 10).
- Cite it next to Session230's pointer in the credit table.

## 6. Priority and literature search

**Repository**, by grep of 01-canon and 05-knowledge; the files read are cited.
- **The explicit quadratic attachment form is KNOWN (M4).** It is HYP-2879 / THM-4097 sections 2-3, with the same
  subset-DP coefficient Q_ab = W[a][b], and THM-4114 section 6, the submodular directed-cut quadratic. THM-4118 has
  the ear-response lattice.
  - **Erratum to this audit.** The audit's first draft said "No earlier explicit statement was found". Its grep
    keywords (attachment, quadratic, response polynomial, x_a x_b, M_ab) missed THM-4097's vocabulary: insertion cut
    x_S, exposed slots, Q_ab, cut field. The coordinator found it through the concurrent Session230.
  - Further antecedents:
    - THM-030 and THM-050: the split products E_a(S + a) B_b(R + b);
    - THM-012b and THM-4094: the insertion decomposition;
    - THM-1900 (B): the neutral patterns;
    - Session189 (`deletion_comparison/PROOF.md`): degree <= 2 two-port responses, obtained by Boolean Mobius
      inversion;
    - THM-4868, and Session189 `cutarc_incoming/S10_NOTE.md`: position-of-v decompositions with two-path covers;
    - THM-4882 and Sessions 197-198: the one-apex grammar with ordered two-path covers.
  - The prime-deletion route of THM-4169 (Schmerl-Trotter) is related context for section 3.1 (m18).
- **Arc-strong connectivity and its lower bounds.**
  - THM-4811 (B) is a direct priority for Lemma 3.2 and the per-lambda counts (M3).
  - THM-4869, HYP-9459 and HYP-9472 study cycle floors under connectivity (m8).
  - THM-4867 and THM-4874 use cut arcs and lambda = 1.
  - No in-repo statement of Hamiltonian-path floors f_k for k >= 2, or of "lambda >= 3 implies first creator", was
    found. Those appear new.
- **Strong fibers as antichains.** No prior statement found. The "antichain" hits concern module trees: Session188 and
  Session189.
- **The Busch bound.** klein (THM-4889, HYP-9481 (c)) defines B(h) and studies where windows end, e(h) in
  {B - 1, B}. "n(h) = B(h)" is new as a statement, but readable off klein's length-1 windows (2.7).

**Literature.**
- Moon 1972 and Busch 2006 (Electron. J. Combin. 13, #N3) give f_1 and are credited.
- For k >= 2, a web search found no exact or asymptotic minimum of the number of Hamiltonian **paths** in k-arc-strong
  tournaments. The known bounds concern Hamiltonian cycles of 2-connected tournaments:
  - Thomassen 1980: alpha^(n/32 - 1), alpha = 6^(1/4);
  - Busch 2006, Corollary 2: beta^(n/32 - 1), beta = 5^(1/3).

  Both are vacuous at the orders here, and repo files (THM-4869's audit) already read them.
- "lambda = delta^0 when n <= 4 delta^0 + 1" is an elementary semi-degree count. Its in-repo statement is THM-4811
  (B). No literature citation is claimed here.
- On this evidence, f_2-f_4 at n <= 10, f_3(11), and "lambda >= 3 implies first creator through order 10, failing at
  11" are not known in the literature. A full literature search was not done.

## 7. Files and reproduction

All in this folder.

| File | Role |
| --- | --- |
| `aud_attach.c`, `aud_attach.out` | item 1: the form against direct DP, every pattern; consequences |
| `aud_gen.c`, `aud_gen.out`, `aud_gen.log`, `aud_S1.txt` ... `aud_S10.txt`, `aud_rest_nonpendant.txt` | items 2-3 and the registration's data: orders 1-10, census, floors, degree tables, fertility, shares |
| `aud_lemma.py`, `aud_lemma.out` | item 4: sharpness at n = 14 and random checks |
| `aud_h385.py`, `aud_h385.out` | the order-9 exception |
| `aud_backarc.py`, `aud_backarc.out` | the one-back-arc formula quoted in m9 |
| `aud_busch.py`, `aud_busch.out` | item 6 and the S(n) comparison with klein |
| `aud_o11.c`, `aud_o11.out`, `aud_o11_lambda3.txt`, `aud_o11_hypo.txt` | P2-b and P2-c: the exhaustive order-11 run, with sample brute force |
| `aud_o11_check.c`, `aud_o11_witness.out`, `aud_o11_fert_root63.out`, `aud_o11_lambda3_classes.out`, `aud_lineage_cmp.py`, `aud_lineage_cmp.out` | P2-a, R4 extremal classes, root-63, P2-d, P2-e. The full rebuilt ledger `aud_lineage_le10.txt` is 22 MB, unabbreviated, so it is not kept here; mode `lineage` regenerates it in about 13 minutes |
| `aud_o11_cls.c`, `aud_o11_cls.out`, `aud_o11_cls_control_order8.out`, `aud_o11_cls_control_order10.out` | class-level genealogy at order 11, and its controls |

Reproduce (gcc -O2 -fopenmp; Python 3; about 25 minutes on 12 cores in total; KL is klein's
`../../path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt`):

    gcc -O2 -o aud_attach aud_attach.c && ./aud_attach 20261108
    gcc -O2 -fopenmp -o aud_gen aud_gen.c && ./aud_gen . 12 SCRATCH            # writes SCRATCH/aud_o10_codes.bin
    python aud_busch.py KL; python aud_lemma.py; python aud_h385.py; python aud_backarc.py
    gcc -O2 -fopenmp -o aud_o11 aud_o11.c && ./aud_o11 SCRATCH/aud_o10_codes.bin aud_S9.txt aud_S10.txt KL . 12 2500
    gcc -O2 -o aud_o11_check aud_o11_check.c
    ./aud_o11_check wit ../witnesses_order11.txt aud_S9.txt aud_S10.txt KL
    ./aud_o11_check fert aud_S10.txt 35572154624 18176436867200 4956527526144; ./aud_o11_check root63 aud_S7.txt aud_S8.txt
    ./aud_o11_check lam3 ../lambda3_witnesses_order11.txt aud_o11_lambda3.txt aud_S9.txt aud_S10.txt
    ./aud_o11_check lineage . SCRATCH/aud_o10_codes.bin aud_lineage_le10.txt && python aud_lineage_cmp.py
    gcc -O2 -fopenmp -o aud_o11_cls aud_o11_cls.c && ./aud_o11_cls SCRATCH/aud_o10_codes.bin aud_S9.txt aud_S10.txt aud_o11_cls.out 12
    # controls: aud_o11_check dump 7|9 F.bin; gcc -DN=7|9 aud_o11_cls.c; run with aud_S6/aud_S7 resp. aud_S8/aud_S9

The 78 MB file of order-10 class codes is not committed. `aud_gen` regenerates it in about 2 minutes.
