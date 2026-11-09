# Repunits in all bases, the Carmichael conjecture, first creators, and the objects behind 7 and 21

opus-2026-11-08-S86.
- Status: PROVED (THM-4940 (1)-(5), (7)), FINITE-EXACT complete at all orders (THM-4940 (6) to alpha-sum 80,
  (8) below 1024 and to 8192 with the audit, (9) base 3 to 5398), pre-registered P1-P5 CONFIRMED and P6-P8 REFUTED,
  HYP-9530 OPEN, KNOWN results merged and credited, analogies typed. Audit PASS WITH CORRECTIONS (MISTAKE-739).
- No new H-value. Chan-Pak is OPEN.

The owner's directive asked for five things:
- merge in repunit primes in all bases;
- investigate prior work on the Carmichael conjecture in depth;
- track, for every odd number except 7 and 21, the structure that first created it as the order grows, and see how
  all tournaments can be generated from these first appearances with creative enough operations;
- relate the objects used in the proofs that 7 and 21 never occur to fundamental concepts and constructions;
- use a pasted "self-dual core / chiral pair" ledger as inspiration.

## 0. Outcome

1. **Odd-cycle packing counts in every base (THM-4940, new).** The odd-cycle formula H = sum alpha_k 2^k is a base-2
   numeral whose digits are packing counts. Its natural generalisation is H_b(T) = I(Omega(T), b), the number of
   b-coloured odd-cycle packings.
   - **The 7-hole is base-free.** 1 + 3b is never a value, in every base b >= 2, at every order. This proves the
     "moving gap at 1 + 3x", observed at n <= 6 in March 2026 (HYP-1147), and corrects its stated range: it fails at
     x = 0 and x = 1.
   - **The 21-hole is a carry obstruction.** 1 + 10b is a value exactly in the bases 4 and 8, via the carried
     vectors (6, 1) and (2, 1).
   - **The carry-free layer.** Below b^2, Spec_b equals 1 + bA, where A is the set of odd-cycle counts of tournaments
     with no two disjoint cycles (Session 16's intercyclic spectrum). A has infinitely many holes: 3, 10, 14, 17, 20,
     21, 27, ...
   - **Exact tables at all orders to alpha-sum 80.** Base 2 and base 3 miss only the alpha-sums 3 and 10 (the values
     7, 21 and 10, 31). Every other base misses more.
   - So 7 and 21 are the two holes of the carry-free layer that base-2 carries provably cannot fill. Every other
     hole a < 2^34 is filled in base 2 (all-odd coverage to 2^35), and every hole is filled if Chan-Pak holds.
   - Base 3 has the same two unfillable holes. It has no other gap up to alpha-sum 80 (pre-registered P1,
     confirmed), and none up to alpha-sum 1799 by explicit constructions (THM-4940 (9)).
2. **Repunits in all bases.**
   - The repunit 111_b is a base-b packing count iff b + 1 is in A (THM-4940 (5)). It fails for b = 2, 9, 13, 16, 19,
     20, 26, ...
   - The base-2 failure is 7 = 111_2. The failure tied to 21's obstruction is 91 = 111_9; 21 = 111_4 itself is a
     base-4 value.
   - Repunit primes failing in their own base: 7, 421, 757, 2971, 4423, 6481, 8191, ...
   - Earlier repunit work (Session218, Session219, S51, S52, Session27) is merged in section 2.
3. **The Carmichael conjecture.**
   - **Already done.** The literature review is kind-pasteur S14's `CARMICHAEL_MAP.md`; the same request was answered
     on 2026-11-07 by Session218, S14, S15 with Session223, and Session219. Section 3 merges the repository's
     Carmichael dictionary.
   - **New: a carry reading.** The binary "Carmichael exceptions" (H-values with a unique strong class: 3, 5, 11, 13
     below 625) are exactly the carry-free singletons with alpha-sum <= 7. The carry-free layer has 26 singletons
     below 1024, the largest 933 (64 below 8192 by the audit's extension). Each with a >= 8 has a carried binary twin
     (THM-4940 (8)).
   - In base 2, carries both fill holes and create twins:
     - only the gaps 7, 21 are known, and Chan-Pak says they are the only ones;
     - only the singletons 3, 5, 11, 13 are known below 625, and HYP-9477 (A) says they are the only ones.
   - The known exceptional values of both kinds lie in the carry-free range a <= 10.
4. **First creators and generation (section 4).** Merged from Session224, kind-pasteur S16 (THM-4922), mac-mini
   (THM-4930, audit K1 applied, MISTAKE-738) and kind-pasteur S17 (THM-4942, audited), plus the construction lanes Sessions 221-229.
   - The pasted ledger is corrected value by value. For example, 5 is born from one self-converse class, 9 and 15 from
     two self-converse classes each, and 35, 39 and 61 at order 7.
   - New observation: same-H elementary moves (vertex deletion or one arc reversal) do not connect fibers to their
     first creators.
     - In the move graph on orders <= 8, 632 of the 6,880 order-8 classes lie in no component with a first creator
       of their value.
     - The H = 9 fiber has exactly two components at every order (proved by the audit), so this is not an artefact of
       the truncation.
     - Generation needs H-changing operations, as every lane in fact uses.
   - The current best answer to "generate every tournament from first creators" is kind-pasteur S17's THM-4942
     (audited PASS WITH CORRECTIONS, MISTAKE-737). Through order 11, every tournament is one of:
     - a first creator;
     - a first creator plus a vertex;
     - an ordinal sum;
     - a substitution;
     - an ordinal sum re-connected by at most two arcs.
5. **The objects in the 7/21 proofs (section 5).** A typed dictionary: proved identities, analogies, and numerology
   forbidden by MISTAKE-529 and MISTAKE-559. Its one new proved entry is the numeral-system reading of item 1.

## 1. Odd-cycle packing counts in every base (THM-4940)

**Definitions.**
- alpha_k(T) is the number of sets of k pairwise vertex-disjoint directed odd cycles, with alpha_0 = 1.
- H_b(T) = I(Omega(T), b) = sum_k alpha_k b^k.
- Writing H_b = 1 + b m, the alpha-sum is m = sum_{k>=1} alpha_k b^(k-1). The alpha-sum in base 2 is (H - 1)/2, which
  is Session 210's projective-carrier count q.

**Why "carries".**
- In base 2 the alpha-sum is alpha_1 + 2 alpha_2 + 4 alpha_3 + ...
- A tournament with no two disjoint cycles contributes only its first digit alpha_1. Disjoint packings add the
  higher digits, which is a carry.
- The order join multiplies H_b. In base 2 this gives Session 210's join law a, a' -> a + a' + 2aa' on alpha-sums.

**Results** (THM-4940):

| | statement | type |
|---|---|---|
| block lemma | for a < b, 1 + ab is a value iff a is in A | PROVED |
| 7-hole | 1 + 3b is never a value (every b >= 2; every real x > 0 except x = 1) | PROVED, all orders |
| 21-hole | 1 + 10b is a value iff b is 4 or 8 | PROVED (b >= 11 by THM-1370; b <= 10 by complete exhaustion) |
| repunits | 111_b is a value iff b + 1 is in A | PROVED |
| gap tables | exact for alpha-sums <= 80, every base <= 90 | FINITE-EXACT, complete |
| only 2, 3 | no base b >= 4 has gap set {3, 10} | PROVED |
| singletons | 26 carry-free singletons below 1024 (64 below 8192, audit); every one with a >= 8 has a carried binary twin, so binary singletons only a <= 7 | FINITE-EXACT, complete |
| base 3 to 5398 | every value = 1 (mod 3) up to 5398 except 10, 31 is attained (already at order <= 9), so Spec_3 is exact there | FINITE-EXACT, constructive (exploratory; audit follow-up reproduced) |

**The decompositions behind 7 and 21.** The vectors are (alpha_1, alpha_2). Listed are those that also satisfy
alpha_2 <= C(alpha_1, 2) and alpha_3 >= 1 => alpha_2 >= 3. `base_spectra.out` lists all seven chain-compatible
vectors at b = 2.

| base | alpha-sum 3 | alpha-sum 10 |
|---|---|---|
| 2 | (3, 0) impossible | (10, 0), (8, 1), (6, 2), (4, 3): all impossible |
| 3 | (3, 0) impossible | (10, 0), (7, 1), (4, 2): all impossible |
| 4 | (3, 0) impossible | (6, 1) REALISED (a strong 6-tournament, H_4 = 41) |
| 5, 6, 7 | (3, 0) impossible | (10 - b, 1) impossible |
| 8 | (3, 0) impossible | (2, 1) REALISED (C3 => C3, H_8 = 81 = 9^2) |
| >= 9 | (3, 0) impossible | (10, 0) impossible |

**Exact gap tables** (`base_spectra_80.out`). Gap_b is the set of alpha-sums m <= 80 not attained.

| b | gaps | values |
|---|---|---|
| 2 | 3, 10 | 7, 21 |
| 3 | 3, 10 | 10, 31 |
| 4 | 3, 14, 27 | 13, 57, 109 |
| 5 | 3, 10, 27 | 16, 51, 136 |
| 6 | 3, 10, 14, 21, 46 | 19, 61, 85, 127, 277 |
| 7 | 3, 10, 14, 29 | 22, 71, 99, 204 |
| 8 | 3, 35, 52 | 25, 281, 417 |
| 9 | 3, 10, 14, 17, 27, 33, 37, 44, ... | |
| >= 81 | the holes of A: 3, 10, 14, 17, 20, 21, 27, 28, 29, 33, 35, 37, 38, 41, 44, 46, 50, 52, 54, 55, 56, 59, 62, 65, 66, 67, 68, 73, 74, 77, 80 | |

**Rescue sets.** For each hole a of A, these are the bases in which 1 + ab is attained anyway.
- a = 3: none.
- a = 10: {4, 8}.
- Every other hole a <= 80: always contains 2 and 3. For example, 14 -> {2, 3, 5, 8, 12} and
  17 -> {2, 3, 4, 5, 6, 7, 8, 11, 15}.
- So, up to a = 80, the carries rescue every hole of the carry-free layer except 3 and 10, in both bases 2 and 3.

**What is new and what is not.**
- **Not new:**
  - the intercyclic spectrum and its classification (Session 16, audited, KNOWN);
  - the observation of a gap at 1 + 3x for n <= 6 (HYP-1147 = HYP-1109, opus-S71e, March 2026);
  - the impossibility of alpha = (1, 1) (HYP-1015);
  - the evaluation I(Omega, -1) (THM-4754).
- **New:**
  - every all-order statement in the table above;
  - the exact base set {4, 8} for the 21-hole;
  - the block lemma and the repunit criterion;
  - the all-order gap tables and rescue sets;
  - the carry-free singletons;
  - an exhaustive confirmation of Session 16's spectrum below 1024 by two generators (sharing one canonical-form
    routine): the intercyclic generator to order 13, and the few-cycle inventory to alpha_1 <= 80. The audit's
    independent enumerator confirms it below 8192.
- The March row said "for ALL x >= 0". That fails at x = 0 (value 1) and x = 1 (I(Omega(C3 => C3), 1) = 4). The
  correct range is every real x > 0 except x = 1. The historical rows are annotated.

## 2. Repunits in all bases

**Merged prior work.**
- **Session218** (`session218_smith_repunit_totient_20261107/`, audited).
  - Repunit-Smith bridges in every base; R_b(uv) = R_b(u) R_(b^u)(v).
  - Every prime p is R_(p-1)(2), so repunit primes are trivially infinite once the base may vary.
  - The Mersenne-Wagstaff nontotient bridge. It explicitly attaches no tournament map: halving the nontotients 14
    and 26 gives 7 and 13, but 13 is an H-value.
- **Session219.** Among binary Mersenne numbers the only multiplicative relation is M_6 = M_2^2 M_3 = 63. It passes
  through 7, which is not a path count.
- **S51 (THM-4748, THM-4749, HYP-9359).**
  - Binary repunits 2^k - 1 need three reversed arcs when k >= 5 (the Fibonacci fan).
  - Phi_3(2) = 7 divides 2^k - 1 iff 3 | k.
  - 63 and 511 are "blocked": every factorisation uses 7 or 21.
- **S52.** Every 111_q is a strong value below 2100 except 7 and 21. Its projective-plane reading is NUMEROLOGY
  (MISTAKE-661).
- **Session27.** N_63(n) = 2(n - 7).
- **THM-4754.** The base-4 repunits 1, 5, 21, 85 appear as coefficients J_k.
- **Historical rows.** HYP-1043 (Phi_3(2^k) breaks at 73 = 111_8), HYP-1015, and MISTAKE-024/050/653 on 63.

**What this session adds** (THM-4940 (5)).
- The repunit question becomes exact once each base is read as its own packing count.
- **The criterion.** 111_b = 1 + b(b + 1) is a base-b packing count iff b + 1 is the odd-cycle count of an intercyclic
  tournament.
- **Failing bases.** b = 2, 9, 13, 16, 19, 20, 26, 27, 28, 32, 34, 36, 37, 40, 43, ... This set has positive upper
  density by Session 16's bound.
- **Repunit primes.** Among bases with 111_b prime, the failures are 7 (b = 2), 421 (20), 757 (27), 2971 (54),
  4423 (66), 6481 (80), 8191 (90), ...
- **Numerology, typed.**
  - The Goormaghtigh numbers split: 31 = 11111_2 = 111_5 is a base-5 value, while 8191 = 111_90 is not a base-90
    value, because 91 is a hole of A.
  - In the two bases with gap set {3, 10}, the second gap is a repunit: 21 = 111_4 and 31 = 111_5.
  - No mechanism; MISTAKE-559.
- **Binary statement.** Within the verified coverage (2^32 by THM-4875 (D); 2^35 by the HYP-9478 addendum), the only
  *odd* repunits with at least three digits in some base that are not H-values are 7 = 111_2 and 21 = 111_4. Even
  repunits, such as 1111_3 = 40, are never H-values (Redei). This is a restatement of Chan-Pak's two exceptions, not a
  new fact.

## 3. The Carmichael conjecture

**Literature.**
- Read kind-pasteur S14's `../root_change_carmichael_chamber_20261107/CARMICHAEL_MAP.md` section 1. It covers:
  - Carmichael 1907 (flawed proof) and 1922 (no counterexample below 10^37);
  - Sierpinski's conjecture, proved by Ford-Konyagin for even k and by Ford, Ann. of Math. 150 (1999) 283-311, for
    all k;
  - lower bounds: Klee 10^400, Masai-Valette 10^10000, Schlafly-Wagon 10^(10^7), Ford 10^(10^10);
  - Ford's dichotomy: either no counterexample, or a positive proportion of totients are counterexamples;
  - the Carmichael-Klee lemma, which forces 2^2 3^2 7^2 43^2 | x;
  - Pomerance's criterion;
  - Erdos propagation.
- This session spot-checked the citations: Ford's Annals paper and arXiv:1104.3264 (*The distribution of totients*).
  It adds nothing to that review.

**The repository's Carmichael dictionary** (sense: the totient conjecture). It transfers proof shapes, not numbers.

| totient side | tournament side | where |
|---|---|---|
| phi | H on strong classes / the O-Z-J source-word readout | THM-4884; THM-4904; Session200 |
| multiplicity A(m) | strong multiplicity m(h) / profiles pf / word multiplicity M_K(N) | THM-4884; Session200; THM-4904 |
| nontotients | the gaps {7, 21} / root-change misses | THM-1370; HYP-9494 |
| twin phi(x) = phi(2x) | the twins O/J and wZ/wJ | Session200; THM-4904 |
| Carmichael-Klee descent | the dagger-fixed chamber | THM-4904 |
| Ford's dichotomy | density zero for root-change singletons | THM-4904 |
| a counterexample | EXISTS: strong singletons 3, 5, 11, 13; profile singletons 3, 11 | THM-4884, THM-4810; Session200 |

- **Other senses.** Korselt's criterion appears as orbit-counting congruences (THM-4746, THM-4884 (1)), and
  Carmichael numbers are Camion pseudoprimes (THM-4746).
- **No map.** Session218 and Session219 record that the analogy supplies no map. Halving the nontotients 14 and 26
  gives 7 and 13, but 13 is an H-value; and exact totient multiplicity is not a map to tournament-class multiplicity.
- **Open.**
  - HYP-9477 (A)-(D);
  - HYP-9494 (A)-(E), whose (D) is an effective threshold that would prove Chan-Pak;
  - HYP-9506;
  - Session200's two-profile conjecture;
  - HYP-9358.

**What this session adds: the carry reading** (THM-4940 (8), HYP-9530 (C)).
- **The structure.**
  - Write h = 1 + 2a.
  - If a <= 7, every strong tournament with H = h is intercyclic. A strong tournament with two disjoint cycles has
    order >= 6, alpha_1 >= 6 and alpha_2 >= 1, so its binary alpha-sum is >= 8.
  - So the binary strong singletons with h <= 15 are exactly the carry-free singletons a = 1, 2, 5, 6, i.e.
    h = 3, 5, 11, 13.
- **The carry-free layer alone** has 26 singletons below 1024 (a = 1, 2, 5, 6, 8, 15, 32, 45, 89, ..., 933), all
  self-converse. In every base b > 933 they are singletons of the first block.
- **Binary carries.** Every carry-free singleton with 8 <= a < 312 gets a carried twin in base 2 (THM-4884's census).
  That is why the binary list below 625 stops at 13. The audit extends the twin statement to every carry-free
  singleton with 8 <= a < 8192.
- **The parallel.** It is with Carmichael's mechanism, where uniqueness of phi(x) = m is destroyed by twins such as
  phi(x) = phi(2x). Here uniqueness is destroyed by carries, i.e. by disjoint packings. This is ANALOGY: a shared
  shape, with no map between the fibers.
- **The known binary exceptions of both problems live in the carry-free range a <= 10** (all of them, if Chan-Pak
  and HYP-9477 (A) hold):
  - the gaps 7, 21 (a = 3, 10);
  - the singletons 3, 5, 11, 13 (a = 1, 2, 5, 6).

**Numerology, typed.** The Carmichael-Klee chain 2, 3, 7, 43 is Sylvester's sequence. Its terms satisfy
s_(n+1) = 111_(s_n - 1), so 7 = 111_2 and 43 = 111_6. 43 is a base-6 packing count (7 is in A). No mechanism; this is
MISTAKE-559 territory. Sylvester's sequence already appears in several repository lanes, none of which attaches a
tournament mechanism:
- the covering-min and Heegner hypotheses HYP-3724, HYP-3730 and HYP-4079;
- THM-2005;
- S33's number-tournament bridges;
- the s90/s116 scripts, including `sigma_phi_uniqueness_s116n.py` with the chain 2, 3, 7, 43.

## 4. First creators and generation: what the concurrent lanes established

The same directive went to Session224, kind-pasteur-2026-11-07-S16, mac-mini-2026-11-08-firstrealizers and
kind-pasteur-2026-11-08-S17. Session230 (checkpoints 31d57faefb to 5e537fc629) audited their grammar. This session merges their results and runs no
order-10 or order-11 genealogy or generation census.

### 4.1 The pasted S/C ledger, corrected

Birth fibers, written h(s,p): s self-converse classes and p converse pairs at the least order n(h). Source: Session224
`census/birth_fibers.csv`, reproduced by kind-pasteur S16's two generators.

| n | new values h(s,p) |
|---|---|
| 1 | 1(1,0) |
| 3 | 3(1,0) |
| 4 | 5(1,0) |
| 5 | 9(2,0) 11(1,0) 13(1,0) 15(2,0) |
| 6 | 17(2,0) 19(0,1) 23(0,2) 25(0,1) 27(0,1) 29(2,1) 31(0,1) 33(1,2) 37(1,2) 41(1,0) 43(0,1) 45(2,0) |
| 7 | 58 values, e.g. 35(0,1) 39(0,1) 47(0,1) 49(1,3) 61(1,0) 85(1,0) 87(1,0) 125(1,0) 189(1,0) |

Against the pasted text:
- **Order 4.** 5 is born from ONE class, the unique strong 4-tournament, which is self-converse. That is category S,
  not C. The converse pair at order 4 (source => C3 and C3 => sink) repeats 3.
- **Order 5.** 9 and 15 are each born from TWO self-converse classes, not from a converse pair. 11 and 13 are single
  self-converse classes. The dichotomy "one S class or one C pair" first fails at h = 9 (Session224, PROVED).
- **Order 6.**
  - The new values are 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45.
  - 35, 39 and 61 are born at order 7. 61 cannot occur at order 6, since max H(6) = 45.
  - New S values do occur: 41 is a single self-converse class, and 17 and 45 are two self-converse classes each.
  - 25 and 27 are single converse pairs (two classes each, not three). The new values with an odd number of classes
    are 33 and 37 (five classes each) and 41 (one class).
- **Repeats.** "11 fragments into two classes at order 6" is padding. The only order-6 classes with H = 11 are
  source => T11 and T11 => sink, since strong 6-tournaments have H >= 15. They are converses of each other.
- **The law overall.** Through order 10 it holds for 124 of 7,663 values (kind-pasteur S16 audit,
  `aud_fibers_le10.txt`). The stable colour is mac-mini's eternal colour (THM-4930, audit K1 applied, MISTAKE-738): whether some
  self-converse tournament of any order has H = h. 7, 21, 19, 23, 35, ... are eternally black.

### 4.2 Generation from first creators, in the merged record

- **Ordinal sums.** Session224, PROVED: every birth fiber is synthesised from minimal strong inventories by ordinal
  sums, b(q) = min over d | q, d > 1, of [s(d) + b(q/d)].
- **Deletion genealogy.** kind-pasteur S16, THM-4922: PROVED above n(h) max H(n(h) - 2), FINITE-EXACT at orders 4-10.
  Every first creator has a first-creator vertex deletion, except the orphan pair creating 63 (the Fibonacci fan F_6
  and its converse).
- **Extensions and the generation grammar.** kind-pasteur-2026-11-08-S17 (THM-4942 / HYP-9532; checkpoints 14487ef208,
  f870c2abf1; audited PASS WITH CORRECTIONS, MISTAKE-737):
  - The order-11 pre-registered predictions R1-R5 are all CONFIRMED.
  - Each of the 36,248 values new at order 11 is a one-vertex extension of an order-10 first creator, so 63 stays the
    only orphan value.
  - H(T + v) is an explicit quadratic form in the attachment pattern.
  - Through order 11 every tournament is a first creator, a first creator plus a vertex, an ordinal sum, a
    substitution, or an ordinal sum re-connected by at most two arcs.
- **An all-order grammar.** Session230 (`session230_creator_grammar_20261109`, COMPLETE, audited). It uses the cited
  Schmerl-Trotter classification of prime (indecomposable) tournaments.
  - **The grammar.** It is creator-rooted, for all orders: proper substitution, prime-preserving one-vertex growth, and
    three explicit critical-family two-vertex lifts. It is not a minimal basis, and it does not prove that every odd
    count occurs.
  - **Birth minima are not structural primitives.** The birth library is redundant under substitution already at
    H = 5. Yet substitution alone is not enough: the prime W7 has the old count 31.
  - **63.** It is an operation-relative orphan: its two first creators are W5 substitutions.
  - **Growth law.** The alternating-hub family obeys H(W_(n+2)) = 4 H(W_n) - n.
  - The all-order five-operation grammar remains OPEN.
- **Constructions.** Sessions 221-229, audited: grafts, O/Z/J gates, mobile triangles, substitution grammars.
  - Finite menus reach O((log X)^k) values; recursive grammars miss 19 and 29.
  - All-odd coverage is known only to finite bounds: 2^32 by THM-4875 (D), and 2^35 by the HYP-9478 addendum (S82).
- **H-preserving deletions.** Session133's core lemma: they only strip singleton strong components.

**This session's observation** (FINITE-EXACT, orders <= 8; `fibre_moves.c`, `fibre_moves_le8.out`).
- **Setup.** Join two classes of equal H when they differ by one vertex deletion or by one arc reversal.
- **Result.** Then 3 of the 56 order-6 classes, 34 of 456 at order 7 and 632 of 6,880 at order 8 are in no component
  containing a first creator of their value. These are components of the move graph on orders <= 8.
- **Example.** The H = 9 fiber has two components: the two strong order-5 creators with their paddings, and C3 => C3
  with its paddings.
  - This holds at every order (audit, AUDIT.md U3).
  - An H-preserving deletion only strips singleton components.
  - A single arc reversal cannot turn a padding of C3, C3 into a padding of a strong 5-class, because any 6-vertex
    strong component has H >= 15.
- **Consequence.** Same-H elementary moves cannot generate a fiber from its first creators. Generation needs
  H-changing operations: ordinal sums (multiplying H_b in every base), extensions with known increments (S17's
  attachment form), grafts. The base-b reading of section 1 sorts these operations:
  - joins multiply;
  - joins of pieces with cycles create carries;
  - first creators of new values are strong through order 11 (THM-4922, THM-4942; in general HYP-9481 (b), OPEN),
    so there they are produced by extensions, not joins.

## 5. The objects in the 7/21 proofs and the constructions they resemble

Sources: THM-002, THM-029, THM-338, THM-343, THM-1370, THM-4094, THM-1862, THM-517, THM-4754, Session 16,
Session 210, THM-4892. The "type" column follows MISTAKE-529 and MISTAKE-559: an arrow must be exhibited for anything
stronger than ANALOGY.

| object | role in the proofs | construction it resembles | type |
|---|---|---|---|
| OCF, H = sum alpha_k 2^k | the arithmetic of every exclusion | a positional numeral in base 2 with packing digits; disjoint packings are carries | PROVED (THM-4940): the 7-hole is base-free, the 21-hole a base-2/3 carry obstruction |
| chain constraint alpha_k >= 1 => alpha_1 >= k; alpha_2 <= C(alpha_1, 2) | kills (1, 1) and the stray 21 cases | f-vector constraints of the flag complex Ind(Omega) | ANALOGY (no Kruskal-Katona theorem used); HYP-1015 |
| three pairwise intersecting odd cycles force a fourth | alpha = (3, 0) is impossible (three odd cycles, pairwise intersecting). alpha_1 = 3 itself occurs: C3 => (strong 4-tournament) has alpha = (3, 2), H = 15 (THM-029 correction) | intersecting families; the intercyclic class (cross-intersecting inversion families, Bang-Jensen-Bessy) | PROVED (THM-343; Session 16). The Helly version is FALSE (THM-029 repair). |
| strong components, order join | H multiplicative; 21 = 3 * 7 recursion | unique factorisation in a free monoid (the strong classes under ordinal sum). H maps it onto the odd numbers without 7 and 21. That monoid is atomic but not factorial (35 * 77 = 49 * 5 * 11; 63^2 = 3^4 * 49). Its irreducibles are the odd primes other than 7, the products 7p (p prime, p not 3 or 7), and 49, 63, 343 | PROVED: THM-1862; THM-4094 (prime, 7p and carry-atom lanes) |
| Moon pancyclicity and the sharp count c_l >= n - l + 1 | lower bounds on alpha_1 of a component | isoperimetric growth; makes every finite exhaustion complete | PROVED (Moon); used for completeness in THM-4940 |
| Busch's floor f(3q) = 3 * 5^(q-1), f(3q+1) = 5^q, f(3q+2) = 9 * 5^(q-1) | 7 lies in (5, 9), 21 in (15, 25) | gaps between consecutive floors | KNOWN (Busch, EJC 13 (2006) N3, Cor. 1, cited in THM-4744; f(5..8) exhaustive in THM-1370); the gap picture is ANALOGY |
| insertion lemma H(T - v) <= H(T) | strong-deletion induction | heredity | PROVED; holds in every base (THM-4940 (1)) |
| Redei parity | H odd | the unique fixed point of an involution (Burnside, necklace congruences) | PROVED: THM-4892 (C2), Grinberg-Stanley; Korselt form THM-4884 (1) |
| projective carrier, \|Q(T)\| = (H - 1)/2 | 7, 21 correspond to 3, 10 | orbit counting; join law a + a' + 2aa' | PROVED (Session 210). It is the base-2 alpha-sum product of THM-4940. The star law is older: paired_prime_collatz_20261003.md, tiling_modular_atoms_20261003.md, Session142 (forbidden 3, 10) |
| hard-core gas at fugacity 2 | cluster identities: 7 is the unique excluded cluster | statistical mechanics | PROVED formally (THM-517) |
| 7 = \|PG(2,2)\| = 111_2, 21 = \|PG(2,4)\| = 111_4 | none | projective planes, repunits | NUMEROLOGY (MISTAKE-559; breaks at 73 = 111_8). THM-4940 (5) places the alpha-sum-10 repunit at 111_9 = 91, not 111_4. |
| Carmichael-Klee chain 2, 3, 7, 43 | none | Sylvester's sequence, s_(n+1) = 111_(s_n - 1) | NUMEROLOGY |
| Fano / topos / Lawvere | none | fixed-point theorems | REFUTED (MISTAKE-529) |

Concurrent credit:
- kind-pasteur S17's THM-4942 note, section 6, "the 7/21 proof objects are the generation operations" (typed
  analogy), was pushed at 19:49, before this session's pre-registration.
- Session230 (COMPLETE, audited) made one proof object an actual operation. An explicit positive
  odd-cycle/even-path kernel computes every one-vertex attachment as H + 2 * (a directed cut), with a minimum-cut
  birth certificate. A 31-word certificate excludes 7 and 21 for every apex-transitive (one-apex) tournament
  without using the global exclusion. That family lies inside the intercyclic layer, where THM-4940 (2) and
  Session 16 give the same exclusion as 3, 10 not in A.

**Reading.** The table's one new proved row is also the clearest answer to "relate the proof objects to fundamental
constructions":
- The two exclusions are statements about a numeral system with constrained digits.
- 7 fails before any carry can occur. Its only digit string is the forbidden first digit 3.
- 21 fails because every carried decomposition of 10 is blocked in base 2 (and in base 3), although bases 4 and 8
  carry it through.

## 6. Checks and reproduction

Standard library Python 3 and gcc. `CLS` is a directory of kind-pasteur's class files.

1. **The class census.**
   - Build it with `gcc -O2 -o classgen ../first_creators_genealogy_20261107/classgen.c`, then `classgen 9 CLS/cls`.
   - Then `gcc -O2 -o alpha_census alpha_census.c` and, for each n, `alpha_census n < CLS/cls_n.txt > CLS/alpha_n.txt`.
   - Packing vectors for all 198,948 classes of orders <= 9; OCF is checked on every class (0 mismatches).
2. **`few_cycle_strong.c M`.** Every strong class with alpha_1 <= M, any order, with packing vectors (192-bit
   canonical codes).
   - M = 44: 901 classes, orders <= 14 (`few_cycle_strong_44.txt`).
   - M = 80: 5,324 classes, orders <= 18 (`few_cycle_strong_80.txt`), about 3 minutes.
   - Both reproduce the census at orders <= 9.
3. **`base_spectra.py`.**
   - `base_spectra.py CLS 19 40` writes `base_spectra.out`, from the census. It includes the decomposition table.
   - `base_spectra.py few_cycle_strong_44.txt 44 60` writes `base_spectra_44.out`.
   - `base_spectra.py few_cycle_strong_80.txt 80 90` writes `base_spectra_80.out`.
4. **`intercyclic_gen.c 13`** writes `intercyclic_strong_le13.txt`. These are the strong intercyclic classes of orders
   3-13; the counts 1, 1, 6, 11, 20, 37, 70 match the census, and OCF holds as H = 1 + 2 alpha_1 to order 13.
   `carry_free_layer.py` writes `carry_free_layer.out`: A below 1024 equals Session 16's formula; 26 singletons.
5. **`single_apex_spectrum.py 16`** writes `single_apex_spectrum.out`: the one-apex values, checked against brute
   force on all words of length <= 9.
6. **`fibre_moves.c 8 CLS 9`** writes `fibre_moves_le8.out`: same-H move components, with the H = 9 listing.
7. **`coverage_products.py few_cycle_strong_80.txt`** writes `coverage_products.out`, the refuted P6-P8 test.
8. **`ternary_attainment.py CLS few_cycle_strong_80.txt`** writes `ternary_attainment.out`, THM-4940 (9).
9. **`audit/run_audit.sh WORKDIR`** reruns the whole independent audit, including the follow-up, in about 2
   minutes. It reproduced every output byte for byte.

Bugs found and fixed during the session:
- The first `few_cycle_strong` run silently truncated 136-bit codes into 128 bits at order 17. The OCF cross-check
  caught it ("alpha_1 mismatch").
- The base-table loader did not read alpha_6, which occurs once, at order 18; it is now recovered exactly from H.
- The first `ternary_attainment.py` truncated inventory vectors to alpha_1..alpha_4. That created 23 fake piece
  values, all above the reported ranges. The audit follow-up found it; it is fixed and rerun, and the base-3
  uncovered count to 10^4 is 1701.

## 7. Pre-registration

- **`PREREGISTRATION.md`** (P1-P5) was pushed in commit 18463fd3cf before the budget-80 analysis. The audit confirmed
  the push order from the reflog. P1 (the ternary two-gap law) and P2-P5 are all CONFIRMED.
- **`PREREGISTRATION-2.md`** (P6-P8) was pushed in ebe3096762 before `coverage_products.py` existed. It predicted that
  explicit products of the alpha_1 <= 80 inventory and the intercyclic families would cover every alpha-sum in
  (80, 10^5] in bases 2-8.
  - **REFUTED**, the binary control included: 18-31% of alpha-sums are uncovered. The construction universe is too
    small; nothing about gaps follows.
  - The in-sample analogue would have predicted this (MISTAKE-739): pieces with alpha_1 <= 20 already leave 14
    alpha-sums of (20, 80] uncovered.
  - Calibration afterwards with the complete order-<= 9 census as pieces led to the constructive range of
    THM-4940 (9). That is a verification by construction, not a test.

## 8. Next obligations

- **Ternary Chan-Pak (HYP-9530 (A)).**
  - Prove Gap_3 = {3, 10}, or find a third gap.
  - A ternary analogue of THM-4094's lane reduction is the natural route. Spec_3 is a multiplicative monoid under
    order joins, with its own atoms.
- **Finite gap sets in every base (HYP-9530 (B)).** Even one base b >= 4 with a proved finite gap set would be
  informative.
- **The carry principle as a theorem.** Prove that every carry-free singleton a >= 8 acquires a carried twin in base 2
  (HYP-9530 (C), second half; the binary half of HYP-9477 (A)).
- **Generation.** THM-4942's five-operation grammar is exact through order 11, and Session230's Schmerl-Trotter
  grammar holds at all orders but is not minimal. The base-b grading (alpha-sums and carries) is offered as the
  bookkeeping for "creative operations":
  - joins multiply;
  - carried joins add higher digits;
  - extensions create new first digits.
  An open question: which of the five operations change only the first digit (the carry-free layer)?
- **A ternary coverage engine.** A base-3 analogue of THM-4875 (D) is needed: one-arc joins of one-apex pieces
  realise every binary value to 2^32. The calibration shows that the order-<= 9 census runs out near alpha-sum 1800.

## 9. Audit

Independent audit (background agent; `audit/AUDIT.md`, scripts and outputs in `audit/`): **PASS WITH CORRECTIONS**.
- **The mathematics is correct.**
  - Proofs (1)-(5), the real-x extension, and the completeness arguments all check out. Moon (1966) Cor. 2.1 is the
    right citation for the sharp count.
  - Every number was reproduced by an independent enumerator with a different canonical form and cycle counter:
    - the inventory (an identical multiset, with order 19 empty);
    - the gap tables, rescue sets, the base set {4, 8} and the witnesses;
    - A and the singletons, extended to 8192;
    - the fibre moves, and Session224's birth fibers.
- **MAJOR (2), both on wording.**
  - (8)'s title and status went beyond its proof (a < 312). The audit's explicit twins close it to 8192.
  - A section 5 row typed "alpha_1 = 3 impossible" PROVED. Only alpha = (3, 0) is impossible.
- **MINOR (15).**
  - the odd-repunit restriction;
  - open conjectures typed as facts;
  - incomplete decomposition lists;
  - the block-lemma and (7) wording;
  - missing credits (ip_hierarchy.py's ternary pair, the 2026-05 reflections, HYP-886/892/1144);
  - the Session218/219 and Sylvester sentences;
  - the stale S17 merge;
  - the fibre-move scope, the floor and monoid rows, the INDEX wording, and the Moon citation.
- **Upgrades adopted.**
  - completeness from Moon's deletion lemma alone, since the next level is empty;
  - carried twins for every carry-free singleton below 8192;
  - the H = 9 fibre has exactly two components at every order;
  - A below 8192.
- **Follow-up audit of the three post-audit commits: PASS WITH CORRECTIONS (minor).**
  - the P6-P8 ordering is clean;
  - (9) was reproduced independently, and the base-3 range is already attained at order <= 9;
  - a vector-truncation bug in `ternary_attainment.py` was fixed;
  - the PREREGISTRATION-2 wording was corrected.
- **Applied.** Every correction is applied; the errors are recorded as MISTAKE-739.
