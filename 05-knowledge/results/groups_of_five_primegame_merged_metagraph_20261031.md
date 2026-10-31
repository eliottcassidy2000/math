# Groups of five, PRIMEGAME over the first ten primes and the n = 5 merged metagraph: the residual lanes (exact step counts of both published lists, the owner-gauge label dictionary, and tested correspondences)

mac-mini-2026-10-31-fivegroups, 2026-10-31.

## 0. The owner's directive (verbatim)

"merge in {15,16,17,23,25}  and other groups of 5 we have studied like the fermat primes, the platonic solids etc. and
more consider how 17/91,78/85,19/51,23/38,29/33,77/29,95/23,77/19,1/17,11/13,13/11,15/2,1/7,55/1 relate to the first
10 primes, 2 through 29, and how their prime factorization relates to the tournament merged metagraph at n=5 where 12
isomorphism classes become 10 nodes, and there are 3 all blue points of symmetry, which relate to {7,17,55} and
{p,p^3,p^2qr} regarding F=U+S, think cellular automata and turing machines"

(opus-2026-10-31-S70 received the same text, ending "think cellular automata and conenct fractran with prior ideas
around the 15 supersingular primes".)

## 1. Status, scope and concurrency

**Three sessions answered this directive on the same day.** This note keeps to what the other two do not cover, and
credits them for everything else.

- **opus-2026-10-31-S70.** Note `primegame_variant_selector_automaton_metagraph_20261031.md`, THM-4852 and HYP-9448,
  audited 12:09 (MISTAKE-697). It covers:
  - the two lists and their equivalence: Guy 1983 = the 15/14 "fast version" A007546, Conway 1987 = the owner's list
    A007547;
  - the phase law, which is Conway 1987 section 5 Figure 2 (HYP-9448 CLOSED as KNOWN);
  - PRIMEGAME as a selector (quasi-)cellular automaton in base 29#;
  - no adjacency bijection between the register-flow graph and the merged metagraph;
  - the count coincidences 4+6 and 3 = 3;
  - the F = U + S defects of the fraction parts;
  - the supersingular (Ogg) primes;
  - a short groups-of-five section.
- **session126.** Packet `session126_primegame_symmetry_20261031/` (checkpoints 5635ffe0f1, 522535a7b6). It covers:
  - the primary-source check: the owner's list is **Conway 1987, Theorem 1**, and the 15/14 list is the variant printed
    by Endrullis, Grabmayer and Hendriks 2009;
  - a divisor-loop proof with instruction words and block lengths;
  - the 7-state, 4-counter machine;
  - an actual radius-one cellular automaton;
  - a six-cell perimeter rule-150 automaton on the n = 5 tiling masks, sending 7, 17, 55 to its fixed points 26, 63, 37;
  - F = U + S balance dynamics;
  - a five-families packet (square-sum onsets, Fermat primes, Platonic solids; the shared 2-group mechanism; an A5
    meeting point).

**Typing.**
- STRUCTURAL = a map that preserves a stated predicate.
- NUMEROLOGY = equal counts or labels, with no such map. This includes coincidences that do not survive the program's
  own symmetries.
- GENERIC = true for every object of the class, so it carries no information here.

**Status of this note.**
- **PROVED:** THM-4853.
- **FINITE-EXACT:**
  - the n = 5 merged metagraph in the owner's gauge;
  - the census of THM-643 for n = 3..6;
  - the pre-registered tests P1-P7.
- **CITED:** the inherited items, credited where they are used.
- **Pre-registration:** `correspondence_tests.py` was committed before its run (commit 3a4164776e).
  - P1, P2, P3, P6 and P7 PASSED, and the "95 and 77" half of P4 PASSED.
  - **P4 (a) and P5 FAILED as stated** (sections 6.4 and 6.3). The failures are reported here, not repaired.
- **No HYP was opened.** Nothing conjectural survived: the step counts are proved, and every proposed correspondence
  is NUMEROLOGY, GENERIC, or an inherited STRUCTURAL item.
- **Audit.** Independently audited: SOUND WITH CORRECTIONS, all applied (section 10, MISTAKE-698).
  - After the audit, the per-phase formula for A007546 turned out to be A. P. Heinz's 2009 OEIS Maple program. It is
    now credited, and THM-4853 proves it.

## 2. Results first

1. **Exact step counts of both published lists (THM-4853).**
   - The 15/14 list (Guy 1983) is OEIS **A007546** ("fast version"). The owner's list (Conway 1987) is **A007547**
     ("slow version"; its orbit is A007542). OEIS itself attributes the programs.
   - Closed forms:
     - s_G(k) = sum_{N=2}^{p_k} [7N + 1 + 2 spf(N) + (N-1-D(N))(6N+2) + 2 sum_{D(N)<d<N} floor(N/d)]. This is A. P.
       Heinz's 2009 formula on A007546.
     - s_O(k) = s_G(k) + sum_{n=2}^{p_k-1} (D(n) - 1).
   - **Verification.**
     - Session: the 34 OEIS data terms of each, and 168 primes below 1000.
     - Audit: all 5000 terms of the A007546 b-file and all 1000 of the A007547 b-file.
   - **Growth.**
     - s_G ~ 2(1-c) p^3, with c = sum_p p^-2 prod_{q<p}(1-1/q) = 0.330098.
     - s_O - s_G ~ (c/2) p^2. This is Eldar's A032742 asymptotic of sum D(n).
   - **Credits.**
     - Heinz 2009 for the A007546 formula.
     - session126 for the per-block lengths of the owner's list (6n+2q+2 and 7n+2q+d+3).
     - S70 and Conway for the reset words.
     - Eldar 2022 for c/2.
   - **New:** the proof of Heinz's formula, the A007547 closed form, the difference law, and the constant 2(1-c).
     These are elementary corollaries once the block costs are known.
2. **The merged metagraph in the owner's own gauge** (`merged_metagraph_n5.py`, replicating the owner's explorer
   `inbox/processed/2026-03-05/new/tournament_n5_v5.html`).
   - **Nodes.** 12 classes merge to 10 nodes: the 2 converse pairs (both pure-black) merge.
   - **Lines.** There are 8 blue lines, all simple. There are 24 black lines, in 12 pairs, because the grid transpose
     moves every black line. Two of the pairs are loops.
   - **Points of symmetry.** #0 = TT5, #3 = 1->C3->1 and #11 = T5.
   - **Label dictionary.** The table in section 5 translates between the owner's labels #0..#11, THM-4729 / S70's node
     numbers and session126's least-mask numbers. My independent graph equals S70's line list exactly under this
     dictionary.
3. **10 primes vs 10 nodes: NUMEROLOGY.**
   - S70's no-bijection result is re-verified: the flow graph has 19 pairs, the metagraph 18 edges.
   - A second natural prime graph, the co-occurrence graph (28 edges), is not isomorphic either. Both failures are
     already decided by the edge counts.
   - Under a one-way reading of "preserves adjacency", the metagraph support does embed as a spanning subgraph of the
     co-occurrence graph (audit C(e)). That is weak: the co-occurrence graph is dense.
4. **{1/7, 1/17, 55/1} vs the 3 points of symmetry: NUMEROLOGY.** The set of trivial-sided fractions:
   - has 3 elements in the owner's list and 2 in the 15/14 list, although the two lists compute the same primes;
   - changes under relabelling of the states;
   - read as tiling masks, 7, 17, 55 lie in the mixed classes #4, #9, #10, not in the points of symmetry. session126
     found the same.
5. **{p, p^3, p^2qr} and F = U + S vs the 3 points of symmetry: NUMEROLOGY.**
   - 7 and 17 are balanced, and 55 = pq has defect -2 (also in S70).
   - Of 16 tournament invariants, only 2*c3 = (0, 2, 10) reproduces F at (TT5, 1->C3->1, T5). Neither S = (0, 1, 7)
     nor U = (0, 1, 3) is reproduced, so the identity F = S + U itself has no counterpart.
   - Calibration: 13 of the 136 exhaustively enumerated targets (0, a, b) are hit (9.6%). Pre-registered **P5
     FAILED** (predicted at least 20%).
6. **The strong points of symmetry are the circulants of order 1, 3, 5 (STRUCTURAL, a corollary of inherited canon).**
   - The strong points are K1, C3 and T5, the centres of the three points of symmetry TT5, 1->C3->1 and T5. Only T5 is
     itself a circulant.
   - These are THM-4849's tournament cycles of period at most 5: the cycle {0}, the 3x+5 cycle {1,2,4}, and T5's two
     3x+23 cycles {5,19,40,20,10} and {7,11,22,28,14}. The last two are already in the cycle101 note.
   - That 3 and 5 are Fermat primes is **NUMEROLOGY**. Order 17 is also rigid (THM-871 / THM-4644), yet none of its
     16 circulant classes is a point of symmetry (audit C(a); THM-4730 (B) excludes every odd order from 7 up).
   - **PRIMEGAME echo (NUMEROLOGY).** PRIMEGAME's 95/23 = (5*19)/23 and its numerator 77 = 7*11 reproduce exactly
     these odd points and the denominator 23 = 2^5 - 3^2.
     - The partial event survives 1/20 of the 17280 label permutations. That is exactly 5%, so **P4 (a) FAILED**
       ("< 5%").
     - The full coincidence survives 1/360.
     - It is a property of Conway's labels, not of the computation.
7. **Cellular automata and Turing machines** (section 7).
   - STRUCTURAL and CITED: FRACTRAN = Minsky machines = generalized Collatz maps (Conway 1972, 1987).
   - The automata are S70's selector (quasi-)automaton (a known formulation per MISTAKE-697: Kari 2012, Kari-Kopra
     2017, Cloney-Goles-Vichniac 1987) and session126's two automata.
   - GENERIC: "PRIMEGAME can simulate the metagraph" holds by Turing completeness and says nothing specific.

## 3. Groups of five: what has been merged before (credited) and what this session adds

The repository has merged these fives many times. The standing verdict, restated by session126's five-families packet
today, is: **equal cardinality, no common carrier, except one shared 2-group mechanism.** The Fermat primes are
exactly the odd primes all of whose circulant tournaments are rigid (THM-4644, with Akin 2018 and THM-871), which
is Gauss-Wantzel constructibility.

| five | strongest status (where) | what links it to this directive |
|---|---|---|
| {15,16,17,23,25}: the first nontrivial square-sum Hamiltonian-path orders (A090461; 25 starts the infinite tail, 26 also works) | FINITE-EXACT prefix; the tail CITED (Gerbicz); THM-4505, THM-4510, THM-4851; `collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`; session126 `five_families/README.md` §1 | Under F = U + S only 17 and 23 are balanced (both prime). 15, 17, 23 occur in PRIMEGAME (15 is the numerator of 15/2; 17 and 23 are state primes); 16 = 2^4 and 25 = 5^2 do not. NUMEROLOGY |
| Fermat primes {3,5,17,257,65537} | five known, completeness OPEN; THM-871, THM-4644, THM-4800, THM-4801, THM-4847, HYP-9410 | Fermat ∩ Ogg = Fermat ∩ {the ten PRIMEGAME primes} = {3,5,17} (`aliquot_66_five_fifteen_20261028/ATLAS.md` §3; S70 §5). 3 and 5 are the orders of the nontrivial strong points of symmetry C3 and T5. This is NUMEROLOGY: Fermat rigidity does not select them, since 17 is a Fermat order and none of its 16 circulant classes is pure-blue (audit C(a)). Uniqueness at 3 and 5 is just the count 2^((q-1)/2) = q - 1 of connection sets, true only for q = 3, 5 |
| the five Platonic solids | CITED classification (Schläfli (p-2)(q-2) < 4); HYP-3771 scope-corrected; THM-436 (the 15 overlapping-triangle pairs of a 5-set = the 15 involutions of A5, realised by C5) | session126 `five_families/A5_FUS_MEETING_POINT.md`: \|A5\| = 60 = p^2qr is balanced (also `supersingular_divisor_shapes_qseries_20261022.md` row C10). No metagraph map. NUMEROLOGY except THM-436 |
| the five 3x+101 cycles = the circulant tournaments on Z_7 | PROVED THM-4849 (C) | the q = 7 rung of the cycle ladder in section 6.4; neither Z_7 circulant is a point of symmetry (P3) |
| the five vertex stars of T5 | PROVED THM-4846 (C) | T5 is the third point of symmetry |
| K5 / Petersen / the 10 positive roots of A4 = the 10 proper divisors of p^2qr | THM-261; THM-4733 (Klein fold of the order-five tiling cube onto the divisors of p^2qr, PROVED); "12 classes = tau(p^2qr)" typed an echo in `natural_number_roles_20261004.md` | the only inherited STRUCTURAL bridge between p^2qr and n = 5. It acts on tilings, not on the 10 merged nodes, and makes no F = U + S claim |
| the five mixed classes at n = 5 | THM-643 census | the 10 merged nodes are 3 pure-blue + 5 mixed + 2 pure-black (merged pairs) |

Earlier merges, with their verdicts:
- `gilbreath_fermat_platonic_20260926.md` §4;
- `collatz_doubling_tower_brocard_20260927.md` §4;
- `number_tournament_bridges_20261014.md` §2;
- `fermat_five_pepin_collatz_20261025/FIVE_LADDER.md` §4;
- `s58_fermat_pepin_collatz_20261026/README.md` §5;
- `fermat_landing_points_pepin_20261026.md` §2;
- `cardinality_bridges_20261018m.md` §6;
- `fivefold_boundary_20261018g.md`.

Corrections on record:
- MISTAKE-549 (five mirrors);
- MISTAKE-608 (THM-4644's converse is Akin 2018);
- MISTAKE-668 and MISTAKE-694 (square-sum priority);
- MISTAKE-673 (Fermat lines);
- MISTAKE-696 (tournament-cycle achirality).

There are six regular 4-polytopes, not five.

**Added here, all small.**
- The F = U + S reading of {15,16,17,23,25} (row 1).
- The observation that the strong points of symmetry are exactly the circulants of order at most 5 (section 6.4).
  That 3 and 5 are Fermat primes is coincidental (row 2).
- The PRIMEGAME echo of T5's cycles (section 6.4).

## 4. PRIMEGAME over 2..29 (exact; mostly inherited)

**Factor table.**
- The full factor table is `primegame_decode.out` section 1, identical to S70 §3 and to session126 §2.
- Every numerator and denominator is squarefree. The ten primes are exactly 2..29.
- The lcm of the denominators is 29# = 6469693230 for both lists.

**Machine.** The one-hot lemma is checked in `primegame_decode.out` section 2: every state has an unconditional exit
listed before the first state-free fraction, and only 55/1 creates a state. The machine table of both lists is in
section 3. It agrees with session126 §3:
- data registers 2, 3, 5, 7;
- control states 11, 13, 17, 19, 23, 29 and "none";
- state digraph arcs (1,1), (1,11), (11,13), (11,29), (13,11), (13,17), (17,1), (17,13), (17,19), (19,11), (19,23),
  (23,19), (29,11).

**The two lists differ only in the state-free rows.**
- 15/14 list: `r2>0 and r7>0: 2->3,5`; `r2>0: 2->3,5`; `55/1`.
- Owner's list: `r2>0: 2->3,5`; `r7>0: erase 7`; `55/1`.

**Run from 2.**
- `primegame_decode.out` section 4 gives the trace with exponent vectors through 2^3, and the emissions.
- 15/14 list: (2,19), (3,69), (5,280), (7,707), (11,2363), (13,3876), (17,8068), (19,11319), (23,19201), (29,36866).
- Owner's list: (2,19), (3,69), (5,281), (7,710), (11,2375), (13,3893), (17,8102), (19,11361), (23,19268), (29,36981).
- The owner's list agrees with session126's replay.

**Residue-class form.** `primegame_decode.out` section 5 gives the exact selection densities on residues mod 29#.
- The first eleven fractions have identical densities in both lists, from 1/91 for 17/91 to 129024/2800733 for 13/11.
- The state-free three are:
  - (138240, 829440, 967680)/2800733 for 15/14, 15/2, 55/1;
  - (967680, 138240, 829440)/2800733 for 15/2, 1/7, 55/1.

  This is the same multiset, permuted. "No earlier fraction applies" means exactly that n is coprime to the six state
  primes, because 1/17, 11/13, 13/11, 77/19, 95/23 and 77/29 have bare state denominators. So that event has density
  q = prod_{11<=p<=29} (1 - 1/p) = 21288960/30808063 and is independent of 2 | n and 7 | n. The three densities are
  therefore q/14, 3q/7 and q/2 in both lists.

**List and label dependence of {1/7, 1/17, 55/1}** (section 6). Prime renaming conjugates the dynamics (session126
§2), so any statement about the numbers 7, 17, 55 is about labels.

## 5. The n = 5 merged metagraph in the owner's gauge (FINITE-EXACT)

**Gauge** (the owner's explorer, replicated line by line):
- vertices labelled 5, 4, 3, 2, 1, with base path 5 -> 4 -> 3 -> 2 -> 1;
- tiles (5,1), (4,1), (3,1), (5,2), (4,2), (5,3) are bits 0..5. A clear bit orients the tile's arc x -> y, from
  the larger label to the smaller (the transitive direction), and a set bit reverses it;
- the grid transpose swaps bits 1 <-> 3 and 2 <-> 5;
- a class label #ci is the rank of its lexicographically least adjacency string.

**Census.** `merged_metagraph_n5.out` reproduces THM-643's census rows for n = 3, 4, 5, 6, giving 2, 3, 10 and 34 merged
nodes.

| owner #ci | scores | H | \|Aut\| | tilings | type | owner display node | THM-4729 table / S70 §4 node | session126 node (least mask) |
|---|---|---|---|---|---|---|---|---|
| #0 TT5 | 0,1,2,3,4 | 1 | 1 | 1 | **pure-blue** | 0 | 0 | 0 (0) |
| #1 + #7 | 0,1,3,3,3 + 1,1,1,3,4 | 3 | 3 | 1+1 | pure-black pair | 1 | 1 | 4 (4) |
| #2 + #5 | 0,2,2,3,3 + 1,1,2,2,4 | 5 | 1 | 5+5 | pure-black pair | 2 | 2 | 2 (2) |
| #3 1->C3->1 | 0,2,2,2,4 | 3 | 3 | 1 | **pure-blue** | 3 | 3 | 8 (16) |
| #4 | 1,1,2,3,3 | 9 | 1 | 9 (1 blue) | mixed | 4 | 6 | 1 (1) |
| #6 | 1,1,2,3,3 | 9 | 1 | 9 (3 blue) | mixed | 5 | 7 | 5 (5) |
| #8 | 1,2,2,2,3 | 11 | 1 | 11 (3 blue) | mixed | 6 | 5 | 3 (3) |
| #9 | 1,2,2,2,3 | 15 | 3 | 5 (1 blue) | mixed | 7 | 4 | 9 (17) |
| #10 | 1,2,2,2,3 | 13 | 1 | 13 (3 blue) | mixed | 8 | 8 | 6 (10) |
| #11 T5 | 2,2,2,2,2 | 15 | 5 | 3 | **pure-blue** | 9 | 9 | 7 (11) |

**Lines** (owner display nodes; `merged_metagraph_n5.out`).
- **Blue:** 0-5, 3-8, 4-6, 5-6, 5-8, 6-9, 7-9, 8-9.
  - Each is a single line, and none is a loop.
  - The blue graph has exactly one cycle, 5-6-9-8.
- **Black:** 1-4, 2-4, 2-5, 2-6, 2-7, 2-8, 4-5, 4-8, 5-7, 6-6, 6-8, 8-8. Each carries exactly two lines.
- **Check against S70.** Under the dictionary, these lists are exactly S70's §4 lists. This is an independent check:
  different code and a different gauge reading.

**Why black lines come in pairs.** This is THM-4731 (3), the sigma-twin theorem (2026-10-20), credited; the argument
is recalled here and was re-checked for n = 3..6 by the audit.
- The transpose sigma maps a line {t, t+1} to {sigma t, sigma t + 1}, between the same merged nodes.
- sigma fixes a black line only if sigma t = t + 1 (all ones). That is impossible, because sigma fixes the
  floor((n-1)/2) >= 1 anti-diagonal tiles.
- sigma fixes every blue line.
- So for every n >= 3, black multiplicities in the merged metagraph are even.

**The three points of symmetry** (THM-4730, pure-blue: every Hamiltonian path palindromic).
- #0 TT5: mask {0}.
- #3 1->C3->1: mask {16}.
- #11 T5: masks {11, 46, 53}. THM-643's 2026-10-20 correction: three lines to three distinct nodes.

## 6. The tested correspondences

### 6.1 Ten primes vs ten nodes: NUMEROLOGY

**Map tried.** Any bijection between the primes and the nodes that preserves adjacency in some natural graph on each
side.

**Test** (P6, PASSED).
- The register-flow graph has 24 arcs and 19 pairs, with no loops. This is S70's graph; we reproduced it.
- The support of the merged metagraph has 18 edges plus 2 loops.
- The co-occurrence graph has 28 edges, with degree sequence (3,3,5,6,6,6,6,7,7,7). The metagraph's degree sequence is
  (1,1,1,3,3,5,5,5,6,6).
- Neither graph is isomorphic to the metagraph support. The edge counts alone already decide this.
- Under a one-way reading of "preserves adjacency", the metagraph support embeds as a spanning subgraph of the dense
  co-occurrence graph (audit C(e)). That embedding is weak evidence of nothing.

**What matches.**
- The split 10 = 4 + 6 (data + control) against 4 base-path arcs + 6 tiles. This is S70's DICTIONARY, a language, not
  a model.
- The splits of the merged nodes do not match it: they are 3 + 5 + 2 and 8 + 2.

**Missing coordinate.** An adjacency on primes with three leaves and two loops. PRIMEGAME has none: its flow graph has
no loops.

### 6.2 {1/7, 1/17, 55/1} vs {TT5, 1->C3->1, T5}: NUMEROLOGY

**Map tried.** Trivial-sided fraction to pure-blue node. A 2+1 split would be natural:
- by side: 1/7 and 1/17 have trivial numerator, 55/1 trivial denominator;
- by fibre size: TT5 and 1->C3->1 have one tiling each (leaves), T5 has three.

**Hostile tests.**
1. **List dependence (PROVED, `primegame_decode.out` section 6).** The 15/14 list has only two trivial-sided fractions,
   1/17 and 55/1, and it emits the same primes.
2. **Label dependence.** The transposition 13 <-> 17 of states turns 1/17 into 1/13 (session126 §2: renaming
   conjugates the dynamics).
3. **Mask reading.** As tiling masks in the owner's gauge, 7 lies in #4 (black), 17 in #9 (blue) and 55 in #10
   (black). All three are mixed classes; session126 found the same. Mask 17's line does reach T5 (17 -- 46).

**Verdict.** The trivial-sided set is not an invariant of the computation.

**First failed implication.** "Three trivial sides" is false for an equivalent program.

### 6.3 {p, p^3, p^2qr} and F = U + S: NUMEROLOGY

**Definitions** (`supersingular_divisor_shapes_qseries_20261022.md` DB1).
- D(N) = {d : 1 < d < N, d | N}.
- F = |D|, S = #squarefree elements of D, U = #prime elements of D.
- F = U + S holds exactly at N = p, p^3, p^2qr (the 2026-09-25 divisor lane).
- U ⊂ S, so F = U + S is an equality of counts, not a disjoint union.

**{7, 17, 55}.** 7 and 17 are balanced. 55 = pq has F = 2 and U + S = 4. So {7, 17, 55} is not the shape set (S70
also).

**Shapes vs points.**
- Pre-registered test P5 (`correspondence_tests.out`) used 16 invariants on (TT5, 1->C3->1, T5), each allowed a
  multiplier in {1/2, 1, 2, 3} and an offset in {-1, 0, 1}.
  - F = (0, 2, 10) is matched by 2*c3 alone.
  - S = (0, 1, 7), U = (0, 1, 3) and tau = (2, 4, 12) are matched by nothing.
- **Calibration.**
  - Of the 136 targets (0, a, b) with a <= b <= 15, enumerated exhaustively (not sampled), 13 are matched (9.6%).
    **P5 FAILED** (predicted at least 20%).
  - So a single match is less cheap than predicted.
  - With four targets tried, at least one match had chance about 0.33 (post-hoc, `posthoc_p4_p5.out`). That figure is
    only indicative: it applies the (0, a, b) rate to tau = (2, 4, 12), which lies outside that family.
  - The decisive point is unaffected: the additive identity F = S + U has no counterpart, because neither summand is
    matched.
- **Pairing the atoms.** K1, C3, T5 have orders 1, 3, 5, while Omega(p), Omega(p^3), Omega(p^2qr) are 1, 3, 4. They
  already differ at the third member.

**Inherited STRUCTURAL item.** THM-4733's Klein fold acts on tilings and divisors, not on these nodes.

### 6.4 Points of symmetry as short tournament cycles: STRUCTURAL (inherited); the PRIMEGAME echo: NUMEROLOGY

**Proposition** (corollary; pre-registered as P3, PASSED). The strong points of symmetry K1, C3, T5 are exactly the
circulant tournaments of order at most 5.
- THM-4730 (A): a strong pure-blue tournament of odd order is circulant.
- The repaired (B) of THM-4730: no odd circulant of order >= 7 is pure-blue.
- Orders 1, 3, 5 each have one circulant class.
- Checked: all circulants of orders 1, 3, 5 are pure-blue, and none of order 7 is.

**Under THM-4849's dictionary** (connection set read as a parity word), they are the tournament cycles of period q at
most 5:

| q | point | cycles of x -> x/2, (3x + D_q)/2, D_q = 2^q - 3^((q-1)/2) | odd points |
|---|---|---|---|
| 1 | K1 | {0}, D = 1 | none |
| 3 | C3 | {4, 2, 1}, D = 5 | 1 |
| 5 | T5 | {40, 20, 10, 5, 19} and {28, 14, 7, 22, 11}, D = 23 | {5, 19} and {7, 11} |
| 7 | (R7, QR7: not points of symmetry) | the five 3x+101 cycles of THM-4849 (C) | — |

The q = 5 row is already in `cycle101_tournament_collatz_dictionary_20261029.md` ("3x+23: the two T5 cycles {5, 19} and
{7, 11}"). That q = 3 and 5 are Fermat primes is NUMEROLOGY here (section 3, row 2).

**The PRIMEGAME echo** (exploratory observation, then pre-registered test P4).
- PRIMEGAME contains 95/23 = (5*19)/23, and 77 = 7*11 as the numerator of 77/29 and 77/19. These are exactly T5's
  odd-point products, and 23 = D_5.
- Among the 4! * 6! = 17280 relabellings of registers and states:
  - "some fraction is (odd-point product)/23" holds for exactly 864 = 1/20. **P4 (a) FAILED** at the boundary, since
    the prediction was "< 5%".
    - The post-hoc anatomy explains the exact value. Each of the three (register * state)/state fractions hits with
      probability 1/60, and the three events are disjoint.
  - "95 and 77 are both numerators" holds for 0.56%. This part of P4 PASSED.
  - The full observed coincidence (one fraction over 23 with one product, another numerator with the other) holds for
    1/360.
- **Verdict.** The coincidence depends on Conway's choice of labels. Renaming the states, which conjugates the dynamics
  exactly, destroys it. NUMEROLOGY, however striking.
- **Missing coordinate.** A reason for the state 23 to carry 2^5 - 3^2. PRIMEGAME's 23 is a control label, and
  95/23 means "add one to r5 and return to state 19".

### 6.5 Cellular automata and Turing machines

**STRUCTURAL (CITED).**
- FRACTRAN programs are Minsky register machines (here 4 counters, 7 control states). They are also generalized Collatz
  maps, residue-class-wise linear mod lcm(denominators) = 29# (Conway 1972, 1987; `collatz_procgen_20260922_barrier_atlas.md`).
- The automaton readings are S70's THM-4852 (3) and session126's two automata. S70's is a selector
  (quasi-)automaton, a known formulation (MISTAKE-697). session126's are:
  - the radius-one automaton on 385,025 symbols;
  - the six-cell perimeter rule 150 on tiling masks.

**GENERIC.** Any computable relation between PRIMEGAME and the metagraph can be compiled into FRACTRAN. Existence of a
simulation is therefore no evidence of a correspondence.

**The one new exact result of this session** is THM-4853. It concerns PRIMEGAME alone, not the metagraph.
- The step counts are an invariant of each list, so they distinguish the two equivalent programs.
- The difference is the sum of D(n) - 1 that the 1/7 loop spends erasing the divisor counter.

## 7. Inheritance pass and concept board

**Inheritance.**
- **Closest proved mechanism:** Conway's divisor loop (session126 §4).
- **Canonical hostile examples:**
  - the input 14, where the two lists differ in one step (105 against 15; session126);
  - configuration 15 on the seed-2 orbit, which is state-free but not a phase entry (session126 on the wording of
    HYP-9448).
- **Corrected near misses** (MISTAKE-698):
  - The first draft labelled the 15/14 list "Conway 1987". This was corrected before publication, after session126's
    primary-source check.
  - The pre-audit THM-4853 typed the per-phase formula for A007546 as new. It is A. P. Heinz's 2009 OEIS Maple program.
- **Least-used sidecar:** the OEIS entries' program fields, Maple included, not only their comment and formula
  fields.

**Concept board.**
1. Program equivalence vs label invariance: decided every owner correspondence here.
2. Merged metagraph = tilings / (S5 x transpose): the owner gauge, and the evenness of black multiplicities.
3. Short tournament cycles: the strong points of symmetry and D_q. The Fermat orders were tested and are numerology.
4. F = U + S as a count identity with no summand-level counterpart.
5. Step counts as list invariants (THM-4853).

The board came out unchanged by section 6: no correspondence survived to move it.

## 8. Reproduction

Run these from `04-computation/experiments/groups_of_five_primegame_20261031/`:

```text
python3 merged_metagraph_n5.py    > merged_metagraph_n5.out    # ~4 s; writes merged_metagraph_n5.json
python3 primegame_decode.py       > primegame_decode.out       # <1 s
python3 correspondence_tests.py   > correspondence_tests.out   # ~40 s; exit 1 because P4(a), P5 failed as stated
python3 posthoc_p4_p5.py          > posthoc_p4_p5.out          # post-hoc anatomy of the two failures
```

## 9. Next obligations

- **HYP-9448** was closed as KNOWN by S70's audit (MISTAKE-697). The phase-entry form is Conway 1987 section 5. No
  action is left.
- **OEIS.** Heinz's 2009 formula already gives A007546. The A007547 closed form and the difference law
  A007547(k) - A007546(k) = sum_{n=2}^{prime(k)-1} (A032742(n) - 1) could be submitted as formulas for A007547. This
  is the owner's decision.

## 10. Audit

**Verdict.** Independent audit (`04-computation/experiments/groups_of_five_primegame_20261031/audit/AUDIT.md`, own code
in `audit_primegame.py`, `audit_metagraph.py` and `audit_correspondences.py`): **SOUND WITH CORRECTIONS**, all applied
(MISTAKE-698).

**Confirmed independently.**
- **THM-4853.**
  - The proof was checked line by line.
  - The block costs were checked in 3540 cases.
  - The phase entries were checked by literal FRACTRAN through 25 primes.
  - The closed forms equal all 5000 terms of the A007546 b-file and all 1000 of the A007547 b-file.
  - The difference law holds on those 1000 terms.
  - The owner's orbit equals the A007542 b-file.
  - c = 0.33009804.
- **Metagraph.**
  - The full class table and the 32 lines.
  - The S70 / THM-4729 dictionary. It is the only colour-preserving isomorphism.
  - session126's numbering.
- **Correspondence tests.**
  - All circulants of orders 1, 3, 5 are pure-blue; none of order 7 is (the audit also checked order 17).
  - The Terras cycles.
  - The counts 864, 96, 48 and 13/136.
  - The residue densities.
  - F = U + S to 20000.
  - The byte-identity of the pre-registered script.

**Corrections applied.**
1. **ERROR.** "The strong points of symmetry are exactly the Fermat-order circulants plus K1" was false, because
   17 fails. It now reads "the circulants of order at most 5".
2. **Typing.** The Fermat reading of C3 and T5 is retyped NUMEROLOGY, and "explains" is withdrawn.
3. **Wording.** "Points of symmetry are circulants" became "strong points of symmetry".
4. **Credit.** The black-line pairing is credited to THM-4731 (3).
5. **Credit.**
   - The OEIS attribution of both lists is credited to OEIS itself.
   - The asymptotic c/2 is credited to Eldar's A032742 formula. Only 2(1-c) is new.
6. **THM-4853 wording.** It now says "phase entry" (configurations produced by 1/17) and "phase-entry form of
   HYP-9448".
7. **Cross-references** for P4(a) and P5 were corrected.
8. **Non-isomorphisms.** These are decided by edge counts, and the one-way embedding is noted.
9. **Calibration.** It is exhaustive, not random, and the post-hoc 0.33 is caveated.
10. **Wording fixes:** the bit gauge, the shared machine, the candidate divisors, "count-level contact".

**Found by the session after the audit** (MISTAKE-698 (2)). The per-phase formula of THM-4853 (a) for A007546 is
A. P. Heinz's Maple program on that OEIS entry (2009). It is now CITED, and THM-4853 states that it proves Heinz's
formula. Neither the session nor the audit had read the OEIS program fields. opus-2026-10-31-S70's audited THM-4852
cites it, and that citation is how it was noticed.
