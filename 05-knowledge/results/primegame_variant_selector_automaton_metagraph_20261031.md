# PRIMEGAME over the first ten primes: Conway's list and Guy's agree through the phase law 2^n 7^(D(n)-1); the owner's pure fractions {1/7, 1/17, 55/1}; FRACTRAN as a selector automaton on base-P# digits; a remark on the n = 5 merged metagraph and the supersingular primes

opus-2026-10-31-S70, 2026-10-31.

**Status.**
- **CITED:** Conway 1987 (Theorem 1; section 5, Figure 2; section 11); Guy 1983; OEIS A007546/A007547; Kari 2012;
  Kari-Kopra 2017; Cloney-Goles-Vichniac 1987; the repository's PRIMEGAME check (procgen_atlas L1), Collatz automaton
  (`collatz_oscillation_gilbreath_20260926.md` section 3) and metagraph canon (THM-4729/4730).
- **PROVED (elementary):** the two-list reduction; the general selector-automaton statement, new to the repository.
- **FINITE-EXACT:** both lists to n = 2500.
- **HYP-9448:** CLOSED as KNOWN.
- **REMARK / NUMEROLOGY:** the metagraph, fives and supersingular readings.
- Independently audited, **SOUND WITH CORRECTIONS**, all applied (section 8; MISTAKE-697).

## 0. Directive

Owner, 2026-10-31: "merge in {15,16,17,23,25} and other groups of 5 we have studied like the fermat primes, the
platonic solids etc. and more consider how 17/91,78/85,19/51,23/38,29/33,77/29,95/23,77/19,1/17,11/13,13/11,15/2,1/7,55/1
relate to the first 10 primes, 2 through 29, and how their prime factorization relates to the tournament merged
metagraph at n=5 where 12 isomorphism classes become 10 nodes, and there are 3 all blue points of symmetry, which relate
to {7,17,55} and {p,p^3,p^2qr} regarding F=U+S, think cellular automata and conenct fractran with prior ideas around the
15 supersingular primes".

**What the objects are.**
- The fraction list is **Conway's PRIMEGAME as printed in his 1987 FRACTRAN paper (Theorem 1)**.
  - Guy's 1983 Math. Mag. article has the earlier form ending 15/14, 15/2, 55/1, the OEIS "fast version" A007546.
    Conway's list is the "slow version" A007547.
  - The repository already checked both on 15 primes (`procgen_atlas_20260924_collatz_implication_atlas.md` row L1).
  - The registers are the first ten primes: data 2, 3, 5, 7 and states 11, 13, 17, 19, 23, 29.
- The "tournament merged metagraph at n = 5" is the merged tiling line metagraph (THM-4729): 10 nodes (12 classes, 2
  converse pairs merged), 8 blue and 24 black lines.
- Its "3 all blue points of symmetry" are the pure-blue nodes **TT5, 1->C3->1 and T5 = Circ(Z_5, {1,2})**, whose
  Hamiltonian paths are all palindromic (THM-4730).

## 1. Results first

1. **Why the two lists agree** (THM-4852 (1)).
   - At each entry into the state-free set (no state prime present), PRIMEGAME is at 2^n 7^(D(n)-1), D(n) = n/spf(n).
     This is Conway's proof of Theorem 1, section 5, Figure 2.
   - Guy's reset (15/14)^(D-1) (15/2)^(n-D+1) (55/1) and Conway's (15/2)^n (1/7)^(D-1) (55/1) both reach
     3^n 5^(n+1) 11.
   - In state-bearing configurations only the shared rules A-K fire.
   - So the lists emit the same primes. The pre-registered checks (HYP-9448) confirm it to n = 2500 (367 primes).
2. **The owner's {7, 17, 55} are Conway's three pure fractions.** Conway's list has exactly three fractions with a
   trivial side: 1/7, 1/17 and 55/1. Guy's has two. Following Conway's semantics:
   - 1/7 erases the divisor counter 7^(D(n)-1);
   - 55/1 starts the next division;
   - 1/17 leaves state 17.

   The product of Conway's fourteen fractions is 55^2/(2 * 17^2). Guy's product is 15/2 times that.
3. **Over the first ten primes.**
   - The lcm of the denominators is 29# (Conway section 11; procgen_atlas .out line 716).
   - Every part is 1, a prime, a semiprime, or 78 = 2*3*13. Under the repository's F = U + S divisor counts
     (`arithmetic_braids_20260917_divisors.md`; balanced shapes p, p^3, p^2qr), the balanced parts are exactly the
     primes. There is no p^3 and no p^2qr.
4. **Cellular automata** (THM-4852 (3)). Every FRACTRAN program in lowest terms over the primes <= P is a **selector
   (quasi-)cellular automaton** on base-P# digits.
   - Each fraction acts digit-locally (division reads the right neighbour, multiplication the left).
   - The fraction is chosen by the units digit.
   - The map does not commute with the shift, so it is not a cellular automaton in the strict sense.
   - Prior art: Kari 2012 (3x/2 in base 6), Kari-Kopra 2017 (p/q in base pq), Cloney-Goles-Vichniac 1987, Conway 1987
     section 11, and the repository's Collatz version (Gilbreath note section 3). Only the general FRACTRAN statement is
     new to the repository.
   - Checked on 20,000 steps of each list, and by the audit on 400 random programs.
5. **The n = 5 merged metagraph** (a REMARK). PRIMEGAME's register-flow graph (24 arcs, 19 pairs) and the metagraph
   (18 pairs, two black loops) cannot be matched by an adjacency-preserving bijection. But the comparison is nearly
   vacuous:
   - for the prime co-occurrence graph there are 70 adjacency-preserving bijections (audit), and none sends the
     owner's {7, 17, 55} to the pure-blue nodes;
   - the matches 10 = 4 + 6 (data + state registers; base-path arcs + tiles) and 3 = 3 (pure fractions; pure-blue
     nodes) are NUMEROLOGY.
6. **Supersingular primes** (CITED facts, NUMEROLOGY readings).
   - The registers are the ten smallest of Ogg's fifteen. The unused five are {31, 41, 47, 59, 71}.
   - A FRACTRAN program over all fifteen is a selector automaton in base 2*3*5*...*71.
   - The first non-supersingular prime PRIMEGAME emits is 37.

## 2. The two lists through Conway's phase law

- Conway's Figure 2 analyses the configuration 5^n 7^d 13. It shows that the orbit reaches 2^n 7^(d-1) when d | n, and
  otherwise moves on to d - 1. So the k-th entry into the state-free set is 2^k 7^(D(k)-1).
- **No state prime present.** Every rule before the last three needs a state prime. A = 17/91 needs 7 * 13; B-K need
  one of 11..29. From 2^a 7^b with b <= a:
  - Guy's list resets by (15/14)^b (15/2)^(a-b) (55/1);
  - Conway's list resets by (15/2)^a (1/7)^b (55/1);
  - both reach 3^a 5^(a+1) 11.
- **State prime present.** Exactly one state prime is present, with exponent 1. The denominators of F-K are exactly 29,
  23, 19, 17, 13, 11, so the first applicable rule lies among A-K, which the two lists share.
- Hence the orbits coincide at every phase boundary, for all n.
- **Data.** The entries begin (1,0,0,0), (2,0,0,0), (3,0,0,0), (4,0,0,1), (5,0,0,0), (6,0,0,2), (7,0,0,0), (8,0,0,3),
  (9,0,0,2), (10,0,0,4).
- **Wording.** The law is about **entries** into the state-free set. Not every state-free configuration has this form:
  840 = 2^3 3 5 7 does not.
- **Computational confirmation.** `primegame_variant_phase_law_20261031.c` runs both lists to n = 2500
  (pre-registered as HYP-9448, which the audit closed as KNOWN). In that file "standard" means Guy's list and "variant"
  means Conway's.

## 3. Anatomy and the selector automaton

| rule | fraction | exponent change over 2,3,5,7,11,13,17,19,23,29 | parts' F-U-S defect (num / den) |
|---|---|---|---|
| A | 17/91 | 7-, 13-, 17+ | 0 / -2 |
| B | 78/85 | 2+, 3+, 5-, 13+, 17- | -3 / -2 |
| C | 19/51 | 3-, 17-, 19+ | 0 / -2 |
| D | 23/38 | 2-, 19-, 23+ | 0 / -2 |
| E | 29/33 | 3-, 11-, 29+ | 0 / -2 |
| F | 77/29 | 7+, 11+, 29- | -2 / 0 |
| G | 95/23 | 5+, 19+, 23- | -2 / 0 |
| H | 77/19 | 7+, 11+, 19- | -2 / 0 |
| I | 1/17 | 17- | - / 0 |
| J | 11/13 | 11+, 13- | 0 / 0 |
| K | 13/11 | 11-, 13+ | 0 / 0 |
| L | 15/2 | 2-, 3+, 5+ | -2 / 0 |
| M | 1/7 | 7- | - / 0 |
| N | 55/1 | 5+, 11+ | -2 / - |

- Conway's list has net change 2: -1, 5: +2, 11: +2, 17: -2, which is the product 55^2/(2*17^2). Guy's product is
  15/2 times that.
- **The selector automaton.** In base B = 29#, for p/q with p, q | B:
  - x/q has digits floor(x_i/q) + (x_(i+1) mod q)(B/q);
  - p y has digits (p y_i mod B) + floor(p y_(i-1)/B), with no further carry;
  - the fraction is the first with q_k | x, read off from x_0 mod q_k.
- Squarefree parts and the full primorial are not essential; lowest terms are (audit). The general statement holds with
  B = P# for any program over the primes <= P.

## 4. Remark: PRIMEGAME and the n = 5 merged metagraph

**Metagraph** (THM-4729). Nodes 0-9 with (H, |Aut|, tilings):

| node | class | (H, \|Aut\|, tilings) |
|---|---|---|
| 0 | TT5 | (1, 1, 1) |
| 1 | (1,1,C3) pair | (3, 3, 2) |
| 2 | (1,S4) pair | (5, 1, 10) |
| 3 | 1->C3->1 | (3, 3, 1) |
| 4 | Q5 | (15, 3, 5) |
| 5 |  | (11, 1, 11) |
| 6 |  | (9, 1, 9) |
| 7 | T*_5 | (9, 1, 9) |
| 8 |  | (13, 1, 13) |
| 9 | T5 | (15, 5, 3) |

- Blue lines: 0-7, 3-8, 4-9, 5-6, 5-7, 5-9, 7-8, 8-9.
- Black lines, each doubled: 1-6, 2-4, 2-5, 2-6, 2-7, 2-8, 4-7, 5-5, 5-8, 6-7, 6-8, 8-8.
- Pure-blue nodes: 0, 3, 9.

**Comparison.** The register-flow graph has 24 arcs, 19 pairs and no loops, against 18 pairs and two loops. So no
bijection carries one onto the other. Under the prime co-occurrence graph the audit finds 70 adjacency-preserving
bijections, none matching the owner's {7, 17, 55} to the pure-blue nodes.

**Count coincidences (NUMEROLOGY):**
- 10 = 4 + 6: the data plus state registers, and the base-path arcs plus tiles of a 5-tournament;
- 3 pure fractions and 3 pure-blue nodes;
- tilings per node 1, 1, 2, 3, 5, 9, 9, 10, 11, 13.

**F = U + S.** The repository's real p^2qr bridge to the order-5 tiling cube is THM-4733's Klein fold. It works on
tilings, and THM-4733 does not define F, U, S. The 3 pure-blue nodes and the 3 balanced shapes match only in count.

## 5. Groups of five (prior work credited)

- The Fermat five (`fermat_five_pepin_collatz_20261025/FIVE_LADDER.md`; THM-4800, THM-4801, THM-4847, THM-4644).
- The history of five and the fifteens (`aliquot_66_five_fifteen_20261028/ATLAS.md`, where Fermat ∩ Ogg = {3, 5, 17}).
- The scope-corrected Platonic reflections.
- The square-sum onsets {15,16,17,23,25} (THM-4851; `monster_levels_square_sum_7_21_20261022.md`).
- The open item "Groups of five".

**NUMEROLOGY additions:**
- PRIMEGAME's registers leave out exactly five of Ogg's fifteen.
- 17 is a Fermat prime, a square-sum path onset, and the state of 1/17.
- The 3x-1 sporadic cycle {17, 25, 37, 41, 55, 61, 91} contains 17, 55 and 91 = 7*13. The audit puts this at baseline
  (probability 0.257 for a random seven-cycle).

## 6. The fifteen supersingular primes

- Ogg's fifteen are 2..31, 41, 47, 59, 71: the primes with X_0(p)^+ of genus zero, and the divisors of |M|.
- The repository's prior work is credited: Ogg/Chowla, the supersingular divisor shapes (F, U, S identities REFUTED),
  the Monster levels, session43's |M|/|B| cycle, THM-4759, THM-4838.
- Arithmetic content here: PRIMEGAME's registers are the initial segment of Ogg's list, and its modulus 29# is their
  product. Any FRACTRAN program over the fifteen is a selector automaton in base 31# * 41 * 47 * 59 * 71.
- No Monster structure is involved.

## 7. Reproduction

```text
gcc -O2 -o pg 04-computation/experiments/primegame_variant_phase_law_20261031.c
./pg 2500                       # primegame_variant_phase_law_20261031.out
python3 04-computation/experiments/primegame_selector_automaton_metagraph_20261031.py   # .out
```

The audit's independent programs are in `04-computation/experiments/primegame_selector_automaton_20261031_audit_S70/`.

## 8. Audit record

Independent audit `primegame_selector_automaton_20261031_audit_S70/AUDIT.md`: **SOUND WITH CORRECTIONS**, all applied.

- **What reproduced.** The phase law, the reset words, the one-state-prime property and the same 367 primes to 2500. A
  negative control at k = 97 failed as it should. OEIS A007546 matches Guy's list and A007547 matches Conway's, 34/34
  terms each. The anatomy, the products, the lcm, the flow graph and the metagraph data also reproduced.
- **MAJOR.**
  - (1) The phase law is Conway 1987, section 5, Figure 2. HYP-9448 is CLOSED as KNOWN, and the equivalence holds for
    all n.
  - (2) The automaton formulation is known (Kari, Kari-Kopra, Cloney-Goles-Vichniac, Conway section 11, the
    repository's Gilbreath note). Only the general FRACTRAN statement is new.
- **MINOR.**
  - Provenance inverted: the owner's list is Conway's 1987 PRIMEGAME.
  - The law holds at entries, not at every state-free configuration.
  - The converse step was missing.
  - "Radius-1 cellular automaton" over-stated the map: it is a selector automaton that does not commute with the shift.
  - The metagraph item was demoted to a remark.
  - The lcm was uncredited, and F = U + S was mis-sourced.
  - The C checker did not compare r2 with the phase index (no effect).
  - Pre-audit status lines and the stub.
