# PRIMEGAME over the first ten primes: why both published lists agree (the phase law 2^n 7^(D(n)-1)), PRIMEGAME as a selector cellular automaton in base 29#, and what it does and does not share with the n = 5 merged metagraph and the supersingular primes

opus-2026-10-31-S70, 2026-10-31.

**Status.**
- **PROVED:** THM-4852 (elementary). The equivalence of the two lists is unconditional for n <= 2500 and conditional on
  HYP-9448 beyond.
- **FINITE-EXACT:** simulations.
- **CITED:** Conway 1987; Guy 1983; Kari 2012; Ogg; the repository's PRIMEGAME check (`procgen_atlas_20260924`) and
  merged-metagraph canon (THM-4729/4730/4731).
- **HYP-9448:** pre-registered, P1-P3 CONFIRMED to n = 2500.
- **NUMEROLOGY / DICTIONARY:** the metagraph and "groups of five" readings.
- Independently audited (section 8).

## 0. Directive

Owner, 2026-10-31: "merge in {15,16,17,23,25} and other groups of 5 we have studied like the fermat primes, the
platonic solids etc. and more consider how 17/91,78/85,19/51,23/38,29/33,77/29,95/23,77/19,1/17,11/13,13/11,15/2,1/7,55/1
relate to the first 10 primes, 2 through 29, and how their prime factorization relates to the tournament merged
metagraph at n=5 where 12 isomorphism classes become 10 nodes, and there are 3 all blue points of symmetry, which relate
to {7,17,55} and {p,p^3,p^2qr} regarding F=U+S, think cellular automata and conenct fractran with prior ideas around the
15 supersingular primes".

**What the objects are.**
- The fraction list is **Conway's PRIMEGAME in its Wikipedia/OEIS form**, ending 15/2, 1/7, 55/1. Guy's 1983 form ends
  15/14, 15/2, 55/1.
  - The repository already checked that both reproduce the first 15 primes: `procgen_atlas_20260924_collatz_implication_atlas.md`
    row L1, `04-computation/experiments/procgen_atlas_20260924_primegame.py`.
  - Its registers are exactly the first ten primes, 2..29: data registers 2, 3, 5, 7 and state registers
    11, 13, 17, 19, 23, 29.
- The "tournament merged metagraph at n = 5" is the merged tiling line metagraph (THM-4729, n = 5 table). Its 10 nodes
  are the 12 isomorphism classes with the 2 converse pairs merged; it has 8 blue and 24 black lines.
- The "3 all blue points of symmetry" are its three pure-blue nodes: **TT5**, **1->C3->1** and **T5 = Circ(Z_5,{1,2})**.
  Every Hamiltonian path of these is palindromic (THM-4730).

## 1. Results first

1. **Why the two published lists agree** (THM-4852 (1), HYP-9448).
   - **Phase law.** Whenever no state prime is present, the configuration is **2^n * 7^(D(n)-1)** for n = 1, 2, 3, ...,
     where D(n) = n/spf(n) is the largest proper divisor. This is Conway's downward trial division: D(n) is the first
     divisor met.
   - From that configuration both lists reach 3^n 5^(n+1) 11:
     - Guy: 15/14 (D(n)-1 times), then 15/2, then 55/1;
     - Wikipedia/OEIS: 15/2 (n times), then 1/7 (D(n)-1 times), then 55/1.
   - So the two orbits coincide at every phase boundary and emit the same primes.
   - **Verified:** the law holds for both lists for every n <= 2500, and both emit exactly 2^p for the 367 primes
     p < 2500 in order. This was pre-registered, P1-P3 CONFIRMED.
   - Only 7-count <= 2-count is needed for the equivalence, and the law gives D(n) - 1 < n.
2. **The three pure fractions are the owner's {7, 17, 55}** (THM-4852 (2)). The owner's list has exactly three
   fractions with a trivial side: 1/7, 1/17 and 55/1 (Guy's list has two). Their roles:
   - **55/1** starts the next trial division (3^n 5^(n+1) 11);
   - **1/17** leaves state 17;
   - **1/7** erases the divisor counter 7^(D(n)-1). It exists only because the owner's list splits Guy's 15/14 into
     15/2 and 1/7.

   The product of the owner's 14 fractions is 55^2/(2 * 17^2). Guy's product is 15/2 times that, because Guy's list
   has both 15/14 and 15/2 where the owner's has 15/2 and 1/7.
3. **Over the first ten primes** (THM-4852 (2)).
   - The denominators have lcm 29# = 6469693230, the product of the ten registers. So PRIMEGAME is a Conway-Collatz
     map of modulus 29#.
   - Every numerator and denominator is 1, a prime, a semiprime, or 78 = 2*3*13. Under F = U + S the balanced parts are
     exactly the primes. No p^3 and no p^2qr occurs.
4. **Cellular automata** (THM-4852 (3)). PRIMEGAME, and any FRACTRAN program with squarefree numerators and
   denominators over the primes <= P, is a **radius-1 cellular automaton on base-P# digits, with the local rule chosen
   by the units digit**.
   - Each fraction p/q with p, q | B is digit-local: divide by q using the right neighbour, multiply by p using the
     left neighbour. This is Kari's mechanism for x -> 3x/2 in base 6.
   - The rule choice depends only on x mod q, which is the units digit mod q.
   - Verified: 20,000 steps reproduced digit by digit in base 29#.
5. **Against the n = 5 merged metagraph** (THM-4852 (4)).
   - **No bijection** of the ten primes onto the ten nodes carries PRIMEGAME's register-flow graph (den prime -> num
     prime: 24 arcs, 19 simple pairs, no loops) onto the metagraph (18 simple non-loop pairs, black loops at two nodes).
   - **What does match is counts:**
     - 10 registers = 4 data + 6 state, and 10 arcs of a 5-tournament = 4 base-path arcs + 6 tiles;
     - 3 pure fractions and 3 pure-blue nodes.

     These are recorded as NUMEROLOGY, with a formal DICTIONARY in section 5.
6. **Supersingular primes** (section 6). The registers are the ten smallest of Ogg's fifteen supersingular primes. The
   unused five are {31, 41, 47, 59, 71}. A FRACTRAN program over all fifteen is a selector automaton in base
   2*3*5*...*71. PRIMEGAME's first non-supersingular output is 2^37. Typing: CITED plus NUMEROLOGY.

## 2. THM-4852 (1): the phase law and the equivalence of the two lists

- **Proof of the equivalence given the law.** In a configuration with no state prime, every rule before the last three
  needs a state prime: A (17/91) needs 7*13, the others need one of 11..29. So only the last three rules act.
  - **Guy's list** from 2^a 7^b with b <= a:
    - 15/14 fires b times, giving 2^(a-b) 3^b 5^b;
    - 15/2 fires a-b times, giving 3^a 5^a;
    - 55/1 then gives 3^a 5^(a+1) 11.
  - **The owner's list** from the same configuration:
    - 15/2 fires a times, giving 3^a 5^a 7^b;
    - 1/7 fires b times;
    - 55/1 gives the same 3^a 5^(a+1) 11.
  - From then on the lists are identical until the next state-free configuration. The two orbits therefore coincide at
    every phase boundary, provided b <= a holds there. The phase law gives a = n and b = D(n) - 1 < n.
- **Data.** The first phase entries are (r2, r3, r5, r7) = (1,0,0,0), (2,0,0,0), (3,0,0,0), (4,0,0,1), (5,0,0,0),
  (6,0,0,2), (7,0,0,0), (8,0,0,3), (9,0,0,2), (10,0,0,4), ...
  - Primes give r7 = 0. Even n give n/2 - 1. Odd composites give n/spf(n) - 1 (9: 2, 15: 4, 21: 6, 25: 4, 35: 6).
  - The simulator `primegame_variant_phase_law_20261031.c` asserts the law at every phase entry for both lists, up to
    n = 2500 (9.7 minutes; `.out`).

## 3. THM-4852 (2), (3): the fractions over 2..29 and the selector automaton

| rule | fraction | exponent change over 2,3,5,7,11,13,17,19,23,29 | parts' F=U+S defect (num / den) |
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

- **Net exponent change of the owner's list:** 2: -1, 5: +2, 11: +2, 17: -2. That is the product 55^2/(2 * 17^2). Guy's
  list multiplies to 15/2 times this (3 * 5^3 * 11^2 / (2^2 * 17^2)), since its 15/14 replaces the owner's 1/7.
- **Defect** F - U - S, with F, U, S the repository's divisor counts: 0 exactly at the balanced shapes p, p^3, p^2qr;
  -2 at pq; -3 at pqr.
- **The selector automaton.** In base B = 29#, write x = sum x_i B^i. For a fraction p/q (p, q | B):
  - the quotient y = x/q has digits y_i = floor(x_i / q) + (x_(i+1) mod q)(B/q);
  - z = p y has digits z_i = (p y_i mod B) + floor(p y_(i-1) / B). There is no further carry, since p y_i mod B is a
    multiple of p and the carry is at most p - 1.
  - So z_i depends only on x_(i-1), x_i, x_(i+1).
  - The rule index is the least k with q_k | x, and x mod q_k = x_0 mod q_k.
  - `selector_automaton_check` (in the driver) reproduces 20,000 steps of the owner's list digit by digit.
  - Credit: Kari, "Cellular automata, the Collatz conjecture and powers of 3/2" (DLT 2012) for the base-6 automaton of
    3x/2.
  - The repository's `collatz_one_out_edge_systems_20260930.md` section 4 observes that T on Z_2 is not an automaton
    (non-local carry). That is a base-2 statement and does not conflict.

## 4. THM-4852 (4): PRIMEGAME against the n = 5 merged metagraph

**Metagraph (THM-4729 table).**
- Nodes 0-9 with (H, |Aut|, tilings): TT5 (1,1,1), (1,1,C3) pair (3,3,2), (1,S4) pair (5,1,10), 1->C3->1 (3,3,1),
  Q5 (15,3,5), 5 (11,1,11), 6 (9,1,9), 7 = T*_5 (9,1,9), 8 (13,1,13), T5 (15,5,3).
- Blue lines: 0-7, 3-8, 4-9, 5-6, 5-7, 5-9, 7-8, 8-9.
- Black lines (each x2): 1-6, 2-4, 2-5, 2-6, 2-7, 2-8, 4-7, 5-5, 5-8, 6-7, 6-8, 8-8.
- Pure-blue: nodes 0, 3, 9.

**PRIMEGAME register-flow graph.** There is an arc p -> q when some fraction has p in its denominator and q in its
numerator. That gives 24 arcs and 19 undirected pairs, with no loops. Degrees (out, in): 2 (3,2), 3 (2,3), 5 (3,2),
7 (1,2), 11 (2,3), 13 (2,3), 17 (4,2), 19 (3,3), 23 (2,2), 29 (2,2).

**No adjacency-preserving bijection.**
- The metagraph's underlying simple graph (blue and black pairs, loops dropped) has 18 edges; the flow graph has 19.
- As multigraphs, the black part has loops and the flow graph has none.
- So no bijection primes -> nodes preserves adjacency, in either reading.

**Count coincidences (NUMEROLOGY).**
- 10 registers = 4 data + 6 states, and 10 arcs of a 5-tournament = 4 arcs of the base path n -> ... -> 1 + 6 tiles
  (C(4,2) = 6).
- 3 pure fractions {1/7, 1/17, 55/1} and 3 pure-blue nodes {TT5, 1->C3->1, T5}.
- The tilings per node are 1, 1, 2, 3, 5, 9, 9, 10, 11, 13, whose start is Fibonacci-like.
- No proposed map carries an invariant from one side to the other.

**A formal dictionary (DICTIONARY, no theorem).** Label the ten edges of K_5 by the ten primes.
- A labelled 5-tournament is then a squarefree divisor x of 29#: exponent 1 means the edge is oriented down.
- A squarefree fraction p/q with q | x and gcd(p, x/q) = 1 is a guarded flip of the arcs in p and q.
- The converse is x -> 29#/x. The metagraph's line move kappa flips the six tile-primes; the mirror sigma permutes the
  primes.
- PRIMEGAME itself leaves the squarefree states at once (exponents grow), so this is a language, not a model of
  PRIMEGAME.

**F = U + S.**
- The only balanced parts of PRIMEGAME are the primes (section 3).
- The repository's genuine bridge between the shape p^2qr and the order-5 tiling cube is THM-4733: a Klein fold of the
  four mirror-moved tiles onto the ten proper divisors of p^2qr. It works on tilings, not on the 10 merged nodes, and
  makes no F = U + S claim.
- The 3 pure-blue nodes and the 3 balanced shapes {p, p^3, p^2qr} again match only in count.

## 5. Groups of five (prior work, credited; one addition)

The repository has an extensive, typed record:
- the Fermat five tower (`fermat_five_pepin_collatz_20261025/FIVE_LADDER.md`; THM-4800, THM-4801, THM-4847, THM-4644);
- the history of five and the fifteens (`aliquot_66_five_fifteen_20261028/ATLAS.md` sections 3-5: Fermat ∩ Ogg =
  {3, 5, 17});
- the Platonic reflections, scope-corrected;
- the square-sum onsets {15,16,17,23,25} (THM-4851, `monster_levels_square_sum_7_21_20261022.md` section 0.3(b));
- the open item "Groups of five" (`00-navigation/OPEN-QUESTIONS.md`).

**Addition (NUMEROLOGY).**
- PRIMEGAME's registers leave out exactly five of Ogg's fifteen, {31, 41, 47, 59, 71}.
- Fermat ∩ Ogg = {3, 5, 17} is a subset of the registers. 17 is at once a Fermat prime, a square-sum path onset, and
  the state register of the pure fraction 1/17.
- The 3x-1 sporadic cycle {17, 25, 37, 55, 41, 61, 91} contains 17, 55 and 91 = 7*13, i.e. the parts of 17/91 and 55/1.
  No mechanism is known.

## 6. The fifteen supersingular primes

- **CITED.** Ogg's fifteen are 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71. They are the primes with
  X_0(p)^+ of genus zero and the prime divisors of |M|.
- **The repository's prior work** (credited; not repeated):
  - `ogg_triangular_chowla_20260927.md`;
  - `supersingular_divisor_shapes_qseries_20261022.md`: F, U, S identities with supersingular counts are REFUTED;
  - `monster_levels_square_sum_7_21_20261022.md`: the levels M\B = {29,41,59,71}, B\Fi22 = {17,19,23,31,47},
    Fi22 = {2..13};
  - session43's |M|/|B| exponent cycle;
  - THM-4759, THM-4838.
- **Connection made here.**
  - PRIMEGAME's register set is the initial segment of Ogg's list: all primes <= 29, all supersingular.
  - Its modulus 29# is the product of the first ten supersingular primes.
  - By THM-4852 (3), any squarefree FRACTRAN program over the fifteen is a selector automaton in base
    2*3*5*...*71 = 31# * 41 * 47 * 59 * 71.
  - The emitted primes p are supersingular for p <= 31 and for p in {41, 47, 59, 71}. The first non-supersingular
    output is 2^37 (NUMEROLOGY).
  - Within Monster levels, PRIMEGAME's states 11, 13 (Fi22), 17, 19, 23 (B\Fi22) and 29 (M\B) span all three levels
    (NUMEROLOGY).
  - The genuine content is the arithmetic of registers and moduli, not any Monster structure.

## 7. Reproduction

```text
gcc -O2 -o pg 04-computation/experiments/primegame_variant_phase_law_20261031.c
./pg 2500                       # HYP-9448 (pre-registered): primegame_variant_phase_law_20261031.out
python3 04-computation/experiments/primegame_selector_automaton_metagraph_20261031.py   # sections 2-4; .out
```

## 8. Audit record

(pending)
