# Audit: THM-4855, THM-4856, results note (mac-mini-2026-11-01-fractranTM)

**Verdict: SOUND WITH CORRECTIONS.** This is an independent audit of commit 1224a73dfc, dated 2026-11-01. None of the proofs contains a mathematical error. What needs fixing is one wrong count, one wrong novelty claim, some overclaims, and some typing and numeric slips.

All code is new: Python standard library plus networkx. Run it from the repo root. Each script writes a matching `.out` file.

| Script | Covers | Run time |
|---|---|---|
| `audit_circulants.py` | THM-4855 | 25 s |
| `audit_fractran.py` | THM-4856 | 3 min |
| `audit_tm25.py` | the Turing machine | <1 s |
| `audit_prime_groups.py` | the prime groups | <1 s |

## Findings

1. **ERROR (numeric, not load-bearing).** THM-4855's status and Verification, and note §2.1, say "131,071 connection sets ... at every odd q from 7 to 35". The correct number is 2^17 − 4 = **131,068**, which is also the sum of the per-q lines in `witness_families.out`.

2. **CREDIT: family B is not new.**
   - THM-4855's title and (3), and note §2.1, call family B new.
   - B's path (0..u−1, q−1, u..q−2) is THM-4730's three-interval path P(u, q−1−u), with a + b = q+1.
   - That is exactly the path (B) uses for "one small element a ≥ 3, d = a−1". The audit checked this on all 491,520 (S, u) pairs.
   - What is new is the failing segment [u, u+1].

3. **OVERCLAIM: "closes THM-4730's open item".**
   - The arc-level proof now covers a subfamily of (B)'s paths: b = 2, a + b = q+1, and descending pairs. That subfamily is sufficient.
   - (B1) for general admissible pairs (a, b) is still unproved. The audit found no counterexample among 437,376 such paths at composite q ≤ 27.
   - Reword accordingly.

4. **CREDIT: wrong citation in the proof of (2).**
   - The proof says "every isomorphism between circulants on Z_q is affine (Turner; Babai)". As stated this is false (take S = ∅). Turner and Babai give only *some* affine isomorphism.
   - What the step needs is Aut(T) ≤ AGL(1, q). For tournaments this follows from Burnside's prime-degree theorem (Alspach, JCTB 15, 1973) together with Sylow.
   - Applying it to anti-automorphisms, i.e. maps Circ(S) → Circ(−S), and the step a = −1, are both correct, including in the Paley case.

5. **OVERCLAIM: typing.**
   - Note §2.3 types the 25-state machine's relation to PRIMEGAME as "STRUCTURAL", and THM-4856 calls them the "Same algorithm". No predicate-preserving map is given.
   - The two methods differ. The machine does ascending trial division and only tries splits with x ≤ y (Lean `leap_10_prime`, `leap_17`). PRIMEGAME descends from the top.
   - Retype as DICTIONARY. Note that DICTIONARY is not defined in §1.
   - §2.4's "No predicate-preserving map exists" should read "none was found".

6. **MINOR.**
   - THM-4856 gives "about 6N^2 steps per target (8.97e6 at N = 1000)". The measured cost is 8.70 N² to 8.97 N² for N from 100 to 1000, so ≈9N². The "n^3.26" exponent is a least-squares fit.
   - The title says "20-25-fraction controller", but the controllers have 26, 25 and 22 fractions.

7. **MINOR.**
   - The "1,741,634 Hamiltonian paths" are paths from vertex 0, with one S per multiplier class. The recount is exact.
   - THM-4855 (4)'s "jumps T^s" is T^{−s} under THM-4849 (A)'s labelling T(z_v) = z_{v−1}. State which convention is meant.

8. **CREDIT (minor).**
   - T5's three paths from vertex 0 are already listed in `UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md` §6.
   - "Defined by divisibility of Monster-group orders" holds only for groups A, B and C (S131). D meets π(M) only in 71, and E not at all.
   - Session132 §1 named reading 3, and its §3.3 proved the thresholds at 31 and 32.
   - The repository's atlases already placed Goldbach and Legendre as generic FRACTRAN Π1 statements: rows GOLD and LEG of `procgen_atlas_20260924_graph.py`, and `PROBLEM_CLASSIFICATION.md` §5. THM-4856's contribution is the explicit, verified lists.

9. **MINOR: gaps in the session's own checks, now filled.**
   - `witness_families.py` accepted any failing segment, with one witness per S.
   - `prime_groups_third_indexing.py` hard-codes the admissible orders.

## Independently confirmed

**THM-4855**
- (1) was re-derived. The literal arc test agrees with it on 1,844,579 words.
- (2) holds exhaustively for prime q ≤ 13, and in 9,000 samples at q = 17, 19 and 23.
- The Z_9 hostile has 66 paths.
- (3) holds for every S with q from 7 to 35, in every case:
  - the case split is a partition;
  - the A2 and D derivations are correct;
  - every admissible s and u has its named segment both sum-correct and literally violated.
- (4): the enumerations for q = 3 and q = 5; |Aut T5| = 5; and the THM-4849 cycles.

**THM-4856**
- Program shape:
  - the three published copies of the lists agree;
  - the counts are 37/30, 36/32 and 33/29;
  - the start integers are correct;
  - the state/register split is derived and control is one-hot.
- Literal integer runs, in lockstep with the exponent vectors:
  - GOLDBACHGAME up to N = 100;
  - FERMATGAME from m = 0..3;
  - LEGENDREGAME up to n = 12.
- Exact accelerated runs reproduce every session step count: 2,991,801,422 for N ≤ 1000, 25,772,326,379 from m = 4, and 138,922,424,635 for n ≤ 120.
- The kernel contract (pre- and post-conditions) holds on every call.
- The numbers tested match the number-theoretic prediction.
- No spurious halts occur.
- The G4 and L3 FAILED reports are accurate. An independently built budget control halts at n = 5.

**The 25-state machine**
- One strongly connected component with 14 self-loops.
- Stage times 45, 132, 301, 596, …, 1,793,334.
- The busiest states are 11, 13, 10 and 7 (88.4%), with slope 3.263.

**Prime groups**
- The table is correct; Hamiltonicity was decided, with a DP cross-check for n ≤ 24.
- The interleaving argument is valid but trivial, because it only covers properties monotone in the index.
- The NUMEROLOGY typing is fair.

## Priority (web search)

- No FRACTRAN program for Goldbach, Legendre or Fermat primes turned up ([bbchallenge Fractran](https://wiki.bbchallenge.org/wiki/Fractran), [Cryptid](https://wiki.bbchallenge.org/wiki/Cryptid), [Esolang](https://esolangs.org/wiki/Fractran), [OEIS A395424](https://oeis.org/A395424)).
- Nothing turned up on palindromic Hamiltonian paths in circulant tournaments.
- In the literature, "sequencing" orders distinct elements ([Costa–Della Fiore, EJC 29(3)](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v29i3p33)). An "S-sequencing" is just a Hamiltonian path of Cay(Z_q, S).
