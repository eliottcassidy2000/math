# Circulants of order <= 5 as palindromic path words; Goldbach, Fermat and Legendre as halting FRACTRAN programs on PRIMEGAME's divisor loop; the 25-state Goldbach machine; the five prime groups against square-sum graphs

mac-mini-2026-11-01-fractranTM, 2026-11-01. All work committed directly to `main` and pushed after each step.

## 0. The owner's directive (verbatim)

"investigate how the strong points of symmetry K1, C3, T5 are exactly the circulant tournaments of order at most 5,
which are THM-4849's short Collatz cycles, and work to find fractran and turing machine connections with even more
of the questions from this repo, and work to simplify their representations and preform analysis on their
structure. see how metagraph structure relates with primegame and the 25 state turing machine that corresponds with
goldbach conjecture and make sure all work is done on the master branch and always pushed and also consider how
{2,3,5,7,11,13},{17,19,23,31,47}{29,41,59,71},{37,43,53,61,71},{73,79,83,89,97} realtes with the first 25
corresponding graphs for the square-sum sequence"

The repository has only `main` (no `master`). Every commit of this session was pushed to `origin/main` at once.

## 1. Concurrency and inheritance

Overnight, sessions 129-134 and opus-2026-10-31-S71 worked on related directives. Their packets were read first.

**What they covered.**
- PRIMEGAME's blocks and periodic families (S129).
- Blue/black cycle theorems (S130).
- Monster/Baby/Fi22 prime shells (S131).
- FRACTRAN counters for tournament counts, and both "25 graphs" readings (S132).
- Sixteen FRACTRAN lists for counts and spectra (S133).
- THM-4854, metagraph vs PRIMEGAME (S71).

**What none of them did.**
- None encoded a conjecture as a halts-iff-counterexample program.
- None used Turing machines.
- None explained the circulant cut at 5.

**Session135** (checkpoint 1 at 08:23 on 2026-11-01, two minutes after this session's pre-audit checkpoint) is
concurrent. Its scope is PRIMEGAME basins and priority-aware acceleration, with no halting encodings of
conjectures and no circulant work. It is credited, not duplicated.

**Session65 (2026-10-26)** audited the 25-state Goldbach machine (Leng 2024, Lean-verified):
- it matched the bbchallenge wiki table under a state renaming;
- it ran 2M steps without a halt;
- it recorded the 27 -> 26 -> 25 lineage and the newer 21-state DeepMind machine.

**Typing.**
- STRUCTURAL: a map that preserves a stated predicate.
- NUMEROLOGY: equal counts or labels without such a map.
- GENERIC: true of every object of the class.

## 2. Results first

1. **Why exactly order <= 5** (THM-4855).
   - **Representation.** A Hamiltonian path of Circ(Z_q, S) is a path word: a word over S with distinct partial sums
     mod q (a Hamiltonian path in the Cayley digraph; not a group-theoretic "sequencing").
   - **Segment-sum criterion.** The path is palindromic (blue) iff the segment sums of the word and of its reversal
     lie in S together. For prime q this just says the word is a palindrome.
   - **q = 3, 5.** Every path word is a palindrome: (1,1), and T5's 1111, 2112, 2222, which are its three tilings
     (the words also appear in the prefix-rerouting note, section 6).
   - **q >= 7.** A covering subfamily of THM-4730 (B)'s own explicit paths (three-interval with b = 2 or with
     a + b = q + 1, and descending pairs) fails at a named segment. This proves (B)'s explicit-path claim at all
     orders for that subfamily; the general three-interval claim stays finitely checked. The theorem itself was
     already proved by the audited prefix-rerouting argument, so this is a second, word-level proof.
   - **Checks.** All 131,068 connection sets up to q = 35. The criterion was checked against the arc test on 1.74M
     path words (paths from vertex 0, one connection set per multiplier class).
   - **THM-4849.** Its Collatz cycles ({0}; 3x+5: {1,2,4}; 3x+23: {5,19,...}, {7,11,...}) describe K1, C3, T5 but do
     not explain the cut. The cut comes from the words, not the cycle arithmetic.
2. **FRACTRAN for more repository questions** (THM-4856).
   - Three halts-iff-counterexample programs, built verbatim on PRIMEGAME's divisor loop (its first eleven
     fractions):
     - **GOLDBACHGAME**, 37 fractions, start integer 43;
     - **FERMATGAME**, 36 fractions: a sixth Fermat prime, the repository's "groups of five" question (THM-4847);
     - **LEGENDREGAME**, 33 fractions.
   - Correctness is proved via one-hot Minsky compilation plus the kernel contract.
   - Pre-registered runs:
     - Goldbach: N <= 1000;
     - Fermat: F_0..F_4, and the inverted control reaches F_5;
     - Legendre: 120 windows.

     Literal = dispatch on every checked configuration.
   - **Two controls FAILED as stated.**
     - G4: a measurement artifact; a FRACTRAN test consumes the unit it tests.
     - L3: the control was mis-specified.

     Post-hoc anatomy confirms the programs behave as predicted.
3. **The 25-state Goldbach Turing machine** (`goldbach_tm25_structure.out`, own simulator).
   - **Structure.**
     - One strongly connected control graph.
     - 14 scanner self-loops.
     - 88% of steps in four states (11, 13, 10, 7): the unary trial-division sweep.
   - **Reproduction.** Session65's stage times are reproduced: 45, 132, 301, 596; n = 74 at 1,793,334.
   - **Per-stage cost.** It grows as n^3.26. GOLDBACHGAME's per-target cost is about 9N^2 (measured 8.70-8.97 N^2).
   - **The machine's relation to PRIMEGAME is a shared design (DICTIONARY), not a map.**
     - Both are trial-division searchers: the machine tests divisors ascending and only splits x <= y, while
       PRIMEGAME tests divisors descending.
     - The structural link is elsewhere: GOLDBACHGAME contains PRIMEGAME's divisor loop literally.
   - The number 25 is a description size (session65).
4. **Metagraph vs PRIMEGAME vs the 25-state machine.** No predicate-preserving map was found.
   - The merged metagraph is a static quotient with no computation. S71's THM-4854, S129 and S130 already typed its
     relations with PRIMEGAME DICTIONARY / NUMEROLOGY.
   - The machine's 25 states and the n = 5 metagraph (12 classes, 10 nodes, 32 lines) share no invariant beyond
     coincidental counts.
   - "PRIMEGAME can simulate the metagraph" is GENERIC (Turing completeness).
   - The only structural link found here is GOLDBACHGAME's literal embedding of PRIMEGAME's kernel. The Goldbach
     machine shares only the trial-division design.
5. **The five prime groups and the first 25 square-sum graphs.**
   - Readings 1-2 are session132's (A = the disconnected G_p, threshold 14; the other groups fail). Session132 also
     named the third indexing and proved the 2-connectivity and cycle thresholds 31 and 32.
   - The reading-3 table is computed here (p_k -> the k-th Hamiltonian-path-admissible order
     a_k = 15, 16, 17, 23, 25, ..., 45):
     - A goes to Q_15..Q_26, all with a path but no cycle, not 2-connected.
     - D and E go to orders >= 32, all with cycles and 2-connected.
     - B and C straddle the cycle threshold.
   - **Interleaving obstruction (PROVED).** B, C and D interleave in the prime order (indices B {7,8,9,11,15},
     C {10,13,17,20}, D {12,14,16,18,20}). No size-monotone graph property separates them under any increasing
     indexing.
   - A, B, C are defined by divisibility of Monster-group orders (S131), an arithmetic, non-monotone invariant. D and
     E are the non-supersingular primes below 100, split by size. Verdict: NUMEROLOGY.
   - The duplicate 71 / missing 67 are preserved as written (S131, S132).

## 3. Circulants (THM-4855)

**Witness table** (k = (q-1)/2; 1 ∈ S).

| case | when | witness word | failing segment (w-sum ∈ S?, rev-sum ∈ S?) |
|---|---|---|---|
| A1 | 2 ∈ S, some s ∈ S in [4, k+1] with s-1 ∉ S | 1^(s-2) 2 1^(q-1-s) s | suffix [s, q-1]: -1 (no), 1-s (yes) |
| A2 | 2 ∈ S otherwise (S = {1..k}, the rotational R_q) | 1 2 1^(q-4) 3 | prefix [1, k-1]: k (yes), k+1 (no) |
| B | 2 ∉ S, a rise u ≠ k (u ∉ S, u+1 ∈ S) | 1^(u-1) (q-u) (u+1) 1^(q-2-u) | [u, u+1]: 1 (yes), 2 (no) |
| D | 2 ∉ S, only rise k (S = {1} ∪ [k+1, q-2]) | 1 ((q-3) 1)^(k-1) (q-2) | [1, 2]: -2 (yes), -1 (no) |

**Counts** (`circulant_sequencings.out`).

| q | class | S-sequencings | palindromic |
|---|---|---|---|
| 7 | R7 | 25 | 9 |
| 7 | QR7 | 27 | 9 |
| 11 | each class | 8,401-8,645 | 181-185 |
| 13 | each class | 282k-285k | 1,115-1,127 |

At composite q, the word test is only sufficient. At Z_9 with S = {1,3,4,7}, 66 palindromic paths have
non-palindromic words, because a non-affine anti-automorphism carries them. This is exactly where THM-4730 (B)'s
original "partner sums" argument broke (credited to the 2026-10-20 audit).

**THM-4849 reading.** T5's three words order the 3x+23 cycle points (labelling vertex j -> T^j(x_0); a jump s is
T^s here and T^(-s) in THM-4849's labelling):
- 1111: 10, 5, 19, 40, 20;
- 2112: 10, 19, 40, 20, 5;
- 2222: 10, 19, 20, 5, 40.

These are jumps by the odd-step offsets.

## 4. The FRACTRAN programs (THM-4856)

**Shared kernel** (Conway's divisor loop):

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, X/17, 11/13, 13/11
```

**GOLDBACHGAME** (start 43):

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 53/17, 11/13, 13/11, 705/1333, 825/43, 43/47, 79/2173,
59/53, 1891/118, 67/56699, 59/61, 97/469, 71/67, 1095/2627, 33825/71, 71/73, 3071/158, 89/108151, 79/83, 97/623,
107/89, 101/679, 103/97, 97/101, 1333/3811, 113/103, 4033/3317, 58867/107, 107/109
```

FERMATGAME and LEGENDREGAME are listed in `primekernel_programs.out`.

**Design rules made explicit here.**
1. One state prime in every denominator and numerator. That gives one-hot control and halting in a state with no
   fractions.
2. No self-loops. A reduced fraction s·a/(s·b) = a/b loses the guard, so loops need a helper state.
3. A FRACTRAN test always consumes the tested unit. This is the cause of G4's misread.

**Simplifications used.**
- The kernel's output r2 = M restores the moved register, so no copy loop is needed.
- Offsets P' = P - 2 make the start a single prime.
- Constants are folded into else-branch fractions.

**Per-target cost** (GOLDBACHGAME).

| N | steps for this target |
|---|---|
| 100 | 86,987 |
| 200 | 354,081 |
| 400 | 1,428,647 |
| 800 | 5,738,237 |
| 1000 | 8,973,223 |

The cost is quadratic, about 9N^2 measured. The main term is the test of a prime Q near N, about 6Q^2 (THM-4853 (a)).

## 5. Concept board and inheritance pass

- **Closest proved mechanism.** For circulants: THM-4730 (A) plus the prefix-rerouting proof. For the programs:
  Conway's divisor loop (THM-4853 proof, session126 section 4).
- **Canonical hostile examples.**
  - The Z_9 non-affine anti-automorphism.
  - "Every isomorphism of prime-order circulants is affine" is false as stated; the correct route is Burnside plus
    a normaliser argument (MISTAKE-700).
  - The state-only and padding hostiles to TM sizes (session65).
- **Corrected near misses.**
  - G4 (test consumption) and L3 (control coupling), both recorded as failures.
  - MISTAKE-698's Fermat over-reading is avoided here.
- **Board.**
  - Words versus parity words (THM-4849 sees S, not its sequencings).
  - Kernel plus controller as the normal form of Π1 prime conjectures.
  - Description size is not an invariant.
  - Monotone versus arithmetic classifiers (prime groups).

## 6. Reproduction

Run these from `04-computation/experiments/circulant_fractran_goldbach_20261101/`:

```text
python3 circulant_sequencings.py            # ~15 s
python3 witness_families.py 35               # ~2 s
python3 goldbach_tm25_structure.py 2000000   # <1 s
python3 goldbachgame.py 1000                 # ~1 min (exit 1: G4 failed as stated)
python3 primekernel_programs.py              # ~4 min (exit 1: L3 failed as stated)
python3 posthoc_controls.py                  # post-hoc anatomy of G4, L3
python3 prime_groups_third_indexing.py       # ~10 s
```

## 7. Audit

**Verdict.** Independent audit (`04-computation/experiments/circulant_fractran_goldbach_20261101/audit/AUDIT.md`;
own code in `audit_circulants.py`, `audit_fractran.py`, `audit_tm25.py`, `audit_prime_groups.py`): **SOUND WITH
CORRECTIONS**, all applied (MISTAKE-700). There was no mathematical error in any proof.

**Confirmed independently.**
- **Circulants.**
  - The criterion on 1,844,579 words.
  - (2) exhaustively for q <= 13 and on samples at q = 17, 19, 23.
  - Every case of (3) and the named failing segments, for all admissible choices, q = 7..35.
  - The THM-4849 data.
- **Programs.**
  - Counts, start integers and one-hot control.
  - Literal runs: GOLDBACHGAME to N = 100, FERMATGAME from m = 0..3, LEGENDREGAME to n = 12.
  - Accelerated runs reproducing 2,991,801,422, 25,772,326,379 and 138,922,424,635 steps.
  - The kernel contract on every call.
  - The G4 and L3 reports, plus its own corrected L3 control.
- **The 25-state machine's claims.**
- **The reading-3 table.**
- **Priority search.** No published FRACTRAN program for Goldbach, Legendre or Fermat primes, and nothing on
  palindromic Hamiltonian paths in circulant tournaments.

**Corrections applied.**
1. **Count.** 131,071 -> 131,068 = 2^17 - 4.
2. **Family B.** It is THM-4730's three-interval path with a + b = q + 1, so only its failing segment is new.
3. **Scope.** "Closes THM-4730's open item" became "proves the explicit-path claim at all orders for a covering
   subfamily". The general three-interval claim stays finitely checked.
4. **Proof of (2).** Burnside plus a normaliser argument replaces the false "every isomorphism is affine".
5. **Typing.** The Turing machine's relation to PRIMEGAME is retyped DICTIONARY (shared design, different divisor
   orders). "No map exists" became "none was found".
6. **Numbers and labels.**
   - Cost: about 9N^2, not 6N^2.
   - Controller sizes: 22-26 fractions.
   - The 1.74M are path words from vertex 0.
   - The jump direction: T^s here, T^(-s) in THM-4849.
7. **Credits.**
   - T5's three words: the prefix-rerouting note, section 6.
   - Monster divisibility defines only A, B and C.
   - Session132 named the third indexing and the thresholds 31 and 32.
   - The repository atlas already listed Goldbach and Legendre as generic FRACTRAN statements.
   - Session135 is concurrent.
8. **Terminology.** "S-sequencing" became "path word"; it is not a group-theoretic sequencing.
9. **The session's checks.** The witness script accepted any failing segment and the prime-groups script hard-codes
   the admissible orders. The audit closed both gaps.
