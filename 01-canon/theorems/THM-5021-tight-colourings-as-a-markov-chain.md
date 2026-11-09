---
id: THM-5021
title: "Tight colourings form a finite automaton (a Markov chain after normalisation). Reading a word letter by letter, the greedy packings of its two colours and the running excess E = (P_A + P_B) - ceil(m/3) form an automaton; a tight word never has prefix excess > 1. (1) Its trim block automaton (three letters per step; 5 live states) has transfer matrix with characteristic polynomial x^2 (x - 2)^3 and a single Jordan block of size 3 at the eigenvalue 2: two independent colour choices per block, plus a nilpotent memory of the slack. (2) The tight colourings of length L = 3k, 3k - 1, 3k - 2 (n = 3k + 2, 3k + 1, 3k) have generating functions 2x/(1 - 2x), 2x(1 - x)/(1 - 2x)^2, 2x(1 - x - x^2)/(1 - 2x)^3 (in x^k), so their numbers are 2^k, (k + 1) 2^(k-1) and (k^2 + 9k - 2) 2^(k-3); the pole order at x = 1/2 is one more than the slack 3k - L. (3) At n = 3k exactly (k - 2) 2^(k-3) tight colourings contain 11011 (a 0-junction) and (k - 2) 2^(k-3) contain 00100, and none contains both."
status: >
  PROVED: (1)-(3).
  - The automaton's correctness and the prefix bound are proved below.
  - The matrices, characteristic polynomials, Jordan forms and generating functions are exact rational computations
    (sympy) on the stated finite automata.
  - The automata agree with brute force for L <= 16 (all tight colourings) and k <= 6 (junctions). The audit agrees
    with its own automaton to k = 40, and by brute force for L <= 27.
  - (2) proves the count (k^2 + 9k - 2) 2^(k-3), which THM-4951 (3c) had only fitted.
  ANALOGY (typed, not used in proofs):
  - the eigenvalue 2 as "a chain of independent trials" (the Euler-product reading of independent prime factors);
  - the Jordan block as bounded structure, the slack;
  - the residue chain mod 3 as the last-digit Markov chains of primes (Uppuluri 1968; Lemke Oliver-Soundararajan
    2016; THM-4473).
  The owner's "Marty's Primes as a Markov Chain" was not identified (see the README).
  KNOWN:
  - The automaton is the greedy-state DP of S87 (ternary_metagraph_20261109/staircase_check.py, part (d)), which
    produced THM-4951 (3c)'s fitted counts. New here: the prefix bound that makes it finite, the trim automaton and
    its spectrum, the exact generating functions and the proof of the counts.
  - Transfer matrices for packings on a path: 07-reflections/the-slope-ladder-of-pascal-and-the-growth-constants-kps.md
    (3-separated subsets, Narayana's recurrence).
  - THM-4952: the ribbon packing polynomial is finite-state.
  - THM-4821/4822: Ulam's and Tao's chains are rank one plus nilpotent.
  AUDIT: independent audit PASS WITH CORRECTIONS. Its one MAJOR finding (M1) is applied: the first version stated
  the 18-state matrix's x^14 (x - 2)^4 and a Jordan block of size 4, which are artefacts of six unreachable states
  and one dead state. Minors m4, m7, m8 and m9 are applied. MISTAKE-745 (session-level).
session: opus-2026-11-09-S88
source: 05-knowledge/results/word_law_markov_20261109/README.md
related:
  - THM-5020
  - THM-4951
  - HYP-9540
  - THM-4952
  - THM-4821
  - THM-4822
  - THM-4473
---

# THM-5021 - tight colourings as a finite automaton

## The automaton

Read a word w = w_0 w_1 ... w_(L-1) left to right, keeping:
- s_A, s_B in {1, 2, 3}: the distance from the next position to the last greedily chosen 1 (resp. 0), capped at 3.
  The value 3 means a new element of that colour may be chosen.
- c: the number of chosen elements so far.
- E = c - ceil(m/3), after m letters.

On reading the letter 1 (resp. 0):
- if s_A = 3 (resp. s_B = 3), choose it and set that distance to 1;
- otherwise increase it by 1, capped at 3;
- increase the other colour's distance by 1, capped at 3.

E then changes by (chosen) - [m = 0 (mod 3)].

**Correctness.** The greedy rule is optimal for each colour separately, so after L letters c = P(A) + P(B), and w is
tight iff E = 0. Also E >= 0 on every prefix: this is Session26's bound for the prefix.

**Prefix bound.** P(A) >= P(A ∩ [0, m-1]) + P(A ∩ [m+2, L-1]), and likewise for B, so
E(L) >= E(m) + ceil(m/3) + ceil((L - m - 2)/3) - ceil(L/3) >= E(m) - 1. A tight word therefore has E <= 1 on every
prefix, and pruning E >= 2 leaves the 18 triples (s_A, s_B, E). The phase m mod 3 enters only through the E update.

**The states that matter.**
- Only 12 of the 18 triples ever occur.
- At block boundaries (m = 0 mod 3) only 6 occur: (3,3,0), (3,1,1), (3,2,1), (1,3,1), (2,3,1) and (3,3,1).
- (3,3,1) is dead: from it every block must be monochromatic, and E never returns to 0. An example of a word reaching
  it is 100111.
- The trim block automaton therefore has 5 states.

## Statement

**(1) Spectrum.** The trim block transfer matrix T (three letters per step, 5 x 5) has characteristic polynomial
x^2 (x - 2)^3. Its Jordan form is one block of size 3 at the eigenvalue 2 and one nilpotent block of size 2.
- The untrimmed 18 x 18 matrix built from all triples has x^14 (x - 2)^4 and a block of size 4 at 2. The extra size
  comes from the dead state (3,3,1), so it is an encoding artefact, not an invariant of the tight colourings.
- T counts words; it is a transfer matrix, not a stochastic matrix. It becomes a Markov chain after normalisation:
  uniform letters, T/8 per block, killed at E = 2.

**(2) Counts.** Write N(L) for the number of tight colourings of length L, both first letters counted. Then
- sum_(k>=1) N(3k) x^k = 2x/(1 - 2x), so N(3k) = 2^k;
- sum_(k>=1) N(3k - 1) x^k = 2x(1 - x)/(1 - 2x)^2, so N(3k - 1) = (k + 1) 2^(k-1);
- sum_(k>=1) N(3k - 2) x^k = 2x(1 - x - x^2)/(1 - 2x)^3, so N(3k - 2) = (k^2 + 9k - 2) 2^(k-3).

The order of the pole at x = 1/2 is 1, 2, 3, one more than the slack 3k - L = 0, 1, 2. The maximal order 3 equals the
size of the Jordan block in (1). Up to complement the counts are half of these, as THM-4951 (3c) lists them.

**(3) Junctions at n = 3k.** Extend the automaton by the last four letters and two flags, one for "11011 has
occurred" and one for "00100 has occurred". With the phase tracked, 103 states occur (the audit's count, m7); 317
occur if any phase is allowed at any step, which is the count markov_junctions.py prints.
- The tight colourings of length 3k - 2 that contain 11011 have generating function x^3/(1 - 2x)^2. So there are
  (k - 2) 2^(k-3) of them for k >= 3.
- The same holds for 00100.
- Those containing both have generating function 0. The flags record only whether each factor occurs, so the
  automaton sees "no tight colouring has both kinds"; THM-5020 (1) proves the stronger J <= 1.

## Proof

- **(1), (2).** `markov_tight.py` builds the 18-state matrices and computes exactly in sympy the generating functions
  x v0 (I - xT)^(-1) T fin, x v0 (I - xT)^(-1) M0 M1 fin and x v0 (I - xT)^(-1) M0 fin. Its printed first function,
  -2/(2x - 1), lacks the factor x; it is sum N(3k) x^(k-1).
- `trim_check.py` finds the used, dead and trim states, and the trim matrix's characteristic polynomial and Jordan
  form (trim_check.out). The audit's `aud_markov_states.py` agrees, and the trim automaton gives the same three
  generating functions.
- The closed forms follow by partial fractions. For example, the coefficient of x^k in
  2x(1 - x - x^2)/(1 - 2x)^3 is 2 [C(k+1, 2) 2^(k-1) - C(k, 2) 2^(k-2) - C(k-1, 2) 2^(k-3)] = (k^2 + 9k - 2) 2^(k-3).
- **(3).** `markov_junctions.py` (markov_junctions.out); brute force agrees for k <= 6, and the audit for k <= 30.

## Remarks (ANALOGY, typed)

- **A chain of independent trials.** The eigenvalue 2 is "one free colour per block of three positions". A tight
  colouring is a chain of independent binary trials, the block colours, as in the Euler product, where each prime
  divides a random integer independently.
  - The trials are not unconstrained: the slack 3k - L is a budget of defects (short blocks, junctions).
  - Each unit of slack adds one power of k, up to the size-3 Jordan block: the chain's nilpotent memory.
  - Structure (the rigid aligned blocks) is overlaid with randomness (the free colours); the defects are where the
    structure bends.
  - The numerators 1, 1 - x, 1 - x - x^2 look like Fibonacci-type truncations. This is recorded as numerology only
    (MISTAKE-559).
- **Primes as a Markov chain.** The phase (position mod 3) drives the automaton, as the last digit drives the
  observed Markov chains of consecutive primes (Uppuluri 1968; Lemke Oliver-Soundararajan 2016); THM-4473 has exact
  Markov laws for Collatz digit chains.
  - The Ulam and Tao chains are rank one plus nilpotent (THM-4821/4822). T is not: T = 2 Pi + N, where the
    eigenvalue-2 part is a Jordan block of size 3, so Pi has rank 3.
  - The common feature is a single nonzero eigenvalue. No mechanism is claimed.
- **Symmetry overlaid with asymmetry.** Reversal acts on the automaton's language, and the blue (palindromic) tight
  colourings are its fixed points.
  - At n = 3k + 2 these are the palindromic block words, 2^ceil(k/2).
  - They grow like the square root of all tight words: symmetric structures are rare but never absent.
