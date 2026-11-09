---
id: THM-5021
title: "Tight colourings form a finite Markov chain. Reading a word letter by letter, the greedy packings of its two colours and the running excess E = (P_A + P_B) - ceil(m/3) form an automaton; a tight word never has prefix excess > 1, so 18 states suffice. (1) The block transfer matrix T (three letters) has characteristic polynomial x^14 (x - 2)^4 and a single Jordan block of size 4 at the eigenvalue 2: two independent colour choices per block, plus nilpotent corrections. (2) The tight colourings of length L = 3k, 3k - 1, 3k - 2 (n = 3k + 2, 3k + 1, 3k) have generating functions 2x/(1 - 2x), 2x(1 - x)/(1 - 2x)^2, 2x(1 - x - x^2)/(1 - 2x)^3 (in x^k), so their numbers are 2^k, (k + 1) 2^(k-1) and (k^2 + 9k - 2) 2^(k-3); the pole order at x = 1/2 is one more than the slack 3k - L. (3) At n = 3k exactly (k - 2) 2^(k-3) tight colourings have a 0-junction and (k - 2) 2^(k-3) a 1-junction (generating function x^3/(1 - 2x)^2 each), and none has both."
status: >
  PROVED: (1)-(3).
  - The automaton's correctness and the prefix bound are proved below.
  - The matrices, characteristic polynomial, Jordan form and generating functions are exact rational computations
    (sympy) on the stated finite automata.
  - The automata agree with brute force for L <= 16 (all tight colourings) and k <= 6 (junctions).
  - (2) proves the count (k^2 + 9k - 2) 2^(k-4) at n = 3k, which THM-4951 (3c) had only fitted.
  ANALOGY (typed, not used in proofs):
  - the eigenvalue 2 as "a chain of independent trials" (the Euler-product reading of independent prime factors);
  - the Jordan block as bounded structure, the slack;
  - the residue chain mod 3 as the last-digit Markov chains of primes (Uppuluri 1968; Lemke Oliver-Soundararajan
    2016; THM-4473).
  The owner's "Marty's Primes as a Markov Chain" was not identified (see the README).
  KNOWN: transfer matrices for packings on a path (07-reflections/the-slope-ladder-of-pascal-and-the-growth-constants-kps.md:
  the 3-separated subsets of a path are counted by Narayana's recurrence); THM-4952 (the ribbon packing polynomial
  is finite-state); THM-4821/4822 (Ulam's and Tao's chains are rank one plus nilpotent).
  AUDIT: PENDING.
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

# THM-5021 - tight colourings as a finite Markov chain

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
tight iff E = 0.

**Prefix bound.** P(A) >= P(A ∩ [0, m-1]) + P(A ∩ [m+2, L-1]), and likewise for B, so
E(L) >= E(m) + ceil(m/3) + ceil((L - m - 2)/3) - ceil(L/3) >= E(m) - 1. A tight word therefore has E <= 1 on every
prefix. Pruning E >= 2 leaves 18 states (s_A, s_B, E), and the phase m mod 3 enters only through the E update.

## Statement

**(1) Spectrum.** Let T be the product of the one-step matrices for the phases 0, 1, 2, an 18 x 18 integer matrix.
Its characteristic polynomial is x^14 (x - 2)^4. Its Jordan form has one block of size 4 for the eigenvalue 2 and
nilpotent blocks of sizes 3, 1, 1, ..., 1.

**(2) Counts.** Write N(L) for the number of tight colourings of length L, both first letters counted. Then
- sum_(k>=1) N(3k) x^k = 2x/(1 - 2x), so N(3k) = 2^k;
- sum_(k>=1) N(3k - 1) x^k = 2x(1 - x)/(1 - 2x)^2, so N(3k - 1) = (k + 1) 2^(k-1);
- sum_(k>=1) N(3k - 2) x^k = 2x(1 - x - x^2)/(1 - 2x)^3, so N(3k - 2) = (k^2 + 9k - 2) 2^(k-3).

The order of the pole at x = 1/2 is 1, 2, 3, which is one more than the slack 3k - L = 0, 1, 2. Up to complement
the counts are half of these, as THM-4951 (3c) lists them.

**(3) Junctions at n = 3k.** In the automaton extended by the last four letters and two junction flags
(317 reachable states), the tight colourings of length 3k - 2 that contain the factor 11011 have generating function
x^3/(1 - 2x)^2. So there are (k - 2) 2^(k-3) of them for k >= 3, and the same holds for 00100. The generating
function of those containing both factors is 0, which is THM-5020 (1)'s "J <= 1" seen by the automaton.

## Proof

- **(1), (2).** `markov_tight.py` builds the matrices and computes the following exactly in sympy: the
  characteristic polynomial, the Jordan form, and the generating functions v0 (I - xT)^(-1) T fin,
  x v0 (I - xT)^(-1) M0 M1 fin and x v0 (I - xT)^(-1) M0 fin. The closed forms follow by partial fractions. For
  example, the coefficient of x^k in 2x(1 - x - x^2)/(1 - 2x)^3 is
  2 [C(k+1, 2) 2^(k-1) - C(k, 2) 2^(k-2) - C(k-1, 2) 2^(k-3)] = (k^2 + 9k - 2) 2^(k-3).
- The automaton counts agree with brute force for every L <= 16 (markov_tight.out), and the closed forms with the
  automaton for k <= 25.
- **(3).** `markov_junctions.py` (markov_junctions.out); brute force agrees for k <= 6.

## Remarks (ANALOGY, typed)

- **A chain of independent trials.** The eigenvalue 2 of T is "one free colour per block of three positions". A tight
  colouring is a chain of independent binary trials, the block colours, as in the Euler product, where each prime
  divides a random integer independently.
  - Here the trials are not unconstrained: the slack 3k - L is a budget of defects (short blocks, junctions).
  - Each unit of slack adds one power of k: the pole order at 1/2. That is the Jordan block, the nilpotent "memory"
    of the chain. Structure (the rigid aligned blocks) is overlaid with randomness (the free colours). The defects are
    the only places where the structure bends.
  - The numerators 1, 1 - x, 1 - x - x^2 are the first Fibonacci-type truncations. This is recorded as numerology
    only (MISTAKE-559).
- **Primes as a Markov chain.** The residue mod 3 of the position (the phase) drives the automaton, as the last digit
  drives the observed Markov chains of consecutive primes (Uppuluri 1968; Lemke Oliver-Soundararajan 2016). The
  repository's THM-4473 records exact Markov laws for Collatz digit chains and the decaying bias for consecutive-prime
  digits.
  - T has a single nonzero eigenvalue and is otherwise nilpotent, like the Ulam and Tao chains (THM-4821/4822), which
    are rank one plus nilpotent.
  - No mechanism is claimed beyond the automaton itself.
- **Symmetry overlaid with asymmetry.** Complementation (swap the colours) and reversal act on the automaton's
  language. The blue (palindromic) tight colourings are the fixed points of reversal.
  - At n = 3k + 2 these are the palindromic block words, 2^ceil(k/2).
  - Their count grows like the square root of the black count: symmetric structures are rare but never absent.
