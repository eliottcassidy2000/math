---
id: THM-4852
title: "PRIMEGAME over the first ten primes: the two published lists (Wikipedia/OEIS ...15/2, 1/7, 55/1 - the owner's; Guy 1983 ...15/14, 15/2, 55/1) agree because every state-free configuration is 2^n 7^(D(n)-1) (verified n <= 2500; HYP-9448), and from there both reach 3^n 5^(n+1) 11; the owner's three pure fractions 1/7, 1/17, 55/1 = {7, 17, 55} erase the divisor counter, leave state 17 and restart the test; PRIMEGAME is a radius-1 cellular automaton on base-29# digits with the rule selected by the units digit; no adjacency-preserving bijection to the n = 5 merged tiling line metagraph exists"
status: PROVED (elementary; the all-n equivalence conditional on HYP-9448, unconditional for n <= 2500) + FINITE-EXACT + CITED (Conway 1987, Guy 1983, Kari 2012, THM-4729); count coincidences NUMEROLOGY; independently audited (see the note's section 8)
session: opus-2026-10-31-S70
date: 2026-10-31
note: 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
script: 04-computation/experiments/primegame_variant_phase_law_20261031.c (+ .out), primegame_selector_automaton_metagraph_20261031.py (+ .out)
related: [HYP-9448, THM-4729, THM-4730, THM-4733, THM-4838, THM-4759]
external:
  - "J. H. Conway, FRACTRAN: a simple universal programming language for arithmetic (1987); R. K. Guy, Conway's prime producing machine, Math. Mag. 56 (1983)"
  - "J. Kari, Cellular automata, the Collatz conjecture and powers of 3/2, DLT 2012 (x -> 3x/2 as a cellular automaton in base 6)"
  - "repository: procgen_atlas_20260924_collatz_implication_atlas.md row L1 (both published lists reproduce the first 15 primes)"
---

# THM-4852 - PRIMEGAME's phase law, its selector automaton, and the n = 5 metagraph

Write PRIMEGAME's states as exponent vectors over the ten primes 2..29: data registers 2, 3, 5, 7 and state registers
11, 13, 17, 19, 23, 29. The owner's list (the Wikipedia/OEIS form) ends 15/2, 1/7, 55/1. Guy's 1983 list ends 15/14,
15/2, 55/1.

**(1) Equivalence through the phase law.** In a configuration with no state prime, only the last three rules can act.
From 2^a 7^b with b <= a:
- Guy's list applies 15/14 b times, 15/2 a-b times and 55/1, reaching 3^a 5^(a+1) 11;
- the owner's list applies 15/2 a times, 1/7 b times and 55/1, reaching the same configuration.

So whenever every state-free configuration satisfies b <= a, the two orbits coincide at every phase boundary and emit
the same powers of 2.
- **Phase law** (HYP-9448): the k-th state-free configuration is 2^k 7^(D(k)-1), with D(k) = k/spf(k) (D(1) = D(p) = 1).
  This is the trial-division reading (Conway; Guy). It holds for both lists for every k <= 2500, and both emit exactly
  the 367 powers 2^p, p < 2500, in order.
- Hence the lists agree unconditionally for n <= 2500, and for all n given the law. The repository had checked the
  first 15 primes (procgen_atlas_20260924, L1).

**(2) Anatomy over 2..29.**
- The lcm of the denominators is 29# = 6469693230 for both lists, so PRIMEGAME is a Conway-Collatz map of modulus 29#.
- Every numerator and denominator is 1, a prime, a semiprime, or 78 = 2*3*13. Under the repository's F = U + S
  (THM-4733's divisor counts) the balanced parts are exactly the primes; no p^3 or p^2qr occurs.
- The owner's list has exactly three pure fractions, 1/7, 1/17 and 55/1 (Guy's has two):
  - 55/1 starts the next trial division at 3^n 5^(n+1) 11;
  - 1/17 leaves state 17;
  - 1/7 erases the counter 7^(D(n)-1).
- The product of the owner's 14 fractions is 55^2/(2 * 17^2). Guy's product is 15/2 times that.

**(3) Selector cellular automaton.** Let B = P# for a prime P. Every FRACTRAN program whose numerators and denominators
are squarefree products of primes <= P is a radius-1 cellular automaton on base-B digit strings, whose local rule is
selected by the units digit.
- For p, q | B, the quotient x/q has digits floor(x_i/q) + (x_(i+1) mod q) B/q.
- Then p y has digits (p y_i mod B) + floor(p y_(i-1)/B), with no further carry.
- The first applicable fraction is decided by x mod q_k = x_0 mod q_k.
- PRIMEGAME is the case P = 29. 20,000 steps of each list were reproduced digit by digit.
- Kari's base-6 automaton for 3x/2 is the one-rule prototype (CITED).

**(4) The n = 5 merged tiling line metagraph** (THM-4729): 10 nodes, 8 blue and 24 black lines, pure-blue nodes TT5,
1->C3->1 and T5. No bijection of the ten primes onto its ten nodes maps PRIMEGAME's register-flow graph (an arc p -> q
when p divides a denominator and q the numerator of the same fraction: 24 arcs, 19 simple pairs, no loops) onto it:
- it has 18 simple non-loop pairs and black loops;
- the simple degree sequences differ: flow 3,3,3,3,3,4,4,5,5,5 against metagraph 1,1,1,3,3,5,5,5,6,6.

The matches 10 = 4 + 6 (data + state registers, and base-path arcs + tiles) and 3 = 3 (pure fractions and pure-blue
nodes) are NUMEROLOGY.

**Novelty.**
- New to the repository: the phase-law explanation of why both published lists agree; the role of 1/7; the selector
  automaton formulation; the negative metagraph test.
- PRIMEGAME, its trial-division semantics and the 15-prime check are prior work.
