# The join route without binary Goldbach

CONDITIONAL.md and notes/26 (cited in wave 42) reduce Chan–Pak to **binary** Goldbach plus a physical
gluing. Binary Goldbach is OPEN. The conflict-graph view (OCF, THM-002) suggests a route where the
arithmetic input is a **theorem**.

**k-fold joins.** I(G_1 ∨ ⋯ ∨ G_k, 2) = Σ I(G_i, 2) − (k − 1). So a k-fold conflict-graph join of prime
witnesses would have count Σ p_i − (k − 1).

**Dials are joins too.** A leg 3 (one triangle, K_1) is a +2 dial and a leg 5 (K_2) is a +4 dial.
CONDITIONAL.md's three rows are therefore the binary join followed by 0, 1 or 2 universal cycles:
- p + q − 1 = J_2(p, q);
- p + q + 1 = J_3(p, q, 3);
- p + q + 3 = J_3(p, q, 5).

## Reduction (PROVED, conditional on the physical hypothesis)

**Hypothesis (J_{≤5}), unproved.** For some B ≥ 23, for every k ≤ 5 and all odd primes p_1, …, p_k ≠ 7
(repetition allowed, 3 and 5 included) with S = Σ p_i − (k − 1) ≥ B, a tournament with count S can be
built from witnesses of the p_i. As in CONDITIONAL.md, the rule may overlap and prepare its inputs.

**Claim.** (J_{≤5}), together with Helfgott's ternary Goldbach theorem (every odd integer > 5 is a sum of
three primes; Helfgott 2013, arXiv:1312.7748, with the Helfgott–Platt verification; CITED) and the
verified spectrum below B, implies Chan–Pak.

*Proof.* Let N ≥ B be the least odd count that is missing, other than 7 and 21. Apply Helfgott to
M = N + 2: M = p_1 + p_2 + p_3.
- **Twos.** Parity forces either no 2 or exactly two of them.
  - In the second case N = q + 2 = 3 + q − 1 is J_2(3, q), with q = N − 2 < N.
  - If q = 7 then N = 9 < B, so q ≠ 7 here.
- **Three odd primes.** Then N = J_3(p_1, p_2, p_3). Replace every leg 7 by the two legs 3 and 5; this
  keeps the join sum, since 3 + 5 − 1 = 7.
  - Three sevens would make N = 19 < B, so at most five legs remain.
  - Every leg is an odd prime ≠ 7 (hence not 21), and every leg is below N.
- **Close.** By minimality every leg has a witness, and (J_{≤5}) produces count N. ∎

The hypothesis must exclude the forbidden joins below B. Examples:
- J_3(3, 3, 3) = 7 (K_3);
- J_3(3, 3, 17) = J_3(5, 5, 13) = 21;
- J_2(3, 19) = J_2(5, 17) = J_2(11, 11) = 21 (SEVEN_TWENTYONE.md).

## What this changes, and what it does not

**What it changes.**
- Every arithmetic input is a theorem, so binary Goldbach can be dropped.
- Avoiding 7 needs no special trick: CONDITIONAL.md chooses c with N + c − 7 divisible by 3, whereas here
  each 7 splits into the dials 3 and 5.

**What it does not.**
- **The physical side is no easier.** Every pair of legs needs all its odd cycles to meet, so each pair
  inherits the overlap obstructions (ONE_VERTEX.md; PRIME_GLUE.md).
- **Inputs must be rebuilt.** The target map shows that, at counts like 27, 57 and 67, no output
  contains induced witnesses for any admissible pair. A successful (J_{≤5}) therefore cannot keep its
  inputs: it is a construction from the numbers p_i (a canonical family), not from given tournaments.
- **What is still missing.** Wave 42 §2 named exactly this ingredient: a family whose values are affine
  in a free parameter.

**Status.** The reduction is PROVED. (J_{≤5}) is OPEN. Chan–Pak is OPEN.
