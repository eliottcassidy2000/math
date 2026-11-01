---
id: THM-4857
title: "The GF(2) automaton of a tournament. For a tournament T with adjacency matrix A, let Phi_T(x) = Ax over GF(2): each vertex takes the parity of the vertices it beats. (1) Isomorphic tournaments, and a tournament and its converse, have conjugate automata, so Phi_T is an invariant of the merged metagraph's nodes. (2) GF(2)-similarity of A classifies the merged nodes exactly for n <= 5 (2, 3, 10 types: the owner's '12 classes become 10 nodes' is the partition by automaton type at n = 5). It fails at n = 6 (24 types for 34 nodes), at n = 7 (51 types for 272 nodes) and for every n >= 9 (counting). (3) A circulant Circ(Z_q, S) gives the additive cellular automaton with neighbourhood S, of rank q - deg gcd(a(x), x^q - 1) where a(x) = sum_{s in S} x^s; for q = 1 mod 4 the all-ones vector lies in the kernel. (4) A Paley tournament's automaton is reversible iff q = 3 mod 8. For q = 7 mod 8 its image is a binary quadratic-residue code of dimension (q+1)/2: the Hamming [7,4,3] code at q = 7 and the Golay [23,12,7] code at q = 23."
status: >
  PRE-AUDIT (mac-mini-2026-11-01-lenses).
  PROVED: (1); the counting half of (2); (3); (4) (elementary finite-field arguments).
  FINITE-EXACT (pre-registered):
  - (1) on all classes n <= 6 with 20 random relabellings each (A1);
  - (2) n = 3, 4, 5, 6 (S1, S2) and n = 7 (S4);
  - (3) on all 26,039 circulant connection sets of prime order <= 31 (A2), and the kernel claim for q = 5..29 (A4);
  - (4) for 13 Paley orders 3..83, with the minimum distances 3 and 7 computed (A3).
  CITED: binary quadratic-residue codes, and the Golay code as the length-23 QR code (MacWilliams-Sloane, The Theory
  of Error-Correcting Codes, ch. 16 and 20).
  EXPLORATORY first, then pre-registered: the n = 5 classification was observed in A5 (coarse dynamical type)
  before S1 tested it with full similarity.
  OPEN: n = 8, where the counting bound is not enough and nothing has been computed.
source: mac-mini-2026-11-01-lenses (owner directive of 2026-11-01: translate between Turing machines, cellular automata, FRACTRAN and the tournament research)
depends_on:
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
related:
  - 01-canon/theorems/THM-4855-circulant-tournaments-as-s-sequencings-segment-sum-criterion-and-explicit-non-palindromic-hamiltonian-paths-for-every-odd-circulant-of-order-at-least-seven.md
  - 01-canon/theorems/THM-4849-circulant-tournaments-are-rational-collatz-cycles-shift-decimation-reversal-two-base-points-and-the-locking-theorem-behind-35-and-175.md
  - 05-knowledge/results/session126_primegame_symmetry_20261031/metagraph_ca/README.md
  - 05-knowledge/results/lenses_tm_ca_fractran_tournaments_20261101.md
scripts:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/tournament_automaton.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/automaton_similarity_types.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/automaton_similarity_n7.py
outputs:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/tournament_automaton.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/automaton_similarity_types.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/automaton_similarity_n7.out
---

# THM-4857 - the GF(2) automaton of a tournament

## Definition

For a tournament T on V with adjacency matrix A (A[u][v] = 1 iff u -> v), the **GF(2) automaton** of T is the linear
map Phi_T : GF(2)^V -> GF(2)^V, x -> Ax. In words, each vertex's new bit is the parity of the bits held by the
vertices it beats.

For a circulant Circ(Z_q, S), (Phi x)_i = sum_{s ∈ S} x_{i+s}. This is the additive (linear) cellular automaton on
the cyclic line Z_q with neighbourhood S.

## Statement

**(1) Invariance.**
- If π : T -> T' is an isomorphism with permutation matrix P, then Phi_T' = P Phi_T P^(-1).
- The converse T^op has matrix A^T, and A and A^T are similar over every field.

So Phi_T up to conjugacy, i.e. the similarity class of A over GF(2), is an invariant of the node of T in the merged
metagraph (classes up to converse; THM-643, THM-4729).

**(2) Classification by similarity type.**

| n | classes | merged nodes | GF(2)-similarity types |
|---|---|---|---|
| 3 | 2 | 2 | 2 |
| 4 | 4 | 3 | 3 |
| 5 | 12 | 10 | 10 |
| 6 | 56 | 34 | 24 |
| 7 | 456 | 272 | 51 |

- **n <= 5.** Two tournaments have similar GF(2) adjacency matrices iff they are isomorphic or converse to each
  other. At n = 5 the ten automaton types are exactly the owner's ten merged nodes.
- **n = 6.** The classification fails. For example, TT6 and a tournament with scores (2,2,2,3,3,3) both have
  adjacency matrices nilpotent of index 6. Seven similarity types contain 2 or 3 merged nodes.
- **n >= 9.** It fails for every such n.
  - An n × n matrix over GF(2) has at most p(n) · 2^n similarity types, where p is the partition function: a
    rational canonical form is a list of invariant factors with degrees summing to n.
  - There are at least A000568(n)/2 merged nodes.
  - At n = 9 this is 15,360 against 95,768. The ratio 2^C(n,2) / (2 n! p(n) 2^n) increases from there on.

**(3) Circulants.** For q odd and a(x) = sum_{s ∈ S} x^s:
- rank_GF(2)(A) = q - deg gcd(a(x), x^q - 1), the standard circulant rank over a field where x^q - 1 is separable;
- if q = 1 (mod 4), then |S| = (q-1)/2 is even, so A𝟙 = 0.

**(4) Paley tournaments** (q prime, q = 3 mod 4, S = quadratic residues; α a primitive q-th root of unity over
GF(2)). The squaring map gives a(α)^2 = a(α^2), and a(α) + a(α^n) = sum_{x ≠ 0} α^x = 1 for any non-residue n.
- **q = 3 (mod 8).** Here 2 is a non-residue, so y = a(α) satisfies y^2 + y + 1 = 0 and never vanishes. The same
  holds at every α^j. Since a(1) = (q-1)/2 is odd, Phi is invertible: the automaton is reversible. Examples: C3,
  P_11, P_19, P_43, P_59, P_67, P_83.
- **q = 7 (mod 8).** Here 2 is a residue, so a(α) ∈ GF(2). Exactly one of the two classes {α^r : r a residue} and
  {α^n : n a non-residue} consists of roots of a. Hence gcd(a, x^q - 1) has degree (q-1)/2, the rank is (q+1)/2, and
  the image of Phi is the cyclic code generated by that gcd. That code is a binary quadratic-residue code
  (MacWilliams-Sloane, ch. 16).
  - q = 7: the Hamming [7,4,3] code.
  - q = 23: the binary Golay [23,12,7] code (computed: dimension 12, minimum distance 7).

## Proof notes

- **(1)** is immediate.
- **(2), small n.** Similarity is decided by the ranks of p(A)^k for every irreducible p over GF(2) and every k
  (elementary divisors). The tables are exhaustive over the 2^C(n-1,2) path tilings of order n, which realise every
  class (Rédei).
- **(2), n >= 9.** The counting argument above.
- **(3)** is standard.
- **(4)** is given in the statement.

## The n = 5 merged metagraph through this lens (`tournament_automaton.out`, A5)

| owner node | type | GF(2) rank | nilpotent depth | periodic-part cycles |
|---|---|---|---|---|
| #0 TT5 | pure-blue | 4 | 5 | (1) (nilpotent) |
| #1/#7 | pure-black | 4 | 2 | 1,1,3,3 |
| #2/#5 | pure-black | 4 | 1 | 1,1,2,4,4,4 |
| #3 1->C3->1 | pure-blue | 3 | 1 | 1,1,3,3 |
| #4 | mixed | 5 | 0 | 1,1,2,7,7,14 |
| #6 | mixed | 5 | 0 | 1,31 (a single 31-cycle: maximal period) |
| #8 | mixed | 5 | 0 | 1,1,5,5,5,5,5,5 |
| #9 | mixed | 5 | 0 | 1,3,7,21 |
| #10 | mixed | 4 | 1 | 1,1,1,1,2,2,4,4 |
| #11 T5 | pure-blue | 4 | 1 | 1,15 |

- The reversible automata (rank 5) are exactly four of the five mixed nodes.
- Node #6's automaton cycles through all 31 nonzero states, so its characteristic polynomial is primitive of
  degree 5.
- No automaton feature common to exactly the three pure-blue nodes was found. Their ranks are 4, 3 and 4, and TT5's
  automaton is nilpotent while T5's has a 15-cycle. So the pure-blue property is NUMEROLOGY under this lens.
- Session126's six-cell rule-150 automaton acts on tilings (64 masks), not on vertices. It is a different
  automaton: it does not descend to the merged nodes (session126).
