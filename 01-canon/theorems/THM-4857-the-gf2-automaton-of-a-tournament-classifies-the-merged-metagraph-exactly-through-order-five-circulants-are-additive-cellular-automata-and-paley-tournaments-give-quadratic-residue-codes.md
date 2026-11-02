---
id: THM-4857
title: "The GF(2) automaton of a tournament. For a tournament T with adjacency matrix A, let Phi_T(x) = Ax over GF(2): each vertex takes the parity of the vertices it beats. (1) Isomorphic tournaments, and a tournament and its converse, have conjugate automata, so Phi_T is an invariant of the merged metagraph's nodes. (2) GF(2)-similarity of A classifies the merged nodes exactly for n <= 5 (2, 3, 10 types; at n = 5 the owner's ten merged nodes are the ten automaton types, a small-n coincidence). It fails for every n >= 6: n = 6 (24 types for 34 nodes) and n = 7 (51 for 272) by computation, n >= 8 by counting. (3) A circulant Circ(Z_q, S) gives the additive cellular automaton with neighbourhood S, of rank q - deg gcd(a(x), x^q - 1) where a(x) = sum_{s in S} x^s; for q = 1 mod 4 the all-ones vector lies in the kernel. (4) KNOWN (THM-481; Kim-Solé, Assmus-Key, T. S. Michael): for prime q = 3 mod 4, a Paley tournament's automaton is reversible iff q = 3 mod 8, and for q = 7 mod 8 its image is a binary quadratic-residue code of dimension (q+1)/2 (Hamming [7,4,3] at q = 7, Golay [23,12,7] at q = 23)."
status: >
  Independently audited 2026-11-01: SOUND WITH CORRECTIONS, all applied
  (04-computation/experiments/circulant_fractran_goldbach_20261101/audit_lenses/AUDIT.md; MISTAKE-701).
  PROVED: (1); the counting half of (2); (3) (standard).
  KNOWN / CITED: (4).
  - THM-481 (June 2026) identifies the Paley tournament's GF(2) code with the binary QR code, with the Golay code at
    q = 23 and full rank at q = 11. It credits Kim-Solé (Des. Codes Cryptogr. 49 (2008)), Assmus-Key and T. S. Michael.
  - The dichotomy holds for every doubly regular tournament of order n = 4t - 1. There A A^T = tI + (t-1)J, so over
    GF(2): A^(-1) = A^T when n = 3 (mod 8), and A^2 = A when n = 7 (mod 8) (audit).
  - Further credits: Sutner's additive automata on graphs (Complex Systems 2, 1988); Taussky-Zassenhaus (1959) for
    A ~ A^T; de Caen (1991) for rank >= (n-1)/2.
  - No prior source was found for the n <= 5 classification.
  FINITE-EXACT (pre-registered):
  - (1) on all classes n <= 6 with 20 random relabellings each (A1);
  - (2) n = 3, 4, 5, 6 (S1, S2) and n = 7 (S4);
  - (3) on all 26,039 circulant connection sets containing 1, of prime order <= 31 (A2), and the kernel claim for
    q = 5..29 (A4);
  - (4) for 13 Paley orders 3..83, with the minimum distances 3 and 7 computed (A3).
  CITED: binary quadratic-residue codes, and the Golay code as the length-23 QR code (MacWilliams-Sloane, The Theory
  of Error-Correcting Codes, ch. 16 and 20).
  EXPLORATORY first, then pre-registered: the n = 5 classification was observed in A5 (coarse dynamical type)
  before S1 tested it with full similarity. S1 could not fail: coarse types <= similarity types <= merged nodes, and
  A1 had already printed 2, 3, 10 coarse types. S2 and S4 were genuine tests.
  No n is left open: n >= 8 follows from the exact count of similarity classes of 8 x 8 matrices over GF(2), 746
  (OEIS A070933), which is below the at least 6880/2 = 3,440 merged nodes at n = 8.
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

(For comparison, the number of all similarity classes of n x n matrices over GF(2) is 2, 6, 14, 34, 74, 166, 350,
746, ... for n = 1, 2, ..., 8 (OEIS A070933).)

- **n <= 5.** Two tournaments have similar GF(2) adjacency matrices iff they are isomorphic or converse to each
  other. At n = 5 the ten automaton types are exactly the owner's ten merged nodes.
- **n = 6.** The classification fails. For example, TT6 and a tournament with scores (2,2,2,3,3,3) both have
  adjacency matrices nilpotent of index 6. Seven similarity types contain 2 or 3 merged nodes.
- **n >= 8.** It fails for every such n.
  - At n = 8 there are 746 similarity classes of 8 × 8 GF(2) matrices in all (A070933), fewer than the
    at least 6880/2 = 3,440 merged nodes.
  - For n >= 9, an n × n matrix over GF(2) has at most p(n) · 2^n similarity types, where p is the partition
    function: a rational canonical form is a list of invariant factors with degrees summing to n.
  - There are at least A000568(n)/2 merged nodes. At n = 9 that is 15,360 types against 95,768 nodes, and the ratio
    2^C(n,2) / (2 n! p(n) 2^n) increases from there on.
  - (At n = 7 counting is not enough: there are 350 similarity classes against 272 nodes, so the computation was
    needed.)

**(3) Circulants.** For q odd and a(x) = sum_{s ∈ S} x^s:
- rank_GF(2)(A) = q - deg gcd(a(x), x^q - 1), the standard circulant rank over any field (no separability is
  needed);
- if q = 1 (mod 4), then |S| = (q-1)/2 is even, so A𝟙 = 0.

**(4) Paley tournaments** (q prime, q = 3 mod 4, S = quadratic residues; α a primitive q-th root of unity over
GF(2)). The squaring map gives a(α)^2 = a(α^2), and a(α) + a(α^n) = sum_{x ≠ 0} α^x = 1 for any non-residue n.
- **q = 3 (mod 8).** Here 2 is a non-residue, so y = a(α) satisfies y^2 + y + 1 = 0 and never vanishes. The same
  holds at every α^j. Since a(1) = (q-1)/2 is odd, Phi is invertible: the automaton is reversible. Examples: C3,
  P_11, P_19, P_43, P_59, P_67, P_83.
- **q = 7 (mod 8).** Here 2 is a residue, so a(α) ∈ GF(2). Exactly one of the two classes {α^r : r a residue} and
  {α^n : n a non-residue} consists of roots of a. Hence gcd(a, x^q - 1) has degree (q-1)/2, the rank is (q+1)/2, and
  the image of Phi (the column space of A) is the cyclic code generated by gcd(b(x), x^q - 1). Here b is the
  non-residue polynomial; since -1 is a non-residue, it is the reciprocal of gcd(a, x^q - 1). That code is the other
  binary quadratic-residue code, with the same parameters [q, (q+1)/2] (THM-481; MacWilliams-Sloane, ch. 16).
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
  automaton is nilpotent while T5's has a 15-cycle. No structural link to pure-blueness was found under this lens.
- Session126's six-cell rule-150 automaton acts on tilings (64 masks), not on vertices. It is a different
  automaton: it does not descend to the merged nodes (session126).
