# Pre-registration — kind-pasteur-2026-10-27-S8 (written BEFORE the confirmatory runs below)

Provenance note. Commit 181cd92f0f was meant to add this file, and its message says so, but a failed shell step left
it out; only the exploratory census was committed then. The file is added in the next commit. No confirmatory script
(P1–P3 below) had been written or run at either time.

## Exploratory runs so far (NOT pre-registered)

These come from `carrier_census.py`, n = 3..7, committed in 181cd92f0f.

- **Objects.**
  - An oriented excess -1 carrier D has underlying graph G with n vertices and n-1 edges.
  - A = Aut(D). Pi is the set of complement pairs, a signed A-set: an automorphism may swap the two ends of a pair.
- **Linear trijectivity.** D is LINEARLY TRIJECTIVE if Pi is isomorphic, as a signed A-set, to the pair set of some
  A-set X with |X| = n-1. Then {T containing D} is Tour(X), A-equivariantly, by an F2-affine isomorphism.
- **Weak trijectivity.** D is WEAKLY TRIJECTIVE if {T containing D} is merely an A-set isomorphic to Tour(X).
- **Census** of oriented carriers, written linear / weak-only / neither:

  | n | linear | weak-only | neither |
  |---|---|---|---|
  | 3 | 3 | 0 | 0 |
  | 4 | 10 | 0 | 0 |
  | 5 | 39 | 2 | 0 |
  | 6 | 167 | 20 | 0 |
  | 7 | 816 | 188 | 5 |

  At n = 7 the 2 uniform stars (|A| = 720) are not tested by the census; they are counted as linear by the tree
  argument.
- **Agreement with a criterion.** In all 1250 tested carriers, "linear" coincided with criterion
  (C): chi(G^K) = |V^K| - |E^K| = 1 for every subgroup K <= A. Equivalently, [V] = [E] + 1 in the Burnside ring B(A).
  The cyclic-subgroup version of (C) agreed everywhere too.
- **Timing of (C).** (C) was formulated from the n <= 6 data before the n = 7 run, but was not written down in
  advance. So the n = 7 agreement is recorded as exploratory.

## Proof status at the time of writing (to be audited)

- **Sufficiency of (C)**, via marks of A x Z2-sets. (C) implies linear trijectivity, with X = V - {r} ≅ E for any
  A-fixed vertex r.
- **Necessity of (C)** is proved for every 2-subgroup: |V^K| - |X^K| is constant along index-2 chains. Necessity for
  odd-order subgroups is OPEN.

## Predictions (each to be run after this file is pushed)

### (P1) Procedural families, out of sample (n up to 11)

- **(a) Bigons.** B_k is two directed s -> t paths of length k plus one isolated vertex, so n = 2k+1, k = 2..5.
  Prediction: each B_k is NOT linearly trijective, IS weakly trijective, and has exactly 2^(C(2k,2)-1) marked types.
- **(b) Alternating even cycles.** C_{2k} + K_1 with alternating orientation, k = 2..5. Prediction: NOT linear, IS weak.
- **(c) Dense cores.** K_4 (transitive) + K_2 + 2K_1 at n = 8 is not linear. The swap of the two isolated vertices
  fixes a subgraph with chi = -1.
- **(d) Positive controls at n = 8, 9.**
  - **Tree sample.** 300 random orientations of random trees. Prediction: the explicit end-edge map is an equivariant
    bijection; each has (C), and is linear.
  - **Functional sample.** 300 random functional orientations: one tree component plus unicyclic components with
    directed cycles, all arcs toward the root or the cycles. Prediction: the explicit root-swap map is an equivariant
    bijection; each has (C), and is linear.

### (P2) "Pairs and triples" (the Maxwell ladder)

- A carrier with jn - C(j+1,2) edges has a complement of exactly C(n-j,2) pairs.
- Definition: D is linearly j-jective if Pi is isomorphic, as a signed A-set, to the pair set of an A-set of size n-j.
- Prediction, over all oriented carriers for j = 2 (n = 4, 5, 6) and j = 3 (n = 5, 6): linear j-jectivity holds iff the
  equivariant Maxwell count |E^K| = j|V^K| - C(j+1,2) holds for every subgroup K <= A.

### (P3) Coherent carrier sets

The coherent carrier set Sigma(T) ⊆ {L, C, S} records which of three structures T contains:
- L: a directed Hamiltonian path;
- C: a directed (n-1)-cycle plus a solo vertex;
- S: a source or a sink.

Predictions, for every tournament type with 4 <= n <= 7:
- Sigma(T) is a pair or a triple whenever n <= 5.
- Sigma(T) = {L} occurs exactly for the non-strong T whose first and last strong components have >= 3 vertices.
- With s = strong types and t = all types, the number of types per class is:

  | Sigma(T) | number of types |
  |---|---|
  | {L,C} | s_n |
  | {L,C,S} | 2 s_(n-1) |
  | {L,S} | 2 t_(n-1) - t_(n-2) - 2 s_(n-1) |
  | {L} | the remainder |

- The {L} counts are n = 6: 1, n = 7: 3, n = 8: 16.
