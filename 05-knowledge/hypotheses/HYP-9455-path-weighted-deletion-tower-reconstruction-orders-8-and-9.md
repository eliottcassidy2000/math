---
id: HYP-9455
title: "Path-weighted deletion tower: every converse-merged tournament node of order n >= 4 is determined, up to a positive scalar, by its Hamiltonian-endpoint-weighted converse-merged deck W(T) = sum_x (h_s(x) + h_e(x)) [T - x] (equivalently by THM-796's primitive face row), and the face matrix R_n has full column rank |M_(n-1)|; pre-registered tests at orders 8 (3,528 nodes) and 9 (97,144 nodes)"
status: OPEN - PRE-REGISTERED (opus-2026-11-04-S72); predictions P1-P5 committed before the confirmatory run at orders 8 and 9
session: opus-2026-11-04-S72
date: 2026-11-04
related: [THM-4863, THM-796, THM-830, THM-4854, HYP-6865]
note: 05-knowledge/results/primegame_deletion_tower_path_weighted_decks_20261104.md
script: 04-computation/experiments/deletion_tower_census_20261104.c
---

# HYP-9455 - the path-weighted deletion tower

## Setting

For a tournament T and a vertex x, let h_s(x) and h_e(x) be the numbers of Hamiltonian paths of T that start and end
at x, and let [T - x] be the converse-merged isomorphism class (node) of the vertex-deleted subtournament. Put

    W(T) = sum over x of (h_s(x) + h_e(x)) e_[T - x],

a nonnegative integer vector indexed by the merged nodes M_(n-1) of order n - 1, with total weight 2 H(T).

**Lemma W** (THM-4863, proved in the note, section 1). THM-796's face matrix
R_n(u, v) = #{tilings t in the fibre of u : the endpoint deletion of t lies in v} satisfies

    R_n(u, .) = W(T) / (2 |Aut T|)    if u = {T} is self-converse,
    R_n(u, .) = W(T) / |Aut T|        if u = {T, T^op} is a converse pair.

So THM-796's primitive face coordinate Rbar_n(u) is the primitive vector of W(T). The census program checks Lemma W
against every fixed-path tiling, through both endpoint faces, at orders 3-7 (no mismatch).

## Statement

**(A) Reconstruction.** For every n >= 4, two converse-merged nodes of order n with proportional vectors W are equal.
**(B) Rank.** For every n >= 3, the face matrix R_n has rank |M_(n-1)| over Q, i.e. full column rank.

Order 3 is the only exception to (A): TT_3 and C_3 both have W proportional to the single order-2 node.

## Evidence before this pre-registration

- (A) holds at orders 4, 5, 6, 7, with 3, 10, 34 and 272 distinct primitive rows. THM-796's table has the same
  counts, and this session's census program recomputes them independently.
- (B) holds at orders 3-7, with ranks 1, 2, 3, 10, 34. These are certified by the rank of the Gram matrix R^T R
  modulo 2^31 - 1.
- The weights are necessary. The unweighted converse-merged deck (the multiset of [T - x]) has collisions:
  - at order 5: 8 decks for 10 nodes;
  - at order 7: 267 decks for 272 nodes.
  It is injective at order 6. The deck's support set is coarser still: 6 distinct values at order 5, 264 at order 7.
- These order 3-7 runs were exploratory: they validated the program and are not part of the test.

## Pre-registered predictions

Confirmatory run: `deletion_tower_census 9` (program `04-computation/experiments/deletion_tower_census_20261104.c`),
first executed at orders 8 and 9 after the commit containing this file.

- **(P1)** At order 8 the primitive rows of the 3,528 merged nodes are pairwise distinct: (A) holds at n = 8.
- **(P2)** At order 9 the primitive rows of the 97,144 merged nodes are pairwise distinct: (A) holds at n = 9.
- **(P3)** R_8 has rank 272 and R_9 has rank 3,528: (B) holds at n = 8, 9.
- **(P4)** The unweighted converse-merged deck is not injective at order 8, nor at order 9. Both orders have at least
  one collision.
- **(P5)** These sanity checks are not predictions:
  - order 8: 6,880 ordinary classes (OEIS A000568), 176 self-converse (A002785), 3,528 merged nodes;
  - order 9: 191,536 ordinary classes, 2,752 self-converse, 97,144 merged nodes;
  - the node tilings sum to 2^21 and 2^28;
  - Lemma W's integrality holds at every node.

A failure of P1 or P2 refutes (A) at that order and will be reported with the colliding nodes. A failure of P3 is
reported only if the modular certificate falls short, in which case the rank over Q is recomputed exactly. A failure
of P4 (injective unweighted decks) is reported as such.

## Why it would matter

- The owner's picture is "one level determines the next": the three order-4 nodes determine the ten order-5 nodes
  (session129, THM-796), and the prime factorization of a number determines the number. In that picture, (A) is the
  unique-factorization statement of the tower: the order n - 1 nodes are the primes, and W is the exponent vector.
- (B) says that each level needs all the nodes of the previous level as independent coordinates. At order 6 these are
  exactly ten, the number of PRIMEGAME's registers.
- Without the Hamiltonian-path weights, the analogous statement for ordinary decks fails. Tournament reconstruction is
  false (Stockmeyer 1977), and the converse-merged deck fails at orders 5 and 7, as found above.
