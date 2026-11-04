---
id: HYP-9455
title: "Path-weighted deletion tower: every converse-merged tournament node of order n >= 4 is determined, up to a positive scalar, by its Hamiltonian-endpoint-weighted converse-merged deck W(T) = sum_x (h_s(x) + h_e(x)) [T - x] (equivalently by THM-796's primitive face row), and the face matrix R_n has full column rank |M_(n-1)|; pre-registered tests at orders 8 (3,528 nodes) and 9 (97,144 nodes)"
status: OPEN (opus-2026-11-04-S72). P1, P2, P3 CONFIRMED at orders 8 and 9 (3,528 and 97,144 distinct rows; ranks 272 and 3,528); P4 REFUTED at order 8 (the unweighted merged deck is injective there) and confirmed at order 9; second pre-registration P6-P7 (order 10) committed before its run
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

## Results of the confirmatory run (2026-11-04, after the pre-registration commit ce33106a23)

The run was `deletion_tower_census 9` (38 s); its output is `04-computation/experiments/deletion_tower_census_20261104.out`.

**(P5) Sanity checks: all pass.**
- Class counts: 6,880 / 176 / 3,528 at order 8, and 191,536 / 2,752 / 97,144 at order 9.
- The tilings sum to 2^21 and 2^28.
- No Lemma W integrality failure occurred.

**(P1) CONFIRMED.** The 3,528 primitive rows at order 8 are pairwise distinct.

**(P2) CONFIRMED.** The 97,144 primitive rows at order 9 are pairwise distinct.

**(P3) CONFIRMED.** The rank certificates are full: rank_(2^31-1)(R^T R) = 272 at order 8 and 3,528 at order 9, so R_8
and R_9 have full column rank over Q.

**(P4) REFUTED at order 8, confirmed at order 9.**
- At order 8, the unweighted converse-merged deck is injective: 3,528 distinct decks. So is its support set.
- At order 9 it has exactly three collision pairs, giving 97,141 distinct decks.

Together with the exploratory orders, the unweighted merged deck fails exactly at the odd orders tested and succeeds
at the even ones:

| order | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| merged nodes | 2 | 3 | 10 | 34 | 272 | 3,528 | 97,144 |
| distinct unweighted merged decks | 1 | 3 | 8 | 34 | 267 | 3,528 | 97,141 |
| collision classes | 1 | 0 | 2 | 0 | 5 | 0 | 3 |

This parity pattern was found after the fact; it was not predicted. The path weights repair every odd-order failure.

## Second pre-registration: order 10 (committed before the order-10 run)

The run uses a separate memory-light program, `04-computation/experiments/deletion_tower_order10_20261104.c`. It
repeats the census through order 9, then generates order 10. Rows are compared by 64-bit hashes, and every hash
collision is resolved by recomputing the exact vectors.

- **(P6)** At order 10, all 4,870,920 merged nodes have pairwise distinct primitive W rows: (A) holds at n = 10.
- **(P7)** At order 10, the unweighted converse-merged deck is injective: there are no collisions, continuing the
  even-order pattern.
- **(P8) Sanity checks, not predictions:**
  - 9,733,056 ordinary classes (OEIS A000568);
  - 8,784 self-converse classes (A002785);
  - 4,870,920 merged nodes.

The rank statement (B) is not tested at order 10: the Gram matrix would be 97,144 x 97,144.
