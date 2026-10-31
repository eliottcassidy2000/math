---
id: THM-4854
title: "PRIMEGAME in four blocks and the merged metagraphs. (A) CITED (Conway 1987 sections 5 and 7; session126; mac-mini; THM-4852; session129's same-day packet with the same four blocks X0-X3) and PROVED where new: the owner's blocks {11/13, 13/11}, {55/1, 1/7, 1/17, 17/91}, {19/51, 77/19, 77/29, 29/33}, {78/85, 23/38, 95/23, 15/2} (the control toggle, the 7-and-5 block with start and exit, the 3-to-7 block, the 2-and-5 exchange) are not the fibres of the per-instruction data-register map, which has six fibres, but touch the distinct register unions none, {5, 7}, {3, 7}, {2, 3, 5}, and session129's height grading separates them intrinsically; the three data blocks pairwise share two flags and all contain the decision flag 17 (the only control prime with three instructions); the flag digraph without the two state-free loops is bipartite with parts {11, 17, 23} and {13, 19, 29, none} (session129's K_(2,3) with pendants), so every closed walk of flag moves avoiding those loops is even; 11/13, 13/11 is the program's only reciprocal pair. (B) FINITE-EXACT: among all ordered products of at most four fractions, the only products equal to 1 or to a fraction are those with the toggle pair inserted; no product with one factor from each data block is a fraction (Conway's macro pairs AB, EF, DG are products that are not single fractions). (C) CITED (the comment of OEIS A068443) with the proof written out: for primes p < q, pq is triangular iff q = 2p +- 1; the program's triangular entries (session129's identities) are 1, 15 = 3 x 5, 55 = 5 x 11, 91 = 7 x 13 and 78 = 6 x 13, its Fibonacci entries 1, 2, 13, 55 (Ming 1989: the triangular Fibonacci numbers are 1, 3, 21, 55; CITED). (D) Remark (definitional; THM-4809 and THM-4729 Proposition K): the merged tiling line metagraph of order n is the quotient by isomorphism, with converse classes merged, of the labelled complement-move graph of THM-4810 restricted to moves along the base path; THM-4732 Lemma 1 and THM-4810's Theorem P are the odd and even cases of one toggle-count parity (THM-4815 (C)); the L = 2 case of THM-4732 Lemma 1 and THM-4810's two-move criterion are two instances of the composition law; the three-move law (HYP-9347 (1)) says the base path is the symmetric difference of three successive black legal paths. FINITE-EXACT: the merged metagraphs at orders 3, 4, 5 (2, 3, 10 nodes; 1, 4, 32 lines; pure-blue nodes both, T_4, and TT_5 / 1 -> C_3 -> 1 / T_5) recomputed independently and agreeing with THM-4729 / S70 line by line. (E) FINITE-EXACT tests of the owner's mappings: the involution pair against the unique order-3 line is a dictionary match and the three data blocks against the three order-4 nodes is not a graph correspondence (a triangle against a star with a loop), both session129's verdicts reproduced; the ten composite entries against the ten order-5 nodes preserve at most 15 of 20 adjacencies over all 10! bijections, against a mean of 15.25 for random graphs of the same size (NUMEROLOGY)."
status: >
  PROVED (new parts only): the fibre check, the pairwise flag intersections, the hub 17, the unique reciprocal pair; the semiprime
  proof is written out for a known OEIS comment.
  CITED: Conway 1987 (sections 5 and 7), session126 section 3, mac-mini section 4, THM-4852 (the data / flag split and the
  table); session129 (05-knowledge/results/session129_primegame_tournament_partition_20261031: the same blocks X0-X3, the
  height grading, the K_(2,3) flag graph, the block products, the merged-node census 2, 3, 10, 34, the order-3 and order-4
  verdicts, the triangular identities); OEIS A068443; THM-4809 and THM-4729 Proposition K (the quotient); THM-4815 (C)
  (the parity link to THM-4732 Lemma 1).
  FINITE-EXACT: (B), the small metagraphs of (D) and the tests of (E): primegame_blocks_metagraph_20261031.py (+ .out).
  CITED: Conway 1987 Theorem 1 and THM-4852 (the program), THM-4853 (step counts), THM-4729 (merged metagraph, closed and
  open tilings), THM-4731 / THM-4732 / HYP-9347 (lines, black detours, the three-move law; the (n-4, 3) theorem for
  n >= 11 in the endpoint-three packets), THM-4810 (the complement-move graph), L. Ming, Fibonacci Quart. 27 (1989)
  98-108 (triangular Fibonacci numbers), session126 and mac-mini-2026-10-31-fivegroups (the correspondences already
  tested: ten primes against ten nodes, the pure fractions against the pure-blue nodes, F = U + S).
  DICTIONARY / NUMEROLOGY / REFUTED: as typed in (E) and in the note's section 5.
  AUDITED: PASS WITH CORRECTIONS (every number reproduced with the auditor's own code; seven majors, all credit and two
  false sentences, applied; fourteen minors applied; MISTAKE-699; report beside the note).
  OPEN: nothing new; the structural proof of the three-move law for n >= 8 remains HYP-9347 (1).
session: opus-2026-10-31-S71
source: 05-knowledge/results/primegame_blocks_register_decomposition_metagraphs_20261031.md
related:
  - THM-4852
  - THM-4853
  - THM-4729
  - THM-4732
  - THM-4810
  - HYP-9347
---

# PRIMEGAME in four blocks and the merged metagraphs

Statements, tables, proofs and the exact tests: the [source note](../../05-knowledge/results/primegame_blocks_register_decomposition_metagraphs_20261031.md).

**(A)** Each fraction of PRIMEGAME changes one flag (consumed control prime to produced control prime, none when absent) and a set of data registers; the owner's blocks are not the fibres of "which data registers move" (that map has six fibres: none for 1/17, 11/13, 13/11; {2}; {3}; {5}; {7}; {2, 3, 5}); their register unions are none, {5, 7}, {3, 7}, {2, 3, 5}. Flag sets {11, 13}, {11, 13, 17}, {11, 17, 19, 29}, {13, 17, 19, 23}; pairwise intersections of the data blocks {11, 17}, {13, 17}, {17, 19}. Bipartition of the flag digraph: parts {11, 17, 23} and {13, 19, 29, none}; each of the twelve non-loop instructions crosses between them. The only two fractions with product 1 are 11/13 and 13/11.

**(B)** Products of 2, 3, 4 fractions equal to 1 or to a fraction: 2, 78, 6, all of the form (toggle pair inserted); none from B2 x B3 x B4.

**(C)** 2pq = k(k+1) with k, k+1 coprime forces {k, k+1} = {2p, q}, hence q = 2p +- 1 (the other split {p, 2q} is impossible for p < q); checked for every semiprime below 10^5 (OEIS A068443, CITED).

**(D)** A tiling is a point of the marked carrier X_n; kappa is tau_P along the base path; sigma is converse followed by i -> n-1-i. Hence the merged metagraph is the quotient of THM-4810's complement-move graph along the base path, black K moves are re-markings along other paths, THM-4732 Lemma 1 ("an even number L of black moves needs n <= 2L+2 and C(n-1,2) even") and THM-4810's Theorem P are the odd and even cases of the toggle-count parity, as THM-4815 (C) records; the L = 2 case of Lemma 1 and THM-4810's Theorem C2 are two instances of the composition law. Small metagraphs: order 3 (T_3 - C_3, one blue line), order 4 (T_4 - S_4 blue, pair - S_4 black twice, S_4 - S_4 blue loop), order 5 (10 nodes, 8 blue and 24 black lines as in THM-4729).

**(E)** As stated in the title; the calibration uses 20 random 20-edge graphs on 10 vertices with the same maximum-preserved-adjacency statistic over all 10! bijections.
