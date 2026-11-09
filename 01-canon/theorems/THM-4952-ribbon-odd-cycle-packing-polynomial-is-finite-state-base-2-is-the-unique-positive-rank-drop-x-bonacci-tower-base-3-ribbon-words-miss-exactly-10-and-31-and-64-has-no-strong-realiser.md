---
id: THM-4952
title: "Ternary Chan-Pak on ribbon words. (A) The odd-cycle packing polynomial I(Omega(T_u), x) of a span-w ribbon word (THM-4880) is a finite-state series: a frontier automaton on fragments (head age, tail age or free, length parity) with closed odd cycles weighted x, 4 / 20 / 112 states for w = 1 / 2 / 3; at x = 2 it reproduces the Hamiltonian-path automaton (OCF). (B) The span-2 series has minimal dimension 14 over Q(x), dropping exactly at x = 0 (to 1) and x = +-2 (to 13): a nonsingular 14 x 14 minor has determinant -x^13 (x - 2)(x + 2), so the Hamiltonian-path base is the unique positive base of lower dimension, and 13 is THM-4880 (G)'s irreducible dimension. Integer representations: rank 14 with entries <= 5 at x = 3 (span 2), rank 70 with entries <= 4 at x = 3 (span 3). (C) The tower 1-^n has I(Omega, x) = f(n) with f(n) = f(n-1) + f(n-2) + (x - 1) f(n-3), f(0) = f(1) = f(2) = 1 (odd-interval packings); tribonacci at x = 2, A122552 at x = 3 with H_3 = (3 * 2^n + c_n)/7, c_n = 4, 1, -5 periodic: exactly 2^n growth. (D) Complete finite census: every v = 1 (mod 3), v <= 10^4, other than 10 and 31 is H_3 of a span-2 ribbon word (fixed point of the pruned state sets at 1,197,559 states, so for all lengths); strong span-2 words miss exactly 10, 31, 64. (E) 64 is the base-3 packing count of no strong tournament (H <= H_3 = 64 < f(9) forces order <= 8; exhaustive), though 64 = H_3(C3 => T) for any T with H_3 = 16: the smallest value other than 1, in any base, that is attained but has no strong realiser (bases 2, 4-7 at orders <= 8 and the block lemma for b >= 8). (F) NEGATIVE: near-join digits do not realise the ternary return forest: of the 97 odd h <= 199 attained at orders <= 8, only 50 admit a root response R_t = h + 1 and 55 admit R_t = h + 2."
status: >
  PROVED (mac-mini-2026-11-08-ternary): (A) the automaton (the frontier method is KNOWN; the tournament automaton
  is new), verified on 360 random words of orders <= 9 against a set-partition brute force and, at x = 2, against
  THM-4880's path automaton; (B) exact: the determinant of a nonsingular minor, the ranks at its roots, and the
  generic rank 14 certified by rank 14 at more integer specialisations than the degree bound of a 15 x 15 minor
  (the functional rows are products of letter matrices, so their entries have degree 2 at depth 4 and 3 at depth 5,
  bounds 30 and 45; the audit checked 65 and 61 points, rankdrop.py now checks 36 at depth 4; the session's first
  certificate used the bound 15 of a single letter matrix, MISTAKE-740);
  (C) the interval argument, with the recurrence also verified symbolically against the automaton for n <= 13 and
  the base-3 closed form from (x - 2)(x^2 + x + 1); (D) by the finite fixed point: the pruned state sets of
  ribbon_levels_mod.c stabilise from length 19 and the pruning is sound by THM-4940 (1) (monotonicity under vertex
  deletion); (E) by Busch's floor (THM-4744, CITED) and the exhaustive packing polynomials of all 6,880 + ... classes
  of orders <= 8 (spec3_small.py). FINITE-EXACT: (F) readiness.py over all tournaments of orders <= 8 and all
  vertices; the span-3 census to 5,000 and the span-2 census to 10^5 (values complete by lengths 10 and 15; the
  state sets were capped, so the fixed point is claimed only at 10^4, where the audit verified it as set equality). KNOWN/CITED:
  OCF (Grinberg-Stanley), THM-4940 (H_b definitions, monotonicity, multiplicativity, the gap tables reproduced here
  at orders <= 8), THM-4865 (c) (tribonacci tower), A122552 and A000213 (OEIS). NOT CLAIMED: ternary Chan-Pak
  (HYP-9530 (A), OPEN) or Chan-Pak (OPEN). Concurrent: opus-2026-11-09-S87's order-10 strong census and lane
  reduction (THM-4950/4951), not used here. AUDIT: PASS WITH CORRECTIONS (independent, adversarial;
  05-knowledge/results/ternary_ribbon_metagraph_audit_20261108.md; one major, the generic-rank certificate,
  repaired as above; minors: the fixed point verified as set equality, the "first composite" claim justified for
  every base, the base-2 strong spectrum cited as finite data, the x = 2 identification spelled out), all applied on
  2026-11-09; the audit reproduced the automaton, the tower, the base-3 letter polynomials, the order-<= 8 census,
  the near-join formula and the readiness counts with its own code.
session: mac-mini-2026-11-08-ternary
source: 05-knowledge/results/ternary_ribbon_metagraph_lines_20261108.md
scripts:
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/ribbon_packing.py
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/packing_intrep.py
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/rankdrop.py
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/ribbon_levels_mod.c (pk3_s2_1e4_levels.txt, pk3_s2_1e5_levels.txt, pk3_s3_5e3_levels.txt, pk2_s2_5e3_levels.txt)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/spec3_small.py (spec3_orders_le_8.txt)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/readiness.py (readiness_n8.txt)
related:
  - THM-4880
  - THM-4881
  - THM-4940
  - HYP-9530
  - THM-4865
  - THM-4744
  - THM-4873
  - HYP-9453
  - HYP-9541
---

# THM-4952 - ternary Chan-Pak on ribbon words

Setting (THM-4880): a span-w ribbon word has letters a_q subset of {1..w}, back arcs q -> q - s (s in a_q), all other
arcs forward; span-2 letters `--`, `1-`, `-2`, `12`. H_x(T) = I(Omega(T), x) = sum_k alpha_k(T) x^k counts packings
of vertex-disjoint odd cycles weighted x per cycle; H_2 = H (OCF).

**(A) Finite-state packing polynomial.** Process the vertices in order. A packing restricted to the prefix is a set
of completed cycles and of vertex-disjoint directed path fragments of unfinished cycles; a fragment's head will be
entered by a back arc from a later vertex, so it lies in the window of the last w vertices; its tail will be left by a
forward arc, which exists unless the later vertex has a back arc to it, so a tail outside the window is free. With the
state = multiset of (head age, tail age or free, parity), the new vertex is skipped, starts a fragment, extends a tail,
extends a head, or does both (merging two fragments, or closing one into a cycle, allowed iff the length is odd,
weight x); heads reaching age w kill the state, tails reaching age w become free; the empty state is initial and
final. Hence I(Omega(T_u), x) = <empty| M_{a_n}(x) ... M_{a_1}(x) |empty> with entries linear in x; reachable states
4, 20, 112 for w = 1, 2, 3.

**(B) Dimension.** For span 2 the functionals u -> I(Omega(T_u v), x) span a Q(x)-space of dimension exactly 14:
the functional matrix of depth <= 4 (341 rows, entries of degree <= 2, so every 15 x 15 minor has degree <= 30) has
rank 14 at 65 integer points (audit; rankdrop.py: 36 points), and the depth-5 matrix (1,365 rows, degree <= 3,
bound 45) at 61 points, so all 15 x 15 minors vanish identically, the depth-4 span is letter-invariant and is the
whole functional space. The rank is 1 at x = 0 and 13 at x = 2 and x = -2 (exact Z-module stabilisation), and a
nonsingular 14 x 14 minor has determinant -x^13 (x - 2)(x + 2), so no other value lowers the rank. The 13 at x = 2
is the irreducible dimension of THM-4880 (G): by OCF the packing and path series agree on every nonempty word, so
their Hankel ranks on nonempty prefixes coincide (13), and the computed rank 13 of the full packing module adds that
the empty-prefix functional (packing value 1, path value 0) is dependent on the others. Integer representations in the echelon or LLL basis of the Z-module (verified on 400 random
words each): span 2 at x = 2, 3, 4, 5: ranks 13, 14, 14, 14 with entries <= 2, 5, 6, 21; span 3 at x = 3: rank 70,
entries <= 4; at x = 2: rank 61 (= THM-4880's path rank 62 minus the empty word), so the drop at the
Hamiltonian-path base is 9 for span 3. In base 3 the span-2 letters have characteristic polynomials x^13 (x - 1),
x^11 (x - 2)(x^2 + x + 1), x^3 (x - 1)(x + 1) q_9(x), x (x - 1)(x + 1) q_11(x) (q_9, q_11 in the note), growth
2, 2.446327, 2.718009 for the constant-letter runs, and no ordered two-letter block has characteristic polynomial
x^13 (x - c) with c != 1 (in base 2 six blocks have c = 3, THM-4881 (B)).

**(C) The tower.** The odd cycles of the tower 1-^n are the odd intervals of length >= 3 of the vertex order, so
packings are sets of disjoint odd intervals and f(n) = f(n-1) + x sum_{j>=1} f(n-2j-1), i.e. f(n) = f(n-1) + f(n-2) +
(x - 1) f(n-3) with f(0) = f(1) = f(2) = 1: 1, 1, 1, x+1, 2x+1, 4x+1, x^2+6x+1, 3x^2+9x+1, 8x^2+12x+1, ... At x = 2
this is tribonacci (THM-4865 (c)); at x = 3 it is 1, 1, 1, 4, 7, 13, 28, 55, 109, 220, 439, 877, 1756, 3511, 7021,
14044 (A122552), and since (x - 2)(x^2 + x + 1) is the characteristic polynomial, H_3(1-^n) = (3 * 2^n + c_n)/7 with
c_n = 4, 1, -5 for n = 0, 1, 2 (mod 3): the tower grows exactly like 2^n in base 3.

**(D) Base-3 census, span 2.** On the rank-14 base-3 representation with the 17-state condensation automaton of
THM-4880 (B), the breadth-first enumeration pruned by H_3 <= 10^4 has 1, 2, 8, 24, 68, 196, 575, 1,778, 5,573, 17,710,
56,013, 172,412, 432,274, 780,552, 1,050,904, 1,171,582, 1,196,745, 1,197,728, 1,197,559 distinct states at lengths
1..19 and 1,197,559 at every later length to 24: a fixed point, so (the pruning being sound by monotonicity) the set of
reduced states of all words with H_3 <= 10^4 is this set, and the value sets are complete for all lengths (the program records cardinalities; the audit verified the equality of the
state sets S_19 = S_20 = S_21 = S_22, and the misses 10 and 31 are in any case excluded for every tournament by
THM-4940 (3)-(4)). Every v = 1 (mod 3) with v <= 10,000 except 10 and 31 is attained (all by length 13), and every
such v except 10, 31, 64 is attained by a strong word. To 10^5: every v = 1 (mod 3), v <= 100,000, except 10 and 31
is attained by length 15 (5,276,750 states at length 15; strong: except also 64); the enumeration was capped at
14,000,000 states from length 17, so no fixed point is claimed at 10^5. Span 3 (rank 70, the 121-state automaton),
H_3 <= 5,000: every v except 10 and 31 by length 10 (2,187,270 states; strong: except also 64), capped at 5,000,000
states at length 11. Pipeline check: the base-2
packing representation reproduces THM-4880's misses 7, 21, 47 below 5,000.

**(E) 64.** For a strong tournament T, H(T) = I(Omega, 2) <= I(Omega, 3) = H_3(T), so H_3(T) = 64 forces H(T) <= 64 <
75 = f(9) and |T| <= 8 (Busch). The packing polynomials of all tournaments of orders <= 8 show that the base-3 values
1 (mod 3) attained at those orders but by no strong tournament are exactly 1 and 64; so 64 has no strong realiser of
any order, while 64 = 4 * 16 = H_3(C3 => T) for any T with H_3(T) = 16. Across bases: at orders <= 8 the attained
values without a strong realiser are {1} in bases 2, 4, 5 and {1, 175}, {1, 120} in bases 6, 7 (audit census; 175 and
120 exceed f(9) and are undecided), and for b >= 8 every attained value below 64 is 1 + ab with a <= 7 in the
intercyclic spectrum (THM-4940 (2)), each such a being the odd-cycle count of a strong intercyclic tournament of order
<= 5; since a strong realiser of a value v has order <= 8 whenever v < 75, 64 is the smallest value other than 1, in
any base, that is attained but has no strong realiser. The same enumeration reproduces THM-4940's small gaps: base 2: 7, 21; base 3: 10, 31;
base 4: 13, 57, 109; base 5: 16, 51, 136, 371 (then only values beyond the order-8 range).

**(F) The return forest is not fed by near-joins.** The near-join of (T, t) has H = h + 2R_t(T) (HYP-9453) and is strong
with T, so the digit 2j of h -> 3h + 2j (j = 1, 2) is near-join-ready at h iff some realiser has a vertex with
R_t = h + j. Over all tournaments of orders <= 8 and all vertices: 97 odd h <= 199 are attained; 50 have R_t = h + 1,
55 have R_t = h + 2, and h = 3, 15, 27, 35, 41, 43, 57, 67, 83, 101, 105, 111, 121, 131, 139, 141, 143, 149, 151, 153,
157, 159, 165, 189, 195, 199 have neither (readiness_n8.txt). Pre-registered P9 predicted no failures for h >= 9: REFUTED.
