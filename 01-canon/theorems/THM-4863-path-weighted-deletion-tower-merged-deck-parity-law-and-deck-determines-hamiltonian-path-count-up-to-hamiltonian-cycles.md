---
id: THM-4863
title: "The path-weighted deletion tower and the deck. (W) THM-796's face row of a converse-merged tournament node is its Hamiltonian-endpoint-weighted converse-merged deck W(T) = sum_x (h_s(x) + h_e(x)) [T - x], divided by 2|Aut T| (self-converse) or |Aut T| (converse pair). (T) The tower is injective at orders 4-10 (4,870,920 nodes at order 10) and has full column rank at orders 3-9. (P) The unweighted converse-merged deck collides at orders 3, 5, 7, 9 and is injective at orders 4, 6, 8, 10; every collision is classified (source lifts of even-order converse same-deck pairs, plus one non-lift class per odd order), and McKay's same-deck census through order 10 is reproduced (even-order pairs converse, odd-order pairs self-converse; largely published: McKay (f), Boudabbous-Ille 2018, Pouzet 1979 / Dumont via arXiv:1910.13914; new as found: the converse relation at orders 6 and 8). (H) If two tournaments have the same converse-merged deck then H(T) = H(T') at even order and H(T) - H(T') = 2(hc(T) - hc(T')) at odd order. (L) Source lift. (J) Join rule for W. (Phi) The line fraction; its order-4 identity is PROVED, and naming the order-3 nodes 11 and 13 gives {11/13, 1, 13/11} (the PRIMEGAME reading is a DICTIONARY). (B) CITED (THM-851 (14)): line classes number (2^m + 2^h)/4."
status: >
  PROVED: (W), (H), (L), (J), the order-4 line-fraction identity and the basic properties of Phi. (H) is not found stated in the
  literature; it is a short corollary of THM-002 and the Pouzet 1977 (I-2.2(a)) / Tutte / Kocay cycle-cover lemma; its determinant
  analogue is Pouzet 1979 Equality (A).
  FINITE-EXACT: (T), (P), and Phi injective at orders 6-8 (exploratory).
  PRE-REGISTERED: HYP-9455, with P1-P3 and P6-P7 confirmed and P4 refuted at order 8.
  CITED:
  - THM-796 (face matrix; injective rows through order 7);
  - THM-781 (fibres HP/Aut);
  - THM-002 (odd-cycle formula; Grinberg-Stanley arXiv:2307.05569);
  - Kelly's lemma;
  - McKay 2022 (Australas. J. Combin. 83; arXiv:2102.01942), the complete census of non-reconstructible tournaments
    through order 13;
  - Stockmeyer 1977 / 2021 (arXiv:2110.06158);
  - session129 (codec, the block assignment);
  - session29 (the order-5 same-deck pair with H 15 / 13);
  - THM-4854 / S71 (blocks, reciprocal pair);
  - THM-851 (13)-(14) ((B), the line-class count), THM-801 section 4 (the face-row reading), session130 (fibre sizes;
    the line counts at order 5);
  - Pouzet 1977 and 1979; Bankoussou-mabiala et al. arXiv:1910.13914; Boudabbous-Ille 2018; Harary-Palmer 1967;
  - the codex-20260821 deck reflection (H and hc not affine in the deck at n = 7; H - 2hc is).
  DICTIONARY / NUMEROLOGY: every PRIMEGAME label, as typed in the note.
  AUDIT: PASS WITH CORRECTIONS, all applied (MISTAKE-704); report 04-computation/experiments/primegame_deletion_tower_20261104_audit_S72/AUDIT.md.
session: opus-2026-11-04-S72
source: 05-knowledge/results/primegame_deletion_tower_path_weighted_decks_20261104.md
related:
  - THM-796
  - THM-002
  - THM-4854
  - THM-4852
  - HYP-9455
verification:
  - 04-computation/experiments/deletion_tower_census_20261104.c (+ .out)
  - 04-computation/experiments/deletion_tower_order10_20261104.c (+ .out)
  - 04-computation/experiments/deletion_tower_line_fraction_20261104.c (+ .out)
  - 04-computation/experiments/primegame_deletion_tower_analysis_20261104.py (+ .out)
  - 04-computation/experiments/deck_parity_source_lift_check_20261104.py (+ .out)
---

# THM-4863 - the path-weighted deletion tower and the deck

**Notation.**
- h_s(x) and h_e(x) count the Hamiltonian paths of a tournament T that start and end at x.
- [X] is the converse-merged class (merged node) of X.
- W(T) = sum over x of (h_s(x) + h_e(x)) e_[T - x].
- hc(T) counts directed Hamiltonian cycles.
- The deck of T is the multiset of [T - x] (merged deck) or of the classes of T - x (ordinary deck).

## (W) Face row = path-weighted deck

THM-796's face matrix is R_n(u, v) = #{t in F_n(u) : d t in F_(n-1)(v)}. It satisfies:
- R_n(u, .) = W(T) / (2|Aut T|) for a self-converse node u = {T};
- R_n(u, .) = W(T) / |Aut T| for a converse pair u = {T, T^op}.

**Proof.**
1. The fibre is HP(T)/Aut(T) (THM-781), and Aut acts freely on Hamiltonian paths.
2. Endpoint deletion of a path P gives a tiling of T - end(P).
3. The fibre of T^op contributes the start-weighted deck, since h_e^(T^op) = h_s^T and [T^op - x] = [T - x].
4. For self-converse T an anti-automorphism equates the start-weighted and end-weighted decks.

**Check.** The matrix agrees with the recomputation from every fixed-path tiling, through both faces, at orders 3-7.

## (T) Injectivity and rank

The primitive rows of W are pairwise distinct at every order 4 <= n <= 10:

| order | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| distinct rows | 3 | 10 | 34 | 272 | 3,528 | 97,144 | 4,870,920 |

At order 3 they are not: TT_3 and C_3 both lie over the single order-2 node.

R_n has rank |M_(n-1)| = 1, 2, 3, 10, 34, 272, 3528 at n = 3, ..., 9. This is certified by the rank of R^T R modulo
2^31 - 1.

## (P) The parity law of the converse-merged deck, and the classification

Collision classes of the unweighted converse-merged deck:

| order | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| collision classes | 1 | 0 | 2 | 0 | 5 | 0 | 3 | 0 |

**Same ordinary deck.** Pairs with the same ordinary deck number 1, 1, 1, 4, 0, 2, 1, 1 at orders 3-10, which is
McKay's census.
- At orders 4, 6, 8 and 10 each pair consists of mutual converses.
- At orders 3, 5 and 9 each pair consists of two self-converse tournaments, with:
  - H = 1 / 3 and hc = 0 / 1 at order 3;
  - H = 15 / 13 and hc = 3 / 2 at order 5;
  - H = 2919 / 2917 and hc = 211 / 210 at order 9.

**Classification of the merged-deck collisions.** Every merged-deck collision whose members have no Hamiltonian cycle
is a source (or sink) lift {1 -> S, 1 -> S^op} of an order-(n - 1) converse same-deck pair. There are 1, 4 and 2 such
lifts at orders 5, 7 and 9, checked with independent brute-force canonical forms. The one remaining collision class at
each odd order is:
- at order 3, {TT_3, C_3};
- at order 5, the self-converse pair;
- at order 7, a strong collision with H 143 / 141 and hc 15 / 14;
- at order 9, Stockmeyer's self-converse pair.

## (L) Source lift

If S and S^op are non-isomorphic converses with the same ordinary deck, then 1 -> S and 1 -> S^op are distinct merged
nodes with the same converse-merged deck.

**Proof.**
- **Same deck.** The cards are [S] = [S^op] and the multiset {[1 -> (S - x)]}, which the ordinary deck of S determines.
- **Distinct nodes.** Removing the unique source distinguishes 1 -> S from 1 -> S^op. A converse would turn the
  singleton source component into a singleton sink component; the two component sequences could agree only if S were
  transitive, which is excluded because transitive tournaments are self-converse.
- **Corollary.** This forces a merged-deck collision at order 11 (from the order-10 converse pair), and at every order
  m + 1 for which order m has such a pair. The even-order injectivity of (P) is observed, not explained.

## (H) The deck determines H up to Hamiltonian cycles

If T and T' have the same converse-merged deck (in particular, the same ordinary deck), then:
- H(T) = H(T') if n is even;
- H(T) - H(T') = 2(hc(T) - hc(T')) if n is odd.

**Proof.**
1. By THM-002, H(T) = sum over collections of pairwise disjoint odd directed cycles of 2^(number of cycles).
2. **Non-spanning collections.** They are counted by sum over [S] of g(S) s([S], T), with g converse-invariant. Each
   s([S], T) with |S| < n is (1/(n - |S|)) sum over x of s([S], T - x) (Kelly), a function of the merged deck.
3. **Spanning collections with c >= 2 cycles.** Their number is (1/c) sum over [A], [B] of
   hc(A) g_(c-1)(B) N([A], [B]). Here N([A], [B]) counts the partitions V = U ⊔ W with T[U] in [A] and T[W] in [B]. It
   equals s([A], T) s([B], T) minus the overlapping pairs, which live in proper subsets and are counted by Kelly.
4. **What remains** is the single spanning cycle, weight 2 hc(T). It exists only when n is odd.

QED. The identity is checked on every same-deck pair and merged-deck collision through order 10.

The total weight of W is 2H. At odd order, therefore, W carries exactly the coordinate the deck misses, and its
endpoint distribution separates the collisions with equal H.

## (J) Join rule

For T = A => B (every vertex of A beats every vertex of B):

    W(T) = H(B) sum over x in A of h_s^A(x) e_[(A - x) => B]  +  H(A) sum over y in B of h_e^B(y) e_[A => (B - y)].

A Hamiltonian path covers A and then B. The case |A| = 1 is the source extension. The rule is checked on all 55 joins
with |A| + |B| <= 6. It accounts for every Fibonacci-valued run of consecutive values in a face row at orders 6-8:
S_4 => TT_k, S_4 => C_3, S_4 => C_3 => 1 and S_4 => 1 => C_3.

## (Phi) The line fraction and PRIMEGAME's first block

**Definition.** Phi(u) = sum over t in F_n(u) of (e_[kappa d t] - e_[d t]), the fraction prod p_[kappa d t] / p_[d t]
over the order-(n - 1) nodes read as primes.

**Basic properties.** Phi sums to 0 over all nodes, and both endpoint faces give the same Phi. This follows from
d kappa = kappa d, d_L sigma = sigma d_H (THM-796 (2.5)), and the invariance of the fibres under sigma.

**Order 4.** Phi(TT_4) = e_(C_3) - e_(TT_3), Phi(pair) = 0, and Phi(S_4) = e_(TT_3) - e_(C_3). Name TT_3 and C_3 by 11
and 13, in either order. Then {Phi(TT_4), Phi(S_4)} = {11/13, 13/11}, PRIMEGAME's first block, and Phi(pair) = 1.

**Beyond order 4 (FINITE-EXACT, exploratory).** Phi takes 9 values on the 10 nodes of order 5, and is injective at
orders 6, 7 and 8.

## (B) Line classes (CITED: THM-851 (14))

The orbits of <kappa, sigma> on X_n number (2^m + 2^h)/4, with m = binom(n-1, 2) and h = floor((n-1)^2/4). The values
are 1, 3, 20, 272, 8320, ... for n = 3, 4, 5, 6, 7, ... This follows by Burnside, since kappa and sigma kappa fix no
tiling (THM-796 section 4). The owner's "1 node at n = 3, 3 nodes at n = 4" are exact as line classes.

## Typing of the PRIMEGAME readings (see the note)

- The first block as the order-3 pair over the single order-2 node, and as the order-4 line fractions: DICTIONARY with
  an exact identity.
- Blocks as order-4 coordinates: session129's codec, CITED.
- Registers as order-5 coordinates, with rank R_6 = 10: DICTIONARY.
- PRIMEGAME fractions as moves between order-6 node codes: NUMEROLOGY (1 incidence; random median 1; P = 0.58).
- Fibonacci values: NUMEROLOGY, with (J) as their only structure.
