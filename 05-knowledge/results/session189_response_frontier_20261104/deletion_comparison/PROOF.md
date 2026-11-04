# Incoming deletion-tower check and exact named-response exclusions

## 1. Result and scope

The exact incoming order-8/9 deletion census is reproducible. Its three reported order-9 unweighted converse-merged deletion-deck collision pairs are **not** full named two-port response collisions, under any bijection of the nine old names and either converse choice on either side. A linear-coefficient multiset already separates every orientation comparison. All named cubic and quartic coefficients are nevertheless retained and independently recounted.

This is a finite exclusion of six merged classes. It is not a reconstruction theorem, does not classify all order-9 named-response pairs, and does not establish a theorem about arbitrary original-module XOR families. The incoming census itself addresses a different observer.

## 2. Exact incoming source

Pinned commit: [e58a09c18bc7536e2388c85df239afd9d10a9c1d](https://github.com/eliottcassidy2000/math/commit/e58a09c18bc7536e2388c85df239afd9d10a9c1d), dated 2026-11-04 16:58:31 UTC.

- [HYP-9455](https://github.com/eliottcassidy2000/math/blob/e58a09c18bc7536e2388c85df239afd9d10a9c1d/05-knowledge/hypotheses/HYP-9455-path-weighted-deletion-tower-reconstruction-orders-8-and-9.md): Git blob bd2b6ba5fa5fc446414c7726ab158451a3d1ee2b.
- [Census code](https://github.com/eliottcassidy2000/math/blob/e58a09c18bc7536e2388c85df239afd9d10a9c1d/04-computation/experiments/deletion_tower_census_20261104.c): Git blob 98953760287aa68650fb91728adea9686ea03e6d.
- [Output](https://github.com/eliottcassidy2000/math/blob/e58a09c18bc7536e2388c85df239afd9d10a9c1d/04-computation/experiments/deletion_tower_census_20261104.out): Git blob 550a9cf686bdf6f130cc5db71fd2bb8e74a285f0.

The local copies preserve these exact Git blobs. The code was compiled without substantive edits and rerun through order 9. Every non-timing output line agrees. This is a source reproduction, not an independent implementation of the whole tournament census or rank calculation.

At the targeted latest-source check on 2026-11-04 around 17:20 UTC, HYP-9455 and the census code had these same blobs. The recent-commit list extended to 92aaef65 at 17:13:28 UTC and contained no newer deletion-tower result. Order 10 remained preregistered. The note and order-10 C paths named in source_metadata.json returned not found at that check; their unavailable content is not used here.

The incoming observer is

    W(T) = sum_v (S_v(T) + E_v(T)) e_[T-v],

where the card class [T-v] is merged by isomorphism and converse. Its total weight is 2H(T). The primitive observer divides W by its coordinate gcd. The incoming face row R is W/(2|Aut T|) for a self-converse node and W/|Aut T| for a converse pair.

The reproduced finite conclusions are:

- Order 8: 6,880 ordinary classes, 176 self-converse classes, 3,528 merged nodes. All primitive W rows and all unweighted merged decks are distinct. Gram rank is 272/272 modulo 2^31-1.
- Order 9: 191,536 ordinary classes, 2,752 self-converse classes, 97,144 merged nodes. All primitive W rows are distinct. There are 97,141 distinct unweighted decks, with exactly three two-node collisions. Gram rank is 3,528/3,528 modulo 2^31-1.

Full rank of R^T R over this finite field certifies full column rank of R over the rationals: the Gram determinant is a nonzero integer modulo the prime, hence nonzero over Q, and a column dependence of R would also be one for its Gram matrix. The converse implication from deficient modular Gram rank is not asserted.

## 3. The observers must not be conflated

W(T)=W(T^op): converse exchanges start and end counts and leaves every converse-merged card class unchanged. In contrast, the named response F_T(a,c), with the fixed edge x->y and a_v=[v->x], c_v=[v->y], retains the old names and the fixed port orientation. Whole converse is not an equality of this named observer in general; reversing all edges exchanges the ports and complements their incidence variables.

Distinct W rows and full column rank do not by themselves prove that equality of F determines W. Such an implication would require a proved extraction of the merged deleted-card classes and their endpoint weights from F. No such extraction is established in this packet. In particular, the shared-name quartics recover symmetrized two-cover capacities on each named vertex deletion, but this alone has not been shown to reconstruct the card classes used by W.

Likewise, equality of the unweighted merged deck is much weaker than equality of F. The following direct computation resolves all three reported examples without relying on an implication between the observers.

## 4. Independent representatives and deck checks

The exact source node pairs are:

1. 20197 and 48562, both H=633 and automorphism order 3, neither self-converse.
2. 86814 and 89262, H=2919 and H=2917 respectively, automorphism order 1, both self-converse.
3. 23898 and 23984, both H=645 and automorphism order 1, neither self-converse.

prepare.py extracts their adjacency codes from the exact census export. representatives.json preserves the six source lines, codes, physical adjacency rows, and both converse choices. There are twelve response graphs; duplicate isomorphism classes for the self-converse pair are intentionally kept.

literal_cards.cpp independently canonicalizes each whole graph by enumerating all 9! permutations, and each of its deleted cards by enumerating all 8! permutations. It separately minimizes over converse by complementing the maximum adjacency code. Thus this check does not reuse the incoming color-refinement canonicalizer. It also computes Hamilton endpoints with a predecessor subset recurrence. For every reported pair it verifies:

- distinct isomorphism-and-converse-merged whole classes;
- identical unweighted merged deleted-card multisets;
- distinct primitive endpoint-weighted merged decks;
- the incoming H, automorphism order, and self-converse status.

The 108 deleted-card records retain the actual named deleted vertex and its two endpoint counts. This confirms the reported finite examples independently of the source's merged-node labels.

## 5. Full named-response exclusion under all relabelings

For each of the twelve graphs, bank.cpp evaluates F(a,c) at every a,c of support size at most two, then applies Boolean Mobius inversion. A Hamilton path uses at most two edges incident to either port, so F has degree at most two in each variable family. These 46^2=2,116 coefficients per graph are the complete response polynomial, including all zeros.

The independent predecessor.cpp implementation recounts all 25,392 actual values, and verify.py independently repeats the inversion and compares every coefficient. There are 7,776 named cubic and 15,552 named quartic coefficients, all retained in bank.json. The independent implementation uses predecessor rather than successor subset transitions.

For each pair and all four converse choices, compare the multiset

    L_c(T) = multiset of [c_v]F_T, over the nine old names v.

Any old-name bijection permutes this multiset. If it differs, no bijection can make the full polynomials equal. Every one of the twelve orientation comparisons has different L_c. For example, for the canonical orientations of the H=633 pair the multisets are

    784, 784, 784, 852, 996, 1074, 1074, 1074, 8688
    756, 818, 818, 818, 1026, 1026, 1026, 1134, 8688.

For the canonical orientations of the H=645 pair they are

    796, 814, 838, 842, 1024, 1034, 1048, 1098, 8784
    802, 808, 830, 850, 1008, 1050, 1070, 1076, 8784.

The remaining pair is already separated by F(0,0)=H, since 2919 differs from 2917. The receipt retains all twelve linear multisets and all differing coefficient-histogram types. Histograms are used only as sufficient inequality witnesses; histogram equality would not be treated as a full named-response collision.

## 6. Newer reversal-cube preregistration

The later [92aaef65 commit](https://github.com/eliottcassidy2000/math/commit/92aaef65f87897b121539028e4fb50e3c957956f) creates a separate preregistration at 05-knowledge/results/reversal_cube_mod4_anchor_20261104/PREREGISTRATION.md. It proposes a mod-4 separator for a singleton child or an unchanged/conversed transitive child, using a named source or sink and a compressed remainder. THM-4864 and HYP-9456 are RESERVED there, and the Outcomes section is empty in the inspected commit.

If proved with the stated ambient and source-relative hypotheses, its unchanged-transitive-child branch could combine with Session189 distance reconstruction to remove an additional residual branch without the odd-child ordered-witness restriction. This is potential relevance, not an imported theorem. The current packet neither validates the proposed weighted model nor uses the preregistered claims in any proof.
