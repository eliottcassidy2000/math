# Two physical continuations of the three-cycle upper branch

**PROVED for every k>=1.** The upper child from `CYCLE_CROSSING.md`, with

    M_k=31+33·2^k,

has explicit strong descendants with counts 3M_k+2 and 3M_k+4. The lower route changes the attachment root and adds one vertex. The upper route changes one old core arc and adds two vertices. Every intermediate tournament is strong, has at most seven Hamiltonian cycles, and has fewer paths than its route's final target. The routes retain physical names throughout.

This is a continuation of one specified family. It is not closure for arbitrary three-cycle inputs or an indefinitely reusable grammar.

## 1. The common parent and the two final cores

The parent core, with physical attachment root r=1, is

    B=(14,52,40,34,45,1).

Use the same exterior as before: w->Z->B, r->w->B-r, where Z is the transitive k-block. Keep every existing name and the same exterior vertices.

**Lower child.** Change the attachment root from old vertex1 to old vertex4. Then append vertex6, which beats exactly old vertices {1,3}. The final core and root are

    D=(78,52,104,34,109,65,10),  root4,
    (H(D),F_4(D),e_4(D))=(95,99,7).

**Upper child.** Reverse the old arc0->2, giving C below. Append vertex6 beaten only by old vertex1, then append vertex7 beaten only by that same old vertex1. The core sequence is

    C=(10,52,41,34,45,1),            (H,F_1,e_1)=(25,27,2),
    E=(10,116,41,34,45,1,61),        (H,F_1,e_1)=(49,51,2),
    F=(10,244,41,34,45,1,61,125),    (H,F_1,e_1)=(97,99,2).

The upper root stays at the original physical vertex1. In both routes each newborn is dominated by Z and w. The exact row restrictions verify every claimed induced parent and every unchanged old incidence. All five displayed cores are strong.

The additional parent profile needed for re-rooting is

    (H(B),F_1(B),F_4(B),e_1(B),e_4(B))=(31,33,34,3,3).

All these finite coordinates have independent literal-permutation/cut and subset-recurrence certificates in `verify_reroot.py`.

## 2. Exact control during the change of root

Perform re-rooting before the lower birth, in this order:

1. Reverse w->4 to4->w, so w now has in-neighbors {1,4}.
2. Reverse1->w to w->1, leaving4 as the sole in-neighbor.
3. Add vertex6 with the incidences specified above.

Every state is strong: the core is strong, w reaches the nonempty Z, Z reaches the core, and at least one core vertex reaches w. The sole-entry zero-response mark is temporarily lost in the two-in-neighbor state and is restored after step2. No unchanged-mark assertion is made for that intermediate.

Here is an exact two-in-neighbor counting law. Let S be a nonempty set of core vertices that beat w, with w beating the other core vertices, and retain w->Z->B. Let K(q,u) count ordered spanning two-path covers of B whose first path ends q and whose nonempty second path starts u. Then

    H = H(B)+2^k sum_(q in S) F_q(B)
        - sum_(q,u in S, q!=u) K(q,u).               (1)

A path starting w traverses all of Z and then B, giving H(B). Otherwise the predecessor of w is one q in S. The core vertices split into an ordered cover counted by F_q, and each Z-name chooses the block before the first core path or the block immediately after w. There are 2^k choices. Precisely when the entire Z-block is before w and the second core path starts at another u in S, the required w->u edge is unavailable. These invalid choices are counted by K(q,u), once each. This proves (1) for every k>=1.

For the actual parent core,

    K(1,4)=19,       K(4,1)=3.

Thus the bridge with S={1,4} has count9+67·2^k. A cycle has its unique representation starting w, followed by all Z and a core path ending at a member of S. Therefore its cycle count is

    hc = sum_(q in S) e_q(B).                       (2)

The bridge has exactly six cycles. Its numerical count is not inferred from just the two final root profiles: the cross-cover correction22 is essential.

## 3. The two all-clock routes

Write x=2^k>=2. The lower route has exact path/cycle counts

| State | H | hc |
|---|---:|---:|
| Parent, root1 | 31+33x | 3 |
| Two-in-neighbor bridge | 9+67x | 6 |
| Same core, root4 | 31+34x | 3 |
| Lower final D, root4 | 95+99x=3M_k+2 | 7 |

Both strict intermediate inequalities hold:

    (95+99x)-(9+67x)=86+32x>0,
    (95+99x)-(31+34x)=64+65x>0.

The upper route has

| State | H | hc |
|---|---:|---:|
| Parent B, root1 | 31+33x | 3 |
| Old-arc flip C | 25+27x | 2 |
| First birth E | 49+51x | 2 |
| Upper final F | 97+99x=3M_k+4 | 2 |

Here the corresponding differences from the final target are72+72x and48+48x. The parent is below either target as well. These are nonovershooting routes relative to their **final targets**, not a claim that every intermediate is below the original parent. In particular the re-rooting bridge usually has more paths than the parent.

The parent has k+7 vertices, the lower child k+8, and the upper child k+9. The first parent/children triples are

    (97;293,295), (163;491,493), (295;887,889).

There are at most seven cycles at every stage, including the bridge. The final lower and upper children have different orders; no direct one-arc sibling patch between them is asserted.

## 4. Reproduction and limits

    python3 verify_reroot.py
    python3 -O verify_reroot.py

This is a standalone Python script with no input files or third-party dependencies. It independently checks six rooted core profiles and both cross-cover constants; 42 actual wrappers at k=1,...,6; all36 named physical transitions, including unchanged old incidences at births; and all seven actual Hamiltonian-cycle sets at k=1. Both modes have identical stdout and receipt hashes. The all-order assertions follow from the bijections above, not the finite tests.

The path-cover wrapper and single-root cycle law are inherited from the preceding fan packets. This continuation makes two changes that those initial interfaces had deliberately avoided: one route changes the physical attachment root, and the other changes an old core arc before adding vertices. Those permissions are mathematical operations in the displayed construction, not a general theorem that every output can be rearmed. Longer branching, arbitrary-target realization, and any universal bounded-cycle grammar remain open. No priority claim against the wider literature or unexamined incoming repository work is made.
