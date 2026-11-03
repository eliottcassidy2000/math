# An indefinitely repeating strong 3× return

**Status.** This proves an all-depth, bounded two-vertex, physical marked-source return for the 0 ternary branch on two explicit strong seed families. It does not provide +2 or +4 returns on those descendants, a closed bank, or a solution of the tournament-spectrum problem. The construction resembles the attachment mechanism in [S162](https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session162_mrb_square_20261103/square_transfer/NOTE.md); no literature-priority or independent-novelty claim is made.

## The moving-guard lemma

Let X be a nonempty tournament with a vertex q having a unique in-neighbor r. Suppose **no** Hamilton path of X−q starts at r. Add two vertices a,b with

    q→a→b→q,   and   a,b→X−{q}.

Call the resulting tournament F(X,q), and take a as the new guard. Then

1. `H(F(X,q))=3H(X)`.
2. The new guard a has the unique in-neighbor q, and no Hamilton path of F(X,q)−a starts at q. Thus the hypothesis repeats indefinitely.
3. If X is strong, so is F(X,q).

**Proof.** Partition the Hamilton paths of F by their first vertex. If a is first, b must follow a because a is the sole in-neighbor of b; then the rest is any Hamilton path P of X. If b is first, a must immediately follow q; conversely any path P of X gives a path of F by placing b first and inserting a immediately after q. If neither starts the path, the only possible predecessors force the consecutive block q,a,b; conversely insert a,b immediately after q in any path P of X.

In the latter two cases, deleting a,b leaves only one potentially invalid old adjacency, q→r, because q beats every old vertex except r. If that adjacency occurred, q could not have an old predecessor, since its sole old in-neighbor r is already its successor. Thus q is the first old vertex, and the remaining old path is a Hamilton path of X−q starting at r, contrary to hypothesis. These are three disjoint bijections with Hamilton paths of X, proving (1).

In F, only q beats a. On deleting a, b beats every old vertex and hence is a source, so any Hamilton path starts at b, not q. This proves (2). For (3), X is strong and the cycle q→a→b→q connects the new vertices to it. ∎

## A fixed observer port

For an old marked port p≠q, let X^a_p add a vertex beating only p and mark it, and let X^aa_p repeat this query. Suppose the lemma's guard hypothesis holds not only in X but also in both X^a_p and X^aa_p, using the same q and r. Then every iteration scales the entire observer `v=(h,x,y)` by 3, keeping p as its observer port and moving the guard separately to the first new vertex.

Indeed the two-vertex attachment commutes with either query: the new attachment vertices beat both query vertices, and each query vertex beats only its own designated port. Thus `F(X,q)^a_p=F(X^a_p,q)` and similarly for the second query, up to the obvious labelling. The lemma triples all three queried H counts. Its propagation proof applies independently to all three query graphs after the first step, and p remains distinct from every new guard. Induction proves the full-state identity for every subsequent step.

## Explicit strong seeds

Use adjacency row masks, where bit j of row i means i→j. The observer port is p=9 in both graphs.

    T2=(1008,4089,8179,997,776,912,816,592,128,256,3065,1017,12283,4095),
       q=12, r=2,   (H,H_a,H_aa)=(2565,4779,8829),
       v(T2)=(2565,1107,2025).

    T4=(996,993,994,999,1999,912,800,576,128,256,15343,1023,11263,3071),
       q=10, r=4,   (H,H_a,H_aa)=(2619,4779,8829),
       v(T4)=(2619,1080,2025).

Both are strong. In each of the three query graphs, direct endpoint/subset DP verifies the stated unique in-neighbor of q and that X−q has zero Hamilton paths starting at r. The two graphs arise by two-vertex strong returns from the 12-vertex sources of [TWO_EARS.md](TWO_EARS.md), themselves continuing the 95 and 97 offset branches.

Consequently, for all k≥0, repeatedly apply F and move the guard to a. The results are strong and have exact observer states

    v(F^k(T2)) = 3^k (2565,1107,2025),
    v(F^k(T4)) = 3^k (2619,1080,2025).

For example the initial exact strong count chains are

    31→95→285→855→2565→7695→23085→…,
    31→97→291→873→2619→7857→23571→… .

The same source equalities transport through every fixed common suffix W in `{b,q2}`. Specifically, use the earlier 31-source U and write `N_W=H(U W)`. If W has an even number of b atoms, then the inherited offset identities and linear source transfers give, for every k≥0,

    H(F^k(T2) W) = 3^(k+3) (3N_W+2),
    H(F^k(T4) W) = 3^(k+3) (3N_W+4).

All these literal sources are strong. Apply the moving-guard iterations first, retaining p, and then attach W at p; this statement does not require a guard to survive an arbitrary suffix.

These chains exhibit a closed 0-branch mechanism for each of the two named source families, not 1- or 2-branch closure.

## Reproduction

Run `g++ -std=c++17 -O3 -Wall -Wextra moving_guard_check.cpp -o moving_guard_check && ./moving_guard_check`, and optionally compile with `-DNDEBUG`. The checker directly verifies both seed tournaments are strong, all six guard hypotheses, the three query counts, and the first 3× attachment and propagated guard on every query graph. Its explicit failure tests remain active under `-DNDEBUG`. The all-depth claim is the bijection and induction above; the finite computation verifies the literal seeds, not infinitely many steps.
