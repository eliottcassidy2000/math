# Unconditional first step and the marked-offset obstruction

**Status.** This extends the [moving-guard lemma](MOVING_GUARD.md): the first two-vertex step has an exact nonnegative defect even without its no-start hypothesis, and every later step using the propagated guard a scales the full observer state by 3. Two explicit strong examples realize scalar offsets +2 and +4. They do not realize the stronger marked t-offset required by the mixed-prefix compiler, and no complete ternary bank or spectrum theorem follows. No claim of novelty against the literature is made.

Let X be a tournament with q having unique in-neighbor r. Let F(X,q) add a,b satisfying `q→a→b→q` and `a,b→X−{q}`. Put

    G(X,q,r) = # Hamilton paths of X−q starting at r.

Without any no-start hypothesis,

    H(F(X,q)) = 3H(X) + 2G(X,q,r).                 (1)

The proof uses the same three path classes as in [MOVING_GUARD.md](MOVING_GUARD.md). Paths starting a are in bijection with H(X). Paths starting b and paths starting neither a nor b each have H(X)+G possibilities: after removing a,b, the only possible invalid old adjacency is q immediately followed by r; since r is q's only old in-neighbor and is now its successor, q must be the first old vertex. Every Hamilton path of X−q starting at r gives exactly one such extra path in each of those two classes, and no other extra paths exist.

In F(X,q), a has unique in-neighbor q. Deleting a leaves b a source, so the new G is zero. Therefore each subsequent attachment using the propagated guard triples H. If X is strong, so is F(X,q) and every later attachment.

## The full observer correction

Let p≠q be a fixed marked observer port, and let `X_0=X`, `X_1=X^a_p`, `X_2=X^aa_p` be the two successive singleton queries defined in [MOVING_GUARD.md](MOVING_GUARD.md). Write `G_i=G(X_i,q,r)`. The guard q still has unique in-neighbor r in all three query graphs, and F commutes with the queries. Therefore, for `v=(h,x,y)`,

    v(F(X,q)) = (3h+2G_0, 3x+G_1−G_0, 3y+G_2−G_1).     (2)

After this first step, all three guard defects are zero: deleting the new guard a leaves b as a source in each query graph. Hence every further step using the propagated guard triples the entire state, retaining p as observer port.

If `|X|≥3`, then

    G_0 ≤ G_1 ≤ G_2.                                  (3)

To see this, consider a path counted by G_i, starting at r in `X_i−q`. The next query vertex t beats only its marked port s and loses to every other old vertex. If the path ends at a vertex other than s, append t. If it ends at s, insert t immediately before s; the old predecessor of s beats t. The start remains r: when s=r, a path with at least two vertices cannot end at s, so the append case applies. Deleting t recovers the original path, giving an injection. The size condition is necessary: for the two-vertex tournament r→q with p=r, `(G_0,G_1,G_2)=(1,0,1)`.

For `j∈{1,2}`, scalar `H(F)=3h+2j` means `G_0=j`. But (3) forces the second observer coordinate in (2) to be at least 3x, whereas the mixed-prefix target `3v+j(2,−1,0)` requires 3x−j. Thus **this attachment with the same fixed observer port cannot implement either marked +2 or +4 branch on any source of order at least 3**, although it can achieve the scalar counts. This does not exclude re-marking, a different gadget, or a different source interface.

## Explicit scalar +2 and +4

Bit j of row i records i→j. Both source and output graphs below are strongly connected, and the stated p is retained.

* `X=(8,1,3,6)`, q=2, r=p=3 has `v=(5,3,4)` and `(G_0,G_1,G_2)=(1,1,1)`. The output is `(8,1,19,6,43,15)` with `v=(17,9,12)=3v(X)+(2,0,0)`.
* `X=(20,1,2,7,14)`, q=3, r=p=4 has `v=(9,10,11)` and `(G_0,G_1,G_2)=(2,2,2)`. The output is `(20,1,2,39,14,87,31)` with `v=(31,30,33)=3v(X)+(4,0,0)`.

Iterating the automatically valid new guard gives strong families with full states `3^k(17,9,12)` and `3^k(31,30,33)` for every k≥0. The examples' output counts 17 and 31 were already known; the point here is their exact first-step physical interface and its correction.

## Why the next attachment does not repeat these offsets

For any strong X, every old vertex of F(X,q) has at least two in-neighbors: q has the old in-neighbor r and new b, and every other old vertex is beaten by a and b. Thus the only possible next guards are new a and b. For a, the new defect G is zero. For b, its unique in-neighbor is a, and

    G(F(X,q),b,a) = H(X) − H_start(X,q),

because, after b is deleted, a may start a Hamilton path precisely by going to any initial vertex of an X path other than q. In the examples these alternative defects are 3 and 6, respectively; the two possible guards have defects `(0,3)` and `(0,6)`.

More generally, by Camion's theorem a strong tournament has a directed Hamilton cycle, so each vertex is the start of at least one Hamilton path. Consequently the b defect is at least `|X|−1`. A strong tournament with positive initial G has at least four vertices (the unique strong three-vertex tournament has G=0), so that defect is at least 3. Thus this particular attachment cannot produce consecutive scalar +2/+4 steps from a strong source, regardless of how the observer port is re-marked. It can continue by the 0 branch using a, or by a larger scalar offset using b; other constructions remain open.

## Reproduction

Run `python3 first_step_check.py` and `python3 -O first_step_check.py`. The checker verifies both named examples, their two next guards, the formulas and monotonicity for all 312 qualifying marked guard cases among labelled tournaments of orders 3 and 4, strongness preservation, and the order-two boundary example. Checks remain active under optimized Python. The general identities and bound follow from the stated bijection, injection and Camion theorem, not from this finite census.
