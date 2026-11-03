# Square-sum-supported matched source chains

## Theorem: induced square-sum backward-edge gadget

For every integer k >= 0 there is a finite increasing list of 7+k positive integers, and three strong tournaments on these same named vertices with the following properties. Their first tournament reverses every square-sum edge except two specified edges; the second reverses every square-sum edge except one; and the third reverses every square-sum edge. Every nonsquare-sum pair points upward in all three. Successively reversing the two specified square-sum edges yields Hamilton-path counts

    1 + 32*2^k + j*(2 + 16*2^k),    j = 0,1,2.

In particular, k=0 gives 33 -> 51 -> 69 and k=3 gives 257 -> 387 -> 517. The middle tournament is literally the same in both steps. This is not a free jump rule for every starting count.

## Explicit 7-vertex witness

Take the ordered vertices

    X = (2,39,119,254,303),    a=321, b=322.

The only square sums among distinct vertices are

    2+322=18^2,   39+322=19^2,   119+322=21^2,
    254+322=24^2, 303+322=25^2,
    2+119=11^2,  2+254=16^2.

Every other pair has nonsquare sum. Initially orient every edge upward except the five edges from b to X. Then reverse 2->119 and 2->254, in that order. The three graphs are exactly G(Z_0), G(Z_1), G(Z_2) in S161 CHAIN, where G(X) is the cycle substitution X->a->b->X and Z_j has j of the indicated edges reversed. Their counts are 33, 51, 69, and they are all strongly connected. Only the last graph is the unrestricted square-sum-reversal tournament on this chosen set: the first two deliberately omit respectively two and one otherwise eligible square-sum reversals. This distinction is essential.

## Arbitrary k

Start with the seven labels above multiplied by q^2, with q a sufficiently large positive integer. Insert k labels y between 303q^2 and 321q^2, each of the form

    y = t^2 - 322q^2,   25q < t < sqrt(643) q.

Then y+b=t^2. Choose them greedily so that each y has nonsquare sum with X, a, and previously chosen y values. We show that this always succeeds when q is sufficiently large (depending on k). The indicated t interval has (sqrt(643)-25)q+O(1) candidates. Each fixed earlier vertex v is less than 322q^2, and a prohibited equality y+v=s^2 becomes

    t^2 - s^2 = 322q^2 - v > 0.

Each solution corresponds to a factorization (t-s)(t+s) of this positive integer, at most its divisor count. The bound can be made uniform and elementary: tau(N) <= C N^(1/4), where C is an absolute constant. Indeed, for primes p>=16, a+1<=2^a<=p^(a/4), while for each of the finitely many primes p<16 the ratio (a+1)/p^(a/4) is bounded uniformly in a. Multiplying those finitely many maxima proves the bound. Here N<=322q^2, so each fixed v excludes at most C 322^(1/4) sqrt(q) candidates, with the same constant even when v depends on q and previous choices. At step i (0<=i<k), the 6+i existing vertices exclude at most (6+k)C 322^(1/4) sqrt(q) candidates; additionally exclude the i previously chosen t values themselves. Since k is fixed and the available interval has positive linear size, all k choices exist for sufficiently large q. The map t -> t^2-322q^2 is strictly increasing on this positive interval, and removing previous t values makes all new labels distinct. All old square/nonsquare statuses survive scaling by q^2. Hence the induced square-sum graph consists exactly of b joined to each vertex in X and padding, plus 2q^2--119q^2 and 2q^2--254q^2.

Orient upward apart from the star at b and initially leave both extra square edges upward. Since the inserted vertices are after X, the internal tournament is Z_j ⊕ T_k. The ordered path-cover identities F_(Z_j⊕T_k)(r)=r^k F_Zj(r) and H(G(Y))=H(Y)+F_Y(2), together with S161's F_Zj(r)=r^5+2j r^3 and H(Z_j)=1+2j, yield the displayed formula. The wrapper makes all three strong. No finite search is used to establish the all-k existence argument.

## Marking a general physical edge

For any tournament on a prescribed total order of n vertices and any marked pair at positions i<j, one may retain every edge orientation and order-preservingly assign distinct positive integer labels with a square sum on that pair. Put the endpoints at n and n^2+n+1, whose sum is (n+1)^2; place all earlier vertices in 1,...,n-1, all vertices between the marked pair in the large gap, and all later vertices above the second endpoint. Relabeling changes no path count or path-cover response. Therefore this applies equally to both orientations of the named apex edge in every S161 carry-fingerprint pair, preserving the exact original and wrapped counts and response polynomials. It guarantees only that the marked pair has square-sum endpoints, not that other backward edges are square-sum or other square-sum edges point backward.

## Checks and scope

Run `python3 check.py` and `python3 -O check.py`. It constructs explicit k=0,...,4 instances, verifies every square-sum edge against the stated exact induced pattern, and independently enumerates directed Hamilton paths by subset/end-point dynamic programming. It checks counts 33,51,69 and 257,387,517, and all intermediate examples. The divisor argument proves existence for all k; the finite search and DP are controls. No claim is made about vertex sets {1,...,n}, absence of square-sum edges to omitted integers, the full tournament spectrum, or preservation of the global square-sum rule by the two preliminary states.
