# Terminal-pair carriers and strong 3-cycle substitution

Status: proved elementary identities, finite controls included. They give a strong-gadget construction and exact formulas, **not** a classification of the tournament Hamilton-path spectrum. No claim of novelty is made.

## A carrier with the 7, 28, 126 counts

A carrier on a finite labeled set is an unordered set partition whose every block is either a singleton, or a set with an unordered pair of distinct terminals. Its exponential generating function is

    F(x) = exp(x + x^2 exp(x)/2).

Indeed, singleton blocks have EGF x; a nonsingleton block has two unordered terminals and an arbitrary set of interiors, hence EGF x^2 exp(x)/2. The counts on 3, 4, 5 labels are respectively 7, 28, 126.

Let T be any tournament. For a singleton block S assign weight 1. For a nonsingleton terminal block (S,{a,b}), consider every linear order of S\{a,b} and both choices of orientation of the terminal pair. Give the block the number of these candidate sequences that are directed paths in T. This weight may be zero or exceed one. Summing over its terminal pairs gives precisely the number h_T(S) of directed Hamilton paths of T[S]: each path has uniquely determined unordered endpoints, endpoint orientation and interior ordering.

Let C_k(T) be the number of ordered k-tuples of nonempty vertex-disjoint directed paths whose union is V(T), including paths of length zero. Consequently

    C_k(T)/k! = sum_{carriers with k blocks} product of block weights.

The weighting is essential: the unweighted carrier counts 7,28,126 are not themselves the Hamilton-path counts. The carrier and this weighted identity give a partition-based route from those numbers to the path-cover data controlling exact strong tournament gadgets.

## A cyclic gluing formula

Let A, B, C be nonempty tournaments and form T by retaining their internal arcs and orienting all cross-arcs A→B, B→C, C→A. Set a_k=C_k(A), b_k=C_k(B), c_k=C_k(C), extending by zero outside their natural range. Then

    H(T) = sum_{k>=1} [
        3 a_k b_k c_k
        + a_{k+1} b_k c_k + a_k b_{k+1} c_k + a_k b_k c_{k+1}
        + a_{k+1} b_{k+1} c_k + a_{k+1} b_k c_{k+1}
        + a_k b_{k+1} c_{k+1} ].

Proof: Cut a directed Hamilton path at each transition between blocks. Its maximal runs are nonempty directed paths internal to the blocks; their block labels must follow A,B,C,A,... . Since all three blocks occur, the number of runs is 3k, 3k+1 or 3k+2 with k>=1. In the first case there are three possible starting blocks and all blocks have k runs. In the second, exactly the starting block has k+1 runs. In the third, exactly the first two of the cyclic run pattern have k+1 runs. Conversely, an ordered path cover for each block, combined with any of these run patterns, uniquely reconstructs a Hamilton path, since every cross-run edge is oriented correctly. This is a bijection, including singleton runs; no overcounting occurs.

Moreover T is strongly connected: for two vertices in distinct blocks follow at most two cross-edges in the cyclic direction; two vertices in the same block are connected by a path through the two other nonempty blocks. Thus the identity is an exact operation producing strong tournaments from arbitrary nonempty inputs.

For A a singleton, put h_B=C_1(B), q_B=C_2(B), and similarly for C. Only k=1 contributes and the formula becomes

    H(T) = (h_B+q_B)(h_C+q_C)+2 h_B h_C.

For C also a singleton this simplifies further to H(T)=3h_B+q_B. This shows that the correct input invariant is the full cover profile, even for questions about Hamilton paths alone.

## Explicit fixed-order flexibility

Both of the following 5-vertex tournaments have H=15, yet have different q=C_2. Begin with the transitive orientation i→j for 0<=i<j<=4. Reverse

    B_1: {0,4}, {1,3};                 H=15, C_2=78;
    B_2: {0,3}, {0,4}, {1,4};          H=15, C_2=90.

The edge sets specify reversed arcs, not undirected tournament edges. Inserting B_i into the 3-cycle together with two singleton blocks produces strong 7-vertex tournaments with 123 and 135 Hamilton paths, respectively. Thus identical input order and H do not determine this output: the two-path-cover coordinate supplies genuine extra freedom. These numbers are exact constructive witnesses, with no implication that other specific values, including the excluded 7 and 21, occur.

For a readable independent formula for q, sum over all vertex permutations and let d be the number of backward adjacent arcs. A two-path cover must cut every backward adjacency and make exactly one cut. Hence

    C_2(T) = (n-1) H(T) + #{permutations with exactly one backward adjacency}.

In particular the second term is not generally zero, and one must not replace C_2 by (n-1)H.

## Reproduction and scope

Run `python3 check.py` and `python3 -O check.py`. Explicit condition checks remain active under -O. The program verifies 7,28,126; independently computes all cover profiles through five vertices by permutation/cut enumeration and H by subset/endpoint DP; checks the gluing identity on several exhaustive small triples; and verifies both 7-vertex examples by subset/endpoint DP. Finite computation is a hostile control, not the proof of the general bijection. These constructions do not establish global exact-spectrum coverage, positive density, or a new exclusion theorem.

The companion [transitive-block packet](../cyclic_blocks/NOTE.md) specializes this identity via C_k(T_n)=k! S(n,k), giving an EGF and explicit families.
