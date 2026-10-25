# Cartesian-tree fibres: exact query equivalence is not full object equivalence

Session56, 25 October 2026. PROVED elementary statements, with finite controls. Classical mechanisms are rederived; no external novelty claim or tournament counting-time lower bound is made.

## 1. An exact quotient for a declared observer

For an array a[0..n-1], use the stable total key (a[i],i): equal values choose the leftmost position. Its Cartesian tree has the minimum key at its root and recursively represents the left and right intervals. Node labels are fixed by inorder position, so an ordered shape determines the index labels.

For every interval [i,j], its range-minimum index equals LCA(i,j) in this tree. Proof: the first recursive root separating the endpoints lies inside the interval and precedes every other member in heap order. If both endpoints lie on one side, recurse on that side. Conversely, all interval-minimum answers reconstruct the root from [0,n-1] and then each subtree recursively. Hence two arrays have exactly the same stable RMQ answers iff their Cartesian trees agree.

This is a sufficient statistic for the DECLARED query family, not for arbitrary questions about a. Arrays (1,0,2) and (2,0,1) have identical RMQ signatures and trees, but the comparison of positions0 and2 reverses. Arrays (0,0,0) and (1,2,3) also give the same leftmost RMQ answers while equality and magnitudes differ. Arithmetic root status, carries, mod-8 guards, and cross-branch relations require their own retained data.

Primary classical reference: Bender and Farach-Colton, *The LCA Problem Revisited*, LATIN2000, pp.88–94, Sections2–5, https://cadmo.ethz.ch/education/lectures/HS18/SAADS/papers/lca.pdf . The first mirror timed out; the ETH-hosted primary PDF was read. Fischer's *Optimal Succinctness for Range Minimum Queries*, https://ae.iti.kit.edu/download/fischer10optimal.pdf , gives the succinct indexing context. No full article text is republished.

## 2. Fibres carry a nonuniform hook weight

Let T be an ordered binary tree with n nodes, and let s_v be the number of nodes in the subtree rooted at v. Exactly

    f(T)=n! / product_v s_v

permutations of0,...,n-1 have Cartesian shape T.

Proof by induction. The root receives the unique smallest value. If its two subtree sizes are l,r, choose which l of the remaining n-1 values belong on the left, then choose valid relative orders in each subtree. Thus f(T)=binom(n-1,l)f(L)f(R). Substitution gives (n-1)!/(product_L s_v product_R s_v)=n!/(n product_L s_v product_R s_v), as required.

Every ordered binary tree occurs. There are C_n=binom(2n,n)/(n+1) such shapes, but their fibres are not equal. Summing over shapes gives

    sum_(|T|=n) product_v (1/s_v)=1.

Equivalently, if W_n is that sum, W_0=1 and W_n=(1/n)sum_(l+r=n-1)W_l W_r=1. The ordinary generating function of weighted shape mass is1/(1-z), whereas the unweighted Catalan generating function satisfies C(z)=1+z C(z)^2. This is an exact example of the repository's support-versus-multiplicity distinction: count shapes for query types; retain fibre weights to count source permutations.

The same distribution appears in the usual random-binary-search-tree shape construction: conditioning on the root rank splits the remaining labels uniformly into the same two sizes. This is a distributional equality of ordered tree shapes, not an identification of arbitrary tournament or arithmetic objects with a random tree.

## 3. A quantitative obstruction to tree-only tournament encoding

There are2^(n(n-1)/2) labelled tournaments, and each isomorphism class has at most n! labelled copies. Therefore the number t_n of classes satisfies

    t_n >= 2^(n(n-1)/2)/n!.

Suppose an injective encoding of all unlabelled tournament classes uses one ordered binary tree with EXACTLY n nodes, inorder positions as its only labels, and a fixed-width b-bit sidecar. Then

    C_n 2^b >= t_n,
    b >= max(0, binom(n,2) - floor(log2(C_n n!))).

The inequality follows by counting encoding states; no assumption about an efficient canonicalization algorithm is needed. Since C_n<=4^n, it forces b>=binom(n,2)-log2(n!)-2n asymptotically, up to integer rounding. A tree plus a small index cannot encode the full quadratic information of all tournaments merely by accelerating its LCA queries.

If extra labels, values, quotient tournaments, multiple trees, or a different node count are allowed, their capacity must be counted separately. This bound is NOT a lower bound on computing the scalar sequence A000568, and says nothing against exact decomposition methods on structured families. Prime-node quotient data is precisely one possible sidecar.

For calibration only, the published OEIS A000568 values give t_7=456 while C_7=429, so tree-only capacity already fails at n=7. The cited values for n<=6 are at most C_n (equality at n=0,1). Thus7 is the first such cardinality failure using that table; no new tournament census is claimed. Source https://oeis.org/A000568 , offset0. The difference27 is a numerical capacity deficit here, not a proved connection to G27, Collatz27, or a finite field.

## 4. What transfers to current work

- RMQ/LCA gives a model of a lossless quotient only after its query family is specified.
- The hook formula shows exactly why keeping only support can lose source multiplicity.
- A strong-module tree may compress a tournament, but prime quotients and child ports carry information that generic LCA alone does not return.
- The capacity lower bound rules out a proposed all-tournament tree-only classifier; it does not rule out faster enumeration or a smaller sufficient statistic for some narrower query family.

## Verification

Run `python3 -B verify.py` and `python3 -O -B verify.py`. All checks use explicit exceptions. Every permutation through n=8 is grouped by Cartesian shape and compared with the hook formula; every resulting shape count equals the Catalan number and weights sum exactly to1. Ternary arrays through length6 test ties and RMQ-equivalence. Capacity ceilings are checked using exact integer arithmetic, not floating logs. These finite controls supplement the proofs.
