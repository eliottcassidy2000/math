# A universal root-deletion SCC divisibility obstruction

Status: PROVED by the endpoint-block partition below; independent audit requested. This is a necessary arithmetic condition for every rooted tournament, not a sufficient state-membership test or scalar-classification theorem.

## 1. Statement

Let (A,a) be any finite nonempty rooted tournament. Write the strong components of A−a in their unique linear order as B_1,...,B_m, so all edges between distinct components point forward. A singleton is permitted as a component, and m=0 when A is the singleton root. Put

    H_i=H(B_i),
    Q_i=sum_(S⊆V(B_i)) H(B_i[S]) H(B_i[V(B_i)−S]),
    g_i=gcd(H_i,Q_i).

As usual H(empty)=1. The quantity Q_i is the ordered two-bin spanning path-cover count with empty bins allowed. Define

    G=product_i g_i,
    L=lcm { product_(i∈S) g_i : S⊆{1,...,m}, |S|≤2 },
    K=G/L.

The empty subset contributes1; thus these definitions also cover m=0,1. Every candidate in the lcm divides G, so K is a positive integer. Then

    K divides H(A).                                  (1)

Equivalently, for any prime p let a_i=v_p(g_i), and let a_(1),a_(2) be the two largest entries, padding with zeros when necessary. Then

    v_p(H(A)) ≥ sum_i a_i − a_(1) − a_(2).           (2)

The right side is nonnegative under the stated padding convention. The same theorem holds for any ordered decomposition into nonempty tournament blocks, without requiring them to be strong. Using SCCs makes it available canonically for every rooted tournament.

## 2. Proof

A tournament's SCC condensation is an acyclic tournament, hence transitive; therefore its strong components have the asserted linear order. This standard elementary SCC fact is the only reason the statement applies to every rooted tournament.

Delete a from a Hamilton path of A. Its prefix P before a and suffix R after a are ordered paths partitioning A−a. P must be empty or end at an inneighbor of a; R must be empty or start at an outneighbor of a. Conversely each such pair produces exactly one Hamilton path by P,a,R.

If R is nonempty, specify the block i of its first vertex; if P is nonempty, specify the block j of its last vertex. These choices partition the path pairs into disjoint cases. They identify a set C of at most two distinct boundary blocks, on which the endpoint conditions may impose arbitrary local restrictions.

Every nonboundary block has an unrestricted local count H_k or Q_k, or the entire case is impossible. Specifically, when i<j, blocks before i lie wholly in P, blocks after j wholly in R, and blocks strictly between i and j split freely into the two paths. When i=j, all other blocks lie wholly in one of the paths. When j<i, blocks through j lie wholly in P and blocks from i lie wholly in R; any intervening nonempty block makes the case impossible. If P or R is empty, there is at most one boundary block and all others belong wholly to the other path. If both are empty, m=0 and H(A)=1.

Once local path pairs are selected, their nonempty pieces concatenate uniquely in block order, because all cross-block edges point forward. Thus a nonzero case has an integer local boundary factor multiplied by a product whose factor for every k∉C is H_k or Q_k. Its count is divisible by

    product_(k∉C) g_k = G / product_(k∈C) g_k.

Since the boundary product divides L, every such count is divisible by K=G/L. Summing the disjoint cases proves(1). Taking p-adic valuations in the lcm gives the sum of the two largest a_i, because its indexing sets have cardinality at most two. This proves(2). In particular repeated boundary indices count only once; the proof never subtracts one block's valuation twice.

## 3. Consequences and limits

For a triangle block H=3,Q=12,g=3; for a singleton H=1,Q=2,g=1. The bound therefore reduces exactly to3^max(t−2,0) for t triangle blocks, recovering DIVISIBILITY.md. Its witnesses show this universal bound can be attained at every t in that subclass.

If H(A)=q is prime, then for every prime p≠q the sum of all but the two largest valuations v_p(g_i) must be0. For p=q that sum must be at most1. In particular at most two SCCs can have g_i divisible by a prime p≠q, and at most three SCCs can have g_i divisible by q; in the latter case the smallest retained q-valuation must be1. These are necessary restrictions on the deletion decomposition for every choice of root.

No condition Q_i divisible by H_i is assumed. The gcd is essential. The obstruction may be trivial when all but at most two g_i equal1, for example when the deletion graph is itself strong. It gives no converse: passing the divisibility test does not establish physical realizability or a desired Hamilton count.
