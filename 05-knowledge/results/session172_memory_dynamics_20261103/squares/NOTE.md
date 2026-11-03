# Exact quadratic Hamilton memory: graph read heads and square-sum arrivals

**Status: proved.** The packet gives an exact one-vertex contextual response theorem, its sharp linear dimension, and an undirected equal-profile read-head example. The square-sum numerical examples are inherited mechanisms, independently verified here. There is no claim of a new square-sum Hamiltonicity range, a polynomial-time counting algorithm, or an arithmetic realization of the equal-profile seeds.

## 1. Inheritance and the precise new target

Session 32's `session32-square-sum/README.md`, §3, already proves the degree-two arrival decomposition into endpoint extensions and joined spanning two-path covers, including the compulsory ears at 23 and 25. We do not claim that deletion argument or those examples as new. S171's `framework/NOTE.md`, §5, gives a general square-zero visited-set boundary response; S171's `information/NOTE.md` supplies tournament seeds with identical complete profiles and different one-vertex contextual counts.

Here the extension is a complete characterization of **all undirected one-vertex contexts**, with a sharp response-rank theorem. A second concrete result realizes the hidden-profile/read-head phenomenon in ordinary undirected graphs. The old square-sum arrivals are a faithful arithmetic specialization, not merely a visual analogy.

## 2. The exact response polynomial

Fix a simple undirected graph G on labelled vertex set [n], n≥1. Count oriented Hamilton paths: opposite vertex orders count separately, except the one-vertex path counts once. Let H(G) be this count. Put:

- E_i(G): the number of oriented Hamilton paths of G ending at i.
- K_ij(G), i≠j: the number of ordered pairs (P,Q) of nonempty oriented paths whose vertex sets partition [n], with P ending at i and Q starting at j. Both paths must use edges of G. Singleton paths are allowed.

Put K_ii=0. Reversing both paths and exchanging them gives K_ij=K_ji.

For S⊆[n], adjoin one new vertex a with exactly S as its neighborhood, and write R_G(S)=H(G+a).

**Theorem 1 (quadratic contextual memory).**

    R_G(S) = 2 Σ_(i∈S) E_i(G) + 2 Σ_({i,j}⊆S) K_ij(G).                 (1)

Consequently two graphs on the same labelled set have identical responses to every such one-vertex context if and only if their E and K records agree. The singleton and pair contexts already recover the complete record:

    E_i = R_G({i})/2,
    K_ij = [R_G({i,j})−R_G({i})−R_G({j})]/2.                         (2)

**Proof.** If a is an endpoint, deleting it leaves a Hamilton path ending at one of its neighbors, or the reverse of such a path. This gives the first term. If a is interior, deleting it gives a unique ordered pair (P,Q) as defined above, with joining ports in S. Conversely each such pair gives a unique Hamilton path by inserting a. Summing over ordered distinct ports gives the second term. Formula (2) recovers every coefficient, proving the equivalence. ∎

Thus an exponentially indexed context table has an exact quadratic description. In Boolean variables x_i indicating membership in S, it is the multilinear polynomial

    R_G(x)=2 Σ_i E_i x_i + 2 Σ_(i<j) K_ij x_i x_j.

An edge-influence interpretation is useful. If ij is absent, then

    H(G+ij)−H(G)=2K_ij(G).

If ij is present, then

    H(G)−H(G−ij)=2K_ij(G).

In either case deleting that edge from a counted Hamilton path gives precisely the two pieces in K_ij. The pair never uses ij internally, since i and j belong to different pieces, so its count is unchanged by toggling ij.

## 3. Sharp linear response dimension

Let M_n be the matrix with rows indexed by all labelled simple n-vertex graphs and columns indexed by all subsets S⊆[n], with entry R_G(S). Work over Q (or R).

**Theorem 2 (exact rank).** For n≥3,

    rank M_n = n + binom(n,2) = n(n+1)/2.                            (3)

Hence this is the minimum dimension of a **linear** state record that preserves every one-vertex Hamilton response on the vector space formally spanned by graph states. The singleton/pair columns form a basis of its column space.

**Proof.** Theorem 1 supplies the upper bound and makes the coefficient-to-evaluation map injective. It remains to span all coefficient directions.

For three distinct labels i,j,k, take the graph consisting of isolated i and a spanning simple path on [n]\{i} with endpoints j,k. It has exactly n−2 edges and two components. A spanning two-path cover must use every edge, so its only nonzero unordered K coordinates are

    K_ij=K_ik=1,

and E=0. Thus graph states realize every star wedge e_ij+e_ik in the quadratic coordinates. For n≥4, choose j,k,l distinct from i. The identity

    2e_ij=(e_ij+e_ik)+(e_ij+e_il)−(e_ik+e_il)

spans every individual quadratic coordinate. For n=3, the three wedges on the triangle form an invertible matrix of determinant ±2. The entire quadratic coordinate space therefore lies in the span of actual graph records with E=0.

Now take a spanning path graph with endpoints i,j. Its endpoint vector is e_i+e_j. Subtract a linear combination of the already-spanned quadratic records to remove its K component. For n≥3, the vectors e_i+e_j span all endpoint coordinates, using the same three-term identity on distinct i,j,k. This supplies n additional independent directions and proves (3). ∎

The coefficients in these spanning arguments are signed rationals, not graph-building operations. The theorem is an exact rank statement, not a construction of a graph with an arbitrary prescribed nonnegative response. It does not forbid nonlinear encoding of a whole record in one enormous integer, and does not say counting the coefficients is easy. The characteristic-zero hypothesis matters to the displayed divisions by two.

## 4. An actual undirected profile-blind memory bit

Use labels 0,...,4 and these edge sets:

    A: 01,02,03,04,12,13,14;
    B: 01,02,03,04,12,13,23.

A is K₂ joined to an independent set of size three. B is K₄ with one pendant vertex. Their sorted degree sequences are respectively (2,2,2,4,4) and (1,3,3,3,4), proving nonisomorphism.

Define P_j as the number of partitions of all vertices into j nonempty oriented paths, where the components are unordered but each nonsingleton path retains its orientation. Both graphs have

    P(z)=12z+60z²+54z³+14z⁴+z⁵.

For u labelled path slots, including empty slots, both therefore have

    F(u)=Σ_j P_j(u)_j=u²(u+1)²(u+2).

Despite this complete profile equality, their endpoint vectors are

    E(A)=(0,0,4,4,4),
    E(B)=(0,2,2,2,6).

Attach a single new leaf to label 1, using the same attachment rule for both states. Then

    H(A+a)=0,       H(B+a)=4.

So even the complete unpointed path-cover profile forgets a bit that an authenticated one-vertex context reads. In contrast to the S171 tournament seed's 17/19 count readout, this undirected read head distinguishes **nonexistence from existence** of a Hamilton path. The bulk graphs themselves each have twelve oriented Hamilton paths.

The fixed label is part of the interface. No claim is made that these graphs occur as square-sum graphs, induced square-sum subgraphs, or arithmetic blocks. Their role is a rigorous ordinary-graph realization of the same contextual-memory phenomenon.

## 5. The square-sum specialization and inherited resurrection

For the square-sum graph G_n on {1,...,n}, the new vertex n+1 has actual neighborhood

    S_(n+1)={s²−(n+1): n+2≤s²≤2n+1}.

Thus (1) yields the exact arithmetic recurrence

    H(G_(n+1))=2 Σ_(i∈S_(n+1)) E_i(G_n)
                +2 Σ_({i,j}⊆S_(n+1)) K_ij(G_n).                    (4)

This applies at every arrival degree, not only degree two. Crucially, it is not a scalar recurrence in H(G_n): the coefficients are contextual boundary data.

Independent DFS and subset-endpoint dynamic programming give the following oriented counts:

    n       15  16  17  18  22  23  24  25
    H(G_n)   2   2   2   0   0   6   0  20.

At 23 the neighbor set is {2,13}; at 25 it is {11,24}. Independent enumeration of disjoint path banks gives

    K_(2,13)(G_22)=3,       K_(11,24)(G_24)=10.

Since both old Hamilton counts vanish, every E coordinate vanishes. Equation (4) gives 6 and 20 entirely from the two-path memory. These are S32's already-known compulsory-ear phenomena, now interpreted as exact quadratic response coefficients.

One certificate at 23 is the joined pair

    (2),
    (13,12,4,21,15,10,6,19,17,8,1,3,22,14,11,5,20,16,9,7,18).

Insert 23 between these paths. Every resulting consecutive sum is a square. The certificate JSON lists all three ordered joining pairs with first port 2, and all ten with first port 11 for the 25-arrival.

## 6. What remains open or requires more memory

The quadratic record is complete for **one** additional vertex. It is not claimed to be closed under repeated extension. Deleting k new vertices from a Hamilton path can leave up to k+1 nonempty old path pieces, so higher contexts can require higher path-cover port data. S171's square-zero boundary response retains those data, with potentially exponential cost.

A concrete failure of iterative closure uses five old vertices. Let U be the edgeless graph and V consist of disjoint edges 01 and 23 and isolated vertex 4. Neither has a spanning cover by at most two paths, so both have E=K=0 and every one-vertex response is zero. Now add two new vertices, each adjacent to every old vertex but not to the other new vertex. For U the result is K_(5,2), which has no Hamilton path because the two bipartition sizes differ by three. For V, a path alternates its three old components with the two new vertices. There are exactly

    3! × 2! × 2² = 48

oriented Hamilton paths: order the three components, order the two new connectors, and orient each old edge. Thus equal one-step contextual memory is **not** sufficient for two-step contexts. This realizes, rather than merely cautions about, the difference between a finite experiment bank and its closure under continuation.

The exact rank in Theorem 2 ranges over all labelled simple graphs and all one-vertex neighborhoods. Restricting to arithmetic square-sum sources or arithmetic-legal arrivals can reduce it drastically; the theorem does not supply a lower bound for that restricted family. A meaningful next target is the rank of explicitly declared arithmetic block/context families, followed by a closure law that survives the chosen gluing operations.

An arithmetic realization of the two equal-profile seeds with an authenticated square-sum read head would be genuinely stronger than the present result. It is a research question, not a consequence of sharing the words “Hamilton path.” No improvement to the known all-n square-sum existence result is claimed.

## 7. Reproduction

Run:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json

The standalone checker imports no repository files. It exhausts all labelled graphs through order five and all their one-vertex neighborhoods; compares permutation/cut coefficient enumeration to an independent endpoint subset DP; computes exact rational ranks; verifies the complete profiles and nonisomorphism of A,B; and checks the square-sum counts and cut-cover certificates independently. All validation uses explicit failures rather than assert statements.

Verified SHA-256:

    check.py: cf250701e222fd5984929797fc9a585dbd6d8f87c434e3bcc16232a91421cd12
    certificate.json: f45a5ec57555b3b8c0b4838e6268940654a9d9a95bfb99eb932d2f2ed6c3924e

Exact ranks for n=1,2,3,4,5 are 1,2,6,10,15. The exceptional n=2 rank is two, so the theorem's n≥3 hypothesis cannot simply be dropped. There are 33,866 checked extension contexts in total.
