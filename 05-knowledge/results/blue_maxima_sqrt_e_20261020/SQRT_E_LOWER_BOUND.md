# A square-root-of-e lower bound for maximum blue mass

Status: **PROVED; INDEPENDENT MATHEMATICAL REVIEW APPROVED**. The only new external existence input is the cited block-design theorem. No matching asymptotic upper bound is claimed.

## Theorem

With the inherited definition `M_k=max_D h(D)` and its identification with maximum blue mass at odd order `2k+1`,

    liminf_(k -> infinity) M_k/k! >= sqrt(e).

The statement is stronger than the uniform-frame mean or second-moment lower bound. It does not assert that the ratio converges, or that `sqrt(e)` is an upper bound.

## 1. A local four-index model

Take the four-index frame with lexicographic masks

    A=2,     b=11.

Its eight signed vertices have external adjacency row masks

    [56,200,35,224,14,194,21,21].

Matching edges are omitted. Every vertex has outdegree3 in the underlying complete four-partite graph with two vertices per part. The reflection swapping the two vertices of each part reverses every edge. Equivalently, the first three indices induce a cyclic tournament, the fourth dominates them, and `b` is the clique on the first three indices.

Randomize this fixed frame by a uniform permutation of its four pairs and independent uniform exchanges of the two vertices in each pair. This produces a valid hollow frame on the same four index pairs. The exchanges do not have to preserve the center's orientation: this is a distribution on **external** frames, which are completed by the usual fixed center afterward.

For any specified signed vertices from distinct pairs:

1. A specified directed edge occurs with probability `1/2`.
2. Two specified directed edges on four distinct pairs occur jointly with probability `1/4`.
3. A specified directed two-edge path on three distinct pairs occurs with probability `5/16`.
4. A specified directed three-edge path on all four pairs occurs with probability `3/16`.

For1, the signed-permutation group is transitive on ordered signed vertices from distinct pairs, and exactly half of such ordered pairs are edges.

For2, exchange both pairs belonging to the first edge. This reverses that edge by the anti-involution and leaves the second edge unchanged. Independently do the same for the second edge. These two measure-preserving transformations give all four orientation combinations equally often.

For3, transitivity reduces the probability to the proportion of all proper directed two-edge paths. Each of the eight vertices has indegree and outdegree3, hence there are `8*3*3=72` directed two-edge paths before excluding paths whose endpoints lie in one pair. Exactly12 are excluded. To see the latter count for any four-index frame, for each unordered pair of pair-indices, precisely one of its endpoints' signed vertices has mixed in/out incidence to the other pair; its reflected vertex also does. Thus each of the six index pairs accounts for two excluded paths. There remain60 proper paths, out of `4*3*2*2^3=192` ordered signed triples. Their ratio is `60/192=5/16`.

For4, write the three required edge signs as `Z_1,Z_2,Z_3` in `{+1,-1}`, with `+1` indicating agreement with the path. All odd products have expectation zero: exchanging all four pairs reverses every edge and preserves the distribution. Adjacent-edge covariance is `4*(5/16)-1=1/4`, while disjoint-edge covariance is zero by2. Expanding `E[(1+Z_1)(1+Z_2)(1+Z_3)]/8` therefore gives `(1+1/4+1/4)/8=3/16`.

The script `verify_local_design.py` checks these probabilities exactly on all 384 signed permutations, without relying on the counting proof.

## 2. Pack almost all index pairs into four-sets

Let `v` be the largest integer at most `k` with `v=4 (mod12)`. For `k>=4`, `0<=k-v<=11`.

**CITED input:** Hanani, Ray-Chaudhuri and Wilson, *On resolvable designs*, Discrete Mathematics3 (1972),343–357, prove the existence of a resolvable `(v,4,1)` BIBD for every `v=4 (mod12)`. Only its underlying block design, not resolvability, is needed.

- Primary publisher record and theorem statement: <https://doi.org/10.1016/0012-365X(72)90091-X>
- Publisher page: <https://www.sciencedirect.com/science/article/pii/0012365X7290091X>

Consequently there is a family `B` of four-element index subsets with every pair among the first `v` indices in exactly one block. Distinct blocks meet in at most one index, and

    q=|B|=v(v-1)/12 = k^2/12 + O(k).

On each block independently sample the local four-index distribution above. For every uncovered index pair, sample its two external bits independently and fairly. This consistently defines an external frame on all `k` indices, because block edge sets are disjoint.

## 3. Exact probability for every candidate word

Let `F` be the collection of all three-element subsets of blocks. Each triple has a unique parent block. For an index permutation `p` and a block `B`, let `X_B(p)` count its triples that occur as three consecutive indices in `p`.

The index Hamiltonian path intersects a four-element block in a linear forest. It can use zero edges, one edge, two disjoint edges, a two-edge path, or a three-edge spanning path. In these five cases `X_B` is respectively `0,0,0,1,2`. The local probabilities from Section1, relative to independently fair edge orientations, are respectively

    1, 1, 1, 5/4, 3/2.

Thus the local correction factor is exactly `1+X_B/4` in every case. Independence between blocks gives, for every allowed polarity assignment,

    Pr(the candidate signed word is directed)
      = 2^(-(k-1)) product_(B in B) (1+X_B(p)/4).

There are `2^(k-1)` assignments with final sign `x`. Summing their probabilities and then averaging over the index permutations proves

    E_D[h(D)]/k! = E_p[product_(B in B) (1+X_B(p)/4)].       (1)

This expectation uses the newly constructed, correlated distribution on frames, not the uniform labelled-frame distribution from the inherited mean theorem. No exceptional or unclean permutations are discarded.

## 4. A nonnegative expansion and a finite lower bound

Since distinct blocks share at most one vertex,

    f=|F|=4q = k^2/3 + O(k).

Each index lies in at most `k-1` triples of `F`: in the full design on `v` vertices it belongs to `(v-1)/3` blocks and three triples within each block.

Expand the product in (1). In a term choosing `r` distinct blocks, expand each selected `X_B` into the sum of its four triple-consecutiveness indicators. All coefficients are nonnegative. Retain only terms whose selected triples are vertex-disjoint. Such triples automatically have distinct parent blocks, so every retained collection is indeed present once in the product expansion.

Let `d_r` denote the number of unordered collections of `r` vertex-disjoint triples in `F`. Each fixed collection has simultaneous consecutiveness probability

    6^r (k-2r)! / k!,

by contracting the triples and choosing their internal orders. Therefore, for every `R<=floor(k/3)`,

    M_k/k! >= E_D[h(D)]/k!
      >= sum_(r=0,...,R) (3/2)^r d_r (k-2r)!/k!.          (2)

Here `d_0=1`. This is an exact finite inequality for any supplied four-block packing.

A completely explicit bound follows by counting disjoint triples greedily. After `j` triples have been chosen, their `3j` vertices meet at most `3j(k-1)` triples. Thus, whenever the factors are nonnegative,

    d_r >= 1/r! product_(j=0,...,r-1) (4q-3j(k-1)).        (3)

If one factor is negative, use the trivial lower bound `d_r>=0` instead. Substituting (3) into (2) gives a finite rational lower bound depending only on `k` and `q=v(v-1)/12`. It is not claimed sharp for small `k`.

For example the `r=1` term alone gives

    M_k/k! >= 1 + 6q/[k(k-1)]
            = 1 + v(v-1)/[2k(k-1)].                     (4)

When `k=4 (mod12)`, this is exactly `M_k >= (3/2) k!`.

## 5. Finish the asymptotic argument

For every fixed `r`, equations (2) and (3), together with `4q=k^2/3+O(k)`, show that the retained `r`-th contribution has lower limit

    (3/2)^r/r! * (k-2r)!/k! * (k^2/3)^r
      -> 1/(2^r r!).

Consequently, for every fixed `R`,

    liminf M_k/k! >= sum_(r=0,...,R) 1/(2^r r!).

Let `R` tend to infinity. The right side increases to `exp(1/2)=sqrt(e)`, proving the theorem.

The argument uses only finite nonnegative expansions before taking limits. It neither needs Poisson convergence nor interchanges an unbounded exponential expectation with a distributional limit.

## 6. Explicit one-sided error bounds

The same argument gives the following uniform estimates:

    M_k/k! >= sqrt(e) (1-101/(8k))       for k>=4,          (5)
    M_k/k! >= sqrt(e) (1-13/(8k))        if k=4 (mod12).    (6)

These may be combined with the elementary lower bound1; they are not asymptotic equalities.

For (5), first suppose `k>=13`. Since `v>=k-11`,

    4q = v(v-1)/3 >= (k^2-23k+132)/3.

Let `a_r=1/(2^r r!)` and `C_r=23r+(9/2)r(r-1)`. Define the actual retained contribution `t_r` to be `(3/2)^r d_r (k-2r)!/k!` when `r<=floor(k/3)`, and zero otherwise.

If `C_r>=k`, then `t_r>=0>=a_r(1-C_r/k)`. If `C_r<k`, then `r<k/23` for `r>=1`, so `r<=floor(k/3)`. All numbers `(23+9j)/k`, `0<=j<r`, are in `[0,1)`, and (3) together with `k!/(k-2r)!<=k^(2r)` gives

    t_r >= a_r product_(j=0,...,r-1) (1-(23+9j)/k)
        >= a_r (1-C_r/k).

The last inequality is the elementary product bound `product(1-u_j)>=1-sum u_j` for `u_j` in `[0,1]`. This also covers `r=0`. Thus the termwise lower estimate holds for every nonnegative integer `r`, including the artificial zero tail. Summing the absolutely convergent right side gives

    sum a_r = sqrt(e),
    sum r a_r = sqrt(e)/2,
    sum r(r-1) a_r = sqrt(e)/4,

and hence `sum a_r C_r=(101/8)sqrt(e)`. This proves (5) for `k>=13`. For `4<=k<=12` its right side is negative, so the assertion is automatic.

For (6), use `v=k` and the exact factorization

    4q-3j(k-1) = (k^2/3)(1-1/k)(1-9j/k).

Set `C'_r=r+(9/2)r(r-1)`. Whenever `C'_r<k`, all product factors are positive and `r<=floor(k/3)`: for `r>=2`, `C'_r>=3r`, and `r=0,1` are immediate. The same product bound now gives `t_r>=a_r(1-C'_r/k)`. When `C'_r>=k` the comparison is trivial by nonnegativity. Summing uses `sum a_r C'_r=(13/8)sqrt(e)` and proves (6).

The exact local verifier also checks the rational termwise comparisons over a bounded range of `k,r`; the inequalities above are their all-order proof.

## Scope and attribution

The proof adapts the correlated-small-block strategy of Adler, Alon and Ross, *On the maximum number of Hamiltonian paths in tournaments*, Random Structures & Algorithms18 (2001),291–296, DOI <https://doi.org/10.1002/rsa.1010>, author-hosted primary text <https://adler.ieor.berkeley.edu/ilans_pubs/hamilt_2001.pdf>. Their ordinary-tournament theorem has a different normalization and cannot itself be cited as a theorem about these hollow frames.

The local four-index distribution, its `5/16` wedge probability, and the square-root-of-e calculation above are proved explicitly for the present objects. No literature-priority claim is made. The block-design existence theorem is an external input, not newly proved here. The inherited equality between `max_D h(D)` and maximum blue mass uses its already documented regular-subset theorem.
