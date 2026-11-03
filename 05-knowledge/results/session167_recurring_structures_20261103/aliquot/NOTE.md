# Divisor lattices and aliquot pseudoforests inside one path-profile fiber

**Status.** Proved compositional embeddings and finite-exact arithmetic controls. The cube, ordinal-sum product, strong wrapper and six-vertex binary seed are inherited from S166. The ternary seed was independently found by this lane and the same-session fences/switching lane; the latter owns its minimal-order and general ternary-product analysis. The new applications here retain complete divisor lattices, their arithmetic weights and finite aliquot trajectory graphs, and extend the construction to cactus graphs. No aliquot convergence theorem, new divisor-balance classification, or identification with the blue/black complement-line metagraph is claimed. No literature-priority claim is made.

## 1. Exact interface being used

Write M_n for the simple graph of isomorphism classes of n-vertex tournaments, adjacent when representatives differ in one arc. Write S_n for its induced strong-tournament subgraph. These are the single-arc-flip graphs, not the repository's blue/black metagraph.

For a tournament T let p_j(T) count partitions of its vertices into j nonempty directed paths, with each path internally ordered but the collection unordered. Put

\[
P_T(z)=\sum_jp_j(T)z^j,\qquad F_T(t)=\sum_jp_j(T)(t)_j.
\]

Thus H(T)=p_1(T)=F_T(1). The inherited ordinal-sum law is

\[
F_{T\oplus U}(t)=F_T(t)F_U(t).
\]

S166 supplies strong six-vertex tournaments B_0,B_1 with masks 18,22 and

\[
P_{B_0}=P_{B_1}=15z+93z^2+138z^3+73z^4+15z^5+z^6.
\]

Masks use lexicographic unordered pairs and bit 1 means i→j for i<j. Their masks differ in exactly bit 2, the arc on {0,3}. They are nonisomorphic. For w∈{0,1}^k define T_w=B_{w_1}⊕⋯⊕B_{w_k}. Then

\[
w\longmapsto[T_w]
\]

is an **induced** embedding of Q_k in M_{6k}. Every state has F_T=F_B^k and H=15^k.

For completeness, the key quotient check is structural. The strong components of T_w are exactly the six-vertex blocks, and their condensation gives their unique order. Thus isomorphisms cannot exchange block positions. A cross-block backward flip merges an interval of strong components and leaves the family. A flip inside one block either destroys its strongness, leaves its class unchanged, or changes that one coordinate. No extra class-level edge appears.

The inherited wrapper W sends M_n inducedly into S_{5n+10}: take two regular tournaments A,B of order 2n+5 and orient A→T→B→A. The outdegree bands are separated sufficiently to recognize all three sets even after one flip. Substituting modules with equal complete profiles preserves the complete profile, so the wrapper preserves equality of profiles across our families. **The wrapped common H need not equal the unwrapped H.**

Source: [S166 Cartesian flip note, §§1–2 and 5](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/cartesian_flip/NOTE.md).

## 2. Divisor-lattice realization, retaining the missing arithmetic marks

Let N=∏_{i=1}^r p_i^{a_i}≥2, with an ordered prime factorization, and k=Ω(N)=Σa_i. For d=∏p_i^{b_i} dividing N, replace each exponent b_i by its unary word

\[
\underbrace{11\cdots1}_{b_i}\underbrace{00\cdots0}_{a_i-b_i}
\]

and concatenate these words to obtain c(d)∈{0,1}^k.

**Theorem 1.** The map d↦[T_{c(d)}] induces an embedding of the undirected divisor Hasse graph of N into M_{6Ω(N)}, and into S_{30Ω(N)+10} after wrapping. The unwrapped states all have H=15^{Ω(N)} and the same complete path-cover profile. The underlying divisor-to-cube map is isometric:

\[
d_H(c(d),c(e))=\sum_i|\nu_{p_i}(d)-\nu_{p_i}(e)|.
\]

**Proof.** Within one unary block, changing b to c changes exactly |b−c| bits. Sum over coordinates. In particular Hamming distance one is exactly multiplication or division by one prime while staying among divisors of N. The preceding induced cube theorem proves the unlabelled-class statement. The wrapper proves the strong statement. ∎

The squarefree divisors are the marked subcube b_i∈{0,1}; the prime divisors are its r atom vertices. These are retained graph strata, not inferred from H. Prime values and the ordered coordinate partition are marks required to recover arithmetic values.

This realization is compositional in an exact counting sense. Its weighted path-cover partition function is

\[
Z_N(x_1,\ldots,x_r;t)
 =\sum_{d\mid N}\Bigl(\prod_i x_i^{\nu_{p_i}(d)}\Bigr)F_{T_{c(d)}}(t)
 =F_B(t)^k\prod_i(1+x_i+\cdots+x_i^{a_i}).                 \tag{1}
\]

The product follows from the independent exponent coordinates and the ordinal-sum profile law. This counts a **family of marked tournaments**, not the Hamilton paths of one tournament whose count is σ(N).

At x_i=p_i and t=1, equation (1) gives

\[
15^{-k}Z_N(p_1,\ldots,p_r;1)=\sigma(N),\qquad
s(N)=15^{-k}Z_N(p_1,\ldots,p_r;1)-N,
\]

where s is the ordinary aliquot sum. At all x_i=1, the normalized value 15^{-k}Z_N(1,…,1;1) is the divisor cardinality τ(N).

Let F count proper nontrivial divisors, S their squarefree members and U their prime members, as in the inherited divisor-balance note. Deleting the bottom and top vertices gives F, intersecting with the marked squarefree subcube gives S, and restricting to the atoms gives U, with the prime endpoint removed when N itself is prime. Hence

\[
F=\prod_i(a_i+1)-2,\quad S=2^r-1-[N\text{ squarefree}],\quad U=r-[N\text{ prime}].
\]

The inherited classification F=S+U iff N=p,p³,p²qr, where p,q,r are pairwise distinct primes is transported literally to these marked strata. It is **not** a new proof of that classification or a consequence of a special value of the tournament path count. Its source is [arithmetic braids, DB1](https://github.com/eliottcassidy2000/math/blob/99181f1707aff5f8323e21f1cfa86527a3413af1/05-knowledge/results/arithmetic_braids_20260917_divisors.md).

Two information losses are now exact:

* N=8,20,30 all have Ω=3 and thus use the same common tournament profile; their exponent grids have 4,6,8 vertices and their divisor defects are 0,−1,−3. A single state's full profile does not determine the marked grid or its balance.
* N=290=2·5·29 and 1729=7·13·19 have the same unweighted marked exponent grid B_3 and the same represented profile. But σ(290)=540 and σ(1729)=2240, giving s(290)=250 and s(1729)=511. The missing data are the prime weights, not an unexplained numerical exception.

## 3. Every finite rooted tree becomes literal one-arc dynamics

**Theorem 2.** Let G be a finite rooted tree with m≥1 edges, directed toward its root. It has an explicit induced representation in M_{6m}, with all profiles equal and H=15^m. Every directed tree step is a literal arc reversal inside one fixed six-vertex block, with no relabelling needed. A strong representation exists in S_{30m+10}, still with one common complete profile.

**Proof.** Give each tree edge its own bit. The code c(v) is the incidence vector of edges on the unique root-to-v path. For any u,v the symmetric difference of those paths is the u-to-v path. Therefore d_H(c(u),c(v))=dist_G(u,v), so the code is injective and induced. A step to the parent clears exactly the bit of its last edge. Because masks 18 and 22 differ in one literal arc, T_{c(v)} changes to T_{c(parent(v))} by reversing that arc in that block. Apply the inherited wrapper for strong states. ∎

This uses a tree's actual edges to choose the coordinates; it does not automatically compute an unknown arithmetic map or prove that its orbits are trees.

For the ordinary aliquot map s(n)=Σ_{d|n,d<n}d, with s(1)=0, the specified starts have these exact trajectories:

```
20 → 22 → 14 → 10 → 8 → 7 → 1 → 0
21 → 11 → 1 → 0
290 → 250 → 218 → 112 → 136 → 134 → 70 → 74 → 40 → 50 → 43 → 1 → 0
1729 → 511 → 81 → 40 → 50 → 43 → 1 → 0
```

Their union is a tree with **24 vertices and 23 edges**, rooted at 0. The certificate assigns one edge bit to each of its 23 edges and gives all 24 codes. Theorem 2 yields a literal induced copy in M_138, every state having H=15^23 and the same full profile, and a strong copy in S_700. The vertex 0 is a halt state; a loop s(0)=0 is not represented as a nontrivial flip.

Both the forward arithmetic transitions and all nonadjacencies are retained. In particular the join at 40 and the join at 1 are actual graph incidences, while the path-profile statistic is deliberately constant. The numeric sizes of the vertices may increase along an aliquot step, e.g. 20→22; it is tree depth, not integer value, that decreases in this representation.

This verifies precisely the four displayed starting values. It gives no assertion about termination of general aliquot sequences, and does not reinterpret other pseudoforests unless their precise edge rule is first identified with this s-map.

## 4. A cactus extension, and why length five is different here

The same-session shared seed consists of strong six-vertex masks 153,177,181. They are pairwise adjacent, nonisomorphic, and share

\[
P_C(z)=29z+131z^2+162z^3+77z^4+15z^5+z^6.
\]

The companion independently checks every pair against all 720 relabellings and checks the profiles by two algorithms. The fences/switching packet supplies the broader ternary-product treatment and seed minimality. The same ordinal-component proof embeds the ternary Hamming graph K_3^{\square d} inducedly in M_{6d}, with common F_C^d and H=29^d.

**Theorem 3.** A finite cactus graph has an induced embedding in some ternary Hamming graph if and only if none of its cycles has length five. Consequently every such cactus graph is an induced subgraph of a single complete-profile fiber of a tournament flip graph, and also of a strong-tournament flip graph after wrapping.

A cactus means a simple graph whose nontrivial blocks are edges or cycles. Disconnected graphs and isolated vertices are allowed.

**Proof of sufficiency.** Edges and triangles are in K_3. Every even cycle C_{2k}, k≥2, embeds inducedly in Q_k by starting at 0 and flipping coordinates 1,…,k,1,…,k cyclically. Q_k is an induced subgraph of K_3^{\square k}.

An induced C_7 in K_3□Q_2, in cyclic order, is

```
000, 010, 110, 111, 101, 201, 200.
```

To lengthen any induced cycle C_l, l≥4, by two in its Cartesian product with K_2, replace consecutive vertices u,v,w by

```
(u,0), (u,1), (v,1), (w,1), (w,0)
```

and put all other cycle vertices in layer 0. The old vertex (v,0) is omitted. Since u and w are nonadjacent in C_l, the displayed path has no chord; the other old adjacencies are unchanged and cross-layer edges only join equal old vertices. Starting at C_7 gives all odd lengths at least seven. These constructions use floor(l/2) coordinates for l≥4, l≠5, and one coordinate for C_3.

Induced Hamming embeddings are closed under gluing at one vertex. If f:G→K_3^a and g:H→K_3^b are induced and their chosen roots map to r,s, send old vertices x to (f(x),s) and new vertices y to (r,g(y)). Across distinct nonroot vertices the distance is at least two, so no extra edge appears. This permits the block-by-block construction of a connected cactus. For disconnected components, pad their dimensions and append duplicated binary component labels. Distinct labels have Hamming distance at least two, again excluding unwanted edges. Pad with a constant coordinate if necessary for an isolated single vertex.

**Proof of necessity.** An induced C_5 cannot occur in any Hamming graph with complete factors, in particular a ternary one. Label each cycle edge by its changed coordinate. Every used coordinate must occur at least twice because the walk closes. The only possible occurrence patterns over five steps are 5 or 3+2. In the first case all five vertices lie in one complete factor, giving chords (and a ternary factor has only three vertices). In the second, the three-change coordinate uses three values and the two-change coordinate uses two, so the cycle is contained in a triangular prism K_3□K_2. Any five of its six vertices contain all three vertices of one triangle; a five-cycle induced subgraph cannot contain such a triangle. This is impossible. A cycle of a cactus is induced, giving the required exclusion. ∎

For a connected nontrivial cactus a constructive dimension bound is

\[
d=b+\sum_{C_3}1+\sum_{C_l,\,l\ge4,\,l\ne5}\lfloor l/2\rfloor,
\]

where b is the number of bridge edges. Its tournament order is 6d before wrapping, 30d+10 afterward; common unwrapped H is 29^d. These are constructive upper bounds, not minimal tournament orders.

This includes finite trees and unicyclic pseudoforests with no five-cycle, hence many finite functional-graph shapes after directions are kept as marks. Loops and parallel directions are not edges of the simple flip graph; an amicable two-cycle has one underlying undirected edge plus two directed marks. A sociable five-cycle is excluded from **this ternary Hamming construction**, not proved absent from all constant-profile tournament metagraph fibers. Resolving that latter question would need a new seed or mechanism.

## 5. Reproduction and precise scope

Run `python3 check.py` and `python3 -O check.py`. Both produce the same deterministic `certificate.json`; checks raise explicit errors and remain active under optimization. The standard-library companion verifies:

* binary and ternary seeds by all 720 relabellings per class/pair, strongness, and two independent path-cover counts;
* four literal 12-vertex ordinal sums by Hamilton subset DP;
* all 1,224 single-arc flips of the eight literal 18-vertex binary cube states;
* divisor enumeration, unary isometry, Hasse adjacency, strata, and arithmetic weights for N=8,20,21,60,290,1729;
* every step of all four aliquot trajectories by trial-divisor and factorized-σ calculations, all 576 ordered tree-distance checks, and the 24-vertex/23-edge certificate;
* explicit induced cycles of every length 3 through 31 except 5, and all 720 ordered five-vertex choices in the triangular prism.

Finite computations certify the named seeds and examples. The universal divisor, tree and cactus claims rest on the proofs above and the inherited induced ordinal-product and wrapper theorems. There is no claim that ordinary positive profiles decide arithmetic dynamics, that every abstract metagraph is realized, or that the ternary Hamming exclusion is a general tournament obstruction.

In `certificate.json`, vertex codes are printed most-significant bit first; `edge_coordinates[i]` is bit `1<<i`, so coordinate zero is the rightmost bit.
