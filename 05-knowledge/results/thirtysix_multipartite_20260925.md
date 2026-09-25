# Two multipartite relations with 36 edges on the same ten letters

Status: **PROVED**, elementary identities and constructions below; **FINITE-EXACT**, the explicitly bounded census and controls. Date: 2026-09-25. These graphs are not asserted to be Hamiltonian-path kernels, square-sum graphs, or Collatz convergence certificates.

## 1. Inheritance and the object being compared

The closest mechanism is the rank-two alternating-form classification in [duck_decoder_20260925.md, section 4](duck_decoder_20260925.md): a symmetric zero-diagonal binary matrix of rank two is a complete graph between at most three nonzero coordinate classes, with a possible isolated zero class. With all degrees odd, only two classes of odd sizes survive. The six-vertex Hamiltonian-kernel census contained both the planar K1,5 and the nonplanar K3,3, so nonplanarity did not identify that obstruction.

The positive input is the same note's ten-letter object `D=Sym²(V)`, `V=F₂²`, and its [prime-divisor realization](duck_primes_20260925.md). Its letters are unordered pairs with repetition from `V`. The order-three linear map `M` fixes zero and rotates the three nonzero colors `a,b,c`, where `a+b=c`. The distinction between virtual color and the numeric value of a divisor remains essential. The [Fibonacci carry note](duck_zeckendorf_20260925.md) also requires the neutral fourth state and a carry frame; a three-color rotation is not itself a ternary odometer.

Live board: ten-letter rotation channels; XOR fibers; multipartite rank; graph orientation and quotient; dyadic replication; ternary tournament substitution. The anchor is an actual graph operator on the inherited letters. The niche is a rank/parity classification. The wildcard is the nine-vertex regular tournament obtained by deleting the fixed letter and filling three holes.

Here `K_(n₁,...,n_r)` means a **simple undirected complete multipartite graph**, with every `n_i>0`: vertices in different parts are adjacent; vertices in one part are not. The edge count is

\[
 E=\sum_{i<j}n_i n_j
   =\frac{N^2-\sum_i n_i^2}{2},\qquad N=\sum_i n_i.
 \tag{1}
\]

Its quotient by parts is K_r. It is not an r-vertex tournament until directions have been chosen. The nonedges are actual missing edges, not a claim about equality of an external observable.

## 2. An exact six-edge trade with an odd boundary

**M1 (PROVED).** On the same set `D` there are two natural C₃-invariant relations with exactly 36 edges:

* `G_O`: join letters in different C₃ orbits;
* `G_X`: join letters with different pair-XOR charges `q({u,v})=u+v`.

The respective partitions are

\[
\begin{array}{c|l}
\text{orbit parts}&\{00\},\ \{aa,bb,cc\},\ \{0a,0b,0c\},\ \{ab,bc,ca\}\\
\text{XOR fibers}&\{00,aa,bb,cc\},\ \{0a,bc\},\ \{0b,ac\},\ \{0c,ab\}.
\end{array}
\]

Thus

\[
G_O=K_{1,3,3,3},\qquad G_X=K_{4,2,2,2},
\qquad E(G_O)=E(G_X)=\frac{100-28}{2}=36.
\tag{2}
\]

The partitions are different but both have the same sum of squared part sizes, 28. The map M preserves the first partition and permutes the last three fibers of the second, so both edge relations are equivariant under the given action.

**M2 (PROVED).** To turn `G_O` into `G_X`, remove the following six edges:

\[
 00-aa,\ 00-bb,\ 00-cc,
 \qquad 0a-bc,\ 0b-ac,\ 0c-ab,
\tag{3}
\]

and add the two triangles on `{0a,0b,0c}` and `{ab,bc,ca}`. This lists every difference between the two defining equivalence relations. They therefore share 30 edges. Their symmetric difference is the disjoint union of a star K1,3 and a triangular prism. Every one of its ten vertices has odd degree: the three repeated nonzero letters have degree one, and the other seven letters have degree three.

In binary edge-chain notation, the trade `H=G_O+G_X` has boundary `∂H=1_D`. Accordingly it changes odd degree at every vertex into even degree at every vertex:

| Object | Degree multiset | Binary adjacency rank | Triangles |
|---|---|---:|---:|
| G_O | one 9, nine 7s | 4 | 54 |
| G_X | four 6s, six 8s | 4 | 56 |
| G_O symmetric-difference G_X | three 1s, seven 3s | not used | two |

The triangle counts follow by summing `n_i n_j n_k` over triples of parts. The differing degrees already prove that the two graphs are not isomorphic. Their equality of edge count and rank does not identify their edge information.

This is a graph operator connecting rotation channels to an additive color observable, with an explicit parity defect. It preserves the ten-letter set, its C₃ action, the edge count, and the rank. It changes adjacency, vertex degrees, triangles, and the partition being quotiented. Recovering both relations requires the letter or equivalent orbit-and-charge sidecar. Forgetting to just a four-vertex quotient loses the part sizes and the trade.

For the virtual prime realization with `N=60`, the letters may be labeled

\[
00\mapsto30,\quad(0a,0b,0c)\mapsto(2,3,5),\quad
(ab,ac,bc)\mapsto(6,10,15),\quad(aa,bb,cc)\mapsto(4,12,20).
\]

In `G_X`, the letters labeled 2 and 3 are adjacent because their virtual colors differ, although 2+3 is not a square. In `G_O` those same letters are nonadjacent. Neither relation is numeric divisor incidence or the square-sum relation. No claim about graceful labels follows from (2).

## 3. The exact rank and parity content of two, three, or four parts

**M3 (PROVED).** For all positive part sizes,

\[
\operatorname{rank}_{\mathbb F_2} A(K_{n_1,\ldots,n_r})
=\begin{cases}r&r\text{ even},\\r-1&r\text{ odd}.
\end{cases}
\tag{4}
\]

Let P be the vertex-by-part indicator matrix. Over F₂,

\[
 A=P(J_r+I_r)P^T.
\]

Because the parts are disjoint and nonempty, P is injective and Pᵀ is surjective. These compositions preserve the rank of `J+I`. A vector in its kernel satisfies `x=(sum x_i)1`; the all-ones vector belongs to the kernel exactly for odd r. This proves (4), including the case r=1 of an edgeless graph.

Each vertex in part i has degree `N−n_i`. Hence all degrees are odd if and only if N is even and all part sizes are odd. Indeed odd N would require every part size even, contradicting their sum. With N even, all odd part sizes force r even. In particular rank two plus odd degree forces exactly two odd parts. Rank two without the degree condition permits three parts as well.

The color mechanism can be written explicitly. For vectors `u,v` in F₂², put

\[
 u\sim v\quad\Longleftrightarrow\quad
 \det(u,v)=1.
\]

All unequal nonzero colors are adjacent, equal colors are nonadjacent, and color zero is isolated. This is the rank-two tripartite relation. The order-three color rotation preserves it. It differs from **inequality of four colors**, whose quotient is K4 and has rank four. The fourth color cannot simply be inserted into a rank-two symplectic rule as if it were another nonzero vector.

For `G_O`, begin with K3,3,3 and the isolated fixed letter, then add its nine incident edges. For `G_X`, the symplectic charge rule first gives K2,2,2 and four isolated zero-charge letters; changing to charge inequality adds the 24 edges from the zero fiber to the other six letters. Both operations raise the rank from two to four, but with different degree-parity outcomes.

On a fixed tripartite graph, cycling the three parts is an automorphism only when their sizes agree. Unequal sizes still admit the quotient color action and a relabeling isomorphism to the permuted part-size presentation. This distinction is needed, for example, for K2,3,6.

## 4. All complete multipartite graphs with 36 edges

For two through four parts, the classification is elementary:

| Number of parts | All unordered positive part sizes giving 36 edges |
|---:|---|
| 2 | (1,36), (2,18), (3,12), (4,9), (6,6) |
| 3 | (2,2,8), (2,3,6) |
| 4 | (1,1,1,11), (1,3,3,3), (2,2,2,4) |

For two parts this is the divisor-pair classification of 36. For three sorted parts `a≤b≤c`, one has `a≤3` and `(a+b)(a+c)=36+a²`. At a=1 the right side is prime; a=2 gives factor pairs `(4,10),(5,8)`; a=3 gives none with both factors at least six. These produce the listed two triples.

For four sorted parts, `a≤2` since `6a²≤36`. When a=1, the lower bound `3b(b+1)≤36` gives b=1,2,3. The equations for `(c,d)` are respectively `(c+2)(d+2)=39`, `(c+3)(d+3)=43`, and `(c+4)(d+4)=49`, yielding only `(1,11)` and `(3,3)`. When a=2, the bound `3b²+6b≤36` forces b=2; then `(c+4)(d+4)=48` gives `(c,d)=(2,4)`.

There are twelve types in total: the ten above, `(1,1,1,3,4)`, and nine singleton parts (K9). There are no six-, seven-, or eight-part types. For a short proof of this extension, write `n_i=1+x_i`, `t=sum x_i`, and `C=sum_(i<j) x_i x_j`. Then

\[
36=\binom r2+(r-1)t+C,\qquad 0\le C\le\binom t2.
\tag{4a}
\]

If at least two x_i are positive, also `C≥t−1`. For r=9, (4a) forces t=0. For r=8, t≤1 cannot supply the required extra eight edges. For r=7, t≤2 supplies at most thirteen of the required fifteen. For r=6, t≤3 supplies at most eighteen of the required twenty-one, while t=4 would require C=1, which is neither zero nor at least three. For r=5, t≤4 supplies at most twenty-two of the required twenty-six. At t=5 one needs C=6, realized only by positive excesses 3 and 2; at t=6 one needs C=2, impossible because a nonzero C is at least five. The t=5 assertion follows from the seven partitions of five: their C values are 0,4,6,7,8,9,10.

The script independently exhausts **all** numbers of positive parts: `E≥N−1` bounds `N≤37`, while `binom(r,2)≤36` bounds `r≤9`; sorted parts are exhausted with an exact incremental edge budget. Each output type is independently recounted and given a rank and planarity certificate.

Exactly two of these twelve graphs are planar: K1,36 and K2,18. The first is a tree. For K2,n, order the n degree-two vertices forward around one pole and backward around the other; this rotation has n faces, so `V−E+F=(n+2)−2n+n=2`. Every remaining output has an explicitly listed K3,3 subgraph. K3,3 cannot be planar because its bipartite face lengths would imply `E≤2V−4=8`, but E=9. The routine also supports a K5 certificate if needed, although all ten nonplanar 36-edge types have a K3,3 witness.

In particular K3,3,3 has **27**, not 36, edges. K6,6 has 36 edges and twelve vertices. K9 has 36 edges and nine vertices. The two graphs of section 2 have 36 edges and ten vertices. The edge count alone does not specify the vertex set or quotient.

## 5. A ternary construction that really gives a nine-vertex tournament

**M4 (PROVED).** Delete the fixed letter `00` from `G_O`. The three remaining rotation channels each have size three, giving K3,3,3 with 27 edges. Fill the three missing within-channel triangles, adding nine edges. The result is K9, with

\[
 27+3\cdot3=36=\binom92=8\cdot9/2.
\tag{5}
\]

Orient each channel by its existing color rotation. Choose a cyclic order of the three channels and orient all crossing edges uniformly in that order. This gives the regular tournament `C₃[C₃,C₃,C₃]`. Each vertex has one outgoing edge within its channel, three toward the next channel, and four incoming edges. The given color rotation is an automorphism; the cyclic order of the channels is an additional specified gauge, not determined by the bare set with its C₃ action.

This has an exact recursion. On words of length d over `{0,1,2}`, orient a pair by the first position where the words differ, using the directed triangle at that position. The resulting tournament T_d has `3^d` vertices and outdegree `(3^d−1)/2`. Equivalently,

\[
 T_{d+1}=C_3[T_d,T_d,T_d],\qquad
 E_{d+1}=3E_d+3(3^d)^2.
\tag{6}
\]

The proof sums the three internal edge sets and the three full crossing blocks; the degree equation follows by the same split. At d=2 this is the 36-arc tournament above. A marked vertex has eight other vertices, partitioned into four in-neighbors and four out-neighbors. That is a precise role for eight at this scale, but it is not an eightfold arithmetic recursion or a map from the square-sum path endpoints 8 and 9.

More generally, an orientation of a complete multipartite graph descends to a tournament on its parts **if and only if every crossing block has a uniform direction**. Necessity is the well-definedness of a single quotient arc; sufficiency is immediate. Arbitrary crossing orientations need the full block matrix as a sidecar. Filling internal parts with tournaments instead gives a tournament substitution `Q[H₁,...,H_r]` and adds `sum binom(n_i,2)` arcs. Thus a full tournament on the ten letters would have 45 arcs, not 36.

## 6. Exact doubling, an eightfold operation, and their limits

**M5 (PROVED).** Replace every vertex by q independent twins, preserving adjacency between fibers. On complete multipartite graphs this sends

\[
 K_{n_1,\ldots,n_r}\longmapsto K_{qn_1,\ldots,qn_r},\qquad
 N\longmapsto qN,\quad E\longmapsto q^2E.
\tag{7}
\]

The part quotient and binary rank (4) are preserved. For q=2 all degrees become even; three iterations give q=8 and multiply the edge count by 64. This is an actual eightfold vertex replication. It preserves neither the original parity of each degree nor planarity: the sequence `K1,1 -> K2,2 -> K4,4` has constant binary rank two but becomes nonplanar at its second doubling. Similarly `K1,1,1 -> K2,2,2 -> K4,4,4` has rank two throughout; the middle graph is the planar octahedral graph and the last has a K3,3 subgraph.

When every part size is even, halving those sizes recovers a unique graph **up to isomorphism within this specified family**. A part of size two has a unique pairing. A part of size `2m≥4` has `(2m−1)!!` pairings, permuted transitively by its full symmetric automorphism group; none is invariant under that entire group. Thus larger parts need a marked fiber partition to recover a particular vertex-level half. No arithmetic successor relation has been transferred by this size-halving operation.

The exact connection established here is therefore an equivariant edge trade between two partitions, followed by explicitly defined replication or tournament completion. The maps preserve listed structural predicates; they have not been shown to preserve square-sum edges, graceful differences, Hamiltonian-kernel realizability, or the guarded route to a Collatz root. Those obligations remain **OPEN**.

## 7. Reproduction and hostile controls

Run:

```text
python 04-computation/experiments/thirtysix_multipartite_20260925.py
python -O 04-computation/experiments/thirtysix_multipartite_20260925.py
```

Both runs produce the same [deterministic output](thirtysix_multipartite_20260925.out). The checks use explicit exceptions, so optimization cannot remove them. The universe is: all 36-edge multipartite types; 325 additional sorted positive part profiles (2 through 7 parts, each size at most four); all 32,768 labeled simple graphs on six vertices; the ten Sym² letters and all 45 possible pairs; ternary tournament depths 1 through 4; and all 64 quotient orientations on four parts with one crossing-edge hostile mutation per orientation.

There are 651 rank-two graphs in the six-vertex graph universe. Exactly 16 also have every degree odd: six labeled stars K1,5 and ten labeled K3,3 graphs. These are **graph counts**, not the earlier 960 and 720 tournament realizations of Hamiltonian kernels. They independently check the rank-two shape theorem without assuming a part partition at input.

The controls include distinct triangle and degree data for the two 36-edge graphs, the nonsquare numerical sum 2+3 on a virtual-color edge, loss of uniform block orientation after one flip, planar and nonplanar graphs of the same rank and edge count, and the rank/edge behavior of exact doubling and eightfold replication. Independent proof audit by the signed-clock lane checked M1–M3, including all parts being nonempty and the complete six-edge trade.
