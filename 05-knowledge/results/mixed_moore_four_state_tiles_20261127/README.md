# Four-state tournament tiles and exact mixed Moore routing

A loopless binary adjacency matrix allows four states on each unordered pair: forward, backward, both, or neither. A tournament uses just the first two. A path-rooted four-state tiling on \(n\) vertices consequently has \(2^{n-1}4^{(n-1)(n-2)/2}\) choices: each consecutive pair needs its forward arc, but may also have the reverse. If no Hamiltonian path exists a matrix cannot be represented by such a path-rooted tiling.

## 1. The exact equation is much more rigid than a tournament

Use ordinary integer multiplication, not Boolean multiplication. For \(i\ne j\), \((A^2+A)_{ij}=1\) means there is exactly one route of length one or two from \(i\) to \(j\). A missing arc must have exactly one intermediate vertex in that direction. An existing arc cannot have any such intermediate vertex. On the diagonal, \(A^2_{ii}\) counts reciprocal neighbors.

**Theorem.** The only binary matrices satisfying \(A^2+A+I=J\) are the \(1\)-vertex zero matrix and a directed triangle, up to relabeling. In fact the diagonal forces both no loops and no reciprocal pairs. Commutation with \(J\) gives common in- and out-degree \(k\), and row sums give \(n=k^2+k+1\). On the zero-sum subspace the polynomial is \(X^2+X+1\), so the two nonreal cube roots of unity have equal multiplicities there. Taking traces gives \(0=k-(n-1)/2\). Together these imply \(k=0,1\). Both cases are realized.

### Let the diagonal count the reciprocal neighbors

A natural variation keeps the off-diagonal identity, allowing \(A^2+A+I=J+T\) with \(T=\operatorname{diag}(t_i)\) and \(t_i=(A^2)_{ii}\) for loopless \(A\).

**Dichotomy.** Every such matrix is either an all-bidirectional star or has common in- and out-degree \(k\) and constant reciprocal degree \(t\). Indeed \(A\) commutes with \(J+T\). The \(i,i\) entry gives \(d_i^+=d_i^-=d_i\), and entry \(i,j\) gives
\[
d_i-d_j=(t_i-t_j)A_{ij}.
\]
If degrees differ, both directions must be present. Three degree classes would then give a forbidden two-step detour between vertices already adjacent. With two classes both of size at least two, two vertices of one class have at least two routes through the other. Thus one class is a singleton and every arc across classes is reciprocal. Within the other class there can be no arc, since the center already supplies the route. This is a star. In the equal-degree case the same identity propagates equality of \(t_i\) along every arc, and two-step reachability makes the digraph strongly connected. The constant parameter case is the familiar diameter-two [mixed Moore equation](https://arxiv.org/abs/1712.06308), including its small degenerate cases:
\[
A^2+A+(1-t)I=J,\qquad n=k^2+k+1-t.
\]

## 2. Where \(9/4\) appears literally

Completing the square gives
\[
(A+\tfrac12 I)^2=J+(t-\tfrac34)I.
\]
On the zero-sum subspace the right side is scalar. **It is exactly \(9/4\) when each vertex has three reciprocal neighbors** (\(t=3\)), as in the Petersen graph and in the [18-vertex Bosák graph](https://www.sciencedirect.com/science/article/pii/S0012365X06005863). The two nonprincipal eigenvalues are then \(1,-2\). This is a spectral parameter, not the exponent of matrix multiplication. The [tensor-value exposition of the \(9/4\) upper bound](https://arxiv.org/abs/2611.09716) and the [November 24 preprint claiming the matching lower bound](https://github.com/openai/math/tree/main/preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026) concern arbitrary matrices; nothing in this graph argument depends on that lower-bound claim.

The spectral arithmetic also gives a compact degree-57 filter. Set \(c=\sqrt{4t-3}\). If the two nonprincipal roots are nonrational, trace gives \(k=(n-1)/2\); at \(k=57\) this is out of range. Otherwise \(c\) is an odd positive integer and the nonprincipal multiplicities are
\[
m_\pm=\tfrac12\left(n-1\ \pm\frac{n-1-2k}{c}\right).
\]
Thus \(c\mid k^2-k-t\), equivalently \(c\mid(2k-3)(2k+1)\). With \(k=57\) and \(0\le t\le57\), this leaves \(c=1,3,5,15\):

| reciprocal degree \(t\) | one-way out-degree | order | nonprincipal eigenvalues | multiplicities |
|---:|---:|---:|---|---|
| 1 | 56 | 3306 | \(0,-1\) | 3248, 57 |
| 3 | 54 | 3304 | \(1,-2\) | 2183, 1120 |
| 7 | 50 | 3300 | \(2,-3\) | 1968, 1331 |
| 57 | 0 | 3250 | \(7,-8\) | 1729, 1520 |

This is an arithmetic sieve, not an existence claim. In particular the regular mixed relaxation on **3250 vertices** gives no new freedom: \(k^2+1\le n\le k^2+k+1\) pins \(k=57\), then \(t=57\) makes every arc reciprocal. Only the irregular star escapes if the diagonal is left unconstrained.

A rank check for proposed matrices uses \(B=2A+I\) mod an odd prime dividing \(c\): then \(B^2=4J\). If the prime does not divide \(n\), the zero-sum restriction is square-zero, so \(\operatorname{rank} B\le 1+\lfloor(n-1)/2\rfloor\). If it divides \(n\), \(2k+1=0\) mod that prime, hence \(B^3=0\), \(B^2\) has rank one, and there is exactly one size-three nilpotent block. This yields \(\operatorname{rank}B\le 2+\lfloor(n-3)/2\rfloor\). For example the 3304 case has rank at most 1652 mod 3, and the 3300 case rank at most 1650 mod 5.

## 3. The exact 3306-vertex control cannot be trimmed into 3250

The first row is realized by the Kautz construction: take ordered pairs \(ab\) of different symbols from an alphabet of size 58 and give \(ab\to bc\) when \(b\ne c\). For distinct vertices \(ab,cd\), there is a direct arc if \(b=c\), and otherwise a unique two-step route through \(bc\). Reciprocal neighbors pair \(ab,ba\). This host has a Hamiltonian cycle because an Euler tour of the complete loopless digraph on the 58 symbols visits each ordered pair once as an edge.

**Induced rigidity.** In any Kautz graph of this form, a nonempty induced subgraph having exactly one route of length 1 or 2 between distinct vertices is a singleton, a directed triangle on three symbols, or a Kautz graph on a smaller alphabet (including the whole graph). The irregular star case can contain at most one leaf, since the host pairs each vertex with just one reciprocal neighbor. If no reciprocal pairs remain, use the exact theorem above. Otherwise reciprocal degree is one and membership is closed under reversal. View the selected word pairs as edges of an undirected alphabet graph. If it contains \(ab\) and \(ac\) with \(b\ne c\), the selected vertices \(ab\) and \(ca\) need the intermediate \(bc\), so the alphabet graph contains \(bc\). Thus every neighborhood is a clique and each connected component is a clique. Two components cannot support selected vertices in the same strongly connected subgraph, proving the classification. At alphabet size 58, the largest proper perfect induced subgraph is on \(57\cdot56=3192\) vertices, so deleting 56 vertices cannot produce the 3250 target.

## 4. Verification is cheaper than general matrix multiplication

For dense binary input, read all entries into adjacency lists and check degrees and the diagonal parameter. On a regular input passing the order gate, \(k^2\le n-1\). Enumerate direct and two-step routes from each vertex, flagging duplicates and missing destinations and counting backtracks separately. The cost is \(O(n^2+nk^2)=O(n^2)\); no general matrix product is necessary. The unconstrained-diagonal test first classifies any nonregular input as a star or rejects it, in the same time. Thus even the numerical coincidence at \(t=3\) has no algorithmic implication for arbitrary dense products.

Run the JavaScript replay for all loopless labelled matrices through order four, all nonempty induced subsets of the 20-vertex Kautz control, and all ordered-pair entries of the 3306-vertex example. The classification and arithmetic do not depend on the finite checks.
