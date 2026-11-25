# A uniform Hamilton-count window for native degree-57 Moore charts

**Conditional on the existence of a degree-57 Moore graph.** This concerns the native \(3250\)-site chart and its marked tile complement, not arbitrary tournaments at that order or larger arithmetic payload encodings. It substantially narrows the [earlier fixed-order bounds](README.md) in both directions.

## Main result

Let \(G\) be a Moore graph of degree 57 and diameter 2 on \(n=3250\) sites. Choose **any** Hamiltonian path of its undirected complement, use it as the displayed order, and orient each pair forward except the edges of \(G\). The native tournament is \(T\); the marked tile mate is \(T^\ast\), flipping all arcs not on that spine. Such a spine exists conditional on \(G\) by Dirac's theorem: \(\overline G\) has minimum degree 3192. Both counts are odd, and
\[
 2^{12734}<\frac{1513!}{177!}\ \le\ H(T),H(T^\ast),
 \qquad
 H(T)<2^{28860},\qquad H(T^\ast)<2^{28863}.
\]
The factorial lower bound is from the Moore spectrum and actual compatible paths, not a heuristic or an approximation to a count. The refined upper bounds are over 1125 and 1122 bits, respectively, below the uniform labelled mean \(3250!/2^{3249}\). They apply uniformly over all graphs with the Moore parameters and all complement spines.

In particular neither native endpoint can directly realize \(E_{61}\pm1\), or even \(E_p\pm1\) for \(2\le p\le6367\) where \(E_p=2^{p-1}(2^p-1)\). The native chart would be useful as a construction or transfer structure, not as a direct witness for targets in those bands. The result supplies neither an attainable height nor a converse to the interval.

## 1. Balanced matching as a Hamilton-path multiplier

The Moore identity \(A^2+A=56I+J\) gives eigenvalue 57 on constants and eigenvalues 7 and \(-8\) off constants. Put \(L,R\) equal to the first and last 1625 positions. The crossing \(G\)-edges have a matching of size at least **1571**.

Indeed if their maximum matching is \(m\), its minimum bipartite vertex cover has size \(m\) by König. Remove that cover, leaving cross-anticomplete \(X\subset L,Y\subset R\), with \(x+y=n-m\) and each size at most 1625. The shifted operator \(A+\tfrac12 I\) gives top eigenvalue \(57.5\) and nontrivial norm \(7.5\). Because the sets are disjoint, shifting leaves their cross-edge count zero. Centering their indicators gives
\[
 x y\ \le\ \left(\frac{7.5}{57.5}\right)^2 (n-x)(n-y)
 =\frac9{529}(nm+xy).
\]
On the other hand \(xy\ge1625(1625-m)\). Combining yields \(m\ge 1625\cdot 520/538>1570\).

Work in \(T\) or in the **converse** of \(T^\ast\), which has the same \(H\). The latter has exactly the original backward pairs plus the spine's consecutive pairs, since the marked complement does not flip spine arcs. In both charts every matched pair supplies a backward step \(R_i\to L_i\). Regard these \(m\) pairs as ordered two-vertex blocks. A transition from block \(i\) to block \(j\) uses \(L_i\to R_j\). Across the cut, at most \(D=58\) transitions are forbidden per starting or ending site: at most 57 graph neighbors and at most one consecutive spine edge. All unmatched \(L\)-sites can be put into one fixed Hamiltonian-path prefix using Rédei on their induced tournament, and all unmatched \(R\)-sites into one fixed suffix. Either is omitted when empty. Their endpoints impose at most \(D\) bad choices at their respective joins for the same reason.

Order the first \(m-r\) blocks greedily with \(r=4D+3=235\) held back, leaving at least \(i-D\) choices when \(i>r\) blocks remain. In the remaining blocked positions use the undirected graph joining two blocks when **both** transitions are permitted. It has minimum degree at least \(r-1-2D=(r+1)/2\). There is a Hamiltonian path between any two prescribed distinct blocks in this undirected graph: add a vertex adjacent just to the prescribed endpoints; in the enlarged graph of order \(r+1\), each nonadjacent pair of original vertices has degree sum at least \(r+1\), so can be joined using the Hamiltonian-cycle closure lemma. After these edges are added there is a Hamiltonian cycle (the original vertices form a clique), and closure preserves existence. For completeness, the closure lemma follows by opening a cycle at a newly added edge \(uv\); on the resulting spanning path starting at \(u\) and ending at \(v\), a pair of consecutive sites can be used to splice it into a cycle without \(uv\), as the two endpoint degrees sum to at least the order. The first and last reserved blocks can be chosen distinct to fit the two exterior joins since each has at least \(r-D\) eligible blocks. Complete each prefix in one such way. Distinct greedy prefixes give distinct actual directed Hamilton paths, so
\[
 H\ge \prod_{i=r+1}^{m}(i-D)\ge 1513!/177!.
\]
Unlike a count of block covers alone, this accounts for all unmatched vertices and every join.

## 2. Concave score majorization gives a stronger upper bound

Adjoin a bidirectional apex joined to everyone. Round trips through that apex identify Hamilton paths with single-cycle permutations on the augmented graph. The permanent counts all cycle covers, so Brégman's row bound is valid. The apex has row sum \(n\); let \(r_i\) denote the augmented row sums of the original sites in spine order \(i=0,\ldots,n-1\).

For any initial segment of size \(s\), the \(r_i\) sum on it is the corresponding sum for \(n-i\), minus the number of backward arcs across its cut. By the Moore spectrum the graph cut has size at most \(65s(n-s)/n\): center a set indicator and use the eigenvalue bound \(-8\) for the induced quadratic form. In the converse of \(T^\ast\) there is exactly one additional backward arc at each proper prefix cut. This is bounded by \((n/(n-1))s(n-s)/n\). Thus use \(c=65\) or \(c=65+n/(n-1)\), respectively.

Define the decreasing reference sequence
\[
 q_i=n-i-c(n-2i-1)/n.
\]
Its initial segment sums are exactly the transitive sums minus \(cs(n-s)/n\), and its full sum is that of \(r\). Hence the decreasing rearrangement of \(r\) majorizes \(q\): its first \(s\) entries sum to at least the \(s\) entries in the given spine order. The piecewise-linear extension of \(f(k)=\log(k!)/k\) is concave. To see this, set \(\Delta_k=(\log(k+1)-f(k))/(k+1)\) and observe that
\(\log(k+1)-f(k)\ge\frac1k\sum_{i=1}^k(1-i/(k+1))=1/2\).
Then \(\log((k+2)/(k+1))<1/(k+1)\le2\Delta_k\), which gives \(\Delta_{k+1}\le\Delta_k\). Majorization and Brégman now give
\[
 \log H\le f(n)+\sum_i f(q_i).
\]
For the sparse and mate cases these expressions in base 2 are less than \(28859.535\) and \(28862.547\). The replay establishes the strict inequalities with scaled integer logarithm enclosures via the atanh series, not a floating-point cutoff. The corresponding lower logarithm of the factorial bound is greater than \(12734.456\).

## 3. What arithmetic is and is not forced

Both counts are odd by Rédei. The exact triangle shift between endpoints is **3134**: in passing from \(T\) to the converse of \(T^\ast\), all internal scores are unchanged, the first vertex loses one win and the last gains one. Original first and last scores are \(n-1-57=3192\) and \(57\); the directed-triangle identity \(c_3=\binom n3-\sum_i\binom{d_i^+}{2}\) gives \(3192-57-1=3134\). Taking a converse does not change the triangle count.

In particular triangles have the same parity at these two endpoints. This is **not** a proof that the endpoint heights agree modulo four. The [Grinberg–Stanley odd-cycle formula](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf) needs parity of **all** odd directed cycles. The small exact control in the replay uses the native Petersen complement spine with counts \(5625,5627\): its triangle counts are \(28,30\), so even that equal-parity input with the same type of positive shift allows different height residues mod four. Likewise no extra universal odd-prime divisibility is asserted here. If some unmarked even-order endpoint is self-converse the usual \(H\equiv1\pmod4\) applies to that class, not to all native charts.

### Further search filter: a cube of legal graph-preserving arc moves

A standardized complement order can be chosen with **no graph edge between positions at distance at most 3**. To see this, minimize the number of such bad pairs over orders. For a site in a bad pair, swap it with a position more than three away. At most \(6\cdot57\) candidates have a vertex adjacent in \(G\) to one of the sites within three of the bad position, and at most \(6\cdot57\) positions lie within three of a neighbor of the vertex being moved. Avoid those and the at most seven nearby positions. Since \(3250>691\), the swap can remove an existing bad pair without introducing another. This contradicts a positive minimum.

Now independently swap or do not swap in the 1625 disjoint consecutive pairs of positions. Any newly consecutive vertices were originally within three positions, so every resulting order is still a complement spine. For each such adjacent swap, **one tournament arc** reverses on the common labelled vertex set, and the backward graph decoded from the new order is unchanged. The result is a literal \(1625\)-cube on labelled tournaments carrying the same hypothetical graph, all of whose heights obey the bounds above. Vertices of the cube need not have distinct heights or remain distinct after passing to unlabelled classes, and there is no deduction of \(+2\) on its edges. This separates the now-automatic supply of legal count-changing *candidates* from the unresolved problem of selecting an exact target.

## Reproduce

Run `python3 degree57-height.py`. The script certifies all displayed binary bounds by rational arithmetic and independently checks the Petersen endpoint counts and triangle-parity control by subset DP. It does not enumerate a degree-57 graph.
