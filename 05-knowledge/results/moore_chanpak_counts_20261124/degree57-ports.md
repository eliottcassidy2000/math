# What a degree-57 Moore input would and would not buy

**Conditional** on the existence of a diameter-two Moore graph of degree 57. Graph Hamilton paths, tournament Hamilton paths and counts of roots of a tiling class are different objects.

## From Moore stars to variable contacts

For any diameter-two Moore graph \(G\) of degree \(d\) and order \(n=d^2+1\), its graph edges and the complete graphs on open neighborhoods \(N_G(v)\) partition all unordered pairs. A graph edge cannot be in a neighborhood (no triangle); a nonedge is in the neighborhood of its unique common neighbor. Consequently a tournament on those labels is specified by orientations of graph edges and one local tournament at each star. This unrestricted representation is tautologically universal; the point of a graph input is to derive useful restrictions on *which* local orientations to choose. The [nonlinear-pencil experiments](../moore57_nonlinear_tournament_pencils_20261124/) concern a different, stronger projective construction assumption.

There is nontrivial count variation even if every star is oriented transitively. Choose independent uniform orders at all stars and fair orientations of the graph edges. For a fixed candidate Hamilton word, group consecutive transitions in each common star into maximal runs of \(\ell\) edges. Treat each graph edge as a singleton source. The probability of the word is the product of \(1/(\ell+1)!\) over the runs, because runs in the same star use disjoint subsets of a uniform order. The number \(K=\sum(\ell-1)\) of same-star consecutive pairs in a uniformly random word has expectation \((d-1)(d-2)/d\), by counting ordered triples in stars. Relative to the independent-edge probability \(2^{-(n-1)}\), the product is between \(d^{-K}\) and 1, since \(\ell\le d-1\). Jensen's inequality gives
\[
 \frac{n!}{2^{n-1}}\,d^{-(d-1)(d-2)/d} \ \le \ \mathbb E H \ \le \ \frac{n!}{2^{n-1}}.
\]
At degree 57 the lower bound exceeds \(10^{8931}\), whereas a shared global order has \(H=1\) in this same locally transitive model. This says nothing about realizing all intermediate integers. It gives one quantitative reason a formula depending on degree alone cannot describe the spectrum of the encodings.

## Two-contact caps as readout ports

Call \(S\) a cap if \(|S\cap N_G(v)|\le2\) for all \(v\). Make \(S\) contiguous in a global order. Put any tournament on \(S\), and orient all other arrows forward in that order. Each local Moore star remains transitive, since its cap subsequence has at most two vertices and any reversal of that pair can be absorbed in the star's order. The outside vertices are ordinal padding, so the count is exactly the chosen tournament's \(H\). This alone does not realize a new height: it identifies a protected interface, not an arithmetic compiler.

Write \(m=|S|\), \(t_v=|N_G(v)\cap S|\). Counting pairs gives
\[
 \binom m2-e_G(S)=\sum_v\binom{t_v}2 \le dm/2,\qquad e_G(S)\le m .
\]
Hence \(m\le d+3\). A [deterministic replay](moore-star-cap-history.py) gives equality in all known cases: a 5-cycle, a 6-cycle, and two 5-cycles for degrees 2, 3 and 7 respectively. At equality both inequalities must be sharp: \(S\) induces a union of cycles, and every vertex sees either zero or two cap vertices.

With \(x\) the indicator of \(S\) and \(y\) the indicator of the vertices seeing two, \(Ax=2y\), and the Moore identity gives \(Ay=(d-1)x/2+(m/2)\mathbf 1-y\). In particular a degree-57 cap attaining 60 would force the following exact equitable partition into \(S,Y\setminus S,V\setminus Y\):
\[
 (60,1650,1540),\qquad
 \begin{pmatrix}2&55&0\\2&27&28\\0&30&27\end{pmatrix}.
\]
This extrapolation is a testable optional hypothesis, not a condition imposed by Moore equality. A search for a maximal cap can begin in the binary adjacency null code: \(x\) is a weight-60 word whose incident counts are all \(\le2\) and with induced degree two. More generally every nonzero null word has support of even size at least 58. Indeed even neighbor counts give \(\sum_v \binom{t_v}2\ge dm/2\). The code dimension is 1520: over \(\mathbb F_2\), \(A(A+I)=J\) and \((A+I)J=0\) give zero semisimple, and the characteristic polynomial is \(x^{1520}(x+1)^{1730}\), by reducing the Moore spectrum \(57,7,-8\) with multiplicities \(1,1729,1520\).

Finally, the known three Moore graphs are vertex-transitive. The [Mačaj–Širáň bound](https://doi.org/10.1016/j.laa.2009.07.018) implies that degree 57 would have at least nine vertex orbits. Uniformity of the root response in known symmetric encodings is not a sound assumption at 57. Different root classes could instead supply different correction channels, provided we can glue them with **nonuniform** contacts beyond fixed-template substitution.
