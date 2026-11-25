# Joint spectral trace and two nontransitive flanks

This sharpens the [weighted cuts and nested-apex result](WEIGHTED-CUT-UPPER-AND-NESTED-APEX.md) in both directions of the tandem. Neither open problem is settled.

## A further 127 bits off the degree-57 upper bound

**Conditional on existence of the degree-57 Moore graph.** For every complement spine its native tournament and marked mate both satisfy
\[
 2^{18927}<H<\tfrac45\,2^{26894}.
\]
The lower bound is unchanged. In particular neither native endpoint can equal \(E_p\pm1\) for \(p\ge13448\), with \(E_p=2^{p-1}(2^p-1)\).

Retain the previous scaled apex adjacency and its published \(w_i,z_i\) (25 bins of size 130; \(w_a=4,z_a=1\)). For a row bounded by \(w_i\), write it as the integral of its threshold indicators over \(0<t<w_i\). By permanent multilinearity and the binary Brégman bound the permanent is at most the product of the integrals of \((d_i(t)!)^{1/d_i(t)}\), with empty support contributing zero, where \(d_i(t)\) is the threshold degree. The row sums are the integrals of \(d_i(t)\). The replay certifies
\[
 (d!)^{1/d}\le(d+5)/e\quad(1\le d\le3251),\qquad
 (d!)^{1/d}\le(d+3.1)/e\quad(1\le d\le59).
\]
It uses \(e<271828183/10^8\), from the factorial sum through degree 13 plus a geometric bound for the tail, and integer powers. Above threshold 1, a nonapex row has at most 59 entries: 57 backward Moore arcs, at most one backward spine arc in the mate's converse, and at most one apex arc. All ordinary forward arcs have weight at most one. The row integral is therefore at most \((r_i+3.1w_i+1.9)/e\), with \((r_a+5w_a)/e\) at the apex. Apply the previous \(z\)-weighted AM–GM to these row bounds.

Let \(G\) be the symmetric matrix with zero diagonal and previous flip cost \(g(i,j)=z_jw_j/w_i-z_iw_i/w_j\) at both entries for \(i<j\). Let \(B\) have columns the normalized bin indicators, \(A\) denote the Moore adjacency, and \(P\) project onto the orthogonal complement of the bin constants. The total Moore-edge cost (not counting the extra spine flips in the mate) is exactly
\[
 \tfrac{130}{2}\operatorname{tr}\big[B^{\mathsf T}(A+I/2)B\,G\big].
\]
The bin compression acts by 57.5 on constants and has norm at most 7.5 on their orthogonal complement. Thus this cost is at most
\[
 \tfrac{130}{2}(57.5S+7.5\|PGP\|_*),\qquad
 S=\tfrac1{25}\sum_{ij}G_{ij}.
\]
Indeed in an eigenbasis of \(PGP\), each diagonal entry of the centered compression has magnitude at most 7.5, and there are no mixed terms between constants and their complement. This bounds all bin costs at once, rather than optimizing each rectangle independently.

The certificate finds exactly one positive centered eigenvalue, strictly less than \(4.6\). Scale \(G\) by \(D=1000\operatorname{lcm}_i(1000w_i)\) and restrict the form to the basis \(e_i-e_{24}\), \(0\le i<24\). Exact fraction-free elimination gives one positive LDL pivot for this form \(X\) and all positive pivots for \(46D(I+J)-10X\), where \(I+J\) is the Gram matrix of the basis. The trace of \(PGP\) is \(-S\), hence its nuclear norm is at most \(S+9.2\). Add \(\sum_i g(i,i+1)\) for the mate's extra cross-bin spine arcs. The two upward-rounded means (including division by \(3251e\)) are \(369770/1000\) and \(369771/1000\); the integer power comparisons after dividing by \(\prod_i z_i\) give the claimed strict upper bound.

## Two nontransitive flanks yield a new literal seam

Keep \(e_U=1\) on the nested-left block \(U\) but now allow a one-apex tournament \(V\) on the right. If \(h_V\equiv1\bmod h_U\), choose outer outcontacts in \(V\) of total starting-path count \(s_V=h_V-(h_V-1)/h_U\), using the [complete endpoint sums](ENDPOINT-AND-SIGNED-DECK.md). The previous endpoint determinant becomes \(h_V-h_U(h_V-s_V)=1\), giving a literal \(+2\) under the same fork arc. For example \(h_U=5\) from word \(100\) is compatible with every right word \(11\,0^t\), \(t\equiv2\bmod4\). Its count is \(3\cdot2^t-1\): there are two paths with the right inner apex at an end, and the interior apex contributions over the ones \(i\) and zeros \(j>i\) are \(2^{j-i-1}\).

For an explicit new height, use inner words \(1100\) and \(1001110\) left and right, with chain sites first and each inner apex last locally. Put the six-site B fork between them. The outer apex outmasks on the three blocks are \(15,7,47\), respectively. Use backward masks \(11524,11588\) on the fork (local arc \(1,3\)). Then \((h_U,e_U,h_V,s_V)=(11,1,45,41)\), so the determinant is \(45-11\cdot4=1\).

Direct enumeration of the two lanes on each block and a nine-state ordinal join give the 20-vertex twin-prime pair
\[
 16{,}029{,}731\ \longrightarrow\ 16{,}029{,}733 .
\]
The replay independently exhausts both earlier transitive variable-tail charts and the previous nested-left/transitive-right charts for lower queries \(16029729,16029731,16029733\). For the latter, the established length lower bound permits exactly the 65,535 strong words through length 17, and the right-tail bounds give short finite intervals. There are no hits. Thus both endpoints enlarge the known height set of these particular families.

Run:

    node trace-layer-upper.js
    node two-nested-flanks.js
