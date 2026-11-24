# A three-state integer grammar gives exponentially many targets from one node

**PROVED.** For every \(k\ge1\), there is one rigid self-converse tournament class on \(128(k+1)\) vertices with exactly \(5^{4(k+1)}\) directed Hamilton paths. From this **same** source, blue tiling lines reach at least \(2^k\) distinct ordinary Hamilton-path counts, and black lines reach at least \(2^k-1\) distinct counts. Each constructed count has an explicit matrix formula from its binary address. The two lists of target counts may overlap. This quantifies the information in placement even after the source class, source count, order and line color are fixed.

## Arbitrary ordered panels have minimum linear state dimension three

Use the cut notation of [the fixed-height transfer](TARGET.md). Write \(z(L)=(H,s,E)^T\), where \(s\) counts full paths starting at the last vertex of the partner fragment \(L\), and \(E\) counts two-slot covers with the marked path ending at that vertex. Between panels the source is an ordinal sum in ascending path order; let \(L\circ R\) denote the corresponding joining of partner fragments. Every cross arc points from \(R\) to \(L\) except the consecutive spine arc. Give the right fragment the matching first-vertex features by taking its reflected converse \(R^*\). Then
\[
 H(L\circ R)=z(R^*)^T J z(L),\qquad J=\operatorname{diag}(1,-1,1).
\]
By the sole-exit law the effect of a singleton on the left state is \(z\mapsto Mz\) with
\[
 M=\begin{pmatrix}1&-1&1\\1&-1&0\\1&0&1\end{pmatrix}.
\]
The three measurements with zero, one and two trailing singletons recover the full state: from \((d_0,d_1,d_2)\) recover \(H=d_0,\ e=(d_0-2d_1+d_2)/2,\ B=(d_2-d_0)/2\), where \(e=s,B=E\) in the left convention. Therefore appending **any fixed panel** acts linearly on the state. Indeed the three measurements after that panel are linear functionals of the original state by cutting *before* the panel.

Every reachable state has \(s\equiv E\pmod2\): the single-step height \(H-s+E\) and the original height are both odd. The three singleton strings of sizes 1, 2, 3 give states \((1,1,1),(1,0,2),(3,1,3)\), whose determinant is 2. They therefore span the entire integral lattice with the displayed parity restriction. In this basis **every panel acts by an integer matrix** because it takes all three basis states to reachable states. Use row coordinates in this basis. The empty initial coordinate is \((-1,-1,1)\) and the final height readout is its scalar product with \((1,1,3)\).

To compile a panel \(Q\), measure the partner height \((K_Q)_{ij}\) on the source word consisting of \(i\) singletons, \(Q\) and \(j\) singletons for \(1\le i,j\le3\). For the empty panel write \(K\) for the same matrix. The panel's integer update is
\[
 D_Q=K_Q K^{-1},\qquad K=\begin{pmatrix}1&3&5\\3&5&9\\5&9&17\end{pmatrix}.
\]
The same formula handles the empty initial state by the cut identity. Successive panels act by their ordered matrix product. Since \(\det K=-4\), any linear interface over the rationals evaluating every left/right gluing must have at least three coordinates. In particular, the tensors in the fixed-height inverse factor into matrices of size three irrespective of the number of panels.

## Encoding a binary choice inside one four-vertex isomorphism class

The unique strong four-vertex tournament has height 5 and trivial automorphism group. Three presentations of this very same class use the following backward tiles; coordinates start at zero. The replay computes all 27 needed measurements on literal tournaments.

| Presentation | Backward tiles in its four-vertex source | Matrix acting on row coordinates |
|---|---|---|
| \(A\) | \((0,3)\) | \(\begin{pmatrix}-1&5&3\\0&8&5\\0&17&8\end{pmatrix}\) |
| \(B\) | \((0,2),(0,3)\) | \(\begin{pmatrix}-1&4&4\\-1&7&7\\-2&14&11\end{pmatrix}\) |
| \(D\) | \((0,3),(1,3)\) | \(\begin{pmatrix}4&14&-3\\8&25&-6\\15&48&-12\end{pmatrix}\) |

Initialize by appending \(AB\). For each of \(k\) successive bits append \(AB\) or \(AD\). All half-words now consist of \(2k+2\) copies of the same ordinary component. Nevertheless their exact boundary vectors differ. To decode, use the distinct projective mod-three row directions \(p=(2,1,1)\), \(q=(1,2,0)\). The prefix gives \(-p\). Both \(p,q\) go to nonzero multiples of \(p\) under \(AB\), and to nonzero multiples of \(q\) under \(AD\). The last symbol is visible in the resulting direction. Both symbol matrices are invertible over the rationals, so invert that step on the integer state and repeat. This distinguishes all \(2^k\) terminal row coordinates and hence all \(2^k\) triples \(z_w\).

## One physical spacer distinguishes the integer outputs

Put a run of \(r=112(k+1)\) singletons between a half-word \(w\) and the reflection of a half-word \(v\). Every pairing gives the **same ordinary source**: \(2k+2\) copies of the strong four-class, then \(r\) singletons, then \(2k+2\) more copies. The linear order of the strong components and their rigidity make the source rigid. The order is \(16(k+1)+r\), the height is the component product in the theorem, and reflection is an anti-isomorphism. Its blue lines have \(w=v\) in this family. If \(w\ne v\) the half-word presentations differ at a named tile, so the line is black. The exact partner count of each such source tiling is
\[
                         h_{w,v}=z_v^T J M^r z_w .
\]

Here is why all diagonal readings are distinct, and why varying \(w\) with one fixed \(v\) also gives distinct readings. The matrix satisfies \(M^TJ=JM\). Its characteristic polynomial \(f(x)=x^3-x^2-x-1\) is irreducible over \(\mathbb Q\), with roots \(\lambda>1\), \(\beta,\bar\beta\) of modulus \(\lambda^{-1/2}\) for the latter two. For any root \(\alpha\), the spectral projector is
\[
 P_\alpha=\frac{M^2+(\alpha-1)M+(\alpha^2-\alpha-1)I}{f'(\alpha)}.
\]
The rank-one symmetric matrix \(JP_\lambda\) gives a form \(c\ell(x)\ell(y)\), where \(\ell M=\lambda\ell\), \(c\ne0\). Its linear functional \(\ell\) is injective on rational vectors: if \(x\ne0\) were rational with \(\ell(x)=0\), it would annihilate \(x,Mx,M^2x\), a basis by irreducibility. Thus the diagonal coefficients \(c\ell(z_w)^2\) are distinct because no two different triples coincide up to sign (their first coordinates are positive). At fixed \(v\) the coefficients \(c\ell(z_w)\ell(z_v)\) are likewise distinct.

For an explicit separation bound let \(m=k+1\), \(U=9\cdot75^{2m}\). Every triple has coordinates bounded by \(U\). Indeed the three panel matrices have maximum absolute row sum at most 75, the initial coordinate has norm three in \(\ell^1\), and each singleton basis state has largest coordinate at most three. For either difference of readings under consideration, its coefficient at \(\alpha^r\) is \(N_\alpha/f'(\alpha)\). Here \(N_\alpha\) is an algebraic integer with \(|N_\alpha|\le150U^2\) at all three conjugates: \(\|M\|_\infty\le3\), \(|\alpha|<2\), and each of the two bilinear numerators is at most \(75U^2\). Since \(N_\lambda\ne0\), its integer norm gives
\[
 |N_\lambda/f'(\lambda)|\ge [17(150U^2)^2]^{-1}.
\]
We used \(|f'(\alpha)|\le17\). Also \(|\prod_\alpha f'(\alpha)|=44\), hence \(1/|f'(\alpha)|<7\). The combined contribution from the two smaller roots is at most \(2100 U^2\lambda^{-r/2}\). The dominant contribution exceeds this whenever
\[
                   \lambda^{3r/2}> 2100\cdot17\cdot150^2 U^6.
\]
Our common spacer suffices: \(\lambda^3>4\), \(U<2^{4+14m}\), \(2100\cdot17\cdot150^2<2^{30}\), and \(112m>54+84m\) for \(m\ge2\). The blue diagonal readings therefore give \(2^k\) different integer counts. Holding \(v\) fixed and taking the other \(2^k-1\) words gives that many different black readings.

This proves exponential diversity of genuine line partners from a single merged node, with an exact inverse for its encoded internal placement when the three-state boundary reading is supplied. The construction enumerates a specified family of integer heights. Coverage of all odd targets remains open. Run `python3 -B panel_monoid.py` for the matrices, the two mod-three directions, a direct orientation check and the integer outputs at \(k=1,3,8\).
