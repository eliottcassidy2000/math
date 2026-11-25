# The missing Moore graph must contain the first 29 square-sum labels

**PROVED, conditional on existence at degree 57.** Let \(Q_n\) have vertices \(1,\ldots,n\), with distinct vertices adjacent if their sum is a square. The largest \(n\) for which \(Q_n\) occurs as an **induced subgraph** of a diameter-two Moore graph is as follows.

| Moore graph | Pentagon | Petersen | Hoffman–Singleton | Every hypothetical degree-57 graph |
|---|---:|---:|---:|---:|
| Degree | 2 | 3 | 7 | 57 |
| Largest square-sum prefix | 3 | 5 | 14 | **29** |

Every vertex of the last graph can carry label 29. There are at least **13,012,905,007,307,981,605,939,200,000 distinct induced copies** of \(Q_{29}\). Restriction gives \(Q_{25}\), including all its square-sum edges and Hamilton paths. At 30 the three labels 6, 19, 30 form a triangle with sums 25, 36, 49. Thus even an ordinary subgraph embedding of \(Q_{30}\) is impossible. The 32-cycle by itself has a different fate below.

The inherited [root branch matchings](../moore_mobius_prime_packets_20261124/README.md) and [exterior Gram equations](../moore_mobius_prime_packets_20261124/RESIDUAL.md) keep the unique common neighbor of a nonadjacent pair. The obstruction at \(Q_{30}\) is the hostile example for transporting the entire square-sum graph. Our construction uses an internal neighbor as a freely selected new centre and its compulsory two-step midpoints as the new labels. A single-centre extension below can depend on several prescribed square sums simultaneously. The [published branch model](https://doi.org/10.3390/axioms15050332) is also available for direct adjacency interpretation.

## The three small cutoffs

For an induced pattern \(F\) with adjacency \(A\) and exterior contact matrix \(N\) (exterior vertices index its rows), the Moore equation gives
\[
              N^tN=(d-1)I+J-A-A^2 .
\]
A negative determinant rules out an embedding. For the pairs \((d,n)=(2,4),(3,6),(7,15)\) the right hand side for \(Q_n\) has determinant \(-1,-32,-931724928\), respectively. The three preceding prefixes have explicit embeddings, listed in label order.

* In the cyclically numbered pentagon use \(0,3,1\).
* In the Petersen model of disjoint two-subsets of \(\{0,1,2,3,4\}\), use \(01,02,23,03,12\).
* In the Hoffman–Singleton model use \(0,24,1,5,43,2,33,35,7,32,21,9,36,20\). Here \(5i+j\) labels a pentagon site and \(25+5i+j\) a pentagram site, with the second indices stepping by one and two respectively. The cross neighbors of \(5i+j\) are \(25+5k+(ik+j\bmod5)\).

The replay checks every induced pair in these embeddings and computes the three integer determinants. The rank test uses the Moore parameters; ordinary triangle avoidance alone would miss all three cutoffs.

## Extending every embedding of a partially built graph

Fix an induced set \(H\) of size \(h\), with \(e\) edges. Add a centre \(c\), with either no direct contact with \(H\), or one prescribed contact \(a\). For each terminal \(u\) in a prescribed set \(T\subseteq H\), attach a new middle vertex between \(u\) and \(c\). All middle vertices must be distinct. The new induced graph has exactly the listed edges. If \(a\) is specified, assume \(a\notin T\) and no terminal is adjacent to \(a\).

Call a pair in \(H\) unresolved if its endpoints are nonadjacent and have no common neighbor *in* \(H\). Let \(P\) count distinct unresolved unordered pairs meeting \(T\), with each pair counted just once. Let \(U(a)\) count unresolved pairs containing \(a\). The following are lower bounds for the number of centre choices extending **any** embedding:
\[
 \begin{array}{ll}
 d^2+1-(d+1)h+2e-(d-2)P, & \text{no direct contact},\\
 d-\deg_H(a)-U(a)-P, & \text{direct contact }a .
 \end{array}
\]

For the first bound discard \(H\) and its exterior neighbors, at most \((d+1)h-2e\) vertices. For the second start among the \(d\) neighbors of \(a\). Remove its old neighbors, then for each other old vertex that could touch a new candidate remove its unique common neighbor with \(a\). Resolved pairs have their midpoint inside \(H\), already removed, and an adjacent pair has none. These steps enforce the exact direct contacts. For each unresolved pair \(u,z\) meeting \(T\), let \(w\) be its unique common neighbor in the ambient graph and forbid \(c\sim w\). In the first case among the \(d\) neighbors of \(w\), the two old endpoints have already been removed. In the anchored case \(w\ne a\), since no terminal touches \(a\); at most one neighbor of \(a\) is also a neighbor of \(w\). Summing these losses gives the bounds.

For each surviving centre take the unique common neighbor of \(c,u\). It is new, since no old contact of \(c\) touches a terminal. If this middle touched another old vertex \(z\), it would be the forbidden witness of \(u,z\). In particular different middle vertices cannot coincide. They are pairwise nonadjacent since they share their neighbor \(c\). This proves the extension including its induced nonedges.

Here is the construction of \(Q_{29}\), written in **removal order**. Each row removes its centre and its listed middle vertices; running upwards restores the graph. Its contacts and terminals are read directly from the square-sum labels.

| Centre | Middle vertices | Minimum extensions when restoring |
|---:|---|---:|
|17|19|16|
|9|16,27|4|
|10|6,26|17|
|21|15,28|20|
|5|4,20|28|
|13|12,23|39|
|14|2,11,22|2193|
|18|7|2635|
|25|24|2967|
|3|1|3079|
|8|none|3192|
|29|none|3250|

The [replay](moore_chart.py) obtains the lower bounds from the displayed formulas. Multiplication counts labelled embeddings injectively since the centres and their unique midpoints recover each stage. Starting from a common color, refinement by current color and neighbor-color multiset gives respectively 4, 12, 28, 29 classes on \(Q_{29}\). Its automorphism group is trivial, so the same lower bound counts distinct image sets.

## Nine independent internal labels force every prefix from 18 through 29

**PROVED.** In any labelled copy of any \(Q_n\) for \(18\le n\le29\), the same nine images
\[
                 1,2,4,10,11,13,16,17,18
\]
determine all other images using only unique common neighbors. They form an independent set. Restore \(15,14,12,7,9,8,3,6,5\) in that order. Each has exactly two neighbors among the labels already restored. This gives the first 18 labels. Now restore 19 through \(n\) in increasing order. For \(19\le j\le24\) its two parents are \(25-j,36-j\); for \(25\le j\le29\) they are \(36-j,49-j\). In every valid embedding the forced images are distinct and have the required induced edges. For a proposed seed, a collision or an unwanted edge rejects the proposal. The preceding positive construction ensures such seeds actually occur. At 30 the parents 6 and 19 already share an edge, so no possible seed can pass.

Nine is optimal for this internal two-neighbor rule at *each* of these twelve sizes. There are 18 edges at size 18, and every following label adds exactly two, so \(Q_n\) has \(2(n-9)\) edges. Growing from \(k\) seeds using at least two earlier incident edges per new vertex consumes at least \(2(n-k)\) distinct edges. Hence \(k\ge9\). In particular each valid labelled \(Q_{25}\) has at most one extension to a labelled \(Q_{29}\) inside a fixed host. The existence of many extendible copies follows from the lower bound; arbitrary copies could fail a new midpoint test.

Those same nine coordinates, with access to the *ambient* midpoint operation, generate all 3250 vertices. Indeed the resulting closure contains \(Q_{15}\), which the three small cutoffs exclude from every smaller Moore graph. It contains a pentagon and hence cannot close to a star. The classification of closures leaves only the whole Moore graph. The nine-site optimality concerns the internal rule on the prescribed square-sum labels; the [ambient percolation problem](TRIPLE-PORT.md) has minimum three or four.

## Induced polygons up to 33 at every root

For a root of degree \(d\), the \(d\) sets of \(d-1\) leaves hanging below its neighbors are independent, with a perfect matching between any two sets. Take a leaf path \(x_1,\ldots,x_j\) that is induced and uses distinct branches. Among the neighbors of \(x_j\), the unused branches provide \(d-j\) candidates. The adjacent \(x_{j-1}\) can give no chord, and the unique common neighbor of \(x_{j-2},x_j\) is the already used \(x_{j-1}\). At most one candidate is disallowed for each of the \(\max(0,j-3)\) more distant earlier vertices. Thus there are at least \(d-j-\max(0,j-3)\) continuations. With two or more leaves the root and the two terminal branch heads close an induced polygon of length \(j+3\), since neither terminal branch reappears. For \(d=57\) this gives every length **5 through 33** through each root, even after prescribing the first branch and its first leaf. In particular the graph contains induced 13-, 25- and 32-cycles; retaining every edge of the 32-label square-sum graph is impossible.

## A construction test after deleting the square-sum pattern

The internal adjacency also determines the entire exterior Jordan form over the field with thirteen elements. This includes data invisible to the real spectrum. Write \(A\) for the adjacency of the deleted induced pattern \(F\), \(B\) for adjacency on its exterior, and \(N\) for exterior-to-pattern incidence. Use \(n=|F|\), \(U=[\mathbf1,N]\) and
\[
    K=\begin{pmatrix}57&\mathbf1^t\\-\mathbf1&-I-A\end{pmatrix},
    \qquad q(t)=t^2+t-56=(t-7)(t+8).
\]
The Moore equation gives \(BU=UK\) and \(q(B)=U\,\mathrm{diag}(1,-I_n)\,U^t\).

**Excision rule.** Suppose \(U\) has independent columns modulo 13. Compute the Jordan blocks of the small matrix \(K\) over the algebraic closure of \(\mathbf F_{13}\). At the two roots \(7,5\) of \(q\), **increase each block length by one**. Copy all other blocks unchanged. Fill the remaining algebraic multiplicities at 7 and 5 with singleton blocks. The characteristic polynomial supplying those multiplicities, even over the integers, is
\[
 \det(tI-B)=(t-7)^{1729-n}(t+8)^{1520-n}\det(tI-K).
\]

Indeed the full-rank factorization of \(q(B)\) shows that its image is precisely the column space of \(U\), on which \(B\) acts by \(K\). On any primary component away from 7 and 5 this polynomial is invertible. At either root \(b\) its image is the image of \(B-bI\), since the other factor is invertible. Acting on that image shortens each nontrivial Jordan block by one and kills the singletons. This proves the rule in the stated direction. For the characteristic formula write \(D\) for full adjacency. The Moore identity gives \((tI-D)^{-1}=((t+1)I+D+J/(t-57))/q(t)\). Jacobi's complementary minor formula, with \(\det(tI-D)=(t-57)(t-7)^{1729}(t+8)^{1520}\), reduces the claim to the scalar Schur complement of \(t-57\) in \(tI-K\).

For every induced \(F\) of order at most 55, each pattern vertex has an exterior neighbor incident to that vertex alone. Indeed at a vertex of pattern-degree \(h\), there are \(57-h\) exterior contacts, and a contact incident to a second pattern vertex consumes a distinct nonneighbor among at most \(n-1-h\) such vertices. There is also an exterior vertex with no pattern contact, since \(58n<3250\). Thus \(U\) contains the rows \((1,0),(1,e_i)\) for all \(i\) and has full column rank in every characteristic. This verifies the rank hypothesis in all the following cases.

| Deleted pattern | All nonsingleton Jordan blocks of exterior adjacency modulo 13 (eigenvalue, size) |
|---|---|
| The 13-cycle | \(({-3},7),({-3},6),(5,2)\) |
| The first 13 square-sum labels | \((7,2),(8,2)\) |
| The 25-cycle | None |
| The first 25 square-sum labels | None |
| The first 29 square-sum labels | None |
| The 32-cycle | \((7,2)\) |

Here is a derivation of the long blocks. If \(F=C_{13}\), identify its coordinate space with \(R=\mathbf F_{13}[y]/(y^{13})\), using the cyclic shift \(S=1+y\). The coordinate sum is evaluation at \(y=0\), and the all-ones vector is \(u=y^{12}\). Then \(A=S+S^{-1}\) and \(L=A-2=y^2/(1+y)\). On the small coordinates the map \(K+3I\) takes
\[
  (a,x)\longmapsto (8a+x(0),-au-Lx).
\]
Its invariant hyperplane \(8a+x(0)=0\) carries multiplication by \(-y^2/(1+y)+y^{12}/8\), a unit times \(y^2\). The nilpotent ranks there give precisely blocks 7 and 6. The remaining one dimension has \(K\)-eigenvalue 5; the excision rule extends that block to size two in the exterior. Thus \(B\) has 1716 singleton blocks at 7 and 1506 singleton blocks at 5 in addition to the three displayed blocks. A cyclic shift was chosen only on the polygon; the ambient graph need have no order-13 automorphism.

For \(C_{25},C_{32}\), the characteristic does not divide the length. The nonconstant Fourier directions diagonalize in the algebraic closure, with eigenvalues \(-1-z-z^{-1}\) of \(K\) for each nontrivial \(z^n=1\). A collision with 5 or 7 requires \(z+z^{-1}\) to equal 7 or 5, respectively. Reduction of powers via \(z^2=tz-1\) gives orders 7 and 14, neither dividing 25 or 32. The constant sector has characteristic \((t-57)(t+3)+n\); at size 25 its roots modulo 13 are \(3,-1\), and at 32 they are 7 and 8. Only the root 7 at size 32 acquires a new block of length two. For the literal patterns \(Q_{13},Q_{25},Q_{29}\) the [exact replay](moore_chart.py) evaluates the characteristic factors of \(K\) and the nullities of their successive powers over \(\mathbf F_{13}\). At \(Q_{25}\) the nonnative irreducibles have degrees 4, 4 and 18, and at \(Q_{29}\) degrees 1, 3 and 26, all with multiplicity one.

This gives a test for a proposed exterior before restoring deleted vertices. For example the polynomial of the 32-cycle exterior allows two singleton eigenvectors where the necessary Jordan form prescribes one chain of length two. The internal square-sum pattern forces different answers from the polygon carrying its Hamilton cycle. Passing the test establishes this invariant; completion still requires the binary Moore adjacency equation. No unrestricted tournament path count is inferred from the graph cycles.
