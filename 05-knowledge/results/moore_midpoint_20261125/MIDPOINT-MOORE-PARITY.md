# The midpoint 2, Moore parity, and the 4 mod 8 split

25 November 2026. The two ratios from the [quadratic-cycle note](../quadratic_cycle_synchronization_20261125/QUADRATIC-CYCLE-SYNCHRONIZATION.md) and the [matrix multiplication paper](https://mathcommunity.github.io/oai-math-revised/papers/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/paper.pdf) are not the two sides of one universal spectral law. There is, however, a shared structural question: which information can a balanced representation forget without creating spurious solutions? For Moore graphs the two-state representation gives an exact parity certificate. The higher residue \(4\bmod8\) gives a different, impossible everywhere-balanced condition in the 400-coclique branch.

## 1. The two midpoints and the order tag

In the paper's Proposition 8.1 the tripod model has perimeter \(2L\), \(L=a+b+c\), and canonical shared-label rate \(D_*=L-(a^{-1}+b^{-1}+c^{-1})^{-1}\). On positive arms the model ratio \(2L/D_*\) belongs to \((2,9/4]\), since \((a+b+c)(a^{-1}+b^{-1}+c^{-1})\ge9\). It is \(9/4\) at equal arms and tends to \(8/4=2\) when one arm collapses with the other two fixed. Thus 2 is a genuine boundary of this model, and \(7/4\), the quadratic period-three root, is on the *other* side. This calculation is about the canonical rate, not every actual spectral price. Neither rate gives a numerical bound on Moore graphs.

Fix an order of the graph's sites. Reverse exactly the graph edges in the transitive tournament of that order. If \(E\) is the order's skew sign matrix and \(K\) the resulting tournament's skew sign matrix, the exact graph decoder is
\[
 C=E\circ K=J-I-2A,\qquad A=(J-I-C)/2 .
\]
Here \(\circ\) is entrywise multiplication; \(A\) is the symmetric graph adjacency. The order is a shared tag. A different tag read against the same unmarked tournament need not give the same graph. In particular, storing a tournament need not store the original automorphism group: any tournament automorphism group has odd order, because an involution would exchange some pair and reverse its arc. In the pair of matrices, graph relabelling moves *both* matrices. This is the same type of synchronization issue as the parameter tag in the previous session, not a reduction to fast matrix multiplication.

For the user's actual tiling model, choose the order along a Hamiltonian path in the *complement* graph, so none of the constitutive spine arcs is reversed. This is always possible for the Moore degrees under discussion: for (kge3) the complement has minimum degree (k^2-kge n/2) (Dirac), and the complement of a pentagon is a pentagon. Thus the degree-57 graph, if found, can be carried by a valid 3250-site marked tiling with exactly 92,625 intransitive tiles. The tile weight and the parity certificate are separate from the Hamiltonian-path count of that tournament.

On each off-diagonal position the tournament signs are \(\pm1\); mod 2 they are identical, while mod 4 the sign determines the decoded graph bit. Mod 8 a different blind spot appears: on an even number of sites, every tournament's \(K\) has \(\det K\equiv1\bmod8\), since \(\det K\) is a Pfaffian square and \(K\bmod2=J-I\) is invertible. The full decoded adjacency parity is informative; the skew determinant modulo eight alone is not.

## 2. A one-bit certificate of the entire radius-two structure

For a simple \(k\)-regular graph on \(n=k^2+1\) sites, put \(B_{ij}=A_{ij}+(A^2)_{ij}\) for \(i\ne j\). This nonnegative integer is the number of length-one and length-two routes. Every row has sum \(k^2=n-1\) over \(j\ne i\). Consequently
\[
 \big[\forall i\ne j,\ B_{ij}\equiv1\bmod2\big]
 \quad\Longleftrightarrow\quad A^2+A=(k-1)I+J .
\]
Indeed the odd numbers in each row are at least one and their mean is one. The right side is the Moore condition. It gives a graph-to-tournament test over \(\mathbb F_2\) after decoding the signed entries modulo four. If only maximum degree \(k\) was supplied, the oddness still forces full coverage in two steps and the Moore vertex bound is saturated at every root; hence every root has degree \(k\).

There is also a robust version for the regular case. In row \(i\), let \(z_i\) be the number of missing routes (\(B_{ij}=0\)), and let \(e_i\) count parity failures. The total surplus of the entries \(B_{ij}>1\) is \(z_i\), so
\[
 z_i\le e_i\le 2z_i,\qquad
 \sum_{j\ne i}|B_{ij}-1|=2z_i.
\]
Thus the binary residual never hides a larger all-even surplus without also exposing its missing partners.

For \(k=57\), the nontrivial eigenvalues are \(7,-8\), as \(x^2+x-56=(x-7)(x+8)\), and \(57=7\cdot8+1\). Completing the square gives
\[
 (2A+I)^2=225I+4J,\qquad C^2=225I+3024J.
\]
In particular \(C^2\equiv I\bmod16\), whereas the order-decoded \(2A+I\) still has the \(4J\) term modulo eight. On the integral sum-zero lattice the \(+7\) projector is \((A+8I)/15\); it reads as \(-A\bmod8\) and \(A+I\bmod7\). Working on that lattice avoids incorrectly splitting off the all-ones vector at the prime 2.

## 3. The exact obstruction to 4 mod 8

The following branch is conditional on a 400-site independent set \(S\) in a degree-57 Moore graph, not merely on existence of the graph. The Hoffman ratio bound gives \(|S|\le400\), and at equality its characteristic vector satisfies \(A{\bf1}_S=8({\bf1}-{\bf1}_S)\). Thus each of the 2850 outside vertices sees exactly eight points of \(S\). By uniqueness of the common neighbour of two nonadjacent sites, the outside vertices are the blocks of a Steiner \(2\)-\((400,8,1)\) design. The well-known projective-space candidate uses the points and lines of \(\mathrm{PG}(3,7)\); a general candidate need not be that design.

The size 400 itself can be read from the 8-contact condition. More generally, for a nonempty set \(U\) where every vertex has a multiple of \(h\) neighbours in \(U\), put \(t_v=|N(v)\cap U|\). Then
\[
 \sum_v \binom{t_v}{2}=\binom{|U|}{2}-e(U),\qquad
 \binom{t_v}{2}\ge (h-1)t_v/2 .
\]
Since \(\sum_v t_v=k|U|\), one gets \(|U|\ge(h-1)k+1\), with equality forcing independence and \(t_v\in\{0,h\}\). For \(h=8,k=57\) this is exactly 400; the 2850 nonzero contacts then exhaust the outside. For \(h=2,4,8\) the lower bounds are \(58,172,400\).

Write \(N\) for the incidence matrix of \(S\) versus its outside blocks. Then \(NN^\top=56I+J\). For any half \(U\subset S\) with 200 points, let \(t=N^\top{\bf1}_U\) and \(w={\bf1}_U-{\bf1}_S/2\), written in the 400 point coordinates. The exact variance is
\[
          \sum_{\text{outside }v}(t_v-4)^2
          =w^\top NN^\top w
          =56\cdot100=5600.
\]
Because \(0\le t_v\le8\), the class \(4\bmod8\) means exactly 4. Hence **at least 350 outside vertices miss \(4\bmod8\)** in every half, and **at least 184 have a count outside \(\{3,4,5\}\)**. The second follows from \(5600\le2850+15\,\#\{|t_v-4|\ge2\}\). In particular an exact mean of four does not even allow every vertex to stay within one of four. These are lower bounds, not claims of attainability. The same incident bits give both the balanced mean and the unavoidable defect, a more literal analogue of the synchronization lesson than identifying the exponents with eigenvalues.

## 4. Reciprocal spreads and the plane trade

Put \(H\) for the graph induced on the 2850 blocks. In the ordered decomposition \(S,\bar S\), the original graph is \(\left(\begin{smallmatrix}0&N\\N^\top&H\end{smallmatrix}\right)\), with \(H\) also denoting adjacency. The cross part of the Moore equation is
\[
                 N(H+I)=J .
\]
Consequently each block \(b\) is in a designated parallel class \(P(b)=\{b\}\cup N_H(b)\) of 50 blocks partitioning all 400 points. The choices must be **reciprocal**: \(b\in P(c)\) iff \(c\in P(b)\). They are distinct; for different \(b,c\), their intersection has size 2 if adjacent, 0 if the original blocks meet, and 1 if the blocks are disjoint and nonadjacent. Conversely, for a Steiner \(2\)-\((400,8,1)\) design, choosing such reciprocal pointed classes and requiring that \(H\) have no triangles or squares suffices. Distinct incident blocks cannot share two points; a triangle through a point or a square with exactly one point would contradict disjointness within a class. The built graph is 57-regular on 3250 sites with no cycles of length 3 or 4, so each radius-two ball exhausts it. Alternatively only the bottom-right parity equation needs testing, since the other two blocks already satisfy the Moore equations.

There is a clean second-level interpretation: a collection of designated classes \(P(b)\) is a resolution of all blocks iff the \(57\) corresponding design blocks \(b\) pairwise meet. The usual dual pencil of the \(57\) blocks through one point supplies one. For any non-pencil linear family of pairwise meeting \(8\)-sets, the degree of any point in the family is at most eight: a block avoiding that point must meet the blocks through it at distinct points. For a block \(b\), summing over its eight points now bounds the size of the family by \(1+8\cdot7=57\). Equality forces a projective plane of order 7: all used points are on eight lines, and the number of point pairs covered forces every pair to be on a line.

Equivalently, the closest possible distinct 400-point independent sets are related by a **57-for-57 projective-plane trade**. If \(S,T\) are maximum independent sets then in the bipartite incidence between \(S\setminus T\) and \(T\setminus S\) each vertex has degree 8, since no vertex of \(S\cap T\) can be its neighbour. If either part has \(m>0\) points then \(m-1\ge8\cdot7\), by counting the distinct points reached in two steps. At equality the induced incidence is precisely that of a projective plane. Conversely any order-seven subplane of the Steiner design exchanges its 57 points inside \(S\) for its 57 line-vertices: any two lines meet in \(S\), so their vertices cannot be adjacent, and none touches \(S\) off the subplane. Thus \(|S\cap T|\le343\), with equality iff such a trade exists.

The small Moore controls follow the same calculation with outside block size \(q\) and \(k=q(q-1)+1\). Petersen has \(q=2\), a 4-coclique and size-three triangle-design swaps. Hoffman–Singleton has \(q=3\), a 15-coclique and size-seven Fano swaps, agreeing with [Hafner's coclique geometry](https://www.math.auckland.ac.nz/~hafner/his.pdf). The replay finds four trades at a Petersen coclique and fifteen at a Hoffman–Singleton coclique, and identifies the latter neighbourhood design with the lines of \(\mathrm{PG}(3,2)\) by its binary point coordinates. The new order-seven trade is the same mechanism at degree 57, not an assertion that an actual trade exists there.

If the degree-57 design really is \(\mathrm{PG}(3,7)\), its 400 planes would produce 400 trades from \(S\). Every designated spread covers the points of any projective plane by exactly one full line and 49 single intersections, since the plane has 57 points and a spread has 50 lines. Hence **one can dualize any such completion**: use planes instead of points for \(S\), and retain the same outside adjacency \(H\). The dual incidence matrix has the same column Gram matrix (two lines meet iff they lie in a common plane, uniquely when distinct), and satisfies the same cross equation. This preserves the Moore identity without needing the completion to preserve any chosen duality. It does not force the resulting graph to be isomorphic to the original.

The whole projective linear group is too symmetric for this completion: it is transitive on skew line pairs, so invariant outside adjacency would select all or none of the 2401 possible skew-line neighbours of a line; the completion needs exactly 49. Both this symmetry breaking and the split variance are meaningful constraints, not proofs of nonexistence. Also, in the [order-11 fixed-pentagon branch](../moore_fixed_order_spectra_20261124/PENTAGON-PORT-ATLAS.md), no 400-coclique can be invariant: it would have to contain \(400\bmod11=4\) fixed vertices, but an independent subset of the fixed pentagon has at most two. In the [free order-13 branch](../moore_fixed_order_spectra_20261124/SINGER-ORDER13-MODULES.md) an invariant 400-coclique is likewise impossible since \(13\nmid400\).

Replay: **python midpoint_moore_check.py**. It tests the parity certificate against genuine and degree-preserving switched small fixtures, the spread/resolution rules and both trade censuses. The degree-57 statements are derived from the displayed identities conditional on their stated hypotheses. The [degree-57 Moore problem remains open](https://www.mdpi.com/2075-1680/15/5/332).
