# Optimal detours at 55, integral margins, and a conditional Petersen spectrum

Continuation of [55, 13, 33](55-13-33-SYNCHRONIZED-REPAIRS.md) and [the nested-core constraint](NESTED-CORES-AND-RIGID-FOURTH.md). The numerical board is the rigid fourteen-neighborhood fixture, with the actual \((2,7)\) contact retained at source \((1,1)\). None of the results completes a fifteenth neighborhood or proves existence of the degree-57 graph.

## 1. Optimal two-edit detours, even with a fixed spine

There is no single-arc move from height 55 to 53 or 57, at any order, by the [previous theorem](55-13-33-SYNCHRONIZED-REPAIRS.md). For two moves define the first-step excursion to be \(|H_{\mathrm{mid}}-55|\). The minima, over all orders, are
\[
\begin{array}{c|ccc}
\text{desired exit}& 53 & 57 & \text{both from the same intermediate}\\\hline
\text{excursion}&4&6&6 .
\end{array}
\]
All three minima can be attained while keeping one marked seven-site Hamiltonian spine. Use the previous convention (lexicographic nonconsecutive pairs, bit 1 for reversal):
- For 53, take mask \(81\) of height 55, flip \(\{1,3\}\) to mask \(113\) of height 59, then flip \(\{0,4\}\) to mask \(117\) of height 53.
- For either exit, take mask \(628\) of height 55. Flip \(\{0,2\}\) to mask \(629\) of height 49. The next flip is \(\{2,4\}\) for 53 (mask 117), or \(\{1,4\}\) for 57 (mask 565). This improves the old common setup at 65. The original mask 81 also admits a common setup at 49, but its 57 exit from there reverses a spine edge; the mask 628 construction needs no such change of marking.

For optimality, the existing component classification says height 55 comes either from a sole seven-vertex nonlinear strong component, or two strong components of heights 5 and 11 of total order at least nine, besides singletons. Cross-component reversals never lower height by the positive odd-cycle formula. Therefore a first step to 51 must be internal. The product case cannot change by 4 internally, and the replay finds no such step in the seven-vertex case. The same product divisibility rules out an internal flip to 59. A cross flip to prime height 59 in the product case must merge both nontrivial components (the other factor could not remain), requiring a strong component of order at least nine, contrary to [Busch's lower bound](https://doi.org/10.37236/1141). Here a cross flip over an interval of components with a nontrivial endpoint merges those involved components: each endpoint can be reached internally, and the new return closes a route through every intermediate component. In the sole nonlinear case, the only remaining cross possibility for 59 is to merge the old seven-site core with exactly one adjacent singleton, since again a strong height-59 component has order at most eight. The replay checks all those cross arcs and finds none.

It follows that **every 55-to-59 first edit acts in the seven-site active block**. Within 220 fixed-spine representatives the replay tests all 21 arcs and finds 330 such first edits. Each has exactly **one** next arc to height 53, and none to 57. In arbitrary singleton padding any second edit down from 59 must again be internal, so both conclusions transfer. A first step that stays at 55 cannot help, by the all-order unit barrier. Oddness and the absence of the neighbours 53 and 57 now give the stated lower bounds. This is an optimality result for two edits, not unrestricted paths of arbitrary length.

## 2. The marked board passes every first-order margin test

Regenerate the options using the [fourth-site board replay](nested_core_check.py). Force the single marker contact \((2,7)\), remove its branch row, and remove candidates hitting any of its three occupied slots. The residual selection has **52 margins in each of four parts** (one exterior branch part, three central coordinates), with **1,741 candidate columns**. Each column has one 1 in each part. A physical completion would select a zero-one solution of \(Bx=\mathbf1\).

**Signed integers.** Deleting the last row in each of the three central parts identifies the equal-total lattice in \(\mathbb Z^{208}\) with \(\mathbb Z^{205}\). The replay uses exact integer column subtractions, ordinary unit pivots, and Euclidean reductions where no unit is immediately available. It obtains 205 unit pivots, so the projected columns span \(\mathbb Z^{205}\), and the full column lattice is exactly the equal-total lattice. In particular there is no further congruence obstruction, even at a composite modulus, from the linear margin equations alone.

**Strictly positive rational weights.** Begin all integer weights at \(\lfloor 2^{280}/35\rfloor\); do 200 passes through the 208 rows, successively rescaling incident weights to row total \(2^{280}\) and rounding. The exact maximum remaining row error is 15, while the smallest weight exceeds \(4\cdot10^{82}\). Restrict to the 205 free rows. Choose 205 independent original columns, which exist by the lattice certificate. Each has Euclidean norm at most two. Hadamard bounds every cofactor of their square matrix by \(2^{204}\), and the nonzero integer determinant has magnitude at least one. Correcting the error using these columns changes any chosen weight by at most \(205\cdot15\cdot 2^{204}\), strictly below the minimum weight. Thus an exact rational solution exists and is strictly positive. The three omitted equations follow from the equal total of all four parts.

**Consequence: all sufficiently large multiplicities are attainable by a nonnegative integer multicover.** Clear denominators of the positive rational solution to get \(Ba=D\mathbf1\) for some positive integer \(D\) and strictly positive integral \(a\). The signed-lattice test supplies integral \(b\) with \(Bb=\mathbf1\). For \(t=qD+r\), \(0\le r<D\), take \(q a+r b\), which is nonnegative when \(q\) is sufficiently large, uniformly over \(r\). This provides genuine eventual integer saturation in the multiplicity, not a claim about the required multiplicity one. At \(t=1\), nonnegative integers would automatically be zero or one, so that small case is exactly the difficult original selection. The new results rule out a linear or additive-modular certificate of failure for the fixed board, not higher-order matching obstructions.

## 3. A 2,700-site shell of any embedded Petersen

Suppose a degree-57 Moore graph contains Petersen. (Any such copy is induced, since the core already has all its internal two-step routes.) There are ten port colors of 54 sites and then 2700 residual sites, each of residual degree 47. Use \(D\) for their adjacency and \(M\) for the residual-to-port incidence. Each residual vertex has one port in each color, and each port has 50 residual neighbours. The blocks of \(M^{\mathsf T}M\) on one color are \(50I\); between distinct colors, they are zero-one matrices of row and column sum 50, because two distinct ports cannot have two common neighbors. The off-diagonal blocks are therefore \(J-W\) for a 4-regular zero-one bipartite matrix \(W\). The operator norm of \(W\) is at most 4 by the row and column bounds. On the 530-dimensional space of port vectors of sum zero within every color,
\[
             M^{\mathsf T}M \succeq(50-9\cdot4)I=14I .
\]
The constant contributes one more rank. The other nine color-constant directions are in the kernel, so \(\operatorname{rank}M=531\).

The radius-two equation gives \(MM^{\mathsf T}=56I+J-D-D^2\), hence \(56-\lambda-\lambda^2\ge14\) for the coupled centered residual modes, or \(-7\le\lambda\le6\). The 530-dimensional image of the centered ports in the residual is invariant under \(D\): it is the image of \(M\) with its constant removed, and \(D\) commutes with \(MM^{\mathsf T}\).

Here is the multiplicity accounting. The three layer constants consume one copy each of \(57,7,-8\), by the quotient with rows \((3,54,0),(1,6,50),(0,10,47)\). On each of the nine nonconstant core modes, its paired port-color mode yields the matrix
\[
             \begin{pmatrix}a&\sqrt{54}\\\sqrt{54}&-1-a\end{pmatrix},
             \qquad a\in\{1,-2\}.
\]
This has trace \(-1\) and determinant \(-56\), consuming nine further copies each of 7 and \(-8\). The joint space of 530 centered ports and 530 coupled residual modes has both eigenvalues of the ambient graph: each spectral projector there has off-diagonal block \(\pm M/15\), so each has rank at least 530, and the two ranks sum to 1060. All residual modes orthogonal to the 530 coupled ones and the residual constant have no port incidence, and retain the ambient eigenvalues. Subtracting from the global multiplicities \(1729,1520\) shows that the spectrum of \(D\) is
\[
    47^1,\quad 7^{1189},\quad (-8)^{980},\quad 530\text{ roots in }[-7,6].
\]
In particular the residual is connected.

**A pentagon window from the coupled packet.** Girth at least five and 47-regularity fix the first four moments of the whole shell. Subtract the three specified eigenvalue terms: the moments of orders 0 through 4 of the remaining packet are
\[
             (530,-530,3710,-9890,53150).
\]
For the fifth moment \(m_5\), sum the two nonnegative polynomials
\[
 (6-x)(x^2+\tfrac{3079}{1129}x-\tfrac{215472}{59837})^2,\qquad
 (x+7)(x^2+x-\tfrac{268}{53})^2
\]
over the packet. Exact rational expansion yields
\[
 -\tfrac{14556970}{53}\le m_5\le-\tfrac{6202311130}{59837}.
\]
The specified eigenvalues contribute \(217215890\) to the fifth trace. A length-five closed walk here must be a pentagon, because any repeated vertex would give a smaller odd closed walk; each pentagon contributes ten walks. Thus this shell has between **21,694,124 and 21,711,223 pentagons**, inclusive. The bounds provide a separate filter on any proposed residual graph in the Petersen branch, and therefore apply also in the nested Hoffman–Singleton branch. They are not assertions that the extremal endpoints are achievable.

## Replay

~~~sh
c++ -O3 -std=c++17 detour55_check.cpp -o /tmp/detour55_check
/tmp/detour55_check
python3 nested_core_check.py --board | python3 marked_margin_certificate.py
python3 shell_moments_check.py
~~~
