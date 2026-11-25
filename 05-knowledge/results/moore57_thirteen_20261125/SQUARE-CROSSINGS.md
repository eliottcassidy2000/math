# What the square-sum crossing words require in a Moore graph

A crossing means an edge between two specified vertex parts, and an internal edge stays in one part. The square-sum parts are the two parities. A Moore graph of degree 57, if it exists, need not have an equitable equal bipartition. When it does, each vertex has 32 internal and 25 crossing neighbors. These are the hypotheses of the corrected theorem. [Renteln's corrigendum](https://ajc.maths.uq.edu.au/pdf/79/ajc_v79_p193.pdf) · [Earlier crossing-block and reflection theorem](../moore_fixed_order_spectra_20261124/SQUARE-REFLECTION-CROSSING.md)

## The first two spanning thresholds retain different pentagon information

**CITED.** Gerbicz gives a Hamilton path for every square-sum graph \(Q_n\) from 25 onward and a Hamilton cycle from 32 onward. The [three earlier path onsets and their ports](../../../session211/square-sum-port-bridge/PROOF.md) include the first path at 15. [Gerbicz's statement and source code link](https://oeis.org/A090461)

**FINITE-EXACT.** The following calculation counts simple unoriented three-edge paths, those extendible to a simple pentagon by adding two edges, and the parity words of the resulting pentagons. The words are read up to rotation and reversal, with 0 even and 1 odd.

| Graph | Three-edge paths | Extendible | Pentagon words | Hamilton paths using edge 12–13 |
|---|---:|---:|---|---:|
| \(Q_{25}\) | 99 | 10 | 2 copies of 00111 | 8 of 10 |
| \(Q_{32}\) | 188 | 20 | 4 copies of 00111 | 352 of 392 |

The unique Hamilton cycle of \(Q_{32}\) has 16 crossings and 16 internal edges. The Hamilton paths of \(Q_{25}\) have 10 through 15 crossings, with respective multiplicities \(1,1,2,3,2,1\). In \(Q_{32}\) the range is 10 through 19, with multiplicities in increasing crossing order given by
\[
  (3,5,11,26,55,93,75,86,31,7).
\]
At either size vertex 13 has internal neighbors 3 and 23 and crossing neighbor 12. Thus a single specified port can be tracked alongside the global parity word. The replay enumerates simple paths and recovers cycles by joining adjacent endpoints and dividing by the number of possible deleted cycle edges. [Exact replay](square_cross.py)

For a three-edge path with endpoints \(x,x+d\), a two-edge square-sum closure requires a vertex \(z\) with \(x+z=v^2\), \(x+d+z=u^2\). Hence the *complete* set of potential closing vertices at all later orders comes from \((u-v)(u+v)=d\), \(z=v^2-x\). Keep only positive labels inside the graph and outside the original path. At gap 13 the unique factorization gives \(z=36-x\). At gap 7 it gives \(z=9-x\). At gap 21 it gives \(z=4-x\) or \(100-x\). In the two windows the pairs (number of three-edge paths, number extendible) are, at gaps \(7,13,21\), respectively \((4,1),(3,0),(3,1)\) and \((8,2),(6,4),(5,1)\).

For example the path \(7,9,16,20\) cannot close at 25 but acquires its required vertex 29. The path \(17,19,6,30\) can never close in a square-sum graph at any larger order: its only candidate is the already used 19. A Hamilton cycle at size 32 can therefore coexist with permanent failure of local pentagon completion. The numbers 7 and 21 above are endpoint differences. An implication for forbidden tournament path counts would require an additional map.

## Four observed colors determine a pentagon inventory to one parameter

**PROVED.** In every diameter-two Moore graph of girth five, each directed nonbacktracking three-edge path lies in exactly one simple pentagon. Its endpoints are distinct and nonadjacent by the girth. Their unique common neighbor is outside the given path. Thus deleting one vertex of a rooted, oriented pentagon records **all** directed three-edge paths, once each at a specified deletion position. Conversely, in any simple graph of girth at least five and minimum degree at least two, this returns all three-edge paths precisely when every connected component has diameter at most two. Indeed a three-edge geodesic cannot close in two edges. An unclosed three-edge path is a hostile witness for a proposed complete Moore reconstruction.

Color each vertex by its side. Let \(F(i_0,\ldots,i_4)\) count rooted oriented pentagons with those five colors. Every four-coordinate marginal is a directed three-edge path count, in cyclic order starting after the deleted coordinate. If two oriented inventories have the same marginals, their difference changes sign whenever one bit changes. Thus their difference has the form \(h(-1)^{i_0+\cdots+i_4}\), already rotation and reversal invariant. The only missing linear measurement is the five-way parity coefficient. This proof applies to any binary coloring of a Moore graph. With an equitable equal split its marginals are explicit. Put \(d(i,j)=a\) for equal colors and \(b\) for distinct colors, with \(N\) vertices on either side. Then
\[
 P(ijkl)=N d(i,j)\bigl(d(j,k)-[i=k]\bigr)\bigl(d(k,l)-[j=l]\bigr).
\]
The subtractions exclude an immediate return. Girth excludes any further coincidence.

**CONDITIONAL on the equitable equal split.** For each binary dihedral word the following gives the number of unoriented pentagons in the degree-57 case. The free integer \(t\) is the number entirely in side 1. To reproduce the system, contribute once to the corresponding four-letter path bin for each of the five starting positions in both directions per pentagon. Substituting the table gives every \(P(ijkl)\); the preceding sign argument proves completeness of the linear determination.

| Side word | Number | Complement word | Number |
|---|---:|---|---:|
| 00000 | 6,094,400 − \(t\) | 11111 | \(t\) |
| 00001 | \(5t\) − 5,486,000 | 01111 | 24,986,000 − \(5t\) |
| 00011 | 25,636,000 − \(5t\) | 00111 | \(5t\) − 4,836,000 |
| 00101 | 21,086,000 − \(5t\) | 01011 | \(5t\) − 9,386,000 |

In particular the numbers of monochromatic, one-against-four, adjacent two-against-three and separated two-against-three pentagons are respectively **6,094,400, 19,500,000, 20,800,000 and 11,700,000**. Their sum is 58,094,400. The known Petersen and Hoffman–Singleton splits substitute \((N,a,b)=(5,2,1),(25,2,5)\) in the same equations. Nonnegativity alone forces \(t=1,5\) respectively. The numbers in the four categories are \((2,0,10,0)\) and \((10,0,250,1000)\). Thus degree 57 needs *all four* categories in large numbers even though both earlier splits omit the one-against-four category.

## A cut constraint at every vertex

Fix a physical vertex \(x\), call its side 0, and let \(p_x\) be the count of unoriented pentagons through \(x\) lying entirely in that side. The 16 oriented counts for the next four sides form a four-bit tensor. Every three-coordinate marginal counts a three-edge path with \(x\) at a specified position; each is \(P(ijkl)/N\). For a root at either endpoint this is successive neighbor choice; at an interior root choose the two adjacent edges first, excluding immediate reversal, then the remaining edge. The same formula follows because \(d\) is symmetric. The tensor therefore has one free alternating coordinate, fixed by the all-zero entry \(2p_x\). In lexicographic order of the remaining bits \(w\), the exact counts are \(b_w+2(-1)^{|w|}p_x\), with
\[
\begin{split}
(b_w)=(&0,30752,30752,-5952,30752,-11552,-5952,31552,\\
       &30752,-11552,-11552,25952,-5952,25952,31552,-6752).
\end{split}
\]

**CONDITIONAL.** In any such degree-57 split, every vertex satisfies
\[
                8613\le p_x\le10061 .
\]
Consequently **all 16 prescribed patterns occur at least 5674 times apiece at every physical vertex**. Each half contains between 2,824,575 and 3,269,825 internal pentagons. A fixed-point-free 13-action preserving the split forces each half total to be a multiple of 13 because an order-13 permutation cannot act nontrivially on a pentagon or fix it pointwise.

For the bound let \(B\) be the adjacency matrix of the internal graph in the half of \(x\). A vector supported on the half and perpendicular to constants extends by zero to a centered vector on the full Moore graph. Its adjacency Rayleigh quotient is unchanged. The Moore identity \(A^2+A=56I+J\) gives \(-8\le B\le7\) on this centered subspace. The moments of \(B\) at \(x\) through degree five are
\[
  1,\quad0,\quad32,\quad0,\quad32\cdot63,\quad 2p_x .
\]
For degree four only tree backtracks can close. For degree five any closed walk is a pentagon, since a backtrack would leave a closed odd walk of length three. Subtract \(32^j/1625\) from the moment of degree \(j\) and call the result \(s_j\); this removes the projection onto the constant vector. The matrices \((7s_{i+j}-s_{i+j+1})_{0\le i,j\le2}\) and \((8s_{i+j}+s_{i+j+1})_{0\le i,j\le2}\) are positive semidefinite by applying the spectral bound to polynomials of degree at most two. Their determinants are
\[
 -\tfrac{2560}{13}(7p_x-70432),\qquad
  \tfrac{2048}{13}(11p_x-94736).
\]
Nonnegativity and integrality give the claimed bounds. Sum the vertex counts and divide by five in each half, then use that the two half totals sum to 6,094,400. These necessary conditions constrain a hypothetical split without presuming any vertex transitivity.

## The exact size of the missing five-phase data at 13

**CONDITIONAL representation theorem.** If a free order-13 action is available, choose a basepoint on each orbit and label its phases in \(\mathbb Z/13\). This calculation needs no equitable split. Forget orbit IDs, recording just the five phases along a pentagon. The four-phase margins can be measured from three-edge paths. The dimension of the residual linear space of real, rotation-, reversal- and simultaneous-phase-translation-invariant inventories after these margins are fixed is **1980**. It measures freedom in the linear observations on arbitrary arrays; graph realizability can impose additional relations.

In Fourier coordinates, only tuples of five *nonzero* frequencies can remain after fixing all four-coordinate marginals. Simultaneous translation further requires their sum to vanish modulo 13. There are \((12^5-12)/13=19140\) such ordered tuples by the finite-character filter. A nonidentity rotation of five positions can fix none, since five equal nonzero frequencies do not sum to zero. A reflection fixes 132 tuples: its two transposed pairs have frequencies \(u,v\ne0\); the fifth frequency \(-2u-2v\) is nonzero except at \(v=-u\). Burnside therefore gives \((19140+5\cdot132)/10=1980\) dihedral orbits. Complexification preserves the dimension of this real linear kernel. The basepoint choices are part of the encoding. Orbit IDs, edge incidences, and the five-way correlations themselves cannot be recovered from the marginal dimension.

Run `python3 -B square_cross.py` in this directory for both square-sum windows, the three Moore systems, the 16 rooted formulas, the two determinants, and the finite Fourier count. [Replay](square_cross.py)
