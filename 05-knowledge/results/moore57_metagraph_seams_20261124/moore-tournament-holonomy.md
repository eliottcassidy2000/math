# A Grassmannian and tournament-switching quotient of Moore transport

Fix a hypothetical degree-57 Moore graph and its rooted matching description: 57 branches, 56 sheets, reciprocal orthogonal permutation transports \(P_{ij}\) between the leaf fibers. On the direct sum of centered leaf spaces the block adjacency \(B\) has zero diagonal blocks and obeys \(B^2+B=56I\). This follows by restriction from the Moore adjacency identity; branch-to-leaf terms vanish on zero-sum label vectors. All conclusions below are conditional on a compatible transport array, not an existence result.

## Fusion geometry with an exact transport meaning

Let \(U\) be an invariant real subspace of the zero-sum label space, of dimension \(d>0\), with an orthonormal frame. On \(U^{57}\) form the orthogonal projection
\[
 \Pi=\frac{B+8I}{15}
\]
onto the eigenvalue-7 space (the eigenvalues of \(B\) are 7 and \(-8\)). Since the diagonal blocks of \(B\) vanish, \(\Pi_{ii}=(8/15)I_d\), so its rank is \(57\cdot8d/15=(152/5)d\). In particular \(5\mid d\). For \(I_i:U\to U^{57}\) the coordinate embedding, use the frames \(E_i=\sqrt{15/8}\,\Pi I_i\) in \(\operatorname{im}\Pi\). They satisfy
\[
 E_i^\mathsf{T}E_i=I_d,\qquad E_i^\mathsf{T}E_j=P_{ij}/8\ (i\ne j),\qquad
 \sum_i E_iE_i^\mathsf{T}=(15/8)I_{\operatorname{im}\Pi}.
\]
Thus the branches give 57 equi-isoclinic \(d\)-planes, all principal cosines \(1/8\), forming a tight fusion frame. The fiber permutation is precisely the orthogonal part of a cross-Gram block, not merely a visual analogy. A change of basis in each branch is exactly the transport gauge. The path ordered product of polar parts changes by its two endpoint gauges, while a closed product changes by conjugation. One can also read the matrix-level cancellation:
\[
 \sum_{k\ne i,j} P_{ij}^{-1}P_{ik}P_{kj}=-I_U
\]
because the off-diagonal \(ij\) block of \(B^2=-B+56I\) is \(-P_{ij}\). The average triangle transport at that edge is therefore \(-I_U/55\); taking determinants term by term does not commute with this averaging.

## Oriented Plücker lines recover a tournament switching class

Take the wedge of the \(d\) frame vectors for each plane. These unit decomposable vectors are an equiangular system of oriented lines, with off-diagonal Gram entries \(8^{-d}\det(P_{ij})\). Write \(\epsilon_{ij}=\det(P_{ij})\in\{\pm1\}\). It is symmetric in \(i,j\) because inverse has the same determinant. Changing frame orientations by \(s_i\in\{\pm1\}\) gives \(\epsilon_{ij}\mapsto s_i s_j\epsilon_{ij}\). Its triangular sign
\[
 H_{ijk}=\epsilon_{ij}\epsilon_{jk}\epsilon_{ki}
\]
is exactly the determinant of triangle holonomy. These triangular signs form a complete invariant of the sign field up to switches: choose one anchor and gauge all its edges positive; every remaining edge is then its triangle sign through the anchor. On each four vertices the product of the four triangle signs is 1. Conversely any assignment with that tetrahedron identity arises from such a sign field.

Choose an order of the 57 branches and multiply the sign \(\epsilon_{ij}\) by the antisymmetric transitive-order sign of the pair. This yields an actual tournament on branch labels, canonically a *switching class relative to that ordered reference*, not an unmarked Moore invariant independently of the reference order. Along a fixed oriented word, products of its edge signs transform only at the endpoints, exactly as in projective determinant switching. Products along closed words are invariant. The same formulas apply to multiple determinant characters simultaneously. As in other switching constructions, the unsigned count of directed Hamilton words is not claimed to be preserved. Moreover a nontrivial pure sign character cannot be a one-dimensional invariant subspace of the centered label permutation module, since \(5\mid d\); the quotient lives in exterior power, not a scalar leaf mode.

## Gewirtz has two independent determinant layers

For the [Gewirtz label graph](https://aeb.win.tue.nl/graphs/Sims-Gewirtz.html), with adjacency \(A\), the centered eigenspaces have dimensions \(35,20\), eigenvalues \(2,-4\) respectively. Each is invariant under all Gewirtz automorphisms. They give separate orthogonal transport geometries and separate sign fields. The ordinary sign of the 56-label permutation is the product of their determinants, because the constant line contributes +1.

The [standalone exact enumeration](moore-gewirtz-sign.py) constructs a 56-element hyperoval orbit in \(\mathrm{PG}(2,4)\), with adjacency by disjointness, and generates the full group with 80,640 permutations. Determinants of the six elementary projective generators are \((+,+)\) on the \((35,20)\) spaces. Frobenius and polarity give \((+,-)\) and \((-,+)\), respectively. Thus these characters recover the entire four-element outer quotient by \(\mathrm{PSL}_3(4)\), not just ordinary parity. The exact derangement counts are
\[
\begin{array}{c|rrrr}
(\det_{35},\det_{20})&(+,+)&(+,-)&(-,+)&(-,-)\\ \hline
\#&8280&10800&9240&12600.
\end{array}
\]
The arithmetic in the replay uses characters, not floating point eigenvectors. The trace of \(P\) on the two eigenspaces is respectively \((\operatorname{tr}(PA)+4\operatorname{tr}(P)-14)/6\) and \((2\operatorname{tr}(P)-\operatorname{tr}(PA)+8)/6\), by the spectral projectors. For even order \(m\), the multiplicity of the eigenvalue \(-1\) is \(m^{-1}\sum_{k=0}^{m-1}(-1)^k\operatorname{tr}(P^k|U)\); its parity gives the determinant. At odd order the determinant is +1.

These data provide a four-valued switching quotient for any putative Gewirtz-based Moore array, even if abstract labels have no rigid preference. All four derangement classes have members, so assigning a determinant pair alone cannot certify a Moore triangle or rule it out.

## Sharp transport forces the ordinary sign field to be flat

The finite replay supplies a particularly short certificate for ordinary permutation parity. With the sorted hyperovals indexed from zero as in the code, take the following rectangle of two label sets:
\[
 B=\{25,27,33,35,37,39,45,47,49,50,53,54\},
\]
\[
 C=\{9,11,14,15,25,26,29,31,32,34,37,39,41,43,44,47,48,51,53,54\}.
\]
For every Gewirtz automorphism \(g\), exact enumeration verifies \(\operatorname{sgn}(g)=(-1)^{|g(B)\cap C|}\). This can be checked directly for all group elements, without trusting how the certificate was found. If \(S\) is sharply transitive, double count incidences from \(B\) to \(C\):
\[
 \sum_{g\in S}|g(B)\cap C|=|B||C|=240 .
\]
Thus \(S\) always contains an even number of odd permutations. Apply this to **every pair of branches** in the Moore array. Write \(s_{ij}=\operatorname{sgn}(P_{ij})\), \(\delta_i=\prod_{k\ne i}s_{ik}\). The product of the signs of its 56 direct and two-step maps is
\[
 s_{ij}\prod_{k\ne i,j}s_{ik}s_{kj}
 =s_{ij}\delta_i\delta_j=1 .
\]
Consequently \(s_{ij}=\delta_i\delta_j\) and the ordinary sign of every closed holonomy is positive. In particular **the entire array, in the usual identity-anchor gauge, is confined to the even subgroup of order 40,320**. This is stronger than certifying some isolated sharp set; it removes one color coordinate globally.

The remaining nontrivial even-subgroup character \(\chi\) (detecting the equal minus/minus determinants versus plus/plus) admits another exact test. Modulo three, \(\chi\), viewed as a map with values \(\pm1\), is a linear functional of the permutation matrix. The [mod-three elimination replay](moore-gw-mod3-target.py) checks this on all 40,320 even group elements and verifies that \(\chi(J_{56})=0\). Equivalently the binary tag \((1-\chi)/2\) predicts residue 1 on \(J_{56}\). This means any sharply transitive set in the even subgroup uses \(1\bmod3\) members of the outer minus/minus coset. For the global sign matrix \(S_{ij}=\chi(P_{ij})\) off diagonal and zero on diagonal, the same pairwise equation gives
\[
 S^2+S=56I\pmod3,\qquad (S-I)^2=0\pmod3.
\]
In other words the rows of \(S-I\) generate a ternary self-orthogonal code (and its rank is at most 28). These congruences do not rule out all candidates. In particular, the remaining character is *not* a linear function of the matrices modulo two; the finite span test also finds obstructions modulo five and seven.

The rectangle source set B even induces six disjoint edges in the label graph, as the replay's internal degrees show. This is a literal placement-sensitive parity count: one sees which image labels land in C, not simply the number of moved labels. In the anchored even subgroup, trying to peel the remaining sign by lifting the character bit to twice that bit modulo four also fails. [The exact replay](mod4-outer.py) uses unit pivots over Z/4 and reduces the remaining even rows over F2; it produces a zero matrix coordinate with nonzero tag 2, disproving this proposed linear guard. It does not preclude a nonlinear test or a different modular lift.
