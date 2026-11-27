# Six-vertex switch detection in skew conference tournaments

2026-11-27. In the [earlier switch analysis](../balanced_pants_switch_fourier_20261127/README.md),
one three-vertex run measured the off-diagonal skew Gram matrix.
When \(A^2=-(n-1)I\), that measurement is identically zero. This note
gives the first nonzero run measurement, its exact variance, and
an arithmetic constraint from the four-vertex pants Pfaffian that
classifies the quartic magnitudes in orders eight and twelve.

Let \(A_{ij}=1\) for \(i\to j\), \(A_{ji}=-A_{ij}\), and \(A_{ii}=0\).
Assume \(A\) is skew conference, meaning \(A^2=-(n-1)I\).
Then \(I+A\) is a skew Hadamard matrix; see
[Nozaki–Suda](https://arxiv.org/abs/1202.5374).
Let \(F_l\) count ordered path covers with a first run of size \(l\)
and all other runs singleton. For an even set \(S\), write
\(\chi_S(\varepsilon)=\prod_{i\in S}\varepsilon_i\) for its switch
character. All Fourier coefficients use uniform expectation over
the switch group.

## 1. Exact first detection

The profile \(F_l\) is constant under switching for every
\(1\le l\le \min(n,5)\). For \(l=6\) and \(n\ge6\), the only
nonconstant Fourier coefficients are
\[
\widehat F_6(S)=-\frac{3(n-6)!}{4}Q_S \quad (|S|=4),
\qquad Q_S=\sum_x\prod_{i\in S}A_{ix}.
\]
The variance is
\[
\mathrm{Var}(F_6)=
\frac{3}{128}((n-6)!)^2n^2(n-1)(n-2)(n-4)>0.
\]
Thus no skew conference tournament of eligible order hides from
a six-run measurement. The two-run profile \((3,3,1,\ldots,1)\)
has the same quartic coefficients multiplied by two, so its
switch variance is four times as large.

**Vanishing at weights two.** Expand the run-edge indicators
\((1+A_{uv}\varepsilon_u\varepsilon_v)/2\). Under summation over
permutations, selected subpaths of even vertex size cancel by
reversal. The only weight-two contributions for \(l\le6\) come
from a two-edge or four-edge simple selected path. The signed
two-edge sum from \(a\) to \(b\) is \((A^2)_{ab}=0\) for \(a\ne b\).

For the four-edge sum start from \((A^4)_{ab}=0\).
A walk \(a,x,y,z,b\) with a repeated vertex either has a zero
factor immediately or lies in one of the equality events
\(y=a,y=b,x=z,z=a,x=b\). The first three event sums are
\(-(n-1)(A^2)_{ab}\); the last two contain
\((A^3)_{aa}\) or \((A^3)_{bb}\), hence are zero.
For clarity, the only pairwise intersections that can give
nonzero terms before summation are
\[
(y=a,x=z),\ (y=b,x=z),\ (y=a,x=b),\ (y=b,z=a),\ (z=a,x=b).
\]
Their respective sums are \(-(A^2)_{ab}\) for the first four
and \((A^2)_{ba}\) for the last. Every triple gives an adjacent
repeat or incompatible endpoints. Inclusion–exclusion therefore
eliminates all the repeated walks with zero total signed weight.
Free placements of any other run vertices contribute the same
multiplicity to each signed simple path, so the cancellations
apply to the Fourier coefficients.

**Weight four.** A six-run weight-four contribution consists
of two disjoint selected two-edge paths. For an endpoint pairing
\(ab|cd\) of \(S\), set \(f_{ab}(x)=A_{ax}A_{xb}\). Sum
\(f_{ab}(x)f_{cd}(y)\) over distinct centers \(x,y\notin S\).
It is the product of separate sums minus \(Q_S\).
By Gram orthogonality the separate sums are
\(-f_{ab}(c)-f_{ab}(d)\) and \(-f_{cd}(a)-f_{cd}(b)\).
Their products cancel over all three endpoint pairings:
the sum is a relabeling-invariant polynomial in the six local
skew signs with switch boundary \(S\). After using \(A_{ij}^2=1\)
only even-degree simple edge monomials of degree at most four
remain. Four edges cannot make every vertex odd, since the two
degree-three vertices would force the other two degrees to be at
least two. Degree two gives a perfect matching, whose coefficient
is zero by swapping the endpoints of one matched edge, reversing
the sign of that monomial while leaving its edge set fixed.
This leaves \(-3Q_S\). The two orders of the components and
two orientations each give 8, multiplied by
\((n-6)!/2^5\) from the other vertices and run indicators.
The \((3,3,1,\ldots,1)\) profile has the same order and orientation
factor but only four run indicators, explaining its factor of two.
For this profile there are no selected four-edge paths.

**Total quartic mass.** Let \(W\) have row indexed by
\(\{a,b\}\) and entries \(W_{ab,x}=A_{ax}A_{bx}\).
Its column Gram is
\[
W^\top W=\frac{n(n-2)}2 I-\frac{n-2}2J .
\]
The diagonal follows from counting row pairs outside one
column index. The off-diagonal follows by squaring
\(\sum_a A_{ax}A_{ay}=0\) and subtracting the \(n-2\) nonzero
squares, then dividing by two. The row Gram has diagonal
\(n-2\); off-diagonal entries for overlapping pairs have absolute
value 1, by subtracting one term from orthogonal rows. For
disjoint pairs the entry is \(Q_S\) for their union. Equate the
squared Frobenius norms, subtracting
\(\binom n2(n-2)^2+\binom n2\,2(n-2)\) for the diagonal and
overlap, respectively. Each \(S\) contributes \(6Q_S^2\), so
\[
 \sum_{|S|=4}Q_S^2=\frac{n^2(n-1)(n-2)(n-4)}{24}.
\]
Parseval gives the claimed variance.

## 2. The four-vertex Pfaffian constrains the six-run readout

Let \(p_S=|\mathrm{Pf}(A[S])|\), either 1 or 3 as in the
[four-cuff partition model](../balanced_pants_switch_fourier_20261127/README.md).
For every \(S\) of size four,
\[
\begin{array}{c|cc}
 & \text{residue of }Q_S\pmod8 & |Q_S|\text{ at most}\\ \hline
p_S=3 & n-4 & n-4\\
p_S=1 & n & n-8 .
\end{array}
\]
For the residue in an arbitrary labelling and switch state,
note that \(Q_S\) just picks up the factor \(\chi_S\) under
switching; \(n\) is a multiple of four for skew conference orders
above two, so the displayed residues are sign invariant.
To see the divisibility directly, normalize a source by switching.
The other rows have equal positive counts on the core from
orthogonality to the source row. Two core rows disagree on
\((n-2)/2\) entries outside their two labels and the source,
by their zero full inner product. The two positive counts
outside their two labels differ by one, so the disagreement
count is odd; this forces \(n\equiv0\bmod4\).

**Exterior column proof.** On each column \(x\notin S\) put
\(y_i=A_{ix}A_{0x}\) with \(0\) denoting the first vertex of \(S\).
The first coordinate becomes 1 and the other three signs
describe the relative column. Write \(m=n-4\) and \(Q=Q_S\).
The count for any of the eight prescribed \(y=(1,y_1,y_2,y_3)\)
is
\[
 \frac{m+\sum_{i<j} c_{ij}y_i y_j+Qy_1y_2y_3}{8},
 \qquad c_{ij}=-\sum_{x\in S}A_{ix}A_{jx}.
\]
This is the three-sign Fourier inversion formula, with the
six lower moments set by full row orthogonality. In the
Pfaffian-three class the restricted rows are orthogonal, so
the middle term is zero. In the Pfaffian-one class switch the
four-set to a transitive representative (or switch its converse
to a transitive one and reverse the linear order). In transitive
order the middle term is
\(-2[(1+y_2)y_1+(y_2-1)y_3]\). It takes values \(4,-4\), each
with both signs of \(y_1y_2y_3\). Nonnegative integer counts
give the two bounds and residues.

**A universal order-eight/order-twelve law.** The number of
Pfaffian-three four-sets is
\[
N_3=\frac{n^2(n-1)(n-2)}{96}.
\]
Indeed the sum of four-by-four principal determinants is the
coefficient of \(t^{n-4}\) in
\((t^2+n-1)^{n/2}\); each determinant is 1 or 9.
At order eight the residues and bounds force
\(|Q_S|=4\) exactly on the 28 Pfaffian-three four-sets.
At order twelve they force \(|Q_S|=4\) on the 330
Pfaffian-one four-sets; this exhausts the universal squared
mass \(5280\), so all 165 Pfaffian-three coefficients vanish.
These are statements about **every** skew conference matrix
of the respective order, not just a Paley example. The same
local Pfaffian that separates the balanced pants orbits thus
selects the support of the first nonconstant run observable,
with the selection reversed between orders eight and twelve.

The companion replay builds the bordered prime Paley examples
([Pantangi](https://arxiv.org/abs/1905.08568)), then independently
counts directed paths by subset DP and checks every Walsh
coefficient through run size six. The proofs above do not
require a classification of skew conference matrices.
