# Extremal six-run cuts reconstruct the tournament and expose Golay

2026-11-27. The [six-run switch profile](../pants_skew_conference_quartic_detection_20261127/README.md)
first detects a skew conference tournament, and a
[rooted slice](../pants_rooted_six_run_reconstruction_20261127/README.md)
recovers its arrows. A still smaller sufficient observation is the
**labels of the \(2n\) minimum cuts**. In order 24 the even part
of this observable has a code interpretation, and for the
bordered Paley tournament it realizes extended binary Golay
on the same 24 switching coordinates.

## 1. A fourth moment locates all \(2n\) minima

Let \(A^2=-(n-1)I\), with \(A\) the signed skew adjacency of a
labelled tournament, \(n\ge8\). A cut is a sign vector
\(\varepsilon\) modulo simultaneous negation, and \(F_6\) counts
ordered covers with one designated six-vertex path and \(n-6\)
singleton runs. Set \(m=n-1\), \(L=n(n-2)\), and
\(T=\sum_i(A\varepsilon)_i^4\). Exact fourth-power expansion gives
\[
 T=n(3m^2-2m)+24\sum_{|S|=4}Q_S\chi_S(\varepsilon),
 \qquad
 F_6=\frac{(n-6)!}{32}\big[(n)_6+n(3m^2-2m)-T\big].
\]
Here \((n)_6=n(n-1)\cdots(n-5)\). To see the first identity,
the constant is the fourth moment of \(m\) independent signs in
each row. Terms with exponents \(3,1\) sum to zero by column
orthogonality. For exponents \(2,1,1\) with odd coordinates
\(k,l\), the squared-entry factor on coordinate \(j\) removes
only row \(j\); summing over \(j\ne k,l\) again gives zero by
column orthogonality. Terms with distinct coordinates are \(24Q_S\).
The second formula follows from the previously established
quartic coefficients and the uniform-group mean \(n!/32\).

Put \(y_i=(A\varepsilon)_i^2-1\). These are nonnegative, since
each row is a sum of the odd number \(m\) of signs, and
\(\sum_i y_i=L\) by skew orthogonality. Then
\(T=n+2L+\sum_i y_i^2\) is maximal exactly when one \(y_i=L\),
as \(\sum_i y_i^2\le L^2\). This is attainable: for row \(r\)
take \(\varepsilon\) equal to the signs of that row off \(r\),
up to global negation. The free bit \(\varepsilon_r\) makes
\(r\) a source or a sink after switching. Write \(E_r\) for
the resulting two-point edge in the quotient Hamming cube.

For distinct rows \(r,s\), their products agree on exactly half
of the \(n-2\) other coordinates, by orthogonality; both of the
omitted coordinates can match using the two free bits. Thus
the minimum distance between \(E_r,E_s\) in the quotient cube
is \((n-2)/2\). In particular the induced graph on the minima
is the matching of the \(n\) edges \(E_r\).
The changing coordinate names the row \(r\), and their common
off-\(r\) signs recover that row up to sign. Fix one row sign;
skew symmetry then fixes every other row's sign. The remaining
two choices are exactly global converse.

This is compression of the response once labelled minima are
available. It is not a claim about the cost of finding the
minima from black-box path-count queries.

## 2. The next response levels and an exact gap

Write \(d\) for quotient Hamming distance to the nearest \(E_r\),
and \(S_2=\sum_i y_i^2\). The largest absolute row sum is
\(M=n-1-2d\). At \(d=0,1,2\) the responses satisfy:

| Distance | \(S_2\) | Number of cuts |
| --- | --- | ---: |
| 0 | \(L^2\) | \(2n\) |
| 1 | \(((n-2)(n-4))^2+32(n-2)\) | \(2n(n-1)\) |
| 2, \(n\ge12\) | \(((n-4)(n-6))^2+144(n-4)+16n\) | \(2n\binom{n-1}{2}\) |

The stated shells around distinct edges do not overlap, by the
edge distance. For \(d=1\) all other row magnitudes are 1 or 3,
and their squared norm forces \((n-2)/2\) threes.
For \(d=2\), flip \(j,k\ne r\) from a nearest edge vector
\(\varepsilon\) which agrees with row \(r\) off \(r\). At rows
\(i\ne r,j,k\) write \(u_j(i)=\varepsilon_j A_{ij}/(\varepsilon_r A_{ir})\)
and similarly \(u_k(i)\). Before flipping, row \(i\) has sum
\(\varepsilon_r A_{ir}\), by row orthogonality. The new magnitude
is five iff \(u_j(i)=u_k(i)=-1\). By column orthogonality the sums of
\(u_j,u_k\) on these rows add to zero: the removed terms at \(k,j\)
have opposite signs. Their product sum is \(-1\) since the
removed term at \(r\) accounts for \(+1\) in the orthogonal
column product. Thus \((n-4)/4\) rows achieve magnitude five.
Rows \(j,k\) can only have magnitude one or three. The fixed
squared norm now forces exactly \(n/4\) threes in total off \(r\).

All higher distances have smaller \(S_2\) than the listed shells.
For \(d\ge2\), \(S_2\le (n-4)(n-6)L\), which is below the \(d=1\)
value by \(8n(n-2)\). For separation from \(d=2\), at \(n=12\)
and \(d\ge3\) use \(M\le5\), giving \(S_2\le24L<3648\).
At \(n\ge16\), states with \(M\le n/2\) have
\(S_2\le(n^2/4-1)L\). For \(n/2\le M\le n-7\), the distance
between distinct edges forces all other row magnitudes at most
\(n-M\). It follows that
\[
 S_2\le g(M)=(M^2-1)^2+((n-M)^2-1)(L-M^2+1).
\]
The cubic \(g(M)\) has second derivative \(2n(6M-2)>0\) here,
so its largest bound on this interval occurs at an endpoint.
At \(x=n-16\ge0\), the \(d=2\) value exceeds \(g(n/2)\) by
\(3x^4/4+57x^3/2+365x^2+1758x+2272\), and exceeds \(g(n-7)\)
by \(8x^3+240x^2+1984x+3072\). Both are positive. This also covers
the bound for \(M\le n/2\), which equals \(g(n/2)\).

The exact gap between the smallest and next smallest \(F_6\) is
\[
 \frac{(n-6)!\,n(n-2)(n-4)}4.
\]
Thus even noisy measurements of the full response still recover
the tournament, provided their uniform absolute error is below
half this gap and one keeps the bottom \(2n\) labels.

At order eight the first two shells exhaust the cube, giving
1176 on 16 cuts and 1272 on 112 cuts. At order twelve the first
three shells give \(T=14652,6972,3900\) on \(24,264,1320\) cuts.
The overall mean and variance of \(T\) are 4092 and 3041280,
from the quartic mass formula. Subtracting the known shells
from the first and second moments forces the remaining 440
cuts to have mean 2364 and zero variance. The four \(F_6\)
values are \(14731200,14904000,14973120,15007680\), respectively.
These histograms are universal among the corresponding skew
conference matrices. The replay counts actual paths independently.

## 3. The even minima give a binary code on the same contacts

For a skew conference matrix of order 24, normalize a source.
Every row of \(I+A\) has even sign parity: the root row is
positive and the other rows have 12 negative entries
(11 in the regular core and one at the root).
The 48 even minimum vectors, retaining both antipodes,
are exactly these rows and their negatives. The other endpoint
of each minimum edge is obtained by toggling its row coordinate.

The binary span \(C\) of the even minima is doubly even
and self-dual. To see this without a size-24 classification,
flip nonroot rows to make the root column positive, which
does not change the span because the all-one word is present.
The normalized Hadamard binary rows have weight 12 except for
the zero root row; any two distinct nonroot rows meet in 6
places. Remove the root row and column to give a square binary
matrix \(B\) of size 23, regarded also as an integer matrix.
Then \(BB^\top=6(I+J)\). The kernel of \(B^\top\bmod 2\) is
totally isotropic for the alternating form \(I+J\bmod2\):
integer lifts \(x,y\) in the kernel have
\(x^\top BB^\top y\equiv0\bmod4\). The form has rank 22 and
a one-dimensional radical, so the kernel has dimension at most
12. Hence \(B\) has binary rank at least 11, and its even
self-orthogonal rows give rank at most 11. Adding the all-one
word with nonzero root coordinate makes dimension 12.
Every generator has weight divisible by four and the generators
are orthogonal, so every codeword has weight divisible by four.

There is an intrinsic weight-four probe in the same quartic
data. A four-set \(S\) is such a codeword iff the product of
its columns is constantly \(+1\) in the normalized \((I+A)\)
*before* normalizing the root column, because \(C=C^\perp\).
Equivalently, with \(p_S=|\mathrm{Pf}(A[S])|\), this is
\[
 p_S=3,\qquad Q_S=20 h_S,
 \quad h_S=\prod_{j\in S-\{i\}}A_{ij}\ \text{for any }i\in S .
\]
For \(p_S=1\) the prior bound \(|Q_S|\le16\) excludes a
constant product outside \(S\). For \(p_S=3\) the local row
products \(h_S\) are all equal: check a source plus triangle
and note that a switch multiplies all four products by the
same factor. Since four signs are multiplied outside \(S\),
their row products agree with the column-product terms in \(Q_S\).
Thus the whole product is constant iff \(Q_S=20h_S\);
the root row makes that constant \(+1\). This identifies
the obstruction to being Golay, not an assertion that every
skew conference matrix of order 24 gives Golay.
There are [other length-24 skew Hadamard codes of minimum weight four](https://link.springer.com/article/10.1007/s10623-008-9173-y).

## 4. Paley order 24: Golay and the projective action

For the bordered Paley matrix on \(\mathbb F_{23}\cup\{\infty\}\),
take \(A_{ij}=\chi(j-i)\) on finite distinct labels, with the
positive root \(\infty\). In bit order \(0,\ldots,22,\infty\)
the finite row at zero is \(0{\tt xfacca0}\). Exact elimination
and enumeration of the even-minimum span give dimension 12,
minimum weight 8, and the weight enumerator
\[
 1+759z^8+2576z^{12}+759z^{16}+z^{24}.
\]
Thus this is extended binary Golay, the classical
[Hadamard-to-Golay construction](https://arxiv.org/abs/1006.2414)
identified here directly from minimum tournament path counts.

Use \(v_x=(1,x)\) and \(v_\infty=(0,-1)\), so
\(A_{xy}=\chi(\det(v_x,v_y))\) off the diagonal.
For \(g\in\mathrm{SL}_2(23)\), let \(g v_x=\lambda_x v_{gx}\);
then \(A_{gx,gy}=\chi(\lambda_x)\chi(\lambda_y)A_{xy}\).
The signed permutation
\(\varepsilon'_{gx}=\chi(\lambda_x)\varepsilon_x\)
preserves the fourth moment and the minimum matching.
Multiplying the arrow formula along a row gives
\(\prod_x\chi(\lambda_x)=+1\), since all off-diagonal row
products of \(A\) are positive. Thus the signed permutation
preserves even minima and acts on their binary lifts as
\(w\mapsto Pw+t\). Since the zero vector is a minimum,
\(t\) belongs to \(C\). Taking spans then shows that the
coordinate permutation \(P\) alone preserves \(C\).
The projective \(\mathrm{PSL}_2(23)\) action on the same
24 contacts is therefore directly visible in the minimum-cut
experiment. The two generators \(x\mapsto x+1\) and
\(x\mapsto-1/x\) have affine shifts 0 and \(0{\tt xfacca0}\)
in the stated convention.

Run the companion Node.js code for direct counts on every
order-eight and order-twelve switch, the Golay weight
enumerator, weight-four obstruction check, and the generator
controls. The code does not enumerate the order-24 switch group.
