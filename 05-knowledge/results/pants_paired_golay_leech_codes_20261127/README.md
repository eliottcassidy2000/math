# Two Golay codes and two Leech lattices from tournament minima

2026-11-27. The [extremal six-run result](../pants_extremal_cuts_golay_20261127/README.md)
produces one minimum edge per vertex. Its source and sink ends
give a naturally **ordered pair** of binary codes, not just the
single code extracted from even Paley minima before. The pair
recovers the switching class, and a simple weight test describes
its image. At order 24 this gives a pair of Golay codes, with an
exact common symmetry and a common Leech-lattice sublattice.

## 1. Two complementary Type II codes

Let \(A\) be a labelled skew conference matrix of order \(n\)
divisible by eight. Retain both antipodal lifts of each minimum
cut. In each of the source and sink halves, span the binary
differences between its sign masks, obtaining \(C_+,C_-\).
These are defined on the same \(n\) coordinates. The source
endpoints are the signed rows of \(I+A\) and the sink endpoints
are the signed rows of \(I-A\). In each half every pair of
distinct nonopposite rows has distance \(n/2\) by orthogonality,
so their masks have fixed parity and the two halves have
opposite parity (their paired endpoints differ in one bit).
Switching translates all cuts by the same mask, so the
difference spans are switch invariant. Converse swaps them.

Normalize a root \(r\) as a source. Let \(h_i\) be the negative
mask of row \(i\) of \(I+A\), and \(\mathbf1\) the all-one word.
Then \(h_r=0\), and taking differences from each half's root
row gives
\[
 C_+=\langle\mathbf1,h_i\rangle,\qquad
 C_-=\langle\mathbf1,h_i+e_i+e_r\rangle .
\]
The relevant Hadamard row distances \(n/2\) give nonroot
generator weight \(n/2\) and mutual intersection \(n/4\).
As \(8\mid n\), both codes are doubly even and self-orthogonal.
They contain the repetition code \(R=\langle\mathbf1\rangle\),
and their sum contains all \(e_i+e_r\), hence equals the even
space \(E\). The dimension formula forces both dimensions
to be \(n/2\) and the intersection to be exactly \(R\).
Thus they are binary doubly even self-dual codes (Type II);
they split the quadratic space \(E/R\) into two complementary
totally singular halves. This derives the code pair from
actual minimum labels, in the setting of the
[classical skew Hadamard codes](https://link.springer.com/article/10.1007/s10623-008-9173-y).

## 2. Inverse classification by weight-two words

Conversely, take two labelled Type II codes \(C,D\) with
intersection \(R\). Let \(\pi\) project \(E/R\) onto \(C/R\)
along \(D/R\). They arise as the ordered pair above **if and
only if** every projection \(\pi(e_i+e_j+R)\) has weight
\(n/2\). This weight condition is independent of the choice
between the two complementary representatives.

To reconstruct, fix \(r\) and lift the first projection of
\(y_i=e_i+e_r\) to the unique \(c_i\in C\) with \(c_i(i)=0\);
put \(c_r=0\), and \(d_i=y_i+c_i\in D\).
The quadratic function \(q(x)=\mathrm{wt}(x)/2\bmod2\) on \(E\)
has binary dot product as polarization. It is zero on \(C,D\),
so \(q(y_i)=1\) gives \(c_i\cdot y_i=1\), and hence \(c_i(r)=1\)
for \(i\ne r\). For \(i,j\ne r\) distinct, the identity
\(y_i\cdot y_j=1\) gives \(c_i\cdot y_j+c_j\cdot y_i=1\),
since the same-code products are zero. The root bits cancel,
leaving \(c_i(j)+c_j(i)=1\). Thus the \(c_i\) are the negative
row masks of a skew sign matrix \(I+A\) with root \(r\)
a source. The distance between rows \(i,j\) is
\(\mathrm{wt}(c_i+c_j)\), the projected weight for \(e_i+e_j\).
Every such distance equals \(n/2\) precisely when
\((I+A)(I+A)^\top=nI\). The previous construction then gives
exactly the two input codes. This reconstructs the labelled
switching class; swapping the two codes gives its converse.

For **two Golay codes** with intersection \(R\), the extra
test is equivalent to saying no octad of one code meets an
octad of the other in seven coordinates. Indeed, a projection
of \(y=e_i+e_j\) is nonzero modulo \(R\), since neither code
contains a weight-two word. Choose for it a representative of
weight at most 12. If that weight is not 12, it is 8.
Because \(c\cdot y=1\), toggling the two coordinates to get
the other component preserves that weight. The two octads
then differ exactly in the pair \(\{i,j\}\). Conversely, any
two opposite octads with seven common coordinates decompose
their weight-two difference and violate the test.

## 3. The Paley pair has exact projective symmetry

For the order-24 bordered Paley matrix on
\(\mathbb F_{23}\cup\{\infty\}\), the source code is Golay by the
[previous span calculation](../pants_extremal_cuts_golay_20261127/README.md).
The map \(x\mapsto-x\) fixing \(\infty\) sends the normalized
core to its converse and exchanges the codes, so both are
Golay. The determinant representation from the same earlier
note shows that projective \(\mathrm{PSL}_2(23)\) preserves
both difference spans and is transitive on the 24 coordinates.

To bound the common coordinate stabilizer, fix \(\infty\).
Reconstruction implies that any such permutation acts on the
unique source-normalized matrix in the switching class, so
its restriction is an automorphism of the Paley core. That
core group is transitive by translation and not doubly
transitive, since it preserves arrows. By
[Burnside's prime-degree theorem](https://arxiv.org/abs/math/0310200)
it is solvable and hence affine: a minimal nontrivial normal
subgroup is abelian, transitive by primitivity, and regular by
faithful abelian transitivity, so is the translation group
of order 23. Exactly the square multipliers are allowed,
giving point stabilizer order \(23\cdot11=253\).
It follows that the ordered common stabilizer has order
\(24\cdot253=6072\) and is exactly \(\mathrm{PSL}_2(23)\).
The nonsquare map already swaps the pair, so the full unordered
stabilizer is \(\mathrm{PGL}_2(23)\), of order 12144.

Fixing the first code, its coordinate automorphism group
[\(M_{24}\)](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M24/)
therefore gives exactly \(244823040/6072=40320\) partners **in the
orbit of this particular second code**. This does not count
all partners satisfying the weight test. The replay checks
that every first-code octad meets \(13,96,228,256,130,32,4\)
opposite octads in \(0,\ldots,6\) points respectively and none
in 7 or 8, as well as both projective generators and the swap.

## 4. Two Leech lattices, one larger-symmetry intersection

Given each Golay code, put its Leech lattice on these same
coordinates. In the standard integral-coordinate convention
([Allcock](https://doi.org/10.1080/10586458.2005.10128938)),
a vector is \(x/\sqrt8\) with all entries of the same parity
\(m\), with \(x_i\equiv m+2c_i\pmod4\) for a codeword \(c\),
and \(\sum_i x_i\equiv4m\pmod8\).
Write the two lattices \(\Lambda_+,\Lambda_-\). Their
intersection \(L\) has only repetition residues:
\[
 x=4k+j\mathbf1,\quad j\in\{0,1,2,3\},\quad
 \sum_i k_i\equiv j\pmod2 .
\]
These are four cosets of \(\sqrt2 D_{24}\), with covolume
\(2^{13}/4=2^{11}\), so the intersection has index 2048
in each Leech lattice. Its 1152 shortest vectors have
squared norm four: in integer coordinates there are 1104
vectors \((\pm4,\pm4,0^{22})\) and 48 vectors
\(\pm(-3,1^{23})\). For \(j=0\) the even-sum condition
excludes a single nonzero coordinate, for \(j=2\) all
coordinates have absolute value at least 2, and for odd
\(j\) the all-unit pattern is excluded by the sum test.
A single entry of magnitude 3 is necessary and sufficient
to reach norm four in the odd cosets.

The lattice sum is \(L^*=\Lambda_++\Lambda_-\), by dualizing
the intersection of the two unimodular lattices. The residue
map modulo repetition has kernel \(L\) in either lattice:
upon adding two residue vectors, carrying the parity adds
at most the all-one word. Thus
\[
 L^*/L\cong(C_+/R)\oplus(C_-/R)\cong(\mathbb Z/2)^{22}.
\]
Its two halves are the two maximal isotropic subgroups of
this decomposition for the lattice discriminant quadratic
form: vectors in each Leech lattice have even norm.
Each lattice recovers its residue code, so the ordered
lattice pair recovers the same labelled switching class.

**The intersection forgets the projective structure.** In
fact \(\mathrm{Aut}(L)=\{\pm1\}\times S_{24}\). For a shortest
vector, count other shortest vectors at absolute inner
product 2. The 48 special vectors have degree 92, the
552 two-coordinate vectors with equal signs degree 176,
and the 552 two-coordinate vectors with opposite signs
degree 180. For the latter two, the 176 come from sharing
exactly one coordinate with another two-coordinate vector;
opposite-sign vectors also meet four special vectors
at absolute inner product 2. A special vector meets 46 other
special vectors and 46 opposite-sign two-coordinate vectors.
Therefore an isometry must preserve the 24 antipodal
special pairs. Choose their representatives
\((\mathbf1-4e_i)/\sqrt8\); their off-diagonal inner products
are all \(+2\), forcing all signs in a signed permutation
of these representatives to agree. The representatives
span the space, so the isometry is a coordinate permutation
times a global sign, and all such maps preserve \(L\).
In contrast the ordered **pair** of Leech lattices retains
exactly \(\{\pm1\}\times\mathrm{PSL}_2(23)\), and the
unordered pair \(\{\pm1\}\times\mathrm{PGL}_2(23)\).

Run the Node.js replay for two-root reconstruction and
projected weights at orders 8, 16, 24, 32, and 48, the
Golay and projective checks, and the common-short-vector graph.
