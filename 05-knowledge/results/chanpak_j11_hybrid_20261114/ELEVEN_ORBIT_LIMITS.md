# Constraints on hybrid exact decoders for the eleven O/Z/J blocks

The following apply to the eleven-block language of wave 20, not to the
different 13-generator port-chain language. They rule out proposed types of
proof, not cofinite first-coordinate support.

## Positive coordinates and the eleven matrices

Use the known O/Z/J state \((H,E,L,h)\):
\[
O=(H,0,2L,h),\quad
Z=(H+L-E,h,2L-h,h),\quad
J=(H+L-E,H-E,H+L,H).
\]
At block boundaries use \((A,E,L)\) with \(A=H-E\). A block \(Jw\) resets
\(h=H\), so it acts on these three coordinates alone. Use the order
\[
w\in(\epsilon,O,OO,OZ,OOO,OOZ,OZO,OZZ,OOOZ,OZOZ,OZZZ).
\]
Expanding the three moves gives the following rows; all entries are
nonnegative:

| \(w\) | First row | Second row | Third row | det |
|---|---|---|---|---:|
| empty | 0,0,1 | 1,0,0 | 1,1,1 | 1 |
| O | 1,0,1 | 0,0,0 | 2,2,2 | 0 |
| OO | 1,0,1 | 0,0,0 | 4,4,4 | 0 |
| OZ | 2,1,3 | 1,1,0 | 3,3,4 | 4 |
| OOO | 1,0,1 | 0,0,0 | 8,8,8 | 0 |
| OOZ | 4,3,5 | 1,1,0 | 7,7,8 | 8 |
| OZO | 3,2,3 | 0,0,0 | 6,6,8 | 0 |
| OZZ | 4,3,7 | 1,1,0 | 5,5,8 | 8 |
| OOOZ | 8,7,9 | 1,1,0 | 15,15,16 | 16 |
| OZOZ | 8,7,11 | 1,1,0 | 11,11,16 | 16 |
| OZZZ | 8,7,15 | 1,1,0 | 9,9,16 | 16 |

Each preserves \(L\ge H=A+E\) on the nonnegative cone (indeed the row for
\(L'-H'\) is nonnegative). Each has \(A'\ge L\), \(L'\ge H+L\) and
\(H'\ge A+L\). In particular *every two-block word* gives \(H''\ge3H\) in
the physical cone. None of the single letters has eigenvalue 2. There is
therefore no positive physical direction with multiplier exactly 2 along
any return word in the eleven-block alphabet.

## Determinant obstruction for integer radix

Use the common cone test vector \(v=(1,0,7/4)\). The respective
\(c_i=\min((M_i v)_A,(M_i v)_L/(7/4))\) are
\[
(11/7,11/4,11/4,40/7,11/4,12,33/4,76/7,95/4,156/7,148/7).
\]
Thus \(M_i v\ge c_i v\) coordinatewise, and for every full-rank letter
\(c_i>|\det M_i|\). For any word \(M\), positivity and iteration
give \(\rho(M)\ge\prod_i c_i\), where repetition of letters is counted.
If \(M\) is full rank this is strictly larger than \(|\det M|\).
An integer Perron root would divide the nonzero determinant by the
rational-root theorem. Consequently **no full-rank word has an integer
Perron root**. Any exact integer-radix ray self-return of this type has
to use one of the four O-ending rank-deficient letters.

## Quantitative two-dimensional sparsity, even with arbitrary prefix

At \(s=37/20\), upper-bound \(1000c_i^{-s}\) respectively by
\[
(434,154,154,40,154,11,21,13,3,4,4).
\]
These are strict exact inequalities by taking 20th powers, and their sum
is 992. Weight any word by \((\prod c_i)^{-s}\); the sum over all lengths
is at most \(1/(1-.992)\). If a starting state dominates \(\gamma v\),
then \(H\le X\) implies \(\prod c_i\le X/\gamma\), so the number of words,
and hence the number of states, is \(O((X/\gamma)^s)\).

This holds even when the first run is arbitrary, not bounded by four.
A one-apex source indexed by nonnegative integer \(x\) has \(L=x+1\);
after the first block both \(A'\) and \(L'\) are at least \(x+1\).
It therefore dominates \(\gamma v\) with \(\gamma=(4/7)(x+1)\).
Sum the preceding estimate over \(x\ge0\), using \(s>1\) and the
eleven choices for the first block. Where no block follows the source,
there are at most \(X+1\) states in \(H,L\le X\). Hence the common natural
hybrid feeds only \(O(X^{37/20})\) projected pairs in that square.
In particular it cannot fill any two-dimensional lattice cone of
positive area. The complementary method has to select a lower-dimensional
section or produce genuinely new source coordinates; it cannot require
blanket integer-pair state completeness. This is compatible with scalar
cofiniteness.

## An arithmetic restriction on prime words

For a proof by induction using ordinal sums, after the already known seeds
\(63,343\), a smallest additional hole must be a prime or \(7p\) for a
prime \(p\). In fact a composite with at least three prime factors has a
proper prime-factor split avoiding \(7,21\), except \(63,343\). To see this,
if there are four or more factors and one prime other than 7, choose it:
the other factor has at least three prime factors and thus is not 7 or 21.
If all factors are 7, split as \(49\) times a product of at least two 7's.
For exactly three factors, if any prime other than 3 and 7 appears, either
it or another non-7 prime gives a split (and \(7\cdot 7\cdot p\) can split
as \(49\cdot p\)). For just 3 and 7, the only unsplittable three-factor
values are \(63,343\); \(147=3\cdot49\). With two prime factors, the
unsplittable possibilities are \(7p\).

There is an exact subshift restriction for these two remaining target
types. If \(E=0\), then \(JO^j\) with \(j=1,2,3\) maps the state to
\((H+L,0,2^j(H+L))\). The normalized state \((1,0,2^j)\) is itself a
physically realized transitive source (append \(j\) O's to the singleton).
Every ensuing suffix realizes a physical multiplier, not just an algebraic
integer multiplier; a nonempty suffix has multiplier at least 3. Each of
\(JO^j\) and \(JOZO\) leaves \(E=0\), with \(L>0\). Thus the two-block
marker
\[
\{JO,JOO,JOOO,JOZO\}\{JO,JOO,JOOO\}
\]
is excluded in every minimal-length physically valid representative of
a prime or \(7p\), where length counts ordinary letters. The reset
height \(H+L>1\) and the suffix factor would both be proper admissible
counts; for a prime this is impossible, and for \(7p\) one factor
would have to be 7, also impossible. If the reset block is terminal,
delete its terminal O's to shorten the word without changing \(H\).
This is not a projection-level residue obstruction and does not imply
the remaining subshift is cofinite.

## Finite spectral control, not an all-word theorem

For every nonempty word of length at most eight, an exact characteristic
polynomial check gives \(p(1),p(-1)<0\). For full rank, the determinant
obstruction proves that the product of the two remaining roots has modulus
\(|\det M|/\rho(M)<1\); the rank-deficient case has zero determinant.
The two sign tests then imply that both non-Perron roots have modulus
less than one (divide the polynomial by its positive Perron factor,
and use its Monic quadratic sign at both endpoints). The all-word
statement is open; it would strengthen the integer-radix barrier beyond
full rank. Run the control with "g++ -O3 eleven_spectral.cpp -o /tmp/spectrum11;
 /tmp/spectrum11 8"; the integer count totals 235,794,768 nonempty words.

## Why a unit-gap repair cannot be appended across a single seam

Distinguish the port-chain even/odd quantities \(e_{\rm p},o_{\rm p}\) from
the coordinate \(E=o_{\rm p}-e_{\rm p}\) above. In the usual single
reversed-seam gluing with a transitive block \(B\), the increment is
\[
H'-H=2(q_o e_{\rm p}+q_e o_{\rm p}),
\]
where \(q_o\ge1\), \(o_{\rm p}\ge1\), and \(q_e,e_{\rm p}\ge0\).
This is the inherited odd-cycle packing join law.

If \(e_{\rm p}=0\), the old port is a sink: otherwise some two-vertex
path from it contributes at least one to \(e_{\rm p}\). Then
\(o_{\rm p}=H\), because deletion of a sink does not change the count,
and an increment of two forces \(H=1\). Otherwise that increment forces
\(e_{\rm p}=1\). Each two-vertex path \(u,w\) from the old port \(u\)
contributes \(H(T-\{u,w\})\ge1\); hence \(u\) has exactly one
outneighbor \(w\), and the complement is transitive since its count is
one. There can be no four-vertex path starting at \(u\). In that
transitive complement, \(w\) can therefore beat no vertex or just its
final vertex: beating any earlier vertex \(x\) gives a path
\(u,w,x,y\) with a later vertex \(y\). The full old tournament is
transitive in the first case or an ordinal sum of a transitive block
and a directed triangle in the second. Thus \(H\in\{1,3\}\).

It follows that a known scalar orbit witness of height \(H\ge5\)
cannot supply height \(H+2\) by appending a transitive single-seam
reader, regardless of re-rooting before that join. This distinguishes
a true adjustable source from a local adjacent pair constructed at an
unrelated baseline.
