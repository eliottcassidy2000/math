# The binary descent obstruction sees only parity of a transitive port

The earlier [BSD component recursion](../bsd_component_recursion_20261013c.md)
uses the row-zero tournament matrix B over F2 and
Phi(B)=B²+B+I. In the specified family of distinct primes all 3 mod 8,
its nullity is the pure 2-Selmer dimension under the cited Monsky
identification. Here we prove an exact block-elimination rule for
**selective** transitive attachments, and compare it to the growing
Hamilton count. No conclusion about elliptic-curve rank or BSD follows.

## Selective attachment, not just ordinal concatenation

Let an appended transitive module V=TT_m beat old labels in C and
lose to the other labels. Put c=1_C, r=1−c as column vectors,
delta=m mod 2, eps=|C| mod 2, and let D_m be the row-zero
matrix of TT_m. The enlarged row-zero matrix is
```
A = [ B + delta diag(r),  r 1_m^T      ]
    [ 1_m c^T,           D_m + eps I ].
```
Extend scalars to F4 and choose a root lambda of lambda²+lambda+1.
The lower-right block of A−lambda I is invertible: D_m is
upper triangular with entries in F2, whereas lambda is not in F2.
Also D_m 1_m=0, so its inverse applied to 1_m is
(eps−lambda)^(-1) 1_m. Eliminating the module gives exactly
```
B − lambda I + delta diag(r) − delta/(eps−lambda) r c^T.
```
This controls the lambda eigenspace, and its Frobenius conjugate
controls the other eigenspace. Since Phi is squarefree, the F2
nullity of Phi(A) is twice the F4 dimension of the displayed
reduced kernel. Thus it depends only on the **parity of m**:
for even m it is the old nullity; for odd m it is the nullity
after one appended vertex with the same cut. This proof retains
the actual port term; dropping it would not be a valid ordinal
cancellation for odd m.

A companion identity makes the nullity compatible with global
converse, even before taking rank. If J is all ones, the row-zero
converse matrix is A+J+(n mod 2)I. Use AJ=0, JA=(n−1)J,
J²=nJ and Phi(X+I)=Phi(X). Expanding gives exactly the same
Phi matrix. These identities are purely finite binary algebra.

## Native separations on the selective rays

For the [order-13 selective ray](SELECTIVE_RAYS.md), the two
base nullities are 4 and 2. With one new vertex they are 2
and 0. Hence the exact equal-H strong ray keeps a two-dimensional
difference for every m. For the prime-derived mask6144 ray the
base nullities are 2 and 6, and one new vertex gives 2 and 6
again. Thus at every length its two equal-H strong witnesses
have different nullities, even though their named H-return
interface is matched. The arithmetic bridge therefore adds a
genuinely independent observer, not an explanation of H.

For the six [Fano metagraph twists](FANO_METAGRAPH.md), masks
0,1,3,7,11,13 have nullities 2,0,0,0,4,2 respectively.
These are exact binary rank calculations. Through the inherited
prime-realization theorem one can realize these finite orientations
in the specified congruent-number descent setting; this does
not identify the elliptic curves or their ranks from the slot
counts. The F4 elimination and F8 port geometry are different
finite-field roles.
