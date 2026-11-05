# A Berggren–Walsh compiler for physical source words

**PROVED:** three exact signed source-word response sums, two exact C2
bilinear sums, and the exact boundary at which each signed aggregate is itself
a physical source response. **NOT a coverage theorem:** signed summation does
not preserve physical membership, scalar support, or graph isomorphism classes.

## 1. Inheritance and the actual bridge

The source is THM-3357,
`01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md`,
especially its signed branch-sum transfer and Section 9 procedure for another
ternary system. Its surrounding THM-4057 Stern–Brocot work distinguishes
ordered words, scalar depth, and genuine inverse-tree coverage. The useful
mechanism is distributivity of matrix multiplication over an indexed word
language, not an identification of tournament sources with Pythagorean triples.

The target is the complete physical source grammar in Session 198
`bilinear-descent/PROOF.md`, Sections 1–3. Let b=(0,0,0,1)^T and use coordinates

    y=(a,f,e,p),
    (D,E,U,d)=(a+2e+p,e,a+f+2e+2p,e+p).

The three physical one-vertex operations act by

    O y=(a+e,a+2f+e+2p,0,e+p),
    Z y=(2a+f+e+p,f,e+p,0),
    J y=(f+p,0,a+e+p,e).

A chronological word w=w1...wk means y_w=G_wk...G_w1 b. Each word is a
literal rooted tournament with one root-start Hamilton path. Words are the
indices of the sums below: equal source responses and isomorphic graphs are
NOT deduplicated.

The map sends an indexed ternary physical-word family to its four parity
Fourier sums. It preserves signed linear observables and, with two independent
word families, signed bilinear observables. It destroys the individual sources,
positivity, and scalar support. The full indexed physical words are the sidecar
needed to recover those data. The cheapest decisive hostile is the length-two
aggregate in Section 3; graph-isomorphism loss is visible already at that depth.

## 2. Three closed signed sums

Define

    S_Z(k)=sum_{|w|=k} (-1)^(number of Z in w) y_w,
    S_J(k)=sum_{|w|=k} (-1)^(number of J in w) y_w,
    S_ZJ(k)=sum_{|w|=k} (-1)^(number of Z or J in w) y_w.

**Theorem 1.** For every integer k>=0,

    S_Z(k)=(0,2k,0,1),
    S_J(k)=(0,3^k-1,0,1),

    S_ZJ(k)=(-2k,0,2k,1)             if k is even,
             (2k-4,2,-2k,1)          if k is odd.

Proof. Distributivity gives S_Z(k)=(O-Z+J)^k b and analogously for the other
two signs. The three transfer matrices, in the displayed coordinate order, are

    K_Z =[[-1,0,0,0],[1,1,1,2],[1,0,0,0],[0,0,2,1]],
    K_J =[[3,0,2,0],[1,3,1,2],[-1,0,0,0],[0,0,0,1]],
    K_ZJ=[[-1,-2,0,-2],[1,1,1,2],[-1,0,-2,-2],[0,0,0,1]].

Each formula equals b at k=0. Multiplication by its displayed matrix gives
the formula at k+1, switching the two parity cases for S_ZJ. This proves all
k by induction. No finite computation is a premise. QED.

Together with S_0(k)=(O+Z+J)^k b, these are all four Fourier transforms for
the two parities (#Z mod2,#J mod2). The sum over one specified parity class
is recovered by the usual four-term sign average. This is an exact indexed
word compiler. In particular the difference of the two Z-parity population
sizes is 1, while their difference of balanced readout sums is 4k+5:

    sum_{|w|=k} (-1)^#Z (D_w+2U_w)=4k+5.             (1)

## 3. The exact physical boundary of the aggregate

**Theorem 2.** The signed aggregate S_Z(k) is an actual nontrivial source
response if and only if k+1 is a power of two. S_J(k) is an actual source
response if and only if k is 0 or 1. S_ZJ(k) is an actual source response if
and only if k=0.

Proof. Translating S_Z and S_J to the original response coordinates gives

    old(S_Z(k))=(1,0,2k+2,1),
    old(S_J(k))=(1,0,3^k+1,1).

An actual tournament with H=1 is transitive: in its ordered strong-component
decomposition every nonsingleton strong component would contribute at least
three Hamilton paths. In a transitive start-one source the root must be first,
and each choice of the vertices assigned to the optional first path gives
one reserved cover, so U=2^(n-1). Conversely the transitive source of order
n has response (1,0,2^(n-1),1). This proves the claimed criterion for S_Z.
For S_J, k=0,1 work. If k>=2 is even, 3^k+1 is 2 mod8 and is greater than2;
if k>=3 is odd, it is 4 mod8 and is greater than4. Neither is a power of two.

For S_ZJ, odd k gives negative D=-2k-3, so is impossible. Even k>=2 gives
old(S_ZJ(k))=(2k+1,2k,2k+2,2k+1), which violates the necessary physical
inequality D-E>=d. The k=0 vector is b. QED.

A particularly clean hostile is k=2: the signed sum of nine actual source
responses is (1,0,6,1), which is not a source response. At k=4 equation (1)
gives 21; it is a signed readout sum, not a construction of a source or graph
with that Hamilton-path count.

There is also a positive exact dyadic compiler. If k=2^m-1, then

    sum_{|w|=2^m-1} (-1)^#Z y_w=O^m b.              (2)

Both sides are exactly the response vector of a physical transitive source.
The equality identifies the signed sum, not its individual terms or their
positive scalar support.

## 4. Pairing with the complete C2 source forms

Let P_N,P_M be the positive C2 matrices in Session 198 `bilinear-descent`,
Section 4, and write F(u,v)=y_u^T P_F y_v. Every summand here is the count of
an actual C2 tournament. Bilinearity and Theorem 1 give, for every r,s>=0,

    sum_{|u|=r, |v|=s} (-1)^(#Z(u)+#Z(v)) N(u,v)
      =12(r+1)(s+1)+1,

    sum_{|u|=r, |v|=s} (-1)^(#Z(u)+#Z(v)) M(u,v)
      =24rs+26r+26s+29.                            (3)

Indeed only the f,p rows survive: N has that 2-by-2 block
[[3,6],[6,13]], and M has [[6,13],[13,29]]. Substituting (2r,1) and
(2s,1) proves the formulas. More generally, any fixed linear or bilinear
source readout admits the same exact transfer, including the C3 matrices.

Equation (3) is an independent exact checksum for indexed source-word
computations. It does not say that its right side belongs to the corresponding
positive image. At dyadic depths on both sides, equation (2) does supply a
physical transitive-pair realization, but those values were already supplied
by that literal pair.

## 5. Why the sign is not a graph invariant

The chronological words OJ and JZ have the same response (5,1,5,1) and give
rooted-isomorphic graphs. Their explicit out-neighbor bitmasks and roots are

    OJ: [6,12,8,1], root3;
    JZ: [10,12,1,4], root2.

The bijection (0,1,2,3)->(0,1,3,2) preserves every arc and the root. But OJ
has Z-parity zero and JZ has Z-parity one. Thus no graph-isomorphism invariant
can equal this word sign on every presentation. Algebraically the overlap
is J O=Z J, where multiplication acts on column vectors from right to left.

The source grammar is not a free tree like the primitive Berggren tree.
The bridge transfers linear word summation, not unique ancestry, Euclidean
inverse branches, disjoint parent domains, or universal scalar coverage.
This identifies both the valid transferred operation and its first failed
extension.

## 6. Independent verification

Run `python3 verify_walsh.py` and `python3 -O verify_walsh.py`. Outputs agree.
The standalone verifier imports no repository code or external package. It
checks the stated recurrences, explicitly sums every word through depth8,
computes direct subset-DP path and reserved-cover inventories for all364 words
through depth5, checks the double C2 identities, tests the dyadic compiler,
and verifies the root-preserving graph isomorphism above. All gates raise
exceptions explicitly and survive optimized Python. Finite checks validate
implementation; the proofs above establish the all-order statements.

The requested all-odd spectrum theorem remains OPEN. The contribution is a
new exact compiler, an independent source-enumeration audit, and a sharp
physical-realizability boundary for a tempting signed-aggregation bridge.
