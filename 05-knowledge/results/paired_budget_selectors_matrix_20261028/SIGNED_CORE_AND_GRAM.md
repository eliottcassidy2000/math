# The aligned-pair lift as a signed core, a Gram matrix, and one orientation bit

**Status: PROVED algebra; FINITE-EXACT full order-nine audit.**
The class bijection and aligned-pair generator were supplied concurrently by
Session89 and are inherited here, not claimed as a new result of this packet:

* `../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md`
* `../session89_consecutive_tournament_pairs_20261028/generation/ALIGNED_GENERATION.md`

These arrived in `61aee2f07a` through `6e5aa7f5fb` during this session and
were read before extending them. The older anchors are THM-024
(`01-canon/theorems/THM-024-involution-anti-aut.md`), THM-644
(`01-canon/theorems/THM-644-gridsym-fiber-law-and-anti-redei.md`), and the
bordered spectral recursion THM-1920
(`01-canon/theorems/THM-1920-the-spectral-insertion-response.md`).

## The exact object, including the orientation convention

The source is an aligned pair `(T,S,r)`: two tournaments on the same `n`
vertices and a distinguished extra vertex attached to `S`. Write `A` for
the skew adjacency matrix of `T`, with `A_ij=1` meaning `i -> j`. Let `R`
be the ordinary disagreement graph of `T,S`, and let the symmetric matrix
`C` have diagonal one and off-diagonal entries `C_ij=2R_ij-1`. Finally let
`r_i=1` if the extra vertex beats `i`, and `-1` otherwise.

The inherited lift `W` has root `0` and vertex pairs `p_i,q_i`, with
`p_i -> q_i`, plus layer `T`, minus layer `T^op`, cross layer `C`, and root
signs `r,-r`. Its anti-involution exchanges `p_i,q_i` and fixes `0`.
Its skew adjacency matrix, in that order, is

\[
\mathcal S=
\begin{pmatrix}
0&r^T&-r^T\\
-r&A&C\\
r&-C&-A
\end{pmatrix}.
\]

There are no ties. The isomorphism action is a common permutation of the
core labels. It is not ordinary simultaneous row/column permutation on an
arbitrary square bit matrix: the fixed-order bit split carries an orientation
cocycle. Session89 gives the exact actions and the order-two hostile proving
that no equivariant bijection can identify those two natural quotient problems.

## A sparse signed-tournament core

Define

\[
L=A+C=I+2B.
\]

For each unordered pair `{i,j}`, precisely one of `B_ij,B_ji` is nonzero,
and that entry is `+1` or `-1`. Thus `B` is an ordinary tournament support
with a sign on each oriented edge. Conversely every such signed tournament
determines the aligned pair of core tournaments:

\[
A=(L-L^T)/2,\qquad C=(L+L^T)/2.
\]

This is a concrete four-choice object per unordered pair: direction and sign.
It does not identify the four choices with the four isomorphism classes of
order-four tournaments; that would change the vertices and the quotient.
The diagonal root mask supplies the remaining `n` bits, so the total is
`2 binom(n,2)+n=n^2` as required.

**Cycle-cover law.** For every such core, `det L = 1 mod 8`.
In the determinant expansion, a nonidentity contributing permutation is a
union of directed cycles in the support of `B`. A two-cycle is impossible
because the support is a tournament. Every nonconstant term consequently
contains at least three factors of two. In particular `L` is invertible,
with odd determinant. The full expansion retains signed vertex-disjoint
cycle covers, not just the number of cycles or the support graph's cycle space.

The first cyclic case has three vertices. An acyclic support gives determinant
one. For a directed triangle with edge-sign product `sigma=+1` or `-1`,

\[
\det L=1+8\sigma=9\quad\hbox{or}\quad -7.
\]

The Eisenstein factors are exact:

\[
\det(I+2P_3)=3(1+2\omega)(1+2\omega^2)=3\cdot3=9,
\]
\[
\det(I-2P_3)=(-1)(1-2\omega)(1-2\omega^2)=-7,
\]

where `P_3` is a cyclic permutation matrix. Thus the same phase field carries
the ROOT lower-charge norm `Q(1+2 omega)=3` and the norm-seven factor
`Q(1-2 omega)=7`. The earlier norm-seven element `2-omega` is an associate
of the conjugate factor: `2-omega = -omega(1-2 omega^2)`.
This is a character decomposition of a signed three-cycle, not a numerical
identification of Collatz dynamics with tournament dynamics.

Exact labelled signed-core counts provide a cheap control:

| Core size | Determinant : number of labelled cores |
|---:|---|
| 3 | `-7:8, 1:48, 9:8` |
| 4 | `-31:192, -15:384, -7:512, 1:1920, 9:512, 17:384, 33:192` |

For comparison, the signed-support switch `B -> D B D`, with diagonal
`D_ii=±1`, preserves these determinants. Such a switch is **not** being
identified with the named-port isomorphism action on `(T,S,r)`.

## The anti-involution gives a rectangular Gram factor

Use the orthonormal `+` and `-` eigenspaces of the anti-involution: root and
`(p_i+q_i)/sqrt(2)` on one side, `(p_i-q_i)/sqrt(2)` on the other. Arc
reversal gives `tau S tau = -S`, so both same-sign blocks vanish. The mixed
block is

\[
K=\begin{pmatrix}\sqrt2\,r^T\\-L^T\end{pmatrix},
\qquad K^TK=L L^T+2rr^T=:G.
\]

Consequently

\[
\boxed{\det(xI-\mathcal S)=x\det(x^2I_n+G).}
\]

This follows by a block determinant (first for `x!=0`, then as a polynomial
identity). It also gives a version of THM-1920's insertion formula in the
smaller folded core:

\[
\det(x^2I+G)=\det H+2r^T\operatorname{adj}(H)r,
\quad H=x^2I+LL^T.
\]

Deleting the root gives determinant `(det L)^2`. The reduction therefore
retains an exact integer square root of this marked even principal minor,
up to its determinant sign. No spectral approximation or anonymous relabeling
is used.

## An exact orientation bit, and its limit

With `L` retained, the Gram correction is `G-LL^T=2rr^T`. It determines `r`
up to the common sign. One pinned sign `r_j` recovers every other sign by

\[
r_i=(G-LL^T)_{ij}/(2r_j).
\]

The bit is necessary: at `n=1`, both root signs have `L=[1]` and `G=[3]`,
but their lifts are the transitive triple and the directed triangle. They
have the same skew spectrum and different isomorphism classes. The claim
requires the **labelled matrix** `G` and the retained core `L`; its spectrum
alone is much weaker. This is a precise setting where a shared frame reduces
an orientation ambiguity to one bit. It is not a universal one-bit repair
for arbitrary compressed graph or arithmetic data.

There is another exact orientation carrier. Solving `Ly=r` gives a kernel
vector `(1,y,y)` of `S`. Clear denominators and divide the common content;
every coordinate of the resulting primitive integer vector is odd, since
`L=I mod 2` and all entries of `r` are odd. Conversely `L` and that marked
projective kernel vector recover `r=Ly/x_0`.

For the transitive core with code zero, empty disagreement graph and all
root signs negative, `L` is upper triangular with diagonal one and every
entry above it equal to `-2`. Its kernel vector is

\[
(1,-3^{n-1},-3^{n-2},\ldots,-1,
      -3^{n-1},-3^{n-2},\ldots,-1).
\]

The ternary growth follows from the exact recurrence `y_i=3y_(i+1)`, not
from a resemblance to a number sequence. This gives a simple hostile to
bounded-coordinate representations even when the determinant is one.

## K4–K3,3 bridge: test the discarded bits on actual strong objects

The inherited K4 probe records three vertex-boundary parities and edge
parity. It separates the nine F64 ports and identifies the nine-prime K3,3
bank, but its two-dimensional kernel survives on arbitrary six-bit graphs.

Fix the transitive four-vertex `T` (code zero), the all-negative root mask,
and take disagreement masks `R=0` and `R=45` in edge order
`01,02,03,12,13,23`. Both have four-probe output zero. Both lifted
nine-vertex tournaments are **strong**, yet they are not isomorphic:

| R | Canonical nine-vertex code | det L | Primitive marked kernel |
|---:|---:|---:|---|
| 0 | 10850734724 | 1 | `(1,-27,-9,-3,-1,-27,-9,-3,-1)` |
| 45 | 275776032 | -31 | `(31,-17,-19,7,-11,-17,-19,7,-11)` |

This moves the old raw-bit collision into actual strong-atom constructions.
The four-probe readout cannot replace the two kernel bits in a general
attachment grammar. It also shows why matching a square-class code or an
anonymous incidence shape is not sufficient to preserve a later spectral
or orientation test.

There is a further exact dyadic boundary. With this fixed transitive core,
the entire four-probe zero fibre is `R in {0,30,45,51}`. Its determinants
are respectively `1,1,-31,-31`. A determinant observation modulo 64 recovers
one of the two lost bits, while modulo 32 recovers neither. Even the exact
determinant leaves two possibilities in each pair. This is a test of a
specific enlarged observation language, not a claim that the bit can be
computed from the old four probes alone.

More generally `det(I+2B) mod 2^d` only needs cycle-cover terms occupying
fewer than `d` vertices. For example modulo 32 it is
`1 + 8 sum_(directed triangles) product(signs)
 - 16 sum_(directed four-cycles) product(signs)`.
This is a natural arithmetic filtration on signed tournament cores. It
supplies a finite-order statistic and does not reconstruct all core labels.

## Exact category and counting sidecars

Aligned pairs with common-label isomorphisms are equivalent to marked lifts
`(W,tau)` with mark-intertwining isomorphisms. Forgetting `tau` gives a
bijection on isomorphism classes by Sylow conjugacy, as proved in Session89,
but it does not preserve automorphism groups:

\[
\operatorname{Aut}(T,S,r)=C_{\operatorname{Aut}(W)}(\tau),\qquad
\#\{\text{anti-involutions}\}=
\frac{|\operatorname{Aut}(W)|}{|C_{\operatorname{Aut}(W)}(\tau)|}.
\]

This is the exact category-theoretic boundary: the marked correspondence is
an equivalence of groupoids; the unmarked forgetful functor is not full.
Weighted orbit counts must keep the centralizer index. A class bijection
alone cannot transport inverse-automorphism weights or named port maps.

`check_lift_gram.py` extends the incoming lift audit from order seven to
**all 2,752 aligned-pair classes with core size four**. It compares their
image against the complete independent **191,536-class order-nine bank**,
filtered by full converse isomorphism, and gets exact equality. The four
core contributions are `1024,352,1024,352`. It checks every named
centralizer, all the Gram block identities, the odd integer kernels, and
the inherited 81/27 transporter control. At one nine-vertex lift the marked
group has order three while the full group has order 81: forgetting the
mark hides 27 anti-involutions. The receipt lists every group profile.

This extends the validated output scope. It does not claim a new asymptotic
algorithm for A000568, and it does not compete with the incoming measured
2.08x core-first generator on a different output contract.
