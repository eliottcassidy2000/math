# A genuine fixed-degree triad: fields, companion dynamics, recurrence streams

**Status:** PROVED elementary algebra below; FINITE-EXACT corroboration at
orders 3, 4, 5. This is a new Session128 instantiation, not an authenticated
recovery of a historical theorem named “triadic equivalence functors.”
See `RECOVERY.md` for the exact inherited sources and their weaker status.

## 1. Fix the source before comparing its representations

Fix one of k=3,4,5 and

    f_k(t)=t^k−t^(k−1)−...−t−1,
    K=Q[α]=Q[t]/(f_k),   A=Z[α].

Session127 certifies irreducibility for these three orders. Nothing below
identifies different k. The same free-module construction works for the
quotient algebra of any monic f with nonzero constant term, but the term
“number field” then requires irreducibility separately.

Use the ordered power basis e=(1,α,...,α^(k−1)), with **column** coefficient
vectors. Multiplication by α is the integral companion matrix C with ones
on the subdiagonal and ones in the last column. Thus

    C_3 = [[0,0,1],[1,0,1],[0,1,1]].

Since α is a unit, C and C^−1 are integral. This is the same exact unit
source as Session127, not the norm scalar or a selected real approximation.

## 2. Three categories, with all morphisms specified

The following are concrete skeletal additive categories, indexed by r≥0.

1. **Field category KFree.** Object r is K^r. Morphisms r→s are all
   s-by-r matrices over K, acting K-linearly.
2. **Companion category Comp.** Object r is Q^(kr) with time operator
   C_r=I_r⊗C. Morphisms r→s are **all** rational matrices H satisfying
   H C_r=C_s H.
3. **Recurrence category Rec.** Let S be the Q-vector space of bilateral
   streams a:Z→Q satisfying

       a_(n+k)=a_n+a_(n+1)+...+a_(n+k−1), for every n∈Z.

   Its time operator is the forward shift E, (Ea)_n=a_(n+1).
   Object r is S^r. Morphisms r→s are **all** Q-linear maps commuting with
   these shifts. The recurrence is reversible, so arbitrary initial values
   a_0,...,a_(k−1) determine exactly one bilateral stream.

These are not categories of arbitrary sequences, arbitrary matrices,
symbolic substitution words, or arbitrary topological dynamical systems.

### First functor and inverse: coefficient realization

Replace each entry b in a K-matrix by the k-by-k rational multiplication
matrix M_b. On objects this realizes K^r in its named power basis. Denote
this functor F:KFree→Comp.

It preserves identity, sum and composition because M_1=I and M_(ab)=M_a M_b.
Every k-by-k rational matrix H commuting with C is M_b for a unique b:
write H e_0=v and let b have coefficient vector v. The cyclic basis
C^j e_0=e_j gives H e_j=C^j v=M_b e_j for 0≤j<k. Apply this argument to each
block of an arbitrary rectangular intertwiner. Hence F has an explicit
inverse: read the first column of each block as a field element. It is
fully faithful and onto the stated objects.

### Second functor and inverse: observable streams

For v∈Q^k define

    (Lv)_n=e_0^T C^n v,  n∈Z.

Cayley–Hamilton gives the recurrence, and L C=E L. This is not “keep one
scalar”: it retains the **entire** scalar stream, or equivalently k initial
samples. Let O have rows e_0^T C^j, 0≤j<k; initial samples are s=O v.

Write v=(c_0,...,c_(k−1)). The explicit inverse is

    c_0=s_0,
    c_(k−j)=s_j−sum_(i=1)^(j−1) s_i,  1≤j≤k−1.                 (1)

Indeed the coefficient update under multiplication by α is
c'_0=c_(k−1), c'_i=c_(i−1)+c_(k−1). Iterating this update gives
s_j=c_(k−j)+sum_(i=1)^(j−1)s_i for 1≤j≤k−1. This proves (1), not merely a
finite numerical observation. Both O and O^−1 are integral, and

    det O = (−1)^((k−1)(k−2)/2).

Thus L is a Q-linear isomorphism and also identifies the integral lattice
A with the integral recurrence streams (all k initial samples integral).
For r copies write L_r componentwise. Define G:Comp→Rec on a morphism H by

    G(H)=L_s H L_r^−1.

This commutes with E because H intertwines C. Conversely every shift-
commuting map has exactly this form. Composition and identity are preserved
by cancellation of L_r^−1 L_r. Its inverse is the displayed conjugation
in reverse, using (1).

Consequently

    KFree ≅ Comp ≅ Rec

are genuine isomorphic concrete categories with these chosen objects.
The round-trip natural transformations are identity maps. The nontrivial
coordinate maps are explicitly natural: G(H)L_r=L_s H for every H, and
field multiplication realization commutes with every K-linear matrix.
No unspecified natural transformation or choice of “similar numbers” is used.

## 3. The exact paper convention and substitution convention

For ascending initial samples s=(s_0,...,s_(k−1)), forward shift has matrix
R with superdiagonal ones and last row all ones. Let J reverse coordinates,
and put B=J O. Then

    O C=R O,       B C=M B,
    M=J R J = first row all ones, subdiagonal ones.

In particular

    O_3=[[1,0,0],[0,0,1],[0,1,1]],
    B_3=[[0,1,1],[0,0,1],[1,0,0]],
    M_3=[[1,1,1],[1,0,0],[0,1,0]].

M_3 is the companion convention reported by the independent attached-paper
audit. This packet uses only the displayed exact matrix, not an unverified
claim from that paper.

The spectral/substitution lane uses another integral basis U. Its jth
column (1≤j≤k) is the power-basis vector of

    ℓ_j(α)=α^(j−1)−sum_(r=0)^(j−2)α^r.

It has diagonal 1, and U M=C U. The convention difference is exact:

    U B=C^−1,          U^−1=B C.                            (2)

For a proof, U B commutes with C and sends e_0 to the last column of U,
which is α^−1. The centralizer argument above proves (2). Thus the
substitution-coordinate realization U^−1 v equals the reversed-stream
state of **α times the field element**, a single time shift. It must not be
silently identified with Bv.

## 4. Transporting the GSMA sparse gate without losing it

The additive sparse lattice Z+Zα is the condition c_2=...=c_(k−1)=0.
Equation (1) makes this exactly

    s_1=...=s_(k−2)=0.

At such a stream the signed binary norm pair is

    p=s_0, q=−s_(k−1),       x=p−qα,
    N_K/Q(x)=F_k(p,q).

In the substitution basis the same sparse element has coordinates

    U^−1 v=(p−q,−q,0,...,0),
    p−qα=(p−q)ℓ_1(α)−qℓ_2(α).

Thus the Session127 shear r=p−q is literally its first substitution-basis
coordinate. At k=3 its norm is the exact transformed binary form

    F_3(r+q,q)=r^3+2r^2q−2q^3.

This is an actual coordinate map, not an identification of different
number fields through a shared integer.

This gate is an observer/predicate on the category's objects, not a K-linear
submodule or a morphism into one unnamed common category. Its transported
version retains the same source element.

For the distinguished inverse-unit orbit x_n=α^−n:

- n=k has coefficients (2,−1,0,...,0), initial stream (2,0,...,0,−1).
- For k=3, n=16 has coefficients (−103,56,0), initial stream (−103,0,56),
  and reversed paper state (56,0,−103).
- The corresponding substitution state U^−1 v is (−47,56,0), as (2)
  predicts. This is the same field element in a different named basis.
  Its negative coordinate means it is **not** a letter-count vector of any
  finite positive word. Matrix equivalence need not preserve the positive
  word cone or word order. After sign normalization, 103−56α has substitution
  coordinates (47,−56,0); this explains the 47 coordinate without identifying
  the Lucas quadratic source with the cubic field.

A pointed version of each category (object plus a chosen vector, morphisms
preserving that vector) retains the source 1 and therefore the whole inverse-
unit orbit. Norm is preserved under change of representation as a determinant
of multiplication, but norm alone is **not** an equivalence: 1 and α^−2 are
distinct norm-one elements for every requested k, and α^−2 fails the sparse
gate. Also 0 and α have equal initial constant observation, so a single
sample has no inverse. The full k-sample observation is the necessary sidecar.

## 5. Exact algebraic feedback-loop solvability over the fixed field

For an r-by-r matrix D over K and its rational realization F(D),

    det_Q F(D) = N_K/Q(det_K D).                              (3)

Proof: over a splitting field, separability splits
K⊗C as one copy of C for each embedding σ:K→C. On these components F(D)
becomes the block diagonal matrix with blocks σ(D). Its determinant is
prod_σ det σ(D)=prod_σ σ(det_K D). This is (3). The same determinant is
obtained after conjugation into recurrence coordinates.

Hence for γ∈K and K-linear feedback D,

    I−γD is invertible  iff  det_K(I−γD)≠0
                        iff  N_K/Q(det_K(I−γD))≠0.

This is exact; it has no assumption |γ| small or spectral radius below one.
For example D=2I, γ=1 gives I−γD=−I, invertible although its geometric
Neumann series does not converge in the usual norm. The inverse is transported
by the same functors. Precisely, this solves the algebraic feedback loop
`x=u+γD x` uniquely. Invertibility does not assert dynamical stability, nor
that arbitrary additive feedback in an evolution can be absorbed by a
conjugacy. The native gain quotient and any moving-frame absorption theorem
require their own hypotheses.

The hypotheses matter. A rational feedback matrix is in the fixed-field
category exactly when it commutes with C_r. D=diag(1,0,...,0) on one companion
block does not commute with C; no K determinant may be assigned to it by
pretending it is K-linear. Real parameters outside K require a stated scalar
extension, not an unqualified number-field norm. Native sign/gain quotients
from the feedback lane are separate maps; retaining a positive scale can
restore amplitudes, but their sign itinerary is not the whole field module.

## 6. Spectral and nonlinear boundaries

Keeping all embeddings gives the usual invertible Vandermonde spectral
coordinate change after scalar extension. Discarding the expanding embedding
to obtain an internal Pisot contraction is a real-linear projection, not an
equivalence of the full real dynamical spaces. On K alone an individual
embedding is injective, but this does not make the projected lattice discrete
or supply a continuous inverse for its real completion.

Likewise an incidence matrix forgets word order and does not recover the
symbolic substitution. Attractor digits, graph addresses, topology and
admissibility are additional data. This triad transports the exact linear
source; it does not by itself classify its projected attractor.

The paper audit also identifies the sine identity sin(3θ)=3sinθ−4sin³θ.
It is a semiconjugacy of circle tripling to the interval cubic. Sine has a
two-point generic fiber; only the explicitly named involution quotient
θ∼π−θ becomes conjugate to the interval system. This is a different quotient
construction from the field/recurrence equivalence and does not merge their
categories automatically.

## 7. Verification and stopping boundary

Run `python3 -B check_triadic.py` and `python3 -B -O check_triadic.py`.
Their outputs agree. The checker uses only the Python standard library and
explicit exceptions, not assertions. It checks the coefficient cube
{−1,0,1}^k, all three matrix conventions, centralizer dimensions, inverse and
composition laws, independent polynomial multiplication versus scalar
recurrence, sparse-gate transport, source powers n=0..24, and 625 rank-two
K-matrix determinant/norm identities per order. The finite source range is
not a new global unit classification. The proof above, not the count of
checks, supplies the general equivalence and determinant identity.
