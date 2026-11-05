---
id: THM-4465
title: "Weighted tournament shear production and contact kernel"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. For an oriented
  nonnegative zero-diagonal matrix M, 4 tr(sym(M) skew(M)^2) is three times
  the weighted cyclic-triple sum minus the weighted transitive-triple sum.
  This is the instantaneous Euler skew-norm production in the stated
  coordinate cone. Block-constant substitution closes on size, internal
  edge mass, and production; nonuniform exterior contacts require an
  additional quadratic boundary kernel. Identical scalar block states can
  give opposite full signs +2/9 and -2/9. No Euler trajectory, blowup,
  arbitrary-weight amplitude bound, or LRC result is asserted.
source: tournament-continuation-20260908
depends_on: []
related:
  - THM-462-the-cubic-spectrum-of-tournaments-is-gap-free
  - THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels
  - THM-1862-order-join-reduction-principle
  - THM-1926-tournament-zeta-euler-product-over-strong-core
  - THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law
  - THM-2013-coordinates-for-the-continuum-cyclic-temperature
  - THM-2016-the-deep-continuum-and-the-reducibility-ceiling
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-4457-euler-sharp-transverse-shear-distance-budget
  - THM-4466-sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes
script: 04-computation/tournament_shear_cycles_20260908.py
output: 05-knowledge/results/tournament_shear_cycles_20260908.out
script_sha256: 9e30bac8f84a4a3403dd4a65ba5c9e06a3fc854160eb3a5e6a88eb353b5bfd7e
output_sha256: c2566b9c0d2ec1e72ffe856e909abf8db58e7d0f8b4d60f951e995cd8f6b717b
semantic_sha256: 44cb170f8d233865187066a7744e27af0131c6ecb88e4b500b6527c4425e858a
hash_basis: raw LF bytes
audit: >
  Independent root review accepted the complete proof, trace and vorticity
  constants, cut equality cases, weighted reducible positive example,
  substitution compiler, and opposite-sign contact kernel witness. The
  115991 explicit gates survive -O. Normal, optimized, and stored LF output
  agree. Exact matrix multiplication independently checks triangle sums,
  and all three-dimensional weighted controls check the curl contraction.
---

# THM-4465 -- weighted tournament shear production and contact kernel

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The root
independently checked the complete proof and accepted it without a
mathematical correction. The unit-score and reducibility bounds below are
explicitly inherited corollaries. The weighted substitution compiler and
the sign-changing boundary witness identify the additional state needed by
this operation. These are instantaneous matrix statements, not an Euler
blowup theorem.

## Inheritance and active concepts

The closest proved operation mechanism is
[THM-2195 / transitive-quotients-exactly-control-universal-substitution-products](THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md):
fixed-block transport needs its exterior incidence response. The least-used
sidecar here is the **weighted exterior quadratic response of one block**.
[THM-1862 / order-join-reduction-principle](THM-1862-order-join-reduction-principle.md)
and [THM-1926 / tournament-zeta-euler-product-over-strong-core](THM-1926-tournament-zeta-euler-product-over-strong-core.md)
localize directed cycles to strong components. The new production functional
also charges transitive triples, so that localization does not eliminate
the interfaces.

The unit-count background is
[THM-462 / cubic-spectrum-of-tournaments-is-gap-free](THM-462-the-cubic-spectrum-of-tournaments-is-gap-free.md),
[THM-2013 / coordinates-for-the-continuum-cyclic-temperature](THM-2013-coordinates-for-the-continuum-cyclic-temperature.md),
and [THM-2016 / deep-continuum-and-reducibility-ceiling](THM-2016-the-deep-continuum-and-the-reducibility-ceiling.md).
Their score-variance and maximum-reducible-cycle results are not new claims
here. The corrected near miss is MISTAKE-220: a sum of strong-component
cycle counts cannot be bounded by its largest summand. The repaired
THM-2016 argument uses discrete convexity. MISTAKE-217 also prevents
confusing a signed coordinate with its absolute value.

The older triangle atom in
[THM-1805 / vandermonde-is-a-signed-tournament-sum-intransitivity-cancels](THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md)
has a different observer. The existing substitution result
[THM-1960 / tournaments-compose-from-regular-seeds-the-spectral-substitution-law](THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md)
requires regular blocks to decouple its skew spectrum; the cubic compiler
below works with arbitrary weighted blocks but different retained data.

The portfolio is anchor: actual tournament interpretation of shear
production; niche: weighted strong-core interfaces; wildcard: failure of a
three-scalar block state to preserve the sign. The updated concept board is:

| Concept | What the calculation changes |
|---|---|
| Cyclic temperature | Unit weights give a familiar score statistic with a new matrix meaning |
| Strong-core cycle localization | Crossing transitive triples retain a negative production cost |
| Substitution operation | Constant contacts close on size, internal edge mass, and production |
| Nonuniform boundary contacts | Those scalars lose the sign; a quadratic boundary kernel restores it |
| Euler skew evolution | The functional is actual instantaneous skew-norm production, with a basis/trajectory sidecar |

Targeted searches of the named old mechanisms and correction matches found
the unit-count overlap above. No global priority claim is made.

## 1. Exact trace and genuine tournament interpretation

Let M be an n-by-n nonnegative real matrix with zero diagonal and with at
most one of M_ij,M_ji nonzero for each pair. A strictly positive magnitude
on exactly one direction gives a weighted tournament; zeros are permitted
as boundary points and do not receive artificial tie orientations.

Define

```
S=(M+M^T)/2,              K=(M-M^T)/2,
F(M)=4*tr(S*K^2).
```

The vertices are coordinate axes; the pairwise observable is the ordered
off-diagonal matrix entry, in a fixed orthonormal coordinate gauge. Write
w_ij=M_ij+M_ji for its nonnegative edge magnitude. For unordered triples
let sigma=3 for a cyclic triple and sigma=-1 for a transitive triple.
If an edge is zero, its product contribution vanishes regardless of its
unused orientation label.

**Proposition 1 (PROVED).**

```
F(M)=tr(M^3)-tr(M^2*M^T)
    =sum_(i<j<k) sigma_ijk*w_ij*w_ik*w_jk.              (1)
```

Indeed expansion and cyclicity of the trace give
8 tr(SK^2)=2 tr(M^3)-2 tr(M^2 M^T). A cyclic triple contributes its product
three times to tr(M^3), while a transitive triple has exactly one directed
two-edge path closed by a forward edge and contributes once to
tr(M^2 M^T). No repeated-index term survives the zero diagonal and absence
of opposite nonzero entries.

For an incompressible Euler velocity gradient along a particle, its skew
part obeys K'=-(SK+KS), because the pressure Hessian is symmetric. Therefore

```
d/dt ||K||_F^2 = 4 tr(SK^2)=F(M)                       (2)
```

at an instant when the gradient belongs to this weighted tournament cone.
The trace identity itself requires no PDE assumption. In dimension three,
with the usual curl convention

```
omega=(M_32-M_23, M_13-M_31, M_21-M_12),
K^2=(omega*omega^T-|omega|^2 I)/4,
```

the zero trace of M gives F(M)=omega^T M omega. Thus (2) agrees with the
usual vorticity stretching term. No higher-dimensional vector curl is
invented: in larger dimension (2) is the skew-matrix norm statement.

General orthogonal changes of basis preserve the trace but can leave the
nonnegative tournament cone. The triangle coloring is a coordinate-cone
description, not a basis-free coloring of an arbitrary gradient. A full
Euler trajectory need not remain in this cone.

## 2. Unit production and the inherited reducibility threshold

For unit weights, with c_3 directed cyclic triples and outdegrees d_i,

```
F(T)=4c_3(T)-binom(n,3)
    =n(n-1)/2-2*sum_i (d_i-(n-1)/2)^2.                 (3)
```

This is the classical cycle/score identity already proved and used by
THM-462/2013: every transitive triple has one source, so their count is
sum_i binom(d_i,2); substituting sum_i d_i=n(n-1)/2 gives (3). The minimum
centered square sum is zero for odd order and n/4 for even order. Regular
odd tournaments exist by cyclic orientation, and deleting one vertex of a
regular odd tournament gives an even near-regular one. Thus the maximum is
n(n-1)/2 at odd regular order and n(n-2)/2 at even near-regular order.

For a directed cut T=A join B, put a=|A|, b=|B|, n=a+b. Every mixed
triple is transitive, so

```
F(T)=F(A)+F(B)-ab(n-2)/2
 <= -[n(a-1)(b-1)+a*1_(a even)+b*1_(b even)]/2.        (4)
```

The bound follows by inserting the odd/even maxima from (3); the algebra
uses a(a-1)+b(b-1)-ab(n-2)=-n(a-1)(b-1).

**Corollary (PROVED; inherited bound, refined equality description).** Every
reducible unit tournament has F<=0. Equality holds exactly for the order
join of a singleton and an odd regular tournament, in either order. Order
two is included by using a one-vertex regular block. Thus F>0 implies
strong connectivity.

For equality, (4) forces one side to be a singleton and the other side odd
regular. Conversely that join gives zero. A regular tournament of order
greater than one is strong: a proper dominating cut would make the average
global outdegree of its dominating block exceed the common global degree.
So the equality description also specifies the entire strong-component
structure. The inequality itself is a direct corollary of the earlier
THM-2016 reducibility ceiling, not a separate new extremal theorem.

Strong connectivity is insufficient. At order six, orient i->j for i<j
except replace 1->6 by 6->1. This tournament is strong and has exactly the
four cyclic triples {1,j,6}, j=2,3,4,5. Its production is 16-20=-4.
The complete labeled census through order five has no nonpositive strong
example of order at least three; this minimality statement is FINITE-EXACT.

## 3. Strong components retain a weighted exterior cost

Every directed cyclic triple lies inside one strong component. Consequently
for arbitrary nonnegative weights on an oriented tournament,

```
F(T)=sum_(strong components B) F(B)
     -sum_(triples meeting multiple components) product(edge magnitudes).
                                                               (5)
```

The second term is nonnegative. This is why the zeta/cycle strong-core
localization does not make F a strong-core-only observable. Discarding a
source, sink, or transitive block can remove exactly the production cost
that determines the sign.

There is a useful closed operation law when exterior contacts are constant
on blocks. Let Q be a quotient tournament, let blocks B_i have sizes n_i,
internal edge masses e_i=sum_(u<v in B_i) w_uv, and productions F_i. For
each quotient pair ij, all edges between the two blocks have magnitude
lambda_ij>=0 and the orientation prescribed by Q. Then

**Proposition 2 (PROVED, exact substitution law).**

```
F(Q[B_1,...,B_q])
 =sum_i F_i
  -sum_(i<j) lambda_ij^2*(n_j e_i+n_i e_j)
  +sum_(i<j<k) sigma_ijk*n_i*n_j*n_k
                         *lambda_ij*lambda_ik*lambda_jk.        (6)
```

To prove it, classify a triple by the number of blocks it meets. One block
gives F_i. Two vertices in B_i and one in B_j always form a transitive
triple, contributing -lambda_ij^2*n_j*e_i. Three different blocks inherit
the quotient type and have n_i n_j n_k choices. This proves every term.

The corresponding state update is

```
n_total=sum_i n_i,
e_total=sum_i e_i+sum_(i<j) lambda_ij*n_i*n_j,
F_total given by (6).                                      (7)
```

Thus (n,e,F) is an exact finite-dimensional state for recursive
substitution with **block-constant contacts**, even when the blocks are
irregular. The regularity premise of THM-1960's spectral splitting is
unnecessary for this different observer.

The unit strong-connectivity consequence does not extend to arbitrary
positive weights, even with interior weights bounded by the exterior
ones. Replace each vertex of C_3 by T_2 with its internal edge weighted
epsilon, retain weight one on the twelve cross-block edges, and adjoin a
dominant vertex with all six incident weights one. Formula (6) gives

```
F(core)=24-12epsilon,      e(core)=12+3epsilon,
F(whole)=12-15epsilon.                                    (8)
```

At epsilon=1/2, the seven-vertex tournament is reducible, every weight is
in [1/2,1], and F=9/2>0. This refutes both an arbitrary-weight extension
and the proposed repair “exterior weights dominate internal weights.”
The strongest surviving statement is the exact weighted cost (5)/(6).

## 4. Nonuniform contacts require the quadratic boundary kernel

Let B be one weighted block and attach a new source or sink whose edge to
vertex u has magnitude b_u>=0. Define

```
Q_B(b)=sum_(u<v in B) w_uv*b_u*b_v.
```

Every new triple is transitive, hence the exact interface law is

```
F(B with the attached source/sink)=F(B)-Q_B(b).             (9)
```

Equivalently Q_B(b)=b^T W_B b/2, where W_B is the symmetric zero-diagonal
matrix of internal edge magnitudes. This matrix is not claimed positive
semidefinite. Its value on nonnegative b is nonnegative, which is the
property needed in (9).

**Proposition 3 (PROVED obstruction, exact smallest block in this model).**
The state (n,e,F), even together with the labeled boundary vector b, does
not determine the sign after nonuniform attachment. Take the fixed cycle
1->2->3->1 and respectively assign internal magnitudes

```
B_1: w_12=2, w_23=1, w_13=1,
B_2: w_12=1, w_23=1, w_13=2.
```

Both have state (3,4,6), the same orientation and the same internal weight
multiset. Attach the same dominant vertex with

```
b=(2/3,4/3,2).
```

Then

```
Q_(B_1)(b)=52/9,         Q_(B_2)(b)=56/9,
F(whole_1)= 2/9,         F(whole_2)=-2/9.                (10)
```

For a one-vertex block there is no internal edge; for a two-vertex block,
the scalar e is its one edge magnitude and Q_B(b)=e b_1 b_2. Thus size
three is the first possible failure for **one attached vertex and this
specific state**. No unrestricted minimality claim is made.

The first failed implication was that equal aggregate production and
edge mass imply equal interaction with the outside. The missing coordinate
is which internal edge receives which product b_u b_v. The repaired
operation state retains Q_B, or W_B with the matching boundary labels.
This is a direct weighted version of the exterior-response discipline
behind THM-2195, rather than a cosmetic use of tournament language.

## 5. Connection contract and amplitude boundary

**Source:** a positively weighted tournament on fixed coordinate axes.
**Target:** a trace-free gradient matrix at one instant. **Map:** place its
edge magnitude in M_ij for an arc i->j. **Preserved:** (1), and therefore
the sign of the actual skew-norm production under (2). **Lost:** spatial
gradient compatibility, pressure Hessian, future evolution, and invariance
of the positive tournament cone along that evolution. **Needed sidecar:**
those PDE data plus the coordinate gauge. The smallest decisive control
is a cyclic versus transitive triple, giving +3abc versus -abc.

**Operation transfer:** weighted tournament block -> its exterior
production response. The map (7) preserves exact production only for
block-constant contacts. Nonuniform attachment destroys the weighted
incidence arrangement; (10) gives opposite signs with identical retained
scalars. Q_B is the precise repair for the one-vertex attachment operation.

The separate amplitude question is routed to
[THM-4466 / sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes](THM-4466-sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes.md).
That companion is not a dependency of this theorem. In particular neither
(1) nor the contact law is an unrestricted amplitude bound for arbitrary
weighted tournaments; the general higher-order bound remains outside the
proved scope here.

## 6. Reproduction and finite scope

Run the standard-library verifier:

```
python3 04-computation/tournament_shear_cycles_20260908.py
python3 -O 04-computation/tournament_shear_cycles_20260908.py
```

The transcript is
[tournament_shear_cycles_20260908.out](../../05-knowledge/results/tournament_shear_cycles_20260908.out).
It checks all 33,867 labeled unit tournaments through order six, all 4,312
weighted matrix instances in its declared order-three/four alphabets,
4,160 exact substitutions, and 1,529 directed-cut pairs. Direct triangle
products are compared with independent matrix multiplication and the
three-dimensional curl contraction. There are **115,991 explicit gates**;
they remain active under -O. Exact controls verify (8), the two opposite
signs in (10), and the strong negative order-six witness. Normal and
optimized output agree with the stored LF transcript.

The semantic SHA256 is
`44cb170f8d233865187066a7744e27af0131c6ecb88e4b500b6527c4425e858a`.
The full symbolic proofs establish the all-order identities; finite
enumeration is not their source or an uncomputed classification.

Independent root review accepted the trace normalization, Euler skew-norm
factor, unit-cut equality classification, weighted [1/2,1] hostile, and the
boundary taxes 52/9 and 56/9 without a mathematical correction.
