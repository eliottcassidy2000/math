---
id: THM-3811
title: "Ramification class unit criterion and nonlinear cubic packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Units on the
  maximal etale open of
  a normal finite completion are controlled exactly by integer relations
  among the deleted ramification classes.  The nonlinear Delone--Faddeev
  packet (a,b,c,d)=(A,C,7A,C^2-3A) is a normal nonmonogenic S3 cubic with a
  squarefree irreducible rational branch: the four Veronese rays glue into
  one ramification curve, with an explicit simple companion.  Its completion
  has S*=k*, Cl(S)=Z^3, the ramification class is primitive, and its etale
  complement is the explicit affine surface
  Spec S[A/D,omega/D], with units k* and Picard group Z^2.  This surface is
  not the THM-3785 pseudo-plane.  Constructing or obstructing a polynomial
  plane atlas remains OPEN; no Jacobian counterexample is claimed.
source: jc_quartic_c3_construct / nonlinear binary-index design lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn-boundary, 2026-08-23).  The
  divisor exact sequence, norm/companion
  torsion criterion, Delone--Faddeev signs, irreducibility and normality,
  squarefree discriminant, rational ramification normalization, four-branch
  vertex gluing, simple companion, and the bounded affine norm-Pell no-go
  were checked independently inside the proof.  The exact companion has a
  direct coefficient Groebner gate for the affine exponent-one norm cell.
  The audit rederived the normal-surface divisor sequence and both directions
  of the norm criterion, checked that monicity plus integral closedness makes
  the A=1 irreducibility specialization valid, audited the codimension-one
  maximal-order argument and rational branch coverage, and replayed all 72
  exact gates in normal and optimized mode against the frozen output and raw
  hashes.  A second independent hostile audit (jc_sparse_direct_search,
  2026-08-23) rederived the complete six-prime Nagata lattice and unit
  kernel, the two A=0 primes and primitive boundary class, both different
  valuations, the affine-complement Bezout inverse, the q/r chart coverage,
  and the Picard quotient.  It replayed all 94 exact gates in normal and
  optimized mode against the frozen output and raw hashes.  No repair was
  found; the plane-atlas clause remains deliberately OPEN.
depends_on: []
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
script: 04-computation/jc2_nonlinear_binary_cubic_ramification_class_thm3811.py
output: 05-knowledge/results/jc2_nonlinear_binary_cubic_ramification_class_thm3811.out
script_sha256: 1b5f6ecf9a928222683ebd71e8d2c1221b5d6496a1b3a18ec8af53c4aafe5c11
output_sha256: 097b40607df719d548131baa584ca8ac9f8057c0a4ca58cc86a1226ca95ac239
semantic_sha256: 68cc58bc48ee278772d8021bdd7af3c1b62d045e20dda0f14eb4862c7d2a25b0
hash_basis: raw LF bytes
---

# THM-3811 -- constant units are a ramification-class problem

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
algebraically closed of characteristic zero.

## 1. The exact deleted-divisor unit criterion

Let `Xbar=Spec S` be a normal integral affine surface and let
`E_1,...,E_r` be distinct prime divisors.  Put

```text
U=Xbar minus union_i E_i,                 B=Gamma(U,O_U).           (1)
```

Removing an additional codimension-two subset does not change `B`, because
`Xbar` is normal.  There is an exact sequence

```text
1 -> S* -> B* --div_E--> direct_sum_i Z E_i
  -> Cl(S) -> Cl(U) -> 0.                                        (2)
```

The middle map sends a unit on `U` to its Weil divisor on `Xbar`, which is
supported on the deleted primes.  Therefore, if `S*=k*`,

```text
B*=k*
iff the classes [E_1],...,[E_r] are Z-linearly independent in Cl(S). (3)
```

This is the exact obstruction behind THM-3808.  Its four ramification rays
lie in `Cl(S)=Z/3`, so relations are unavoidable.  A minimal escape glues
the four local rays into one global prime `E`; then `(3)` becomes

```text
B*=k* iff [E] has infinite order,             provided S*=k*.      (4)
```

Affineness of `U` is a separate requirement.  Sequence `(2)` computes its
global units whether or not `U` is affine; a counterexample factor surface
needs both `(3)` and affineness.

There is a useful norm form of the same gate.  Suppose `S` is finite cubic
over `R=k[A,C]`, one irreducible squarefree branch is `Gamma=V(delta)`, and
its generic pullback is

```text
pi^*Gamma=2E+F,                                                     (5)
```

with `E` ramified and `F` the simple companion.  Then `[E]` has finite
order if and only if there are `n>=1` and `u in S minus {0}` such that

```text
div_Xbar(u)=nE.                                                     (6)
```

Equivalently, `(6)` is the simultaneous norm/companion packet

```text
Norm_(K/k(A,C))(u)=kappa delta^n,        kappa in k*,
v_F(u)=0.                                                       (7)
```

Indeed a principal relation `n[E]=0` has an effective divisor `nE`, so its
principalizer has no poles and belongs to the normal ring `S`.  Pushing its
divisor to `A2` gives `(7)`.  Conversely, the norm equation together with
the companion valuation says that every height-one zero of `u` is on `E`,
and the total valuation is `n`.  Thus `(7)`, not a discriminant computation
alone, is the cheapest decisive constant-unit test.

## 2. The first nonlinear packet

Let `R=k[A,C]` and give the free module

```text
S=R*1 direct_sum R*omega direct_sum R*theta                         (8)
```

the Delone--Faddeev multiplication for

```text
(a,b,c,d)=(A,C,7A,C^2-3A):                                        (9)

omega^2=-7A^2+C omega-A theta,
omega theta=3A^2-AC^2,
theta^2=3AC-C^3+(C^2-3A)omega-7A theta.                           (10)
```

Its binary index form is

```text
I(x,y)=-(A x^3+C x^2y+7Axy^2+(C^2-3A)y^3).                       (11)
```

Thus every index value lies in `(A,C)` and no global generator exists.  The
only common zero of the four coefficients is the origin, where

```text
S/(A,C)S=k[omega,theta]/(omega,theta)^2.                            (12)
```

The trace discriminant is

```text
Delta = A(C+5A)(4C+19A)(3C-17A)
      + C^2(162A^3+126A^2C-4C^3)-27A^2C^4.                        (13)
```

Exact differentiation gives

```text
gcd(Delta,Delta_A,Delta_C)=1,                                     (14)
```

so the affine branch is reduced.  Unlike THM-3808, `(13)` is not a union of
radial lines.

The characteristic polynomial of `omega` is

```text
f(T)=T^3-CT^2+7A^2T+3A^3-A^2C^2.                                 (15)
```

It is irreducible over `k(A,C)`.  Indeed, reducibility of a cubic would give
a root in `k(A,C)`.  Since the polynomial is monic, that root is integral
over the integrally closed UFD `k[A,C]`, hence belongs to `k[A,C]`; setting
`A=1` would then give a root `v(C) in k[C]` of

```text
v^3-Cv^2+7v+3-C^2=0                                               (16)
```

cannot have a finite pole: its cubic term would have strictly largest pole
order.  Hence `v` is a polynomial.  Degree comparison leaves only
`v=C+beta`; substitution gives

```text
(beta-1)C^2+(2beta^2+7)C+(beta^3+7beta+3),                         (17)
```

which cannot vanish because `beta=1` makes the linear coefficient `9`.
Thus `(15)` is irreducible.  The free ring `(8)` is a domain, and the same
squarefree-discriminant DVR index argument as in THM-3808 gives `R1`; it is
free over `R`, hence `S2`.  Therefore `S` is normal.  Its nonsquare reduced
discriminant and simple inertia give generic Galois group `S3`.

## 3. Four rays glue into one rational ramification curve

On `A!=0`, put

```text
q=omega/A,                 G(q)=q^3+7q+3.                          (18)
```

Equation `(15)` becomes the rational chart

```text
A G(q)=C(C+q^2).                                                   (19)
```

At a double root, `f=f_T=0`.  Solving these two equations gives

```text
R(q)=(q-3)(q+1)(q+2)=q^3-7q-6,

A(q)=-2q^2R(q)/(3q^2+7)^2,
C(q)=-qR(q)/(3q^2+7),
omega(q)=A(q)q,
theta(q)=A(q)(q^2-7)/2.                                          (20)
```

Direct substitution verifies `(10)`, `(15)`, `f_T=0`, and
`Delta(A(q),C(q))=0`.  Conversely, at every point of `V(Delta)` with
`A!=0` and a nontriple cubic, the unique double root recovers `q`, and
`(20)` follows.  No branch component lies in `A=0`, because

```text
Delta(0,C)=-4C^5.                                                  (21)
```

Hence the closure of `(20)` is the whole branch curve.  It is irreducible
and rational; squarefreeness in `(14)` then makes `Delta` its irreducible
equation.

Most importantly,

```text
q=-2,-1,0,3                                                     (22)
```

are four distinct normalization points mapping to

```text
(A,C,omega,theta)=(0,0,0,0).                                     (23)
```

The degree-four initial form of `(13)` is exactly the four-line
discriminant of THM-3808.  Thus the nonlinear term `C^2` performs the
minimal geometric repair: it glues the four radial ramification rays into
one global curve while retaining them as four local branches at the
square-zero vertex.

The simple companion is equally explicit.  On the branch, the roots of
`(15)` are

```text
T=Aq  (double),                 T_comp=C-2Aq  (simple).             (24)
```

The corresponding `theta` is recovered from `(10)`.  Therefore the generic
pullback is exactly `(5)`, with one rational ramification curve and one
residue-degree-one companion.  The packet passes the nonmonogenic,
squarefree, `S3`, and companion gates of THM-3801 without the Veronese
arrangement symmetry.

## 4. First exact norm-Pell obstruction

Write an arbitrary element uniquely as

```text
u=p+q_1 omega+q_2 theta.                                           (25)
```

For the universal Delone--Faddeev law, its norm is

```text
a^2d q_1^3+2abd q_1^2q_2-ac^2 q_1^2q_2-2acd q_1q_2^2
+acp q_1^2-ad^2q_2^3+3adp q_1q_2+b^2d q_1q_2^2
-bcp q_1q_2+bdp q_2^2+bp^2q_1-cp^2q_2+p^3.                       (26)
```

The exact companion substitutes `(9)`, lets each of `p,q_1,q_2` be an
arbitrary affine polynomial in `(A,C)`, and equates `(26)` to
`kappa Delta` with `kappa!=0`.  The resulting 38 coefficient equations in
ten coefficient variables, plus the inverse for `kappa`, have reduced
Groebner basis

```text
[1].                                                               (27)
```

Thus no exponent-one principalizer in `(7)` has affine coordinate profiles.
This is a decomposability obstruction, not merely a linear span test.

There is also an all-exponent degree floor.  If
`deg p,deg q_1,deg q_2<=D`, then `(26)` has total degree at most

```text
3D+5,                                                              (28)
```

whereas `Delta^n` has degree `6n`.  Therefore an order-`n` torsion witness
must have

```text
D>=2n-1;                                                           (29)
```

for `n>=2`, while `(27)` strengthens the `n=1` floor to `D>=2`.

## 5. The class gate closes exactly

Put `Xbar=Spec S` and set

```text
D=f'(omega)=C omega-3A theta-14A^2.                              (30)
```

The apparent four-ray unit problem can be computed without guessing a
principalizer.  Introduce the root chart

```text
T=k[A,C,q]/(A G(q)-C(C+q^2)),       G(q)=q^3+7q+3,                (31)
omega=Aq,                 theta=Cq-A(q^2+7).                     (32)
```

The map `S -> T` is injective and becomes an isomorphism after inverting
`A`.  Indeed `q=omega/A` there, and `(31)` is exactly the characteristic
cubic `(15)`.  The surface `T` is smooth.  For its equation `F`,

```text
F_A=G,                  F_C=-2C-q^2,
F_q=A(3q^2+7)-2Cq.                                           (33)
```

If the first two derivatives and `F` vanished, then `C=-q^2/2` and
`F=q^4/4`; hence `q=0`, contradicting `G(0)=3`.

Let `alpha_1,alpha_2,alpha_3` be the three distinct roots of `G`.
After inverting `A` and `G`, `(31)` becomes the UFD

```text
k[C,q,G^-1,C^-1,(C+q^2)^-1].                                  (34)
```

The six height-one primes deleted by inverting `G` are

```text
P_j=(q-alpha_j,C),             Q_j=(q-alpha_j,C+q^2).           (35)
```

Their complete Nagata relations are

```text
P_j+Q_j=div(q-alpha_j),        sum_j P_j=div(C),
sum_j Q_j=div(C+q^2).                                           (36)
```

The last relation follows from the first four.  A `4 x 4` minor of the
relation matrix is a unit, so there is no hidden torsion and

```text
Cl(S_A)=Z^2.                                                     (37)
```

The same six valuations compute the units.  A unit of `(34)` has the form

```text
c product_j(q-alpha_j)^{n_j} C^u(C+q^2)^v.
```

Zero valuation at all `P_j,Q_j` forces

```text
n_1=n_2=n_3=-u=-v,
```

so it is a scalar times a power of
`C(C+q^2)/G=A`.  Consequently

```text
S_A^*=k^* A^Z.                                                   (38)
```

There are exactly two height-one primes over `A=0`:

```text
P_0=(A,omega),
P_1=(A,omega-C,theta).                                          (39)
```

Their residue fields have degrees two and one over `k(C)`, respectively.
Because `Delta(0,C)=-4C^5`, the cover is etale at both generic points, and

```text
div_Xbar(A)=P_0+P_1.                                             (40)
```

The localization sequence for `S -> S_A`, together with `(37),(38)`, now
has boundary map `A |-> (1,1)`.  It follows first that `S^*=k^*`, and then
that

```text
0 -> Z[P_1] -> Cl(S) -> Cl(S_A) -> 0,
Cl(S)=Z^3,                    [P_0]=-[P_1],                       (41)
```

where `[P_1]` is a primitive infinite-order class.  The displayed splitting
of the free abelian groups is noncanonical; the primitivity assertion is
canonical.

The ramification class is exactly this missing companion class.  On `(31)`,

```text
D=A J,                    J=A(3q^2+7)-2Cq.                       (42)
```

At the generic point of `P_0`, write `C=-q^2`; then `J=2q^3` is a unit.
At `P_1`, `(30)` restricts to `C^2`, also a unit.  Along the simple
ramification curve, `D` has order one.  The norm identity

```text
Norm(D)=-A^2 Delta                                                (43)
```

shows that there is no further codimension-one support.  Therefore

```text
div_Xbar(D)=E+P_0,
div_Xbar(D/A)=E-P_1,
[E]=[P_1].                                                       (44)
```

Thus the class in `(4)` is primitive, not merely nontorsion.  In particular
the norm/companion equation `(7)` has no solution in any degree; the bounded
Groebner gate `(27)` was the first finite shadow of this all-degree fact.

## 6. The etale complement is affine

Let

```text
U=Xbar minus E,             h=A/D,             k_0=omega/D.       (45)
```

Equation `(44)` shows that `h` is regular on `U`.  The same root chart gives
`v_{P_0}(omega)=v_{P_0}(D)=1`; at `P_1` both `omega` and `D` are units.
Hence `k_0` is also regular on `U`.  Put

```text
B=S[h,k_0] subset Frac(S).                                       (46)
```

Besides `hD=A` and `k_0D=omega`, these two functions satisfy the crucial
Bezout identity

```text
C k_0-(3theta+14A)h=1.                                          (47)
```

Let `Z=Spec B -> Xbar` be induced by `S subset B`.  It has no point over
`E`: if `D=0`, the first two graph relations force `A=omega=0`.  On `E`,
`Delta(0,C)=-4C^5` then forces `C=0`; the square-zero fibre `(12)` has
`theta=0` at its only geometric point, contradicting `(47)`.  Thus the image
of `Z` lies in `U`.

Conversely, `(45)` defines a morphism `U -> Z`.  Its composite with
`Z -> U` is the identity on the `S`-coordinates.  The reverse composite
fixes `S,h,k_0` on the dense open `D!=0`, hence fixes the domain `B`
identically.  Therefore

```text
U = Spec S[A/D,omega/D]                                         (48)
```

is affine.  Combining `(2),(41),(44)` gives the exact invariants

```text
U^*=k^*,                       Pic(U)=Cl(U)=Z^2.                  (49)
```

Here `U` is regular: the two root charts below cover it by loci on which the
projection to `A2` has a unit Jacobian.

## 7. Why the one-root chart lied, and the remaining source gate

On `T`, let `L=(A,C)`.  Directly from `(33)` and `(20)`,

```text
div_T(J)=L+E_tilde.                                               (50)
```

This does **not** identify `U` with `D_T(J)`.  The chart `(31)` contracts
`L` to the square-zero vertex and entirely misses the companion divisor
`P_1`.  Downstairs, `(44)` says that the same rational function has a pole
on `P_1`; this is the missing valuation that resolves the apparent
nonconstant-unit paradox.

The second root chart is obtained from `r=(omega-C)/A`:

```text
T_1=k[A,C,r]/(
 C^2r+A(2Cr^2-C^2+7C)+A^2(r^3+7r+3)),                           (51)
omega=C+Ar,                 theta=-Cr-A(r^2+7).                  (52)
```

Its projection Jacobian is

```text
C^2+4ACr+A^2(3r^2+7)=D,                                        (53)
```

which restricts to `C^2` on `P_1`.  Thus

```text
Spec T[J^-1]       and       Spec T_1[D^-1]                      (54)
```

are the ramified-root and companion-root affine charts covering `U`; their
intersection is the `A!=0` root chart.  This is the required companion
sidecar that the single `q`-coordinate destroyed.  More precisely, `q` has
its only missing codimension-one valuation at `P_1`, while `r` has its only
missing codimension-one valuation at `P_0`; inverting the displayed
Jacobians removes their common exceptional vertex and ramification curve.

The map `U -> A2_(A,C)` is etale of generic degree three.  Its image omits
the square-zero origin and the two nonzero triple-root values.  On the branch
normalization the latter are exactly

```text
q^2=7/3,                   C=3Aq;                                (55)
```

the four values `q=-2,-1,0,3` all map to the origin.  Every other branch
point retains the simple companion `(24)`.

Equation `(49)` rules out both `U~=A2` and
`U~=Y_(THM-3785)`, since the latter has Picard group `Z/3`.  It does not rule
out a nontrivial dominant etale morphism

```text
psi:A2_(x,y) -> U.                                                (56)
```

Indeed `psi^*:Pic(U)->Pic(A2)=0` may kill the whole free lattice.  If a
polynomial map `(56)` with the needed codimension-one coverage is found,
then `(A o psi,C o psi)` is a planar Keller map whose finite completion is
the cubic algebra `(8)`.  Conversely, an obstruction to every such atlas
would close this construction lane without changing the positive affine
surface theorem.

```text
OPEN-3811:
  construct a polynomial dominant etale plane atlas (56), with the
  codimension-one coverage needed by the completion, or prove that none
  exists.                                                         (57)
```

Thus the class, constant-unit, and affineness gates are closed.  The only
remaining counterexample gate is the source-plane realization.  **QED for
`(2)--(55)`; `(56),(57)` remain OPEN.**
