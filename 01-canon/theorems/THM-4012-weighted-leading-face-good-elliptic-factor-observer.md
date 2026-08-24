---
id: THM-4012
title: "Weighted leading-face good elliptic-factor observer"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY GEOMETRIC-AUDITED. In THM-3992's
  normalized reduced 2:3 cell,
  the highest wt(p,y)=(2,3) face gives, after one common base change, a
  rational component and the face curve H_M(P,SP)=1, while the target has
  good j=0 reduction. Under the explicitly named face-stable specialization
  hypothesis, a finite source-to-target fibre map forces the Jacobians of
  the normalized positive-genus face components to have a nonzero Hom to
  the target elliptic curve. Singleton faces fail this observer. The
  nonresonant weight-six p^3+y^2 face is exactly a j=0 Hom survivor, but on
  the forced b=d=0 max-weight-six branch its attachment point is provably
  nontorsion (orders 6 and 9 after good reduction at 11 and 17), contradicting
  the six-attachment isogeny invoice. An exact weighted-boundary
  normalization proves face-stability for this max-weight-six polynomial, so
  that live branch is excluded without an extra stable-model hypothesis and
  the actual weighted degree is unconditionally at least seven. Every
  weight-eight p^4+p*y^2 face with both coefficients nonzero is a Bolza
  genus-two curve whose Jacobian is isogenous to E_8000^2; an explicit CM
  mismatch gives Hom(E_8000,E_0)=0, so it fails the observer. This is a
  conditional stable-specialization obstruction outside the exact max-six
  model proved here, not a proof that every weighted closure has the required
  stable model and not a JC(2) result.
source: root + generic_fibre_residual / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (primary exact certificate + independent SymPy-free arithmetic,
  orbit, boundary-normalization, and geometric audit, 2026-08-24). The exact
  SymPy
  certificate checks the common source/target scaling for M=2..12, all 36
  singleton monomials of weights 2..18, the live weight-six coefficient,
  boundary normalization, attachment, two-prime nontorsion certificate, and
  a self-contained Bolza quotient,
  j-invariant, degree-two CM endomorphism, and Hom mismatch. Normal and
  optimized streams byte-match across 206 primary and 58 independent gates.
  The independent audit reconstructs the max-six proper/flat closure, finite
  simultaneous normalization, unique elliptic component, six A_35 rational
  chains, and two-prime torsion separator. The general Hom theorem remains
  deliberately conditional on face-stable specialization; only the exact
  nonresonant max-six family is proved face-stable here.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
related:
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-3999-companion-divisor-boundary-endpoint-and-class-ledger
  - THM-4016-sharp-five-by-five-elliptic-attachment-nontorsion
script: 04-computation/jc2_weighted_face_good_factor_thm4012.py
output: 05-knowledge/results/jc2_weighted_face_good_factor_thm4012.out
script_sha256: 2728ce7f8319e8e2584ff3db75e21403c5c5aeae968c9febd6a5f6b763331de8
output_sha256: 8567b9c68367df5a3581149786405080c3a15c46f5d3fd1eb62710f00f58753b
independent_audit_script: 04-computation/jc2_weighted_face_good_factor_thm4012_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_weighted_face_good_factor_thm4012_independent_audit.out
independent_audit_script_sha256: 06e0f2182a17568a546cccaa59b1e32b47eb298503d4496fc6985770a3fa100e
independent_audit_output_sha256: 7acc2638a6269255cbf27cea86bfe2f73a18434a539e82db50cd857c6a4109ea
independent_audit_semantic_sha256: 0d100f08e59226afa995266af65d35f74667bed5d698ea01cedbeec7c21d0bfc
hash_basis: raw LF bytes
---

# THM-4012 -- weighted leading faces must carry the good elliptic factor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY GEOMETRIC-AUDITED.**

This theorem isolates a genuine necessary observer and three exact
low-weight outcomes. Its central implication is conditional on the
`face-stable` hypothesis in Section 2. That condition is not decorative: a
rational-looking singular special fibre can acquire a positive-genus tail
after base change and stable reduction. Section 4.3 proves the gate separately
for the exact nonresonant max-weight-six polynomial; the general conditional
observer and this exact boundary calculation were independently audited.

## 1. Exact weighted degeneration

Work over an algebraically closed field `k` of characteristic zero. Use the
normalized reduced `(2,3)` cell and notation of THM-3992:

```text
s=xt,                  p=s^2+t,                  y=sp,
u=s^2/t,

G=gamma*u+H(p,y),      H=lambda*p+R(p,y),
a*gamma!=0,            R in (p^2,y).                         (1)
```

The word **leading** below always means the highest weighted face of the
**entire polynomial** `H`, with

```text
wt(p)=2,               wt(y)=3,
M=max{2i+3j:[p^i y^j]H!=0},
H_M(p,y)=sum_(2i+3j=M) c_ij p^i y^j.                     (2)
```

It does not mean the first currently displayed residual row. Since
`lambda!=0`, one has `M>=2`.

Use the target value as parameter, `q=G`. On the source generic fibre,
`t=p-s^2` turns `(1)` into

```text
s^2(q+gamma-H(p,sp))=p(q-H(p,sp)).                       (3)
```

The target generic fibre is

```text
E_q: C^2=A^3-(3a^2/4)A+q-a^3/4.                         (4)
```

If `(A,C)` is a Keller pair, the finite function-field inclusion restricts
to a finite nonconstant morphism between smooth projective generic fibres

```text
phi:C_q -> E_q.                                          (5)
```

Make the common ramified base change

```text
q=rho^(-6M),       s=rho^-6 S,       p=rho^-12 P,
y=rho^-18 SP.                                             (6)
```

Put

```text
H_rho(S,P)=rho^(6M)H(rho^-12 P,rho^-18 SP).
```

Every monomial of weight `w` occurs with coefficient `rho^(6(M-w))`, so
`H_rho` is integral and `H_rho mod rho=H_M(P,SP)`. Multiplying `(3)` by
`rho^(6M+12)` gives the exact family

```text
(S^2-P)(1-H_rho)+gamma*rho^(6M)S^2=0.                   (7)
```

Hence its affine central equation is

```text
(S^2-P)(1-H_M(P,SP))=0.                                 (8)
```

There is also an exact proper weighted closure. In
`P_O(1,2,1)_[S:P:Z]`, where `O=k[[rho]]`, define

```text
Hbar_rho=
 sum c_ij rho^(6(M-2i-3j)) S^j P^(i+j) Z^(M-2i-3j).    (9)
```

Every term has weighted degree `M`, and `(7)` closes to the degree-`M+2`
divisor

```text
(S^2-P)(Z^M-Hbar_rho)+gamma*rho^(6M)S^2 Z^M=0.         (10)
```

Its displayed central support is

```text
C_rat: P=S^2,
D_M:   Z^M=H_M(P,SP).                                   (11)
```

The normalization of `C_rat` is `P1`. The second curve is the weighted face
observer.

On the same base change, put

```text
A=rho^(-2M)X,                  C=rho^(-3M)Y.             (12)
```

Equation `(4)` becomes

```text
Y^2=X^3+1-(3a^2/4)rho^(4M)X-(a^3/4)rho^(6M),           (13)
```

so the target has actual good reduction

```text
E_0:Y^2=X^3+1,                  j(E_0)=0.                (14)
```

## 2. The face-stable hypothesis and the Hom observer

The weighted equation alone is not declared to be the stable model. Say that
`H` is **face-stable at weight M** if, after a finite extension of
`k((rho))`, there is a proper regular semistable model dominating the
normalization of `(10)` with the following positive-genus inventory:

1. for each positive-genus normalized irreducible component of `D_M`, its
   strict transform occurs exactly once, and these strict transforms are
   **precisely** the positive-genus irreducible components of the geometric
   special fibre;
2. every remaining special-fibre component, including all components introduced
   at the weighted boundary and in resolving the total space, is rational;
3. the special fibre is connected, with its actual multiplicities retained.

A sufficient local sidecar is that all remaining total-space and boundary
singularities become toroidal after the chosen base change: toroidal
resolution inserts only rational chains. This sufficient sidecar is not
claimed automatically for an arbitrary collection of lower faces.

At every displayed **finite transverse attachment**, however, the needed
local sidecar is automatic. Such a point has `S!=0`. Put

```text
U=S^2-P,                 V=(1-H_rho)/S^2.                (15a)
```

Transversality says that `U,V,rho` are regular local parameters in the
ambient threefold. Dividing `(7)` by the unit `S^2` gives the exact completed
local equation

```text
U*V=-gamma*rho^(6M).                                    (15b)
```

Thus the attachment has toroidal thickness `6M`. Its regular resolution
inserts a rational chain joining the strict transforms of the two displayed
components. In particular, the six nonresonant weight-six attachments used
below really do transmit equality of specialized map values; no global
address census is being inferred.

Write `D_nu_tilde` for the positive-genus normalized face components and
define

```text
B_M=product_nu Jac(D_nu_tilde).                          (15)
```

The empty product is the zero abelian variety.

### The observer

**THEOREM (under the displayed face-stable hypothesis).** If `(5)` exists
and `H` is face-stable at `M`, then

```text
Hom_k(B_M,E_0) != 0.                                    (16)
```

Equivalently,

```text
Hom_k(B_M,E_0)=0                                        (17)
```

is a no-go certificate for that face-stable lane.

### Direct curve-model proof

Base-change `(5)` to the DVR in the face-stable definition. The target
family `(13)` is a smooth proper elliptic scheme. Resolve the rational map
from the proper regular source model to this target by point blowups in the
special fibre. Every new exceptional component is rational.

If `(16)` failed, every normalized positive-genus face component would map
constantly to `E_0`: a nonconstant curve map to `E_0` induces a nonzero map
on Jacobians. Every rational component also maps constantly by
Riemann--Hurwitz. Constants agree along the connected special fibre, so the
whole special map is constant.

Let `L` be relatively ample on the elliptic scheme and write the full special
fibre as `sum e_i C_i`, retaining all positive multiplicities. Constancy of
relative intersection degree gives

```text
deg(phi_eta^*L)
 =sum_i e_i deg((phi|C_i)^*(L|E_0))=0.                  (18)
```

But `(5)` is finite and nonconstant, so

```text
deg(phi_eta^*L)=deg(phi)*deg(L|E_q)>0,                  (19)
```

a contradiction. This proof does not assume exactness of a Neron-model
quotient functor. It works over every algebraically closed characteristic-zero
`k` and is insensitive to fibre multiplicities created by later blowups.

The observer is only necessary. A nonzero Hom in `(16)` does not construct a
curve map, satisfy attachment compatibility, or produce a Keller pair.

## 3. Singleton top faces are closed in the face-stable universe

Suppose

```text
H_M(p,y)=c p^i y^j,             c!=0,
M=2i+3j.                                                (20)
```

On `Z=1`, the face curve is

```text
c S^j P^(i+j)=1.                                        (21)
```

Put `g=gcd(j,i+j)`. Over the algebraically closed field, `(21)` has exactly
`g` irreducible components. Each is a translate of a one-dimensional torus
when both exponents occur; in the endpoint cases it is an affine line or the
same torus. Every projective normalization is therefore `P1`.

The intersections with `C_rat` are completely nonresonant. Substitution
`P=S^2` gives

```text
c S^M=1,                                                (22)
```

so there are exactly `M` simple attachment points. For

```text
F_0=S^2-P,             F_1=1-cS^jP^(i+j),
```

the gradient determinant on `P=S^2` is

```text
det(dF_0,dF_1)=-M c S^(M-1) !=0.                        (23)
```

Thus the face abelian variety `(15)` is zero. By `(17)`, every singleton top
face is impossible whenever its weighted degeneration is face-stable. This
includes pure `p`, pure `y`, and mixed singleton monomials. THM-4008 remains
stronger on the exact pure-`p` lane because there the full semistable model
and its rational resolution chains were proved rather than assumed.

## 4. Weight six is the exact survivor boundary

The weight-six face has only two possible monomials:

```text
H_6(p,y)=epsilon*p^3+kappa*y^2.                         (24)
```

If either coefficient vanishes, Section 3 applies. Assume

```text
epsilon*kappa!=0.                                       (25)
```

On `D_6`, put `T=SP`. Since `P` cannot vanish there, this is birational with

```text
E_6: kappa*T^2=1-epsilon*P^3.                           (26)
```

Its smooth projective normalization is elliptic. Choose `u,v in k*` with
`epsilon=-u^3` and `kappa=v^2`; then `P=X/u,T=Y/v` sends `(26)` to `(14)`.
Consequently

```text
Jac(D_6_tilde)=E_6 isomorphic to E_0,
Hom(Jac(D_6_tilde),E_0)!=0.                             (27)
```

The observer is therefore satisfied exactly, not violated.

The finite attachments are nonresonant precisely when

```text
epsilon+kappa!=0.                                       (28)
```

Indeed, on `C_rat` they obey

```text
(epsilon+kappa)S^6=1,
det(dF_0,dF_6)=-6(epsilon+kappa)S^5,                    (29)
```

giving six transverse points. At `epsilon+kappa=0`, the finite attachments
move to the weighted boundary, and face-stability must be re-audited there.
Thus `(27)` is the exact boundary where the Hom observer by itself stops.

### 4.1 The six-attachment torsion invoice

There is a stronger condition once face-stability and the six transverse
attachments are both used. In the weight-six face, `E_6` is the only
positive-genus special component. Therefore the degree argument of Section 2
forces the specialized restriction

```text
psi:E_6 -> E_0
```

to be nonconstant. After a target translation, `psi` is an isogeny. The
rational component maps constantly. Each of the six exact local chains
`(15b)` also maps constantly, so all six attachment points on `E_6` have one
common image under `psi`.

They form one orbit under

```text
sigma(P,T)=(zeta_3 P,-T).                               (29a)
```

This is an order-six elliptic automorphism fixing infinity. In the
Eisenstein endomorphism ring, `sigma=-zeta_3`, hence

```text
sigma-1=-(zeta_3+1)=zeta_3^2                            (29b)
```

is a unit. If `P_0` is one attachment, equality of the images of `P_0` and
`sigma(P_0)` gives

```text
psi((sigma-1)P_0)=0.                                    (29c)
```

The isogeny kernel is finite and `(29b)` is invertible, so

```text
P_0 is torsion on E_6.                                  (29d)
```

This conclusion is **conditional on face-stability**, because without the
positive-genus inventory a hidden component could carry the nonconstant
specialized map while `E_6` mapped constantly. The transmission through the
six displayed finite nodes is not an extra assumption: it is the exact
toroidal local calculation `(15a)--(15b)`.

### 4.2 The forced live attachment is nontorsion

For the exact `b=d=0`, max-weight-six branch determined below, it is cleaner
to use normalized residual coefficients

```text
epsilon_tilde=2752/(135A5^3),
e_tilde=-5696/(45A5^3).                                 (29e)
```

Multiplying both by `gamma` gives the raw coefficients in `(24)` and does not
change the following ratios. Under the isomorphism `(26)->(14)`, an
attachment point has

```text
X^3=-epsilon_tilde/(epsilon_tilde+e_tilde)=43/224,
Y^2= e_tilde/(epsilon_tilde+e_tilde)=267/224.            (29f)
```

Let

```text
K=Q(alpha,beta),       alpha^3=43/224,       beta^2=267/224,
P_0=(alpha,beta) in E_0(K).                              (29g)
```

The cubic is irreducible by the rational-cube test, and `267/224` is not a
rational square. A cubic field has no quadratic subfield, so the two
extensions are linearly disjoint and `(29g)` presents a degree-six
compositum; there is no hidden mixed relation between `alpha` and `beta`.

The denominators are units at `11` and `17`. The assignments

```text
mod 11: alpha->2, beta->3,       P_0->(2,3),
mod 17: alpha->7, beta->2,       P_0->(7,2)              (29h)
```

respect the defining equations and therefore select primes of `K` above
those rational primes. The displayed roots are simple modulo both primes, so
these are residue-degree-one primes of the localized order. The curve
`Y^2=X^3+1` has good reduction at both.
Exact addition gives

```text
over F_11: 2P=(0,1), 3P=(10,0), 6P=O,   so ord(P)=6;
over F_17: 3P=(0,16),              9P=O, so ord(P)=9.    (29i)
```

If `P_0` had torsion order `N`, good-reduction injectivity on prime-to-residue
characteristic torsion would give

```text
N/6 is a power of 11,             N/9 is a power of 17. (29j)
```

But `6*11^r=9*17^s` is impossible already at the `2`-adic valuation.
Therefore

```text
P_0 is nontorsion.                                      (29k)
```

Equations `(29d)` and `(29k)` contradict each other. Thus the forced exact
max-weight-six branch is excluded **whenever its degeneration is
face-stable**. This closes the live branch beyond the Hom observer, but it
does not by itself remove the global face-stability hypothesis. The next
subsection supplies that missing model for this exact face.

### 4.3 Exact boundary normalization proves face-stability at max weight six

Assume now that the **entire** polynomial has maximum weight six. Write it,
with unrelated coefficient names, as

```text
H=lambda*p+c*y+b*p^2+h*p*y+epsilon*p^3+kappa*y^2,       (29l)
```

where `epsilon*kappa*(epsilon+kappa)!=0`. Use the proper weighted closure
`(10)`. Its equation is primitive modulo `rho`, so the hypersurface is flat
over `O`. It avoids the only orbifold point `[S:P:Z]=[0:1:0]`, because its
equation there has value `epsilon!=0`.

Every remaining boundary point lies in the smooth chart `S=1`. Put

```text
r=P/S^2,                         z=Z/S.
```

The exact family there is

```text
F=(1-r)(z^6-[epsilon*r^3+kappa*r^2
      +h*rho^6*r^2*z+b*rho^12*r^2*z^2
      +c*rho^18*r*z^3+lambda*rho^24*r*z^4])
      +gamma*rho^36*z^6.                                 (29m)
```

On `rho=0`, the boundary points away from `r=z=0` are `r=1` and
`r=-kappa/epsilon`. The `r`-derivatives of `F` at these points are

```text
epsilon+kappa,
-kappa^2(epsilon+kappa)/epsilon^2,                       (29n)
```

respectively, so both are smooth off the named resonance.

The remaining point `r=z=0` is the persistent two-branch `A_5` boundary
singularity of the projective closure, not a vertical degeneration. Put
`r=z^3w` and divide the strict transform by `z^6`. One obtains

```text
Phi=(1-z^3w)(1-kappa*w^2-c*rho^18*w
      -z(h*rho^6*w^2+lambda*rho^24*w)
      -b*rho^12*z^2*w^2-epsilon*z^3*w^3)
      +gamma*rho^36.                                     (29o)
```

At `z=0`, this is the quadratic

```text
1+gamma*rho^36-c*rho^18*w-kappa*w^2,                    (29p)
```

whose discriminant is

```text
c^2*rho^36+4kappa(1+gamma*rho^36)=4kappa mod rho.        (29q)
```

It is a unit. Hensel factorization therefore gives two distinct smooth
horizontal branches through the completed normalization. Equivalently, the
weighted blow-up chart `r=z^3w` separates the two persistent branches, and
their completed local rings are regular `O[[z]]` rings. Since the total
space is excellent, its normalization is finite; because `(10)` is proper,
this simultaneous normalization is proper as well. It introduces no vertical
component and no positive-genus tail.

On `Z!=0`, the normalized central support is the smooth rational component,
the smooth elliptic curve `(26)`, and their six transverse intersections.
Sections `(29n)--(29q)` cover all of `Z=0`, while the orbifold point is
absent. Thus the only remaining total-space singularities are the six exact
local models

```text
U*V=-gamma*rho^36                                      (29r)
```

from `(15b)`. Their regular resolutions insert rational chains. The result
is a proper regular semistable model whose positive-genus inventory consists
of exactly one strict transform, namely `E_6`. Hence the exact nonresonant
max-weight-six family `(29l)` is face-stable.

Combining this model with `(29d)--(29k)` removes the condition from the live
max-weight-six branch:

```text
b=d=0 and total wt_(2,3)(H)=6 is impossible.             (29s)
```

This unconditional statement is still restricted to the **exact** maximum
six polynomial. It does not say that the currently known rows prevent terms
of weight seven or higher.

## 5. Weight eight is a Bolza no-go

The weight-eight face is

```text
H_8(p,y)=c*p^4+d*p*y^2.                                 (30)
```

If `c=0` or `d=0`, this is a singleton and Section 3 closes it. Assume
`cd!=0`. The face equation is

```text
cP^4+dS^2P^3=1.                                         (31)
```

Put `T=SP` and then `Y=PT`. Since `P!=0` on `(31)`, these are birational
changes, and `(31)` becomes

```text
dY^2=P(1-cP^4).                                         (32)
```

Choose `u,v in k*` with `c=-u^4` and `d=v^2/u`, and set
`P=x/u,Y=y/v`. The smooth projective normalization is the genus-two Bolza
curve

```text
C_B: y^2=x(x^4+1).                                      (33)
```

The quintic is squarefree in characteristic zero. Its finite attachments are

```text
(c+d)S^8=1,
det(dF_0,dF_8)=-8(c+d)S^7.                              (34)
```

Thus `c+d!=0` gives eight transverse points. When `c+d=0`, the same Bolza
normalization and Hom computation below remain valid, but the attachments
collide at the weighted boundary; face-stability is then an additional,
unproved boundary assertion.

### 5.1 Two explicit elliptic quotients

On `(33)`, consider

```text
iota_1(x,y)=(1/x,y/x^3).                                (35)
```

The invariants

```text
u=x+1/x,                 v=y(x+1)/x^2                   (36)
```

give the degree-two quotient

```text
E:v^2=(u+2)(u^2-2)=u^3+2u^2-2u-4.                      (37)
```

For this integral Weierstrass equation,

```text
c_4=160,                 Delta=512,              j=8000. (38)
```

A conjugate involution is

```text
iota_2(x,y)=(-1/x, i*y/x^3).                            (39)
```

Indeed, conjugate `(35)` by `(x,y)->(-i*x,b*y)` with `b^2=-i`.
For the basis

```text
omega_0=dx/y,                    omega_1=x dx/y
```

of `H^0(C_B,Omega^1)`, one computes

```text
iota_1^*(omega_0)=-omega_1,     iota_1^*(omega_1)=-omega_0,
iota_2^*(omega_0)=-i*omega_1,   iota_2^*(omega_1)= i*omega_0.   (40)
```

The two invariant quotient differentials are

```text
omega_0-omega_1,                 omega_0-i*omega_1,      (41)
```

and are linearly independent. Therefore the product of the two quotient
pushforwards has full-rank differential and is an isogeny

```text
Jac(C_B)  ~  E^2.                                        (42)
```

### 5.2 Self-contained CM mismatch

Shift `x=u+2` in `(37)` to write

```text
E:y^2=x^3-4x^2+2x.                                      (43)
```

Quotienting by `(0,0)` gives the explicit degree-two isogeny

```text
X=x-4+2/x,
Y=y(2-x^2)/x^2,
E':Y^2=X^3+8X^2+8X.                                    (44)
```

The change `x=-X/2,y=tau*Y`, where `tau^2=-1/8`, identifies `E'` with
`E`. If `alpha` is the resulting degree-two endomorphism of `E`, direct
substitution into the duplication law gives

```text
alpha^2=[-2].                                            (45)
```

Since an elliptic endomorphism of degree two cannot be an integer
multiplication, `(45)` and the characteristic-zero classification of
elliptic endomorphism algebras give

```text
End^0(E)=Q(sqrt(-2)).                                    (46)
```

On the target `(14)`, the automorphism

```text
beta(X,Y)=(zeta_3 X,Y)
```

satisfies `beta^2+beta+1=0`, so

```text
End^0(E_0)=Q(sqrt(-3)).                                  (47)
```

A nonzero Hom between elliptic curves in characteristic zero is an isogeny
and identifies their rational endomorphism algebras. The two quadratic fields
in `(46)--(47)` differ because `2/3` is not a square in `Q`. Hence

```text
Hom(E,E_0)=0,
Hom(Jac(C_B),E_0)=0.                                    (48)
```

By the observer, every face-stable weight-eight face is impossible. This
includes the two-term Bolza face and both singleton endpoints. For
`c+d=0`, this conclusion still requires the face-stable boundary hypothesis;
equation `(34)` alone does not supply it.

### 5.3 Primary-source corroboration, not a dependency

**CITED / PRIMARY CORROBORATION ONLY.** Enric Florit and Benjamin Smith,
[*An atlas of the Richelot isogeny graph*, arXiv:2101.00917v2,
Section 4.14](https://arxiv.org/html/2101.00917v2#S4.SS14), identify the
Type-VI normal form `y^2=x(x^4+1)` and state that its special Type-`Sigma`
Richelot neighbour is `E^2` with `j(E)=8000`. Their Section 2.1 gives the
general degenerate Richelot product formula. No table, finite-characteristic
specialization, or CM assertion from that paper is used in `(35)--(48)`;
the split and Hom mismatch above are independent exact derivations.

## 6. Sharp combination with the live third row

Now restrict to THM-4007's old minimal-support seam

```text
b=d=0,               A5=a^5,              Rtilde=R/gamma. (49)
```

The already proved rows give

```text
[p^3]Rtilde=2752/(135A5^3),

[p^4]Rtilde+(6/(7A5))[y^2]Rtilde
    =-11392/(105A5^4).                                  (50)
```

THM-4016 studies a different specialization of this affine line: the sharp
fixed-support cancellation has nonzero `p^4`, attachment ratios
`(X^3,Y^2)=(43/84,127/84)`, and only a conditional stable-model application.
The exact `M=6` branch below instead sets `p^4=0` and has ratios
`(43/224,267/224)`. Neither arithmetic certificate may be substituted for
the other.

In particular, the **total** polynomial weighted degree of `H` is at least
six. Assume, as an exact additional hypothesis, that it is at most eight.
There are then only three maxima.

1. If `M=6`, both `[p^2y]Rtilde` and `[p^4]Rtilde` vanish, as does every
   weight-eight coefficient. Equation `(50)` forces

   ```text
   [y^2]Rtilde=-5696/(45A5^3).                           (51)
   ```

   Thus the unique top face is

   ```text
   H_6/gamma
    =(2752/135A5^3)p^3-(5696/45A5^3)y^2,                (52)
   ```

   and its coefficient sum is

   ```text
   -14336/(135A5^3)!=0.                                 (53)
   ```

   This is exactly the nonresonant `j=0` survivor of Section 4. Multiplying
   both coefficients by `gamma` gives the raw face; no normalization is
   silently mixed. The boundary audit `(29l)--(29r)` proves face-stability,
   and the attachment calculation `(29e)--(29k)` then excludes this specific
   survivor unconditionally within the exact max-weight-six lane.

2. If `M=7`, the only weight-seven monomial is `p^2y`. Its coefficient is
   nonzero by definition of `M`, so the top face is singleton and fails the
   observer under face-stability.

3. If `M=8`, the only top monomials are `p^4` and `p*y^2`. A singleton
   endpoint fails Section 3, while a two-term face fails the Bolza Hom test
   `(48)`, again under face-stability.

The row lock and Sections 4.2--4.3 first give the unconditional low-face
exclusion

```text
b=d=0 ==> M!=6.                                         (54)
```

Since THM-4007 already forces `M>=6` on this seam, `(54)` is the
unconditional new floor `M>=7`.

If the actual highest face at `M=7` or `M=8` is face-stable, the singleton
and Bolza observers then give

```text
b=d=0 and face-stable at the actual highest face
    ==> total wt_(2,3)(H)>=9.                            (55)
```

Most importantly, the rows currently computed in THM-4007 do not prove an
upper bound on `M`: an actual polynomial may have later terms of weight at
least nine. Equations `(54)--(55)` concern the actual total polynomial, not a
finite truncation. The max-weight-six transfer is now unconditional, but
without face-stability at weights seven and eight, `(55)` is not an
unconditional floor for every formal residual.

## 7. Information-loss and mistakes audit

The observer deliberately forgets lower weighted faces and the attachment
map. Its typed ledger is

```text
source:       exact polynomial H and a hypothetical finite fibre map;
target:       normalized components of the highest face D_M;
map:          weighted degeneration (6)--(11);
preserved:    positive-genus face Jacobian factors under face-stability;
lost:         lower faces, boundary tails, node graph, attachment images;
sidecar:      semistable boundary resolution plus attachment specialization;
test:         Hom(B_M,E_0).                              (56)
```

The following recorded failure modes are explicitly avoided.

- `MISTAKE-244`: a weighted face is not called a coordinate or a terminating
  descent. Here it supplies only the conditional observer `(16)`.
- `MISTAKE-455`: a homogenized or cleared equation is not automatically a
  saturated model at the exceptional divisor. This is exactly why
  face-stability is a named hypothesis.
- `MISTAKE-453`: a target row-coordinate or elliptic quotient does not count
  the whole source boundary. No such transfer is made here.
- `MISTAKE-481`: the displayed attachment points are not declared a census of
  all source addresses. They are used only for the stated transversality
  control.
- `MISTAKE-483`: the finite replay universe is reported exactly below and is
  not substituted for the symbolic all-`M` proof.
- `MISTAKE-479`: the primary checker uses explicit runtime gates, not Python
  `assert`; optimized mode retains every check.

## 8. Exact certificate and present status

The primary certificate verifies:

1. source scaling, weighted homogenization, central factorization, and target
   good reduction for every `M=2,...,12` in a dense universal coefficient
   polynomial;
2. all `36` singleton monomials of weights `2,...,18`, including component
   counts, attachment counts, and gradient determinants;
3. the weight-six elliptic change, live raw/normalized conversion, exact
   nonresonance `(53)`, the complete `S=1` boundary/Hensel ledger
   `(29m)--(29q)`, attachment coordinates `(29f)`, and reductions of exact
   orders six and nine at `11` and `17`;
4. the weight-eight Bolza change, both involutions, spanning quotient
   differentials, `c_4=160`, `Delta=512`, `j=8000`, the explicit
   `alpha^2=[-2]` endomorphism, and the CM Hom mismatch; and
5. the THM-4007 affine lock and total-weight-at-most-eight ledger.

Reproduction:

```bash
python3 04-computation/jc2_weighted_face_good_factor_thm4012.py
python3 -O 04-computation/jc2_weighted_face_good_factor_thm4012.py
python3 -B 04-computation/jc2_weighted_face_good_factor_thm4012_independent_audit.py
python3 -B -O 04-computation/jc2_weighted_face_good_factor_thm4012_independent_audit.py
```

Both streams byte-match the frozen output and end with
`ALL THM-4012 EXACT CHECKS PASSED` after `206` runtime gates.

**Proved content:** the conditional face-stable Hom observer, singleton
closure in that universe, exact weight-six Hom survival, unconditional
forced max-six attachment no-go and floor `M>=7`, conditional weight-eight
Bolza no-go, and conditional floor `(55)`.

**Not proved:**

1. face-stability for an arbitrary highest face or arbitrary lower-face
   deformation;
2. that nonresonant finite attachments alone imply face-stability;
3. face-stability at weights seven or eight, or at any face beyond the exact
   nonresonant max-weight-six family proved in Section 4.3;
4. that the live residual has total weighted degree at most eight;
5. exclusion of faces of weight at least nine, the full reduced `(2,3)` cell,
   or `JC(2)`.
