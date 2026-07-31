# The exact square prefix kills every one-pole odd response component

**Status:** PROVED as THM-2745 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.  This strengthens THM-2741 from geometrically integral members to arbitrary
reducible or nonreduced members of the same highest-odd response family.
It uses the polynomial exact-square-prefix sidecar, which the abstract
response curve deliberately forgets.  It does not edit or promote THM-2741.

## 1. Exact conclusion and scope

Use the full chosen-sheet split degree-22 response family

```text
F_23=Phi_22+sum_j a_j h^(22-j)Phi_j-lambda h^23,
G_24=Psi_22+sum_j a_j h^(22-j)Psi_j-W h^24
```

in `P(1,2,3,4)_[h,d,q,s]`.  Suppose at least one odd `a_j` is nonzero,
and let `j` be the largest such index.  Put

```text
r=22-j,                         g=gcd(r,6).
```

The proved conclusion is

```text
no irreducible component of the reduced response curve can carry a
physical trajectory from a split polynomial exact-square prefix.       (1)
```

Since a map from the reduced source line kills nilpotents, `(1)` also
excludes a trajectory into a nonreduced response member.  Together with
THM-2725, the precise degree-22 residual is the all-even, zero-first-flux
boundary.  This is not a derivation of the degree-22 normal form for every
Keller pair and is not `JC(2)` or `DC(2)`.

## 2. The complete `h=0` divisor

Remove harmless nonzero constants from the top degree-22 observables.  Write

```text
Phi_22 = unit*q*P,
Psi_22 = unit*Q,
R_22   = unit*q*T,
```

where

```text
P=-84d^2q^4s+3dq^6+280dq^2s^3-21q^4s^2-84s^5,

Q=-224d^3q^6+3360d^2q^4s^2-336dq^6s-3360dq^2s^4
  +3q^8+560q^4s^3+224s^6,

T=-84d^3q^4s+9d^2q^6+280d^2q^2s^3-105dq^4s^2
  -84ds^5+3q^6s+70q^2s^4.                            (2)
```

If `q=0`, then `Q=224s^6`, so the only point is

```text
P_infty=[d:q:s]=[1:0:0].                              (3)
```

On `q!=0`, pass to the `q=1` index chart and then to the coarse
`mu_3`-invariants

```text
z=s^3,                         rho=d/s^2.
```

The two equations become

```text
Pbar=3rho-21+z(-84rho^2+280rho-84),

Qbar=3+z(-336rho+560)
       +z^2(-224rho^3+3360rho^2-3360rho+224).         (4)
```

Their exact lexicographic basis is

```text
200070rho+930804137984z^4-669330178048z^3
 +53434945536z^2-1142572032z-355401,

p5(z)=20141047808z^5-14386462720z^4+1089822720z^3
      -21288960z^2-35910z+81.                         (5)
```

The polynomial `p5` is squarefree.  The ideals obtained by adjoining `z`,
`rho`, or the Jacobian of `(Pbar,Qbar)` are all the unit ideal.  Thus `(5)`
gives exactly five distinct smooth coarse points, and at each of them

```text
q!=0,                 s!=0,                 d!=0,
ord(h)=1.                                                    (6)
```

After removing the nonzero factor `s`, the response at these points is

```text
Tbar=3+z(9rho^2-105rho+70)
       +z^2(-84rho^3+280rho^2-84rho).                 (7)
```

The ideal `(Pbar,Qbar,Tbar)` is the unit ideal.  Hence `R_22` is nonzero
at all five points.  The affine functions therefore have poles

```text
pole_order(q/h^3)=3,                 pole_order(R_25/h^25)=25. (8)
```

At `(3)`, THM-2741's exact highest-odd Newton face gives `3g` coarse
normalization branches, each with

```text
ord(h)=6/g,
pole_order(R_aff)=(150-6r)/g >= 8.                    (9)
```

The full hyperplane invoice is consequently

```text
deg div_0(h)=5*1+(3g)*(6/g)=23.                       (10)
```

This both matches weighted Bezout and proves that no infinity point or
branch is missing.

## 3. A physical image component has exactly one infinity point

Let `Y` be the normalization of the reduced irreducible component containing
the generic image of a physical source map.  The map is nonconstant because
the third-flux identity is

```text
U R_source'=kappa,                         kappa!=0.   (11)
```

Thus it extends to a finite surjective morphism

```text
gamma:P1_x -> Y.                                      (12)
```

The exact gcd `gcd(qP,Q)=1` also shows that the two projective response
equations have no common surface factor: any such factor would restrict at
`h=0` to a positive-degree common factor of the top forms (or would itself
be divisible by `h`).  Hence the nonconstant image closure is indeed one of
the reduced curve components.

Every projective curve component meets `h=0`: a positive-dimensional
projective component cannot be contained in the affine chart `h!=0`, and no
component lies inside `h=0` because the top forms in `(2)` are coprime.
Thus `Y` contains at least one of the normalization points in `(6)` or
`(9)`.

On the source, the rational-primitive classification in THM-2723 is

```text
U in C* and R_source is nonconstant affine-linear,
```

or, after translating `x`,

```text
U=u_0 x^m,
R_source=s_0+s_1 x^(1-m),             m>=2.           (13)
```

In either case `R_source` has exactly one pole on `P1`.  Since every target
infinity point in `(8)`--`(9)` is an `R_aff` pole and `(12)` is surjective,
different target infinity points have disjoint nonempty source fibres.
Therefore `Y` contains at most one infinity normalization point.  It
contains exactly one.

Let `e>=1` be the ramification index over that point.  Pullback multiplies
the target response pole order by `e`.  Equations `(8)`--`(9)` show that the
source pole order is at least `8e>1`.  The constant-`U` alternative in
`(13)` has pole order exactly one, so it is impossible.  The unique source
pole is therefore the finite root `x=0` of the nonconstant polynomial `U`.

This finite-location conclusion is load-bearing.  At a finite source point,
all coefficients of the original exact prefix

```text
P=H^2+L,
H=U^2 z_source^2+B z_source+C,
L=A z_source+E                                      (14)
```

are regular.  One may not assert the same at source infinity.

## 4. The two coordinates forgotten by the response curve

On the split sheet use

```text
w=U z_source+B/(2U),
H=w^2+d,                         L=q w-s.              (15)
```

Pull the weighted projective coordinates back to the source DVR above its
unique pole, passing if necessary to the finite local index cover on which
the weight-one scale `h` is a uniformizing monomial.  If `z_source` denotes
the original polynomial variable, write coefficientwise

```text
Wtilde=h*w=(hU)z_source+beta,             beta=hB/(2U).
```

Since affine `d,q,s` equal `d/h^2,q/h^3,s/h^4`, respectively, `(15)` gives
the exact identities

```text
h^2 H=Wtilde^2+d,
h^4 L=q Wtilde-s.                                      (16)
```

At the finite source pole, `H` and `L` in `(14)` are regular coefficientwise
in `z_source`, while `U` vanishes.  Thus `hU` is regular and has positive
valuation.  The constant coefficient of the first identity in `(16)` is

```text
beta^2+d=h^2 C.
```

It shows that `beta^2` is DVR-regular and hence that `beta` itself is
regular.  Let `omega_0` be its residue.  At any one of the five smooth
infinity points, taking constant coefficients and residues in both
identities gives

```text
omega_0^2+d_0=0,                  q_0 omega_0-s_0=0.  (17)
```

These are not extra equations imposed on an arbitrary abstract response
curve.  They are necessary lift equations for a physical polynomial
exact-square prefix.  They retain precisely the `H,L` coordinates discarded
by the two-flux quotient.

## 5. The five smooth points are incompatible with `(17)`

At a smooth point, `(6)` makes `q_0,d_0,s_0` nonzero.  Substitute

```text
d=-omega_0^2,                         s=q omega_0
```

into the two top equations in `(2)`.  Exact factorization gives

```text
P=-8omega_0^2q^5(56omega_0^3+3q),

Q=q^6(7168omega_0^6+896omega_0^3q+3q^2).            (18)
```

After removing the known nonzero factors, the resultant in `q` is

```text
Res_q(56omega_0^3+3q,
      7168omega_0^6+896omega_0^3q+3q^2)
 =-76608omega_0^6.                                   (19)
```

But `d_0=-omega_0^2!=0`, so `(19)` cannot vanish.  No physical image
component can contain one of the five degree-one infinity points.

## 6. The `G2` point is incompatible with `(17)`

At `P_infty`, the projective residues of `q` and `s` are both zero, so one
must not substitute them directly into `(17)`.  Let `e` be the source
ramification index over the chosen THM-2741 normalization branch.  On the
pulled-back `d=1` index cover, write the first nonzero terms as

```text
h=h_* t^(6e/g)+...,
q=q_* t^(re/g)+...,
s=s_* t^(re/g)+...,                         q_*s_*!=0.
```

The constant coefficient of the first identity in `(16)` still gives

```text
omega_0^2=-1.                                          (20)
```

The constant coefficient of the second identity is

```text
q beta-s=h^4 E.
```

Here `E` is regular because the source pole is finite.  Since
`r<=21<24`, the right side has valuation at least `24e/g`, strictly larger
than the common valuation `re/g` of `q` and `s`.  Comparing their first
nonzero coefficients therefore gives

```text
s_*=q_* omega_0.                                      (21)
```

Every THM-2741 branch has nonzero `q_*` and `s_*`, and its second-flux
tangent is the `G2` arrangement

```text
(q_*^2-s_*^2)(q_*^4-14q_*^2s_*^2+s_*^4)=0.           (22)
```

Substituting `(20)`--`(21)` into the left side of `(22)` gives

```text
q_*^6(1-omega_0^2)(1-14omega_0^2+omega_0^4)
 =32q_*^6!=0.                                         (23)
```

Thus none of the `3g` singular normalization branches can lift to the
polynomial exact-square prefix.  The obstruction is independent of which
odd seed is highest.

Sections 5 and 6 exhaust the unique infinity point required in Section 3,
proving `(1)`.

## 7. Reducible and nonreduced members

No integrality of the whole response member entered the proof.  The closure
of the generic image of `P1` is one reduced irreducible component, and its
normalization is the `Y` used above.  A morphism from a reduced scheme kills
nilpotents in the target coordinate ring, so it factors through the target
reduction.  Embedded or generically nonreduced structure therefore cannot
alter the source function `R_source`, create another infinity fibre, or
remove either leading identity in `(16)`.

The local infinity schemes themselves are reduced: the five points are
transverse, while THM-2741's equations have simple `G2` tangent lines and
reduced factors `h^r-cs^6` in characteristic zero.  This is a useful hostile
check, although factorization through the reduction already suffices.

## 8. Reproduction and exact residual

Run

```bash
python3 04-computation/jc2_degree22_highest_odd_componentwise_closure_20260728.py
python3 -O 04-computation/jc2_degree22_highest_odd_componentwise_closure_20260728.py
```

and compare both transcripts with

```text
05-knowledge/results/jc2_degree22_highest_odd_componentwise_closure_20260728.out.
```

The companion reconstructs the Faber observables, the five-point quotient
Groebner basis, all transversality and nonvanishing ideals, every highest-odd
response pole, and both exact-prefix contradictions.  It uses explicit
exceptions and no optimization-sensitive `assert` statements.

An independent hostile audit reconstructed the Faber coefficients by
generalized multinomial extraction, repeated the five-point elimination and
response checks, and specifically rejected the invalid use of `q=s=0`
projective residues at `P_infty`.  The repaired DVR comparison in Section 6
then passed its weighted-index-cover, ramification, reducibility, and
nilpotent-target audits.

The LF-normalized hashes after normal/optimized byte-identical replay are

```text
script  863923b4886f6ec7dee44e28e152109acc51fcfd778d550d391ee02fc577b48a
output  619d2e49acdf60f5e2aa3ba9076957e2e38e9078363c833f7484f38311df3bb9
```

The exact surviving chart after combining THM-2745 with THM-2725 is

```text
split polynomial exact-square prefix;
reduced degree 22;
all odd Faber coefficients zero;
chosen-sheet first flux lambda=0.                     (23)
```

The broader split branch, other reduced degrees, and `JC(2)`/`DC(2)` remain
open.
