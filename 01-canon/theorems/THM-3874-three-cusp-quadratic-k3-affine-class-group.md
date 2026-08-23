---
id: THM-3874
title: "Three-cusp quadratic K3 and its affine class group"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  THM-3854 three-cusp quadratic resolvent has an elliptic K3 model with
  fibres I6, I3, I4, I2, and III*.  Its Mordell--Weil group is C2 and its
  Neron--Severi discriminant has absolute value 72.  An explicit
  boundary-plus-ADE quotient, not the discriminant alone, gives affine class
  group Z and scalar units.  Hence the affine resolvent has no connected
  codimension-one-unramified cyclic cubic cover, and no normal finite-flat
  cubic or transitive quartic over the plane can have the three-cusp
  quintic as its sole simple branch divisor.  This closes those C3-resolvent
  lanes, not JC(2).
source: jc_quartic_c3_construct / post-THM-3869 elliptic-K3 lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_sparse_direct_search and root,
  2026-08-23).  The
  audit independently rederived all five Kodaira fibres, the K3/Mordell--Weil
  packet, each local UV model and finite-exceptional label, the complete
  boundary exact sequence, the primitive 3d+2e torsion-section relation, the
  unit norm argument, and the Kummer/S3 sole-branch implication.  The root
  strengthening audit separately reconstructed the full 20-by-21 boundary
  presentation, including every killed generator and the doubled E7 weight;
  its Smith form has twenty unit invariant factors and one free column.
  Normal and optimized runs byte-match the frozen output.
depends_on:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
  - THM-3869-three-cusp-square-residual-cardano-line-ramification
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
  - THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate
  - THM-3865-one-place-inverse-discriminant-resolvent-class-group
script: 04-computation/jc2_three_cusp_elliptic_k3_class_group_thm3874.py
output: 05-knowledge/results/jc2_three_cusp_elliptic_k3_class_group_thm3874.out
script_sha256: dd90fbbba3d86aee3a58df5de1d41ee9fd31ffe5527ce35ea70b9e2cee7e6da1
output_sha256: 0bfe2125e864b668bf965e9605cc27990cafbd7394315698781eca2cba234672
semantic_sha256: b135b50909a3ddb57590882e77ce5f7e94e4c034adc5b8b660be01d2b691d09c
hash_basis: raw LF bytes
---

# THM-3874 -- the affine boundary kills the K3 three-primary packet

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Retain
the irreducible three-cusp quintic

```text
Delta=81x^5+90x^4+25x^3+30x^2y^2+30xy^2-y^4+8y^2          (1)
```

from THM-3854 and its affine quadratic resolvent

```text
Q=Spec A,             A=k[x,y,w]/(w^2-Delta).                (2)
```

This is a normal surface.  Indeed `Delta` is irreducible by THM-3854, the
double cover is regular at the generic point of its branch divisor, and its
singular locus lies over the six isolated singular points of `Delta`.
The hypersurface is `S2`, so codimension-one regularity gives normality.

Then

```text
A^*=k^*,                         Cl(A)=Z.                    (3)
```

In particular

```text
Cl(A)[3]=0,                                                     (4)
```

and `Q` has no nontrivial connected cyclic cubic cover which is unramified
at every codimension-one point.  Consequently there is no normal finite-flat
cubic domain over `k[x,y]` whose field-discriminant divisor is exactly the
reduced divisor `V(Delta)`.

The same obstruction excludes every transitive quartic field over
`k(x,y)` whose sole tame field-discriminant divisor is the same reduced
`V(Delta)` with simple transposition inertia.

The mechanism is load-bearing.  The projective K3 has nonzero three-primary
discriminant data; one cannot infer `(4)` from its Neron--Severi determinant.
The affine boundary and the finite ADE exceptional curves leave two classes
`d,e`, and the two-torsion section imposes the primitive relation

```text
3d+2e=0.                                                     (5)
```

Thus the actual quotient is

```text
Cl(A)=Z d direct_sum Z e / <3d+2e> ~= Z,                    (6)
```

which is stronger than the needed absence of three-torsion.

## 1. The elliptic K3 model

Put

```text
a=x+1,                 F=15a^2-15a+4,
L=9a-5,                M=9a-4,                b=1-a.         (7)
```

THM-3869 proves the polarization identity

```text
F^2-a^3L^2=b^3M^2.                                         (8)
```

After multiplying `w` by a scalar square root of `-1`, the generic fibre of
`(2)` over `k(a)` is the monic even quartic

```text
W^2=y^4-2Fy^2+b^3M^2.                                      (9)
```

Its two points at infinity are rational.  Taking the point with
`W/y^2=1` as origin, the exact birational transformation

```text
X=2(W+y^2-F),                         Y=2Xy                 (10)
```

gives

```text
E: Y^2=X(X^2+4FX+4a^3L^2).                                 (11)
```

Indeed `(10)` substituted into `(11)` is `(9)` times
`-8(W+y^2-F)`.  Conversely

```text
y=Y/(2X),                         W=X/2-y^2+F               (12)
```

on the dense chart `X!=0`.  The second quartic point at infinity is the
visible two-torsion section

```text
T=(0,0).                                                     (13)
```

The smooth minimal models of the two generic-fibre presentations are
birational K3 surfaces and hence isomorphic.  We use the elliptic model
`(11)` to compute the projective divisor lattice of `Q`.

## 2. Complete fibre packet

For `(11)`, the standard Weierstrass invariants are

```text
c_4=-64(243a^5-1170a^4+1875a^3-1380a^2+480a-64),           (14)

Delta_E=-4096a^6(a-1)^3L^4M^2.                             (15)
```

The value of `c_4` is nonzero at each finite zero in `(15)`:

```text
c_4(0)=4096,               c_4(1)=1024,
c_4(5/9)=16384/729,        c_4(4/9)=4096/729.                (16)
```

Thus the four finite fibres are multiplicative.  At infinity put `u=1/a`
and make the K3 scaling `X=u^-4 X_infty`, `Y=u^-6Y_infty`.  The transformed
orders are

```text
ord_u(a_2,a_4,c_4,Delta_E)=(2,3,3,9).                       (17)
```

The equation is minimal because the `a_4` order is three.  Tate's
characteristic-zero table gives the exact configuration

```text
a=0          I_6,
a=1          I_3,
L=0          I_4,
M=0          I_2,
a=infinity   III*.                                          (18)
```

Their Euler numbers sum to

```text
6+3+4+2+9=24.                                               (19)
```

The minimal surface is therefore a K3.  The root rank is

```text
5+2+3+1+7=18.                                               (20)
```

Shioda--Tate and the bound `rho<=20` for a characteristic-zero K3 force

```text
rho=20,                         rank MW(E/k(a))=0.           (21)
```

The torsion is exactly `C2`.  At the additive `III*` fibre the torsion of
`E(k(a))` injects into the component group: the formal group and the
additive identity component have no torsion in characteristic zero.  That
component group is `C2`, while `(13)` supplies a nonzero element of order
two.  Hence

```text
MW(E/k(a))=C2.                                               (22)
```

The two infinity sections are disjoint.  The height-zero formula then also
pins down every component address of `T`: it meets component three of `I6`,
component two of `I4`, the nonidentity component of `III*`, and the identity
components of `I3` and `I2`.  The respective nonzero local contributions
are `3/2,1,3/2`, whose sum is `4=2 chi(O_X)`.

The trivial lattice is

```text
Triv=U direct_sum A5 direct_sum A2 direct_sum A3
          direct_sum A1 direct_sum E7.                       (23)
```

It has absolute discriminant

```text
6*3*4*2*2=288                                               (24)
```

and index `|MW_tor|=2` in `NS`.  Therefore

```text
|disc NS|=288/2^2=72.                                       (25)
```

In particular `NS^*/NS` still has three-primary part `(C3)^2`, inherited
from `A5` and `A2`.  This is the hostile against replacing the affine
calculation by a determinant slogan.

## 3. Which fibre components are affine exceptional curves

Let `X` denote the smooth elliptic K3.  Label the nonidentity components of
the `I6` fibre by the `A5` chain

```text
alpha_1--alpha_2--alpha_3--alpha_4--alpha_5                 (26)
```

and those of `I4` by

```text
beta_1--beta_2--beta_3.                                     (27)
```

The six singular points of the affine surface `(2)` are precisely the
double-cover rational double points over the three cusps and three nodes of
`Delta`.

The fibre addresses make their resolution components explicit.

* At `a=0`, the two cusps `y=2,-2` give two `A2` chains.  In `(26)` these
  are `(alpha_1,alpha_2)` and `(alpha_4,alpha_5)`.  The middle class
  `d=alpha_3` is the sole nonexceptional root component.
* At `a=1`, the cusp `y=0` supplies the full `A2` root lattice of the `I3`
  fibre.
* At `L=0`, the two nodes `y^2=8/27` give the two `A1` exceptional curves
  `beta_1,beta_3`.  The middle class `e=beta_2` survives.
* At `M=0`, the node `y=0` supplies the full `A1` root lattice of the `I2`
  fibre.

These labels can also be read directly from the special quartics.  At
`a=0` the two graph components meet at `y=2,-2`, giving the two length-two
resolution chains in an `I6` cycle.  At `L=0` the two graph components meet
at the two roots of `y^2-8/27`, giving the two end curves in an `I4` cycle.
The remaining two fibres have one singular address each and no surviving
nonidentity root component.

The affine chart removes both quartic infinity sections `O,T`.  Compactifying
the base removes the complete `III*` fibre.  Thus, in divisor-class terms,
the boundary kills

```text
O, T, U, E7,                                                 (28)
```

where `U` is generated by `O` and the fibre class.  Together with the finite
exceptional curves just listed, these are exactly the prime divisors killed
when the smooth projective resolution is descended to the normal affine
surface `Q`.

## 4. The boundary quotient

For a normal surface, strict transform gives the standard exact quotient

```text
Cl(A)=NS(X)/N,                                               (29)
```

where `N` is generated by the projective boundary and all exceptional
curves over affine singularities.  Apply the list in Section 3 first inside
`Triv`.  Every basis vector dies except

```text
d=alpha_3,                         e=beta_2.                 (30)
```

Therefore the quotient of `Triv` is initially `Z d direct_sum Z e`.

It remains to impose the unique index-two overlattice class furnished by
the section `T`.  That section meets component three of `I6`, component two
of `I4`, and the nonidentity component of `III*`.  Its Shioda correction is

```text
omega_3(A5)+omega_2(A3)+omega(E7).                           (31)
```

The two needed fundamental-weight identities are

```text
2omega_3(A5)=alpha_1+2alpha_2+3alpha_3+2alpha_4+alpha_5,
2omega_2(A3)=beta_1+2beta_2+beta_3.                          (32)
```

After `(28)` and the exceptional curves are killed, twice the torsion-section
relation becomes exactly

```text
3d+2e=0.                                                     (33)
```

There is no further relation: `NS/Triv` has order exactly two and is
generated by this section.  Since `(3,2)` is primitive,

```text
Z^2/<(3,2)> ~= Z.                                           (34)
```

Equations `(29)-(34)` prove the class-group statement in `(3)` and locate
precisely where the two projective three-primary classes disappear.

As a presentation-level hostile check, order the generators as

```text
O,F; alpha_1,...,alpha_5; A2_1,A2_2;
beta_1,beta_2,beta_3; A1; E7_1,...,E7_7; T.                (35)
```

Take the doubled torsion-section relation

```text
2T-2O-4F+2omega_3(A5)+2omega_2(A3)+2omega(E7)=0            (36)
```

together with the nineteen coordinate relations killing exactly the
boundary and affine exceptional generators listed above.  The resulting
`20`-by-`21` integer matrix has rank `20` and Smith diagonal

```text
1,1,...,1                                                    (20 times). (37)
```

Thus the full quotient has one free generator and no hidden torsion.  This
independently rules out a quotient loss in the shorter calculation `(33)`.

## 5. Units

Every element of `A` has a unique form

```text
r+s w,                         r,s in k[x,y].                (38)
```

If it is a unit, its quadratic norm is a unit of the base:

```text
r^2-s^2 Delta in k^*.                                       (39)
```

If `s!=0` and `r=0`, the norm is already nonconstant.  If both are nonzero,
the two nonconstant terms in `(39)` have total degrees

```text
2 deg r,                         2 deg s+5,                  (40)
```

because `Delta` has odd total degree five.  These degrees have opposite
parity and cannot cancel.  Thus `(39)` cannot be constant.  Hence `s=0`,
and then `r^2 in k^*` forces `r in k^*`.  This proves `A^*=k^*`.

## 6. No unramified cyclic layer and no degree-three or degree-four sole branch

Fix any integer `n>=2`.  Because `k` is algebraically closed, it contains
the `n`th roots of unity and every scalar is an `n`th power.  A cyclic
degree-`n` extension of `Frac(A)` unramified at all height-one valuations has
a Kummer representative `gamma` with

```text
div(gamma)=nD.                                               (41)
```

The class of `D` lies in `Cl(A)[n]`.  By `(3)` it is principal; after
dividing `gamma` by an `n`th power, one obtains a unit.  The unit is scalar
and hence again an `n`th power.  Thus the Kummer class is trivial.  There is
no nontrivial connected cyclic extension of any degree which is unramified
at every height-one valuation of `A`.

Finally suppose a normal finite-flat cubic domain over `k[x,y]` had sole
field-discriminant divisor `V(Delta)` with multiplicity one.  Its generic
Galois group would be `S3`: the cubic is transitive, and the odd irreducible
discriminant is nonsquare.  The quadratic subfield of its Galois closure is
exactly `Frac(A)`.  At the generic point of `Delta`, transposition inertia is
absorbed by this quadratic layer, so the remaining `C3` layer is unramified;
away from `Delta` there is no inertia by hypothesis.  This is the forbidden
cover of the preceding paragraph.  Therefore no such normal cubic exists.

More generally, let `M/k(x,y)` be a finite normal closure with transitive
permutation group `G`, assume that every codimension-one inertia group is
generated by a transposition and that `V(Delta)` is the sole branch divisor,
and put

```text
H=G intersection A_n.                                       (42)
```

The sign quadratic subfield is `k(x,y)(sqrt(Delta))=Frac(A)`.  After this
base change, every transposition inertia disappears because its intersection
with `H` is trivial.  Therefore every quotient of `H` is unramified at all
height-one valuations of `A`.  If `H` has a nontrivial cyclic quotient, the
preceding Kummer argument contradicts `(3)` and the scalar-unit result.

For degrees three and four this closes every transitive possibility.  In
degree three the nonsquare discriminant forces `G=S3`.  In degree four the
transitive groups are `S4,A4,D4,V4,C4`; the nonsquare discriminant excludes
`A4,V4`, while `C4` contains no transposition and cannot realize the stated
simple inertia.  For the remaining cases

```text
G=S3:  H=C3,
G=S4:  H=A4 and H/V4=C3,
G=D4:  H=V4 and H has a C2 quotient.                       (43)
```

Thus no transitive cubic or quartic field has sole simple branch `Delta`.
The mechanism stops sharply at the next full symmetric group: `A5` is
perfect.  This explains why the natural `S5` completion in THM-3854 survives
the sign/cyclic-quotient test even though it fails for its separate monodromy
and atlas reasons.  Higher transitive groups with perfect `H` are outside
this corollary.

There is a useful presentation-independent strengthening.  Suppose an
irreducible depressed cubic has

```text
f=T^3-3PT+2Q,             P^3-Q^2=Delta H_0^2,
gcd(H_0,Delta)=1.                                           (44)
```

For a prime `p|H_0`, the field-discriminant exponent has even parity because
passing from the displayed order to the maximal order changes a discriminant
by a square.  In a tame cubic extension the only ramified generic
decomposition types are `(2,1)`, of discriminant exponent one, and total
ramification `(3)`, of exponent two.  Hence every `p|H_0` is either unramified
or totally ramified in the normalized cubic field.  They cannot all be
unramified, because that would give the sole-branch cubic just excluded.
Thus at least one irreducible factor of `H_0` necessarily survives as a totally
ramified cubic payment.  THM-3869 realizes this boundary with
`H_0=(9x+4)^2` at the displayed-order level and total ramification over the
line `9x+4=0` after normalization.

THM-3869 shows the equality boundary sharply: its explicit nonmonogenic
normal cubic pays the extra line `9x+4`, and its cyclic layer ramifies over
both primes above that line.  THM-3874 proves that some extra divisorial debt
is unavoidable for every cubic completion of this branch, not just for the
displayed Cardano representative.

## 7. Scope and reproduction

This theorem closes every normal degree-three and degree-four sole-simple-
branch lane for the exact THM-3854 three-cusp curve.  It does not exclude
higher-degree completions, different branch curves, nonnormal orders before
normalization, or any general planar Keller map.  In particular `JC(2)`
remains open.

Run

```bash
python3 04-computation/jc2_three_cusp_elliptic_k3_class_group_thm3874.py
python3 -O 04-computation/jc2_three_cusp_elliptic_k3_class_group_thm3874.py
```

and compare both streams byte-for-byte with
`05-knowledge/results/jc2_three_cusp_elliptic_k3_class_group_thm3874.out`.
The companion uses exact rational polynomial and integer-lattice arithmetic
and has no inactive `assert` gates.
