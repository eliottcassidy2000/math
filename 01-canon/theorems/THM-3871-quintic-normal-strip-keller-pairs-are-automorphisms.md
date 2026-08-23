---
id: THM-3871
title: "Quintic normal-strip Keller pairs are automorphisms"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, a polynomial Keller pair in k[s,z] whose two
  transverse z-degrees are at most five is a polynomial automorphism.  After
  a constant target change, the (5,1) row shears to THM-3867 and the apparent
  (5,2), (5,3), and (5,4) rows are empty.  The (5,4) proof uses two conserved
  polynomials and simultaneous regularity of both original arm coefficients.
  Arbitrary polynomial Keller pairs, rational or infinite normal series, and
  the planar Jacobian conjecture remain open.
source: root / quintic_normal_scout / planar JC quintic strip lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE PROOF/REPLAY AUDIT PASSED.  The independent companion
  reconstructs all ten buckets, both constant target charts, every Kummer and
  zero stratum, the (5,0) contradiction, the (5,1) shear, all depressed
  integrations and conserved rows, every finite DVR channel, and both
  constant-scale infinity channels.  Its 16,006 optimization-safe gates
  verify the resultants 441 and 3171942400000, three rational hostiles, the
  arms-only false positive, and arbitrary characteristic-zero base fields.
  Normal and optimized executions of both companions byte-match their frozen
  LF outputs.
depends_on:
  - THM-3867-quartic-normal-strip-keller-pairs-are-automorphisms
related:
  - THM-3861-cubic-normal-strip-keller-pairs-are-automorphisms
  - THM-3856-quadratic-normal-strip-keller-pairs-are-automorphisms
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3846-formal-arm-darboux-lift-and-algebraization-gate
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
script: 04-computation/jc2_quintic_normal_strip_keller_thm3871.py
output: 05-knowledge/results/jc2_quintic_normal_strip_keller_thm3871.out
script_sha256: c110f2d6a1a9e8ae240cb73a6f3ccd13ca132e7793b9a33a85c66a79b27fd83c
output_sha256: f6552fb01d1bcf590b38320953b23b895928466971f3f50049d50fb69ea1107f
semantic_sha256: fbdb1b1b7e784350f14f41a88aa5e890112a2f248cd2b9fade51288e50e8ce85
independent_script: 04-computation/jc2_quintic_normal_strip_keller_independent_audit_thm3871.py
independent_output: 05-knowledge/results/jc2_quintic_normal_strip_keller_independent_audit_thm3871.out
independent_script_sha256: 79adf022910ec505ffb8bb8094955a76f34a413657075a96bbfc8099304feb61
independent_output_sha256: 8da278fd6756731468e899e88d4687e61fe9bf3c166c812e9acbb59dac452d1d
independent_semantic_sha256: 0aeeb394abbde44ccb1ef7a3fe3051a802d3d00125f6eaf2c56cea904f88fd09
hash_basis: raw LF bytes
---

# THM-3871 -- quintic transverse depth is still triangular

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be a
field of characteristic zero.  Suppose

```text
A,C in k[s,z],                 deg_z A,deg_z C <=5,
J_(z,s)(A,C)=A_z C_s-A_s C_z=lambda in k*.                 (1)
```

Then `(A,C)` is a polynomial automorphism.  More precisely, after a constant
target `SL_2(k)` change, every genuinely quintic pair has `deg_z A=5` and lies
in one of the following rows.

1. In row `(5,1)`, there are `rho_5,...,rho_1,beta in k`, with
   `rho_5 beta!=0`, a polynomial `b in k[s]`, and `rho_0 in k`, such that

   ```text
   C=b+beta z,
   A=rho_5 C^5+rho_4 C^4+rho_3 C^3+rho_2 C^2+rho_1 C
     -(lambda/beta)s+rho_0.                                (2)
   ```

2. The apparent rows `(5,2)`, `(5,3)`, and `(5,4)` are empty in `k[s,z]`.

The new phenomenon at degree five occurs in `(5,4)`: simultaneous arm
cancellation by itself has a genuine one-parameter family.  The last two
nonconstant Jacobian buckets give two conserved polynomials, and those
conservation laws have no common point with that arm-cancellation family.

## 1. The ten exact coefficient buckets

Write

```text
A=a+alpha z+u z^2+p z^3+r z^4+w z^5,
C=b+beta z+v z^2+q z^3+t z^4+ell z^5.                      (3)
```

For `m=0,...,9`, the coefficient of `z^m` in the Jacobian is

```text
E_m=sum_(i+j=m+1) (i a_i b_j'-j a_i' b_j).                 (4)
```

In expanded form the ten buckets are

```text
E_0 = alpha b'-a'beta,                                     (5)

E_1 = alpha beta'-alpha'beta+2u b'-2a'v,                   (6)

E_2 = alpha v'-2alpha'v+2u beta'-u'beta
      +3p b'-3a'q,                                         (7)

E_3 = alpha q'-3alpha'q+2u v'-2u'v+3p beta'-p'beta
      +4r b'-4a't,                                         (8)

E_4 = alpha t'-4alpha't+2u q'-3u'q+3p v'-2p'v
      +4r beta'-r'beta+5w b'-5a'ell,                       (9)

E_5 = alpha ell'-5alpha'ell+2u t'-4u't+3p q'-3p'q
      +4r v'-2r'v+5w beta'-w'beta,                        (10)

E_6 = 2u ell'-5u'ell+3p t'-4p't+4r q'-3r'q
      +5w v'-2w'v,                                        (11)

E_7 = 3p ell'-5p'ell+4r t'-4r't+5w q'-3w'q,              (12)

E_8 = 4r ell'-5r'ell+5w t'-4w't,                         (13)

E_9 = 5(w ell'-w'ell).                                    (14)
```

Thus `E_0=lambda` and `E_1=...=E_9=0`.  The count is ten, not nine:
before normalization the Jacobian can have `z`-degree nine.

If `(w,ell)=(0,0)`, THM-3867 applies.  Otherwise `(14)` says that the
nonzero pair `(w,ell)` has constant direction over `k`, because the constant
field of `k(s)` is `k`.  A constant target `SL_2(k)` change therefore puts

```text
ell=0,                         w!=0.                         (15)
```

Let `j=deg_z C<=4`, and let `c_j` be its nonzero top coefficient.  The
highest surviving bucket is

```text
5w c_j'-j w'c_j=0.                                         (16)
```

Since `gcd(5,j)=1` for `j=1,2,3,4`, prime valuations in the UFD `k[s]`
give, without adjoining scalar roots,

```text
w=R h^5,                      c_j=Q h^j,
R,Q in k*,                    0!=h in k[s].                (17)
```

Put `y=hz`.  For arbitrary profiles `F,G`, the Euler terms cancel:

```text
J_(z,s)(F(hz,s),G(hz,s))=h J_(y,s)(F,G).                  (18)
```

Thus the rational profiles in `k(s)[y]` have bracket

```text
delta=lambda/h.                                            (19)
```

A polynomial target shear `A->A-P(C)` can cancel the fifth transverse
coefficient only if

```text
j deg(P)=5.                                                (20)
```

Among `j=1,2,3,4`, this occurs only for `j=1`.  This divisibility check is
the sharp degree-drop boundary; the other three rows require new equations.

## 2. The zero and linear rows

If `j=0`, then `C=b(s)`.  The `z^4` coefficient of the Jacobian is
`5w b'=0`, contradicting `(1)`.  Thus `(5,0)` is empty.

If `j=1`, write the top coefficient as `beta=Qh`.  Equation `(17)` gives

```text
w=(R/Q^5) beta^5.                                         (21)
```

The target shear

```text
A_1=A-(R/Q^5)C^5                                          (22)
```

preserves the Jacobian and lowers `deg_z A_1` to at most four.  THM-3867
then makes `(A_1,C)` an automorphism.  Its `(4,1)` classification, including
all lower-degree edges, gives `(2)`.  Composing with the inverse target shear
shows that `(A,C)` is an automorphism.

## 3. The apparent `(5,2)` row

Assume `j=2`.  After `(17)--(19)`, write the quadratic profile as
`Qy^2+Vy+b`.  Set

```text
g=V/(2Q),                        x=y+g.                    (23)
```

The next bucket makes the `x^4` coefficient of the quintic constant.  Hence

```text
A_*=R x^5+D x^4+U x^3+M x^2+Lx+a,
C_*=Q x^2+b,                      D in k.                  (24)
```

The five remaining buckets are

```text
5R b'-2Q U'=0,
4D b'-2Q M'=0,
3U b'-2Q L'=0,
2M b'-2Q a'=0,
L b'=delta.                                                (25)
```

Put `K=5R/(2Q)`.  Exhaustive integration gives

```text
U=Kb+E,
M=(2D/Q)b+F,
L=(3K/(4Q))b^2+(3E/(2Q))b+G,
a=(D/Q^2)b^2+(F/Q)b+A_0,             E,F,G,A_0 in k.      (26)
```

Suppose `h` is nonconstant and fix an irreducible `pi|h`.  If `b` is
regular at `pi`, `(25)` makes `delta` regular, contrary to `(19)`.  If `b`
has pole order `n>0`, its quadratic term in `L` is uniquely leading, so

```text
ord_pi(delta)=-(3n+1).                                     (27)
```

The original `C` arm is `C_*(g)=Qg^2+b`.  Its regularity forces `g` to have
pole order `n/2` and forces `b=-Qg^2` at leading order.  The order `5n/2`
part of the original `A` arm is then

```text
R g^5+K b g^3+(3K/(4Q))b^2 g
  =(R-KQ+3KQ/4)g^5=(3R/8)g^5!=0.                          (28)
```

Thus `h` cannot be nonconstant.  If `h` is a unit, then `b` is a polynomial.
For `deg b=d>0`, the last product in `(25)` has degree `3d-1` with nonzero
leading coefficient; if `d=0`, it vanishes.  Neither case equals a nonzero
constant.  The `(5,2)` row is empty, including `b=0` and every constant
stratum.

## 4. The apparent `(5,3)` row

Assume `j=3`.  Depress the cubic by `x=y+V/(3Q)`.  Again the next bucket
makes `D` constant, and the profiles become

```text
A_*=R x^5+D x^4+U x^3+M x^2+Lx+a,
C_*=Q x^3+Bx+b.                                            (29)
```

Put

```text
K=5R/(3Q),                        H=4D/(3Q).                (30)
```

The six remaining buckets, from high to low, are

```text
5R B'-3Q U'=0,
5R b'+4D B'-3Q M'=0,
4D b'+3U B'-U'B-3Q L'=0,
3U b'+2M B'-M'B-3Q a'=0,
L B'-L'B+2M b'=0,
L b'-a'B=delta.                                            (31)
```

The first four integrate to

```text
U=KB+E,
M=HB+Kb+F,
L=(K/(3Q))B^2+(E/Q)B+Hb+G,
a=(H/(6Q))B^2+(2K/(3Q))Bb+(2F/(3Q))B+(E/Q)b+A_0.          (32)
```

The penultimate bucket is an exact derivative:

```text
-(1/(9Q^2)) P'=0,                                         (33)

P=(5R/3)B^3-12DQBb-18FQ^2b-9GQ^2B-15QRb^2.              (34)
```

The constant bucket is exactly

```text
9Q^2 delta=
 12DQ b b'-4D B^2 B'-6FQ B B'+9GQ^2 b'
 -5R B^2 b'-10R B b B'.                                  (35)
```

At a finite prime, if only one of `B,b` has a pole, `(34)` has a unique
leading term.  If their pole orders are `m,n>0`, respectively, the three
possible leading orders are

```text
3m,                         m+n,                         2n. (36)
```

The only possible tie not dominated by the third term is

```text
n=3m/2,                       B_0^3=9Q b_0^2.              (37)
```

Here `B_0,b_0` denote leading residues.  The last two terms of `(35)` have
the common pole order `7m/2+1`, and their leading coefficient is a nonzero
multiple of

```text
(5n+10m)R B_0^2b_0.                                      (38)
```

Thus a prime of `h` can be paid only in this balanced channel.  Regularity
of the original `C` arm forces `g` to have pole order `m/2`.  Define the
nonzero leading ratios

```text
X=B/(Qg^2),                       Z=b/(Qg^3).               (39)
```

Equations `(37)` and `C_*(g)=0` give

```text
X^3=9Z^2,                         1+X+Z=0.                 (40)
```

The normalized leading `A` arm is

```text
A_*(g)/(Rg^5)=
1+(5/3)(X+Z)+(5/9)X^2+(10/9)XZ.                           (41)
```

If it vanished, `(40)` would give the two equations

```text
5X^2+10X+6=0,
X^3-9(1+X)^2=0.                                           (42)
```

Their exact resultant is

```text
Res_X(5X^2+10X+6, X^3-9(1+X)^2)=441!=0.                  (43)
```

This integer remains nonzero over every characteristic-zero residue field;
no algebraic closure or square/cube root is being assumed.  Hence the two
arms cannot both be regular, and `h` is a unit.

For constant `h`, `B,b` are polynomials.  The same comparison at infinity
forces `(37)` whenever either is nonconstant.  The leading part of `(35)`
then has degree `7m/2-1` and coefficient
`-(5n+10m)R B_0^2b_0`, so it cannot be a nonzero constant.  If `B,b` are
both constant, `(35)` is zero.  This closes every zero and constant stratum
of `(5,3)`.

## 5. The depressed `(5,4)` packet

Assume finally `j=4`.  Equations `(17)--(19)` give top terms `Rh^5,Qh^4`.
In the scaled variable, write the next coefficient of `C` as `V`.  Put

```text
g=V/(4Q),                         x=y+g.                    (44)
```

The next bucket `5RV'-4QP'=0` makes the translated quartic coefficient of
`A` constant.  Thus

```text
A_*=R x^5+D x^4+U x^3+M x^2+Lx+a,
C_*=Q x^4+B x^2+Nx+b,               D in k.                (45)
```

The seven remaining buckets are exactly

```text
F_6=5R B'-4Q U'=0,
F_5=5R N'+4D B'-4Q M'=0,
F_4=5R b'+4D N'+3U B'-2U'B-4Q L'=0,
F_3=4D b'+3U N'-U'N+2M B'-2M'B-4Q a'=0,
F_2=3U b'+2M N'-M'N+L B'-2L'B=0,
F_1=L N'-L'N+2M b'-2a'B=0,
F_0=L b'-a'N=delta.                                       (46)
```

Set

```text
K=5R/(4Q),                         H=D/Q.                   (47)
```

The first four rows integrate exhaustively to

```text
U=KB+E,
M=HB+KN+F,
L=(K/(8Q))B^2+(3E/(4Q))B+HN+Kb+G,
a=(K/(4Q))BN+(F/(2Q))B+(3E/(4Q))N+Hb+A_0.                (48)
```

The next two rows are exact derivatives of two conserved polynomials:

```text
F_2=-(1/(32Q^2)) mathcal P',
F_1=-(1/(32Q^2)) mathcal T',                              (49)

mathcal W=5RB+12EQ,

mathcal P=
 mathcal W(B^2-8Qb)-20QRN^2-64FQ^2N-32GQ^2B,             (50)

mathcal T=
 N(15RB^2+24EQB-40QRb-32GQ^2)
 +16FQ(B^2-4Qb).                                          (51)
```

Finally the constant row is the one-form

```text
32Q^2 delta=Omega,

Omega=(5RB^2+24EQB+40QRb+32GQ^2)b'
      -(24EQN+10RBN)N'
      -(16FQN+10RN^2)B'.                                 (52)
```

Equations `(48)--(52)` are the complete rational `(5,4)` packet.  In
particular, no division by `mathcal W`, `B`, `N`, or either conserved value
has occurred.

## 6. Nonconstant Kummer scale in `(5,4)`

Suppose `h` is nonconstant and fix an irreducible `pi|h`.  Derivatives of a
function with a pole of order `d` have pole order `d+1`, because `pi` is
separable in characteristic zero.  A function regular at `pi` has regular
derivative.

### 6.1 `mathcal W` has a pole

If `mathcal W` has pole order `m>0`, so does `B`.  Let `n,z` be the pole
orders of `N,b`, with zero meaning regular.  The possible leading orders in
`mathcal P` and `mathcal T` are

```text
mathcal P:       3m, m+z, 2n,
mathcal T:       n+2m, n+z, 2m, z.                         (53)
```

A direct maximum comparison forces

```text
z=2m,                         0<=n<=3m/2.                  (54)
```

If `0<n<3m/2`, the leading `mathcal P` equation gives
`b_0=B_0^2/(8Q)`.  The leading `mathcal T` coefficient is then

```text
N_0(15RB_0^2-40QRb_0)=10RN_0B_0^2!=0,                    (55)
```

a contradiction.

If `n=0`, the same `mathcal P` relation holds, while `mathcal T` fixes the
regular residue

```text
N_bar=-4FQ/(5R).                                           (56)
```

Regularity of `C_*(g)` forces `g` to have pole order `m/2`.  Put

```text
X=B/(Qg^2),                         Z=b/(Qg^4).              (57)
```

At leading order,

```text
Z=X^2/8,                         1+X+Z=0.                  (58)
```

The normalized `A`-arm residual is

```text
A_*(g)/(Rg^5)=(5X^2-8)/32.                                (59)
```

It cannot vanish on `(58)`: substituting `X^2=-8X-8` would force
`X=-6/5`, which does not satisfy `(58)`.

It remains to take the balanced value `n=3m/2`.  With

```text
X=B/(Qg^2),            Y=N/(Qg^3),            Z=b/(Qg^4), (60)
```

the two conserved leading equations and the `C` arm are

```text
X^3-8XZ-4Y^2=0,
Y(3X^2-8Z)=0,
1+X+Y+Z=0.                                               (61)
```

Here `Y!=0`, so

```text
Z=3X^2/8,                         Y^2=-X^3/2.               (62)
```

The normalized `A` arm is

```text
A_*(g)/(Rg^5)=(5X^2+10XY-8)/32.                           (63)
```

If `(63)` vanished, elimination of `Y` from `(61)--(63)` would give

```text
p(X)=15X^3+20X^2+40X+32=0,
q(X)=50X^5+25X^4-80X^2+64=0.                              (64)
```

The exact resultant is

```text
Res_X(p,q)=3171942400000!=0.                              (65)
```

Again this is a scalar-root-free contradiction over every characteristic-
zero residue field.

### 6.2 `mathcal W` is a unit at `pi`

Then `B` is regular.  If both `N,b` had poles, `mathcal P` would force
`pole(b)=2 pole(N)`, after which the `Nb` term of `mathcal T` would be
uniquely leading.  If only one had a pole, `mathcal P` would have a unique
leading `N^2` or `mathcal W b` term.  Thus `B,N,b` and `(52)` are regular,
contrary to the pole of `lambda/h`.

### 6.3 `mathcal W` has a finite zero

Suppose `ord_pi(mathcal W)=r_0>0`.  The same two conserved equations first
force `N` to be regular: if `N,b` both had poles, the `Nb` term of
`mathcal T` would be unique, and if only `N` had a pole, its square in
`mathcal P` would be unique.  If `b` has pole order `z>0`, conservation of
`mathcal P` also gives `z<=r_0`.

In `(52)`, the term `40QR b b'` is then uniquely leading, so `delta` has
pole order `2z+1`.  Regularity of `C_*(g)` forces `g` to have pole order
`z/4` and forces `b=-Qg^4` at leading order.  The leading `A` arm is

```text
R g^5+K b g=(R-KQ)g^5=-R g^5/4!=0.                        (66)
```

If `b` is regular, `(52)` is regular and gives the earlier contradiction.

### 6.4 `mathcal W` is identically zero

Then `B` is constant.  The conserved polynomial `mathcal P` is a quadratic
equation over `k` for `N`.  Since `k` is relatively algebraically closed in
`k(s)`, it forces `N in k`.  The remaining possibilities for `b` are exactly
the regular or pole alternatives in Section 6.3, and `(66)` closes the pole
alternative.  This includes every zero value of `E,F,G`, either conserved
constant, `B`, or `N`.

All behaviors of `mathcal W` at a prime of nonconstant `h` are impossible.
Therefore

```text
h in k*.                                                    (67)
```

## 7. Constant Kummer scale in `(5,4)`

Now the depression is polynomial and `B,N,b in k[s]`.  Suppose first that
`m=deg B>0`.  Comparing degrees in the two conserved polynomials repeats
`(54)--(56)`: either `deg N=0`, or `0<deg N<3m/2` is contradictory, or
`deg N=3m/2` is balanced.

In the constant-`N` channel, the leading coefficient of `Omega` is

```text
(5/2)(Rm/Q) B_0^4                                         (68)
```

at degree `4m-1`.  In the balanced channel, the leading relations are

```text
b_0=3B_0^2/(8Q),                   N_0^2=-B_0^3/(2Q),       (69)
```

and the leading coefficient of `Omega` is

```text
(55/2)(Rm/Q) B_0^4                                         (70)
```

at the same degree `4m-1`.  Both are nonzero in characteristic zero and
neither can equal the constant `32Q^2 delta`.

It remains to let `B` be constant.  If `mathcal W!=0`, the two conserved
polynomials force `N,b` to be constant.  More explicitly, eliminate `b`
from `mathcal P=constant` and substitute in `mathcal T`: the result is a
cubic polynomial over `k` in `N` whose leading coefficient is
`100QR^2/mathcal W!=0`.  Hence `N in k`, and `mathcal P` then gives
`b in k`.  This also covers cancellation among the apparent `Nb`, `b`, and
`N^2` leaders.  Then `Omega=0`.

If `mathcal W=0`, `mathcal P` first forces `N` to be constant.  A
nonconstant polynomial `b` leaves the term `40QRbb'` of degree
`2 deg(b)-1`; constant `b` again gives `Omega=0`.  Thus constant scale is
impossible as well, and the `(5,4)` row is empty.

## 8. Sharp rational hostiles and the apparent modulus

The denominator obstructions above are sharp.  They do not assert that the
rational depressed packets are empty.

### 8.1 A `(5,2)` hostile over `Q`

Take

```text
C=s^14 z^2+2s^6z,

A=s^35z^5+5s^27z^4+(15/2)s^19z^3+(5/2)s^11z^2
  -(5/8)s^3z+3/(8s^5).                                   (71)
```

Direct differentiation gives `J_(z,s)(A,C)=15/4`.  The `C` arm vanishes;
the only nonpolynomial coefficient is the forced `A` arm `3/(8s^5)` from
`(28)`.

### 8.2 A `(5,3)` hostile over an algebraic constant extension

Choose constants `X,Z` satisfying

```text
X^3=9Z^2,                         1+X+Z=0.                 (72)
```

Set `h=s^8`, `g=s^-1`, `B=Xs^-2`, `b=Zs^-3`, and take
`R=Q=1`, `D=E=F=G=0` in `(32)`.  Then `P=0` and

```text
delta=(35/9)X^2Z s^-8.                                    (73)
```

After returning to `z`, every nonarm coefficient is polynomial, the `C`
arm is zero, and the sole bad coefficient is

```text
s^-5(-6+5X^2+10XZ)/9,                                     (74)
```

which is nonzero by `(43)`.

### 8.3 A `(5,4)` hostile over an algebraic constant extension

Choose constants

```text
Z=3X^2/8,                 Y^2=-X^3/2,
1+X+Y+Z=0.                                                (75)
```

Equivalently `X` is a nonzero root of

```text
9X^4+80X^3+112X^2+128X+64.                               (76)
```

Set `h=s^9`, `g=s^-1`, `B=Xs^-2`, `N=Ys^-3`, `b=Zs^-4`, and
take `R=Q=1`, `D=E=F=G=0` in `(48)`.  Both conserved polynomials vanish
and

```text
delta=-(55/32)X^4 s^-9.                                  (77)
```

Again every nonarm coefficient is polynomial and the `C` arm vanishes.  The
sole bad coefficient is

```text
s^-5(5X^2+10XY-8)/32,                                    (78)
```

nonzero by `(65)`.

### 8.4 Why both arms are needed after all

At the raw leading-arm level, the equations are only

```text
1+X+Y+Z=0,                  5X^2+10XY-8=0.                (79)
```

They define a genuine one-parameter cancellation family.  For example,

```text
(X,Y,Z)=(1,3/10,-23/10)                                   (80)
```

cancels both arms.  But it gives the nonzero conserved leading residuals

```text
X^3-8XZ-4Y^2=476/25,
Y(3X^2-8Z)=321/50.                                        (81)
```

Thus degree five does create a new apparent modulus, but the penultimate two
Jacobian buckets destroy it.  Any proof using only one arm or only the top
Kummer equation would miss this mechanism.

## 9. Russell and planar-Jacobian boundary

This closes the next bounded polynomial cell after THM-3867: a canonical
Russell normal expansion that is polynomial in the transverse parameter
through degree five cannot realize the arm self-identification forced by
THM-3843.  The first unclassified bounded polynomial cell is degree six.

This theorem does not truncate an arbitrary element of the Russell
surface, whose canonical normal expansion can be infinite.  It does not
exclude rational higher-normal lifts such as THM-3860, arbitrary polynomial
Keller maps after a nonlinear source change, or the planar Jacobian
conjecture.  Those problems remain open.  **QED.**

## 10. Reproduction

```bash
python3 -B 04-computation/jc2_quintic_normal_strip_keller_thm3871.py
python3 -B -O 04-computation/jc2_quintic_normal_strip_keller_thm3871.py
python3 -B 04-computation/jc2_quintic_normal_strip_keller_independent_audit_thm3871.py
python3 -B -O 04-computation/jc2_quintic_normal_strip_keller_independent_audit_thm3871.py
```

Each primary execution must byte-match
`05-knowledge/results/jc2_quintic_normal_strip_keller_thm3871.out`; each
independent execution must byte-match
`05-knowledge/results/jc2_quintic_normal_strip_keller_independent_audit_thm3871.out`.
