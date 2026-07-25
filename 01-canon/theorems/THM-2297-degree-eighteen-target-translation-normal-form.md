---
id: THM-2297
title: "Degree-eighteen target-translation normal form and sparse-plane closure"
status: >
  PROVED + VERIFIED-EXACT. The centered coordinate y=9s-2alpha in
  THM-2262 is the exact invariant of the legal target translation
  P -> P+2alpha/9. In the translated Faber basis alpha vanishes and the
  remaining degree-eighteen spectrum depends on four explicit invariants
  (B,C,D,W) of weights (2,3,4,5). The spectral cubic has weight six and its
  branch discriminant has weight twelve. Every stratum on which at most one
  of B,C,D,W is nonzero is empty: the B- and D-axes have squarefree branch
  discriminant, the C-axis normalizes to a smooth genus-one cubic, the W-axis
  to a genus-four superelliptic cover, and the invariant origin is a monomial
  cusp killed by the exact degree-eighteen polynomial Faber sidecar. More
  strongly, the whole plane B=D=0 is empty: off its axes its normalization is
  a connected trigonal cover of genus at least two, uniformly through all
  critical-value collisions. Hence every survivor uses at least two
  normalized parameters and at least one of B,D is nonzero. This does not
  close the residual multi-parameter singular locus or prove JC(2).
source: codex-2026-07-25-degree18-target-translation-gauge
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
related:
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2285-centered-grid-footprint-and-generic-keller-lines
script: 04-computation/jc2_degree18_target_translation_normal_form_thm2297.py
output: 05-knowledge/results/jc2_degree18_target_translation_normal_form_thm2297.out
script_sha256: 2d8e42ccb3b3d4d1e698c22bc92a0909b72da191a013a330b4a6d61842b9703b
output_sha256: 476d04c5efbae1d06c77b37b44d632d868bf5f7202e0795611dbdc08467d6b9e
hash_basis: working-tree bytes (LF)
---

# THM-2297 -- the degree-eighteen spectrum has four invariant coordinates

THM-2262 reduces the genuine nonsplit, polynomial exact-square-prefix,
degree-eighteen branch to a singular trigonal spectral curve. Its centered
coordinate

```text
y=9s-2alpha
```

was introduced there as the coefficient that makes the first flux linear in
`Z=T^2`. It has a stronger meaning: it is the exact coordinate left after a
legal translation of the first Keller target. This removes one spurious
parameter, exposes a weighted cone, and makes the first singular strata
small enough to close.

## 1. Inherited degree-eighteen state

Use the notation and scope of THM-2262:

```text
P=w^4+2d w^2+q w+(d^2-s),       q^2=T!=0,

Q=E_18+alpha E_14+beta E_10+gamma E_6+delta E_2,

u=dT,                            Z=T^2,

Phi_Q=0,                         Psi_Q=Psi0 in C,

A(2T F'+F T')=2kappa T,          F=R_Q/q,
                                                    (1)
```

where `kappa!=0`, the prime is differentiation in `x`, and the quadratic
deck sends `q` to `-q`. The wall `y=0` is already empty by THM-2262, so every
trajectory considered below has `y!=0` in `C(x)`.

For `j=1,...,5`, write

```text
mathcal E_j(P)=E_(4j-2)(P).
```

By definition this is the polynomial part at fibre infinity of
`P^(j-1/2)`.

## 2. The finite Faber translation law

For a constant `c in C`, put

```text
P_c=P+c.
```

This preserves monicity, the exact-square prefix, the nonsplit deck, and the
Keller bracket:

```text
[P_c,Q]=[P,Q]=kappa.                                (2)
```

It replaces `s` by

```text
s_c=s-c                                             (3)
```

and leaves `d,q,T,u,Z,A,V` unchanged.

The binomial series gives the exact finite identity

```text
mathcal E_j(P_c)
 =sum_(k=0)^(j-1) binom(j-1/2,k)c^k mathcal E_(j-k)(P).
                                                            (4)
```

Indeed,

```text
(P+c)^(j-1/2)
 =sum_(k>=0) binom(j-1/2,k)c^k P^(j-1/2-k).
```

When `k>=j`, the exponent is at most `-1/2`; since `P=w^4+O(w^2)`,
that summand has only strictly negative powers of `w`. Taking polynomial
parts therefore truncates the series exactly at `j-1`. There is no analytic
convergence issue.

To rewrite the fixed polynomial `Q` in the basis attached to `P_c`, apply
(4) with `P=P_c-c`. The coefficient of `E_14(P_c)` is

```text
alpha_c=alpha-(9/2)c.                               (5)
```

Choose

```text
c=2alpha/9.                                         (6)
```

Then `alpha_c=0`, and

```text
y=9s-2alpha=9s_c.                                   (7)
```

Thus the apparently ad hoc centered variable is precisely the translated
constant coefficient of the quartic spectrum.

## 3. Four translation invariants

Define

```text
B
 =beta-7alpha^2/18,

C
 =gamma-5alpha beta/9+35alpha^3/243,

D
 =delta-alpha gamma/3+5alpha^2 beta/54-35alpha^4/1944,

W
 =Psi0+4alpha delta/9-2alpha^2 gamma/27
       +10alpha^3 beta/729-14alpha^5/6561.           (8)
```

The whole reduced mate has the exact identity

```text
Q
 =E_18(P_c)+B E_10(P_c)+C E_6(P_c)+D E_2(P_c).
                                                            (9)
```

This is an equality of polynomials, not merely of Laurent observables.

The first and third fluxes are unchanged. The second flux changes only by
the constant

```text
Psi0-W
 =2alpha(7alpha^4-45alpha^2 beta
          +243alpha gamma-1458delta)/6561,           (10)
```

so in the translated coordinates its constant value is exactly `W`.

Substituting (8) in every retained equation of THM-2262 gives

```text
G(u,y;alpha,beta,gamma,delta,Psi0)
 =G(u,y;0,B,C,D,W),                                  (11)

N_2(alpha,beta,gamma,delta)
 =N_2(0,B,C,D),                                      (12)

N_3(alpha,beta,gamma,delta)
 =N_3(0,B,C,D).                                      (13)
```

Thus the `Z` elimination, trigonal equation, third flux, and Keller one-form
all survive the normalization. The exact companion verifies (4) as whole
polynomial identities for all five Faber seeds and verifies (9)--(13)
independently through the Laurent recurrence.

## 4. The weighted spectral cone

After setting `alpha=0`, the cubic of THM-2262 becomes

```text
G_0
 =-5878656Wy
  -26040609u^3+49601160Bu^2+1607445u^2y^2
  -20995200B^2u-2857680Buy^2-52907904Du-138915uy^4
  +777600B^2y^2+33592320BD-5598720BCy+78120By^4
  +1959552Dy^2-435456Cy^3+1127y^6.                  (14)
```

Assign weights

```text
wt(y,u,B,C,D,W)=(1,2,2,3,4,5).                      (15)
```

Every monomial in (14) has weight six. Equivalently,

```text
G_0(rho^2u,rho y;
    rho^2B,rho^3C,rho^4D,rho^5W)
 =rho^6 G_0(u,y;B,C,D,W).                            (16)
```

Let

```text
Delta_0(y;B,C,D,W)=Disc_u G_0.                       (17)
```

Then

```text
deg_y Delta_0=12,

Delta_0(rho y;rho^2B,rho^3C,rho^4D,rho^5W)
 =rho^12 Delta_0(y;B,C,D,W).                         (18)
```

The leading `y^12` coefficient is the nonzero constant

```text
-153384762202971019112448.
```

Consequently the repeated-branch resultant

```text
Disc_y Delta_0                                      (19)
```

is weighted homogeneous of weight `12*11=132` in `(B,C,D,W)`. The spectral
singular locus is therefore a cone whose projectivization lies in

```text
P(2,3,4,5).                                         (20)
```

Equation (20) is a spectral covariance statement. It is not permission to
discard the scaling of the Keller one-form or polynomial sidecar in a later
argument.

## 5. The one-sparse strata and the C--W plane are empty

We first close all cases in which at most one of `(B,C,D,W)` is nonzero,
then close the whole coordinate plane joining the `C`- and `W`-axes.

### 5.1 The B- and D-axes are squarefree

Weighted covariance lets one test a nonzero coordinate at a unit
representative. On the `B`-axis, exact factorization gives

```text
Delta_0(y;1,0,0,0)
 =-56684737689882624(7y^2+90)

   *(386561y^10+5402250y^8-52787700y^6
     -178605000y^4+7348320000y^2-23619600000).
                                                            (21)
```

Its gcd with its derivative is `1`. On the `D`-axis,

```text
Delta_0(y;0,0,1,0)
 =-19442865027629740032

   *(7889y^12+7429968y^8+5475968064y^4+793437161472),
                                                            (22)
```

and again the gcd with the derivative is `1`.

Thus every nonzero point on either axis has squarefree branch
discriminant. THM-2262's genus-four argument gives a deck contradiction.
No `B`-only or `D`-only survivor exists.

### 5.2 The C-axis normalizes to an elliptic curve

Suppose

```text
B=D=W=0,                   C!=0.                     (23)
```

Put

```text
v=u/y^2,                   z=1/y,

L(v)=1127-138915v+1607445v^2-26040609v^3.           (24)
```

The discriminant of `L` is

```text
Disc L=-153384762202971019112448!=0.                 (25)
```

Equation (14), divided by `y^6`, is exactly

```text
435456 C z^3=L(v).                                  (26)
```

Its projective closure is a smooth plane cubic. Indeed, a singular point
would have `z=0` and would be a repeated projective root of the squarefree
binary cubic homogenizing `L`. Hence its genus is one.

The rational Keller trajectory gives a map from `P^1` to this smooth cubic,
so it is constant. Therefore `z,y,v,u` are constant. Equation (12) then
makes `N_2` constant; the first-flux formula (11) of THM-2262 makes
`Z=T^2` constant. The nonzero square root `q` then becomes constant over the
algebraically closed constant field, contradicting the genuine deck. The
`C`-axis is empty.

For reference, the nonsquarefree affine branch discriminant is

```text
Delta_0(y;0,1,0,0)
 =-136100055193408180224 y^6
   (1127y^6+149688y^3+25509168).                     (27)
```

The factor `y^6` records the singular affine presentation at the excluded
center; it does not make the normalization rational.

### 5.3 The W-axis has genus four

Suppose

```text
B=C=D=0,                   W!=0.                     (28)
```

With the same `v,z`, equation (14) becomes

```text
5878656 W z^5=L(v).                                  (29)
```

The degree-five map of its normalization to the `v`-line is totally
ramified over the three simple roots of `L` and over infinity. The pole order
of `L` at infinity is three, coprime to five, so these are exactly the four
branch points. Riemann--Hurwitz gives

```text
2g-2=5(-2)+4(5-1)=6,

g=4.                                                 (30)
```

Again a rational trajectory is constant, forcing `Z` and then `q` constant,
contrary to the deck. Thus the `W`-axis is empty.

Its branch discriminant,

```text
Delta_0(y;0,0,0,1)
 =-136100055193408180224 y^2
   (1127y^10+2020788y^5+4649045868),                 (31)
```

shows why the raw squarefree test alone missed this genus-four
normalization.

### 5.4 The whole C--W plane has positive-genus normalization

The two axes conceal a stronger statement. Set

```text
B=D=0,                         C W!=0.               (31a)
```

With `v=u/y^2` and `z=1/y`, equation (14) is

```text
L(v)=435456 C z^3+5878656 W z^5=:R(z).              (31b)
```

This cubic in `v` is irreducible over `C(z)`. Indeed, reducibility would give
a rational root `v=f(z)`. A finite pole of `f` cannot occur because `R` is a
polynomial. If `f` has pole order `m` at infinity, then `L(f)` has pole order
`3m`, whereas `R` has pole order five. This is impossible. Thus the
normalization is a connected degree-three cover of the `z`-line.

The two critical values of `L` are the roots of

```text
27 tau^2+68992 tau+226193408,                        (31c)
```

whose discriminant is

```text
-19668992000!=0.                                    (31d)
```

They are distinct and nonzero. On the other side,

```text
R'(z)=z^2(1306368 C+29393280 W z^2).                (31e)
```

The critical point `z=0` maps to zero, hence to neither critical value of
`L`. The other two critical points are simple and have opposite nonzero
critical values. Let `c` be the number of collisions between these two
values and the two critical values of `L`; then `0<=c<=2`.

Without a collision, the two critical values of `L`, each with five simple
preimages under `R`, contribute ten simple ramification points. At a
collision both `L-L(v_0)` and `R-R(z_0)` vanish to exact order two. The
affine curve has an ordinary node there, and its two normalization branches
are unramified over `z`; the collision therefore removes two from the
ramification count. At infinity the pole orders are three and five. Their
coprimality gives one point of ramification index three, contributing two.
Riemann--Hurwitz on the connected degree-three cover gives

```text
2g-2=-6+(10-2c)+2=6-2c,

g=4-c>=2.                                            (31f)
```

A rational Keller trajectory into this normalization is consequently
constant. Thus `z,y,v,u` are constant; covariance (12) makes `N_2` constant,
and the first-flux formula (11) of THM-2262 makes `Z=T^2`
constant. The nonzero square root `q` then becomes constant over the
algebraically closed constant field, contradicting the genuine deck. Sections
5.2 and 5.3 close the two axes, and the next section closes their
intersection. Hence the entire plane `B=D=0` is empty.

### 5.5 The invariant origin is a sidecar-killed cusp

It remains to exclude

```text
B=C=D=W=0.                                           (32)
```

Since `y!=0`, (14) gives

```text
L(v)=0,                         v=u/y^2 in C.         (33)
```

The first flux gives

```text
Z=T^2=k(v)y^3,

k(v)=-2(6561v^2-810v+5)/729.                         (34)
```

The exact resultant

```text
Res_v(L,6561v^2-810v+5)
 =22220417174433988608!=0                            (35)
```

shows that `k(v)!=0`. Absorbing constant squares and cubes, the rational cusp
has a parameter `h in C(x)^*` with

```text
y=h^2,              T=t_0h^3,              t_0!=0,

u=vh^4,             s=h^2/9,                d=(v/t_0)h.
                                                            (36)
```

After substituting (34), the numerator `N_3` of the third flux is

```text
N_3
 =H(v)y^6,

H(v)
 =-4(6561v^2-810v+5)(19683v^2+4374v-125)/243.       (37)
```

The exact resultant

```text
Res_v(L,numerator(H))
 =-1305960654728555144789507775025586141921280
 !=0,                                                (38)
```

so

```text
F=f_0h^9,                         f_0!=0.             (39)
```

Substitute (36) and (39) in the Keller one-form from (1). After cancelling
the nonzero `t_0h^3`,

```text
A(h^9)' in C*.                                      (40)
```

Apply the rational-primitive lemma used in THM-2262. If `A` is constant,
then `h^9` is affine-linear, impossible because a nonconstant affine
polynomial has a simple zero. Otherwise, after translating `X=x-xi`,

```text
A=a_0X^m,                       m>=2,

h^9=c_0+c_1X^(1-m),             c_1!=0.              (41)
```

If `c_0!=0`, the numerator of the right side has nonzero simple roots,
again impossible for a ninth power. Hence for some integer `ell>=1`,

```text
c_0=0,             m=9ell+1,             h=h_0X^-ell.
                                                            (42)
```

The flux quotient alone leaves this cusp. The polynomial Faber sidecar
removes it. At the invariant origin, equation (38) of THM-2262 reads

```text
E_18-R_5
 =-21T/1024[
    -12 mathcal L^3+24 mathcal L^2s
    +mathcal L(12Td-36s^2)
    +T^2-48Tds+48s^3],                               (43)
```

where `mathcal L=Az+E` is polynomial. On (36), the coefficient of `h^6` in
the part independent of `mathcal L` is

```text
J(v)=-2(6561v^2+1134v-19)/729.                       (44)
```

It is nonzero at every root of `L`, because

```text
Res_v(L,6561v^2+1134v-19)
 =422187926314245783552!=0.                          (45)
```

At `X=0`, the term `T J(v)h^6` in (43) has valuation `-9ell`. The terms
containing respectively one, two, or three powers of the polynomial
`mathcal L` have valuation at least

```text
-7ell,              -5ell,              -3ell.      (46)
```

Therefore the `-9ell` pole cannot cancel. But both `E_18=Q` and the
truncation `R_5` are polynomials. This contradiction closes the origin.

## 6. Exact residual after normalization

Every degree-eighteen survivor in the scope of THM-2262 must now satisfy all
of the following:

```text
alpha=0 after the legal translation P -> P+2alpha/9;

at least two of B,C,D,W are nonzero;

at least one of B,D is nonzero;

Disc_y(Disc_u G_0)=0;

the third flux N_3 and Keller one-form remain active;

the whole-polynomial Faber sidecar remains active.                   (47)
```

The first three lines remove a translation orbit, all one-sparse weighted
strata, and the residual `C`--`W` coordinate plane. The spectral covariance
lowers the parameter geometry from five affine constants to a weighted
four-coordinate cone, or a three-dimensional weighted projective search
before imposing the singular hypersurface.

This is not a closure of the remaining two-, three-, or four-coordinate
singular locus, on which `(B,D)!=(0,0)`. In particular, a repeated branch
value does not identify a specific component of the normalization, and the
spectral scaling alone is not a target-preserving quotient of the Keller
one-form.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_target_translation_normal_form_thm2297.py
python3 -O 04-computation/jc2_degree18_target_translation_normal_form_thm2297.py
```

Both runs are byte-identical to the stored output. The companion verifies:

- the whole-polynomial Faber translation bank in degrees `2,6,10,14,18`;
- the invariant coefficients (8) and the exact identity (9);
- covariance of the first flux, shifted second flux, third flux, `N_2`,
  `N_3`, and the complete spectral cubic;
- the weights, the 32-term branch discriminant, and its degree;
- exact axis factorizations and gcd degrees;
- the infinity-cubic discriminant;
- the vanishing raw `C`--`W` branch resultant, the two critical-value
  polynomial, and its nonzero discriminant used in the genus lower bound; and
- all three nonzero resultants (35), (38), and (45).

The smooth-cubic, trigonal collision, superelliptic genus,
rational-primitive, and sidecar valuation arguments are the mathematical
proof above rather than delegated computer conclusions. The theorem remains
scoped to the genuine nonsplit, polynomial exact-square-prefix, reduced
degree-eighteen branch. It does not close other terminal branches or prove
the planar Jacobian conjecture.
