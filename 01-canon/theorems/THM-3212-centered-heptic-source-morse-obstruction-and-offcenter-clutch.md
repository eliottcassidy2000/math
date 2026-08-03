---
id: THM-3212
title: "Heptic source critical obstruction and complete constant clutch"
status: >
  PROVED + VERIFIED-EXACT.  For every one of the four unmarked S7 covers in
  THM-3123 and every polynomial E_0, the canonical centered source polynomial
  P=V^2 z^4+A_src z+E_0 has five forced, pairwise distinct, reduced Morse
  critical points.  Hence it cannot be one coordinate of a polynomial Keller
  pair.  The obstruction is owned by the squarefree degree-five divisor
  g=ST: both V and A_src vanish on g, while A_src' is a unit there.  In the
  off-center family P=(Vz^2+Bz+C_0)^2+A_src z+E_0, a unit B replaces each
  whole g-fibre by one candidate point, and those five candidates disappear
  exactly when B E_0'-A_src' C_0 is also a unit modulo g.  The control
  B=1,C_0=0,E_0=x passes this local gate but has exactly 52 reduced Morse
  critical points away from g on every cover.  More generally, every constant
  clutch P=(Vz^2+bz+c)^2+A_src z+kx has a critical point: finite-fibre walls
  are explicit, and after they pass the sole S-escape wall removes at most
  three units from a saturated resultant of degree 52 (k!=0) or 28 (k=0).
  Generic members have respectively 52 or 28 reduced Morse points.  Thus the
  entire constant family has no Keller mate and cannot enter the differential-
  owner filtration; this does not prove JC(2).
audit: >
  A self-contained exact companion reconstructs both THM-3123 accessory
  algebras, checks D+A_0=SE^2, the balanced response identities, all
  squarefree/disjointness gates for S,E,T, and the unit resultant of A_src'
  against g.  It then checks the universal centered Hessian calculation and
  the off-center one-jet clutch symbolically.  Ordinary, optimized, and stored
  transcripts are byte-identical; the script has no assert node.  A second
  self-contained quotient-algebra companion independently derives the global
  critical resultant, saturates both passports, and proves the degree-52
  residuals squarefree and disjoint from V by good reductions.  It also
  audits the nearby control; normal/-O/stored transcripts agree exactly.
  A third exact companion derives the universal three-parameter resultant,
  divides all thirteen parameter coefficients by the passport boundary,
  checks both local derivative spectra, proves the S-escape jet PRS has gcd
  one, and verifies degree-52/28 generic and one-root escape controls in both
  accessory fields.  Normal/-O/stored transcripts agree and have no assert.
source: root/frontier-synthesis-cont-2026-08-02
depends_on:
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
related:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2822-sextic-response-centered-lift-mod-three-faber-obstruction
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
  - THM-3172-shear-invariant-differential-owner-filtration-and-transverse-recurrence
script: 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
output: 05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out
script_sha256: 03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448
output_sha256: 729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85
global_script: 04-computation/jc_heptic_offcenter_global_critical_audit_thm3212.py
global_output: 05-knowledge/results/jc_heptic_offcenter_global_critical_audit_thm3212.out
global_script_sha256: d85f861b23ad30337924faff229830bd254774c16d8397f4d77d264fdb840d2e
global_output_sha256: 1a99b53c295d85f90fa4b230fbf1b8e866d97853f0d5b51b43c412e510b9fc60
parameter_script: 04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py
parameter_output: 05-knowledge/results/jc_heptic_constant_parameter_discriminant_audit_thm3212.out
parameter_script_sha256: adeb3c548d5fe3966eefc7ec4eeadfd1410a62356eca8a6c387e39cbe8fc6122
parameter_output_sha256: d170cf2212848ef76722579a40b65506bedf6e65a031012ca06c27dcd1ef77bb
hash_basis: LF-normalized bytes
---

# THM-3212 -- heptic source critical obstruction and complete constant clutch

**PROVED + VERIFIED-EXACT.**

## 1. Scope and inherited response data

Work over an algebraically closed field of characteristic zero.  THM-3123
classifies two heptic pole passports,

```text
(4,1,1,1),                         (3,2,1,1),              (1)
```

and obtains respectively one and three unmarked `S_7` covers.  For either
passport write

```text
D=x^a(x-1)^b(x^2-ux+v),
T=x(x-1)(x^2-ux+v),

E={a(x-1)(x^2-ux+v)+b x(x^2-ux+v)
     +x(x-1)(2x-u)}/7.                                  (2)
```

The accessory algebra, the linear factor `S`, and the nonzero constant `A_0`
are as follows:

```text
(4,1,1,1):
 q=100u^3+244u^2+237u+44,
 v=(8u^2+9u+8)/7,
 S=x+5(u+1)/7,
 A_0=80v^2(u+1)/343;

(3,2,1,1):
 q=75u^3-89u^2-31u+61,
 v=(24u^2-16u-16)/21,
 S=x+(5u-4)/7,
 A_0=9v^2(5u-4)/343.                                  (3)
```

All identities below are taken in the corresponding reduced algebra
`k[u]/(q)`.  THM-3123 proves

```text
D+A_0=SE^2.                                             (4)
```

Put

```text
C=-7A_0,                 M=SET,
F=SE^2/D,                V=4SDT^2/C^2,
G=CE/(2DT),              A_src=VG=2SET/C.               (5)
```

Then

```text
F=VG^2,                  2VG'+V'G=2,
T_F=A_src^2/V=F.                                          (6)
```

Thus `(5)` is the canonical polynomial coefficient suggested by the balanced
response, not an arbitrary choice of source lift.  What THM-3123 did **not**
decide was whether this abstract response enters a two-dimensional polynomial
Keller chart.

## 2. The degree-five owner divisor

Define

```text
g=ST.                                                    (7)
```

For both accessory algebras, `g` is squarefree of degree five and

```text
V/g=4DT/C^2,                     A_src/g=2E/C.            (8)
```

Moreover

```text
gcd(g,A_src')=1.                                         (9)
```

Here is the complete gate behind `(9)`.  The discriminants of `g,T,E` and the
three pairwise resultants of `S,T,E` are units in `k[u]/(q)`.  Hence the five
roots of `ST` are simple, and none is a root of `E`.  At a root `alpha` of
`g`, differentiating `A_src=2gE/C` gives

```text
A_src'(alpha)=2g'(alpha)E(alpha)/C !=0.                  (10)
```

This proves `(9)` and also identifies the missing coordinate: the obstruction
is not merely common support of `V` and `A_src`; it is that their common
degree-five divisor has a transverse `A_src` one-jet.

## 3. Five forced Morse critical points

For an arbitrary polynomial `E_0(x)` consider the centered source family

```text
P_(E_0)(x,z)=V(x)^2 z^4+A_src(x)z+E_0(x).                (11)
```

If `alpha` is any of the five roots of `g`, then `V(alpha)=A_src(alpha)=0`.
Consequently

```text
P_z(alpha,z)=0,
P_x(alpha,z)=A_src'(alpha)z+E_0'(alpha).                 (12)
```

Equation `(10)` makes the second expression affine with nonzero slope.  Thus
there is exactly one critical point over each root of `g`, namely

```text
(alpha, -E_0'(alpha)/A_src'(alpha)).                     (13)
```

The five points are distinct because their `x`-coordinates are distinct.  At
each point,

```text
P_zz=0,                    P_xz=A_src'(alpha),
det Hess(P)=-A_src'(alpha)^2 !=0.                         (14)
```

They are therefore reduced Morse critical points for **every** choice of
`E_0`.  If a polynomial `Q` satisfied `Jac(P_(E_0),Q)=kappa in k*`, the
gradient of `P_(E_0)` could vanish nowhere.  The points `(13)` contradict
that necessary condition.  Hence no member of `(11)` can be a coordinate of
a polynomial Keller pair.

In the Faber notation of THM-3133, this centered lane has

```text
d_F=0,                          s_F=-E_0.                 (15)
```

Allowing arbitrary `E_0` therefore does not repair the obstruction.

## 4. The exact off-center clutch

Now retain the same response coefficients `V,A_src` but allow

```text
P=(Vz^2+Bz+C_0)^2+A_src z+E_0,                           (16)
```

where `B,C_0,E_0 in k[x]`.  At a root `alpha` of `g`, write
`H=Bz+C_0`.  Since `V=A_src=0` there,

```text
P_z(alpha,z)=2B(alpha)(B(alpha)z+C_0(alpha)).             (17)
```

Suppose `B` is a unit modulo `g`.  The unique zero of `(17)` is

```text
z_alpha=-C_0(alpha)/B(alpha).                            (18)
```

At this point `H=0`, so the square in `(16)` contributes nothing to `P_x`
and

```text
P_x(alpha,z_alpha)
 ={B E_0'-A_src'C_0}(alpha)/B(alpha).                    (19)
```

Therefore, inside the unit-`B` lane, all five forced `g`-fibre critical
points disappear **if and only if**

```text
Delta_B=B E_0'-A_src'C_0
```

is a unit modulo `g`.  This is a source, target, map, and sidecar statement:
the centered source `(11)` is the source; `(16)` is the off-center target;
adjoining `B,C_0` is the map; `T_F=A_src^2/V=F` is preserved; the centered
coordinates `d_F,s_F` are changed; and the cheapest decisive sidecar is the
pair of residue classes `(B,Delta_B)` in `k[x]/(g)`.

The smallest local control is

```text
B=1,                    C_0=0,                    E_0=x. (20)
```

Then `Delta_B=1`; at every root of `g`, `P_z=2z` and its only zero `z=0`
has `P_x=1`.  Its Faber coordinates are

```text
d_F=-1/(4V),                     s_F=G/2-x.              (21)
```

Thus `(20)` genuinely removes the five-point local obstruction but moves the
candidate off the centered Faber lane.  It may have critical points away from
`g`, need not satisfy the remaining Faber fluxes, and supplies neither a
polynomial mate nor global chart entry.  This is a precise next test for the
owner filtration of THM-3172, not a proof of `JC(2)` or `DC(2)`.

## 5. The global off-center obstruction

The test in `(20)` is now exact.  On `V!=0`, set `y=Vz`.  Its source becomes

```text
P_tilde=y^2(y+1)^2/V^2+Gy+x.                             (22)
```

After multiplying by units, the gradient ideal is generated by

```text
R_1=4y^3+6y^2+2y+VA_src,
R_2=V^3+V^2y-V'y^2(y+1).                                 (23)
```

Writing `W=VA_src` and `d=V'`, an independent symbolic elimination gives

```text
Res_y(R_1,R_2)
 =64V^9-96V^8+32V^7+48V^6Wd-16V^6W-8V^6d
  -24V^5Wd+8V^5d-4V^4Wd
  +12V^3W^2d^2-4V^3Wd^2+W^3d^3
 =V^3K.                                                   (24)
```

Exact quotient-algebra saturation yields

```text
K=S^3T^8x^9 H_(4111),
K=S^3T^8x^6(x-1)^3 H_(3211),                             (25)
```

respectively.  Both residual polynomials `H` have degree `52` and are
coprime to `V`.  Their successive gcd degrees while saturating by `g` are

```text
(5,5,5,4,4,4,4,4,1,1,1,1,1,1,1,1,1),
(5,5,5,4,4,4,4,4,2,2,2,1,1,1).                          (26)
```

Good reductions `(p,u)=(113,85)` and `(101,64)` preserve degree `52` and
give `gcd(H,H')=gcd(H,V)=1`.  A nonconstant characteristic-zero gcd could be
made monic and could not reduce to a unit at either good place, so these
controls prove squarefreeness and disjointness in the two accessory fields.
The residual coefficient digests are

```text
4111: 693ed1fc7e01eb82f54ec029b5a92fc1a5f6dacf9634368c335717ef3de015bc,
3211: 9728fcc417f5a432c6233f8e2d8813695bf96bfa42554dbd9ff1a441ce2d545a.
                                                               (27)
```

Since `rad(V)=rad(g)`, `(25)--(27)` account for the whole `V!=0` locus.  The
leading y-coefficient of `R_1` is the constant `4`, so elimination creates no
root at infinity.  Each simple resultant root supports one transverse common
zero of `(23)`; transversality of the gradient ideal is precisely a nonzero
Hessian determinant.  Section 4 removes every critical point over `g`.
Therefore `(20)` has **exactly 52 distinct reduced Morse critical points** in
each accessory algebra, hence on all four unmarked `S_7` covers.  It cannot
have a constant-Jacobian mate.

The cheapest adjacent control

```text
B=1,                         C_0=1,                    E_0=x (28)
```

has local clutch `Delta_B=1-A_src'`, a unit modulo `g`.  Nevertheless the
same boundary factor `(25)` remains, and its residual is again squarefree of
degree `52` in both algebras.  The coefficient digests are

```text
b4941e0a8b568d0a0637fb92ad4c5e84547ac839f75930c86c7a0fd77ca1ec2c,
6ab367a3715faae6bcc79fce031885f499c6f7644e8139e8fa9e14494133616b. (29)
```

Thus changing the local clutch did not change the global failure mechanism.

## 6. Lawful Faber and owner-filtration handoff

Neither failed control supplies the input required by THM-3172's `B_1`
owner filtration.  That construction needs a complete marked finite
separable inverse pair (equivalently its inverse minimal polynomial and
marked branch); here only `P` exists, and its critical points exclude every
Keller mate `Q`.  Feeding either control directly into `B_1` would therefore
be a type error.

Formal Faber necessary-condition data remain lawful.  For constant
`B=b,C_0=c,E_0=kx`,

```text
T_F=F,                 d_F=c-b^2/(4V),
s_F=bG/2-kx,           T_Fd_F=cF-(b^2/4)G^2.             (30)
```

Since `A_src=2M/C`, the remaining flux identity
`A_src K_Q=lambda M` forces

```text
K_Q=lambda C/2.                                             (31)
```

For the normalized nonsplit polynomial exact-square-prefix `R>=7` bank this
is the formal condition
`F H_Q=lambda C/2`, together with `phi_Q=0` and constant `psi_Q`.  THM-3133
already excludes this response at the common simple `S`-zero of `V` and `M`.
This is a necessary-condition elimination only; it creates no owner or chart.

## 7. Full constant-clutch parameter obstruction

Work over an algebraically closed field `K_0` of characteristic zero.  For
either cubic accessory algebra `K_i` of THM-3123 and every embedding
`K_i -> K_0`, retain the response polynomials `V,A=A_src` and put `g=ST`.
These base changes include all four unmarked `S_7` covers.  For arbitrary
**constants** `b,c,k in K_0`, consider

```text
P_(b,c,k)(x,z)=(V(x)z^2+bz+c)^2+A(x)z+kx.                 (C1)
```

Then `P_(b,c,k)` has a critical point for every triple `(b,c,k)`.  Therefore
there is no polynomial `Q` with

```text
Jac(P_(b,c,k),Q) in K_0^*.                                  (C2)
```

More precisely:

1. On a nonempty Zariski-open subset with `k!=0`, `P_(b,c,k)` has exactly
   `52` distinct reduced Morse critical points, all away from `g`.
2. On a nonempty Zariski-open subset of the divisor `k=0`, it has exactly
   `28` distinct reduced Morse critical points, all away from `g`.
3. When the finite `g`-fibre clutch passes but the unique `S`-escape wall is
   hit, at most three units of saturated resultant degree can escape through
   `S`.  At least `49` units remain away from `g` if `k!=0`, and at least
   `25` remain if `k=0`.

The last two lower bounds count resultant/intersection multiplicity, not
necessarily distinct critical points.  They are used only to force existence.

### 7.1 Universal critical resultant

On `V!=0` put

```text
y=Vz,                      L=y^2+by+cV,                  (C3)
```

and write `d=V'`.  Using `2VA'-AV'=2V`, the gradient equations are equivalent,
up to units, to

```text
R_1=2L(2y+b)+VA,
R_2=kV^3+V^2y-bd yL.                                     (C4)
```

Exact symbolic elimination gives

```text
Res_y(R_1,R_2)=V^3 K_(b,c,k),                            (C5)
```

where

```text
K = 64k^3V^6+(-96bk^2+64ck)V^5+32b^2ck^2V^4d
    +(32b^2k-32bc)V^4+48bk^2V^4Ad-16V^4A
    +(-8b^4k^2-32b^3ck+32b^2c^2)V^3d
    +(-24b^2k+16bc)V^3Ad+(8b^5k-8b^4c)V^2d
    +16b^3ckV^2Ad^2-4b^3V^2Ad+12b^2kV^2A^2d^2
    -4b^5kVAd^2+2b^4cA^2d^3+b^3A^3d^3.                 (C6)
```

All thirteen parameter-monomial coefficients in `(C6)` are exactly divisible
by the passport-dependent boundary factor

```text
B_4111=S^3T^8x^9,
B_3211=S^3T^8x^6(x-1)^3.                                (C7)
```

Thus

```text
K_(b,c,k)=B_i H_(i;b,c,k),         H_i in K_i[b,c,k][x]. (C8)
```

The divisibility follows intrinsically from the response identity and the
orders of `V,A,V'`; the exact companion also divides every one of the thirteen
coefficients independently.  The local order table is

```text
passport     root       ord(V)       ord(B_i)
4111         S             1             3
4111         x             6            17
4111         x-1,Q_2       3             8
3211         S             1             3
3211         x             5            14
3211         x-1           4            11
3211         Q_2           3             8.             (C9a)
```

At a `T`-root of `V`-order `m>=3`, every coefficient in `(C6)` has order at
least `3m-1`.  At `S`, the apparent order-two terms cancel by
`2VA'-AV'=2V`, leaving order at least three.  These are precisely the orders
in the last column.  In both passports

```text
deg(V)=16,                deg(A)=8,                deg(B_i)=44. (C9)
```

If `k!=0`, the unique degree-96 term is `64k^3V^6`, so

```text
deg_x H_i=52.                                             (C10)
```

If `k=0`, the unique top term is `-16V^4A`, of degree `72`, so

```text
deg_x H_i=28.                                             (C11)
```

These degree statements are uniform on their stated strata, not generic
interpolations.

### 7.2 Complete finite-fibre clutch atlas

At a root `alpha` of `g`, let `m=ord_alpha(V)`.  Since `A` has a simple zero,
the response identity gives

```text
A'(alpha)=2/(2-m).                                       (C12)
```

The multiplicity lists of `V` therefore give the exact spectra

```text
(4,1,1,1):  2, -1/2, -2, -2, -2,
(3,2,1,1):  2, -2/3, -1, -2, -2.                        (C13)
```

If `b=0`, then at every root `alpha` of `g`,

```text
P_z(alpha,z)=0,
P_x(alpha,z)=2cV'(alpha)z^2+A'(alpha)z+k.                (C14)
```

The second polynomial is nonconstant because `A'(alpha)!=0`, so it has a root
over the algebraically closed ground field.  Thus `b=0` is always obstructed.

Suppose `b!=0`.  The unique zero of `P_z(alpha,z)` is `z=-c/b`, and

```text
P_x(alpha,-c/b)=(bk-cA'(alpha))/b.                       (C15)
```

Consequently the finite `g`-fibre clutch passes exactly when the following
product is nonzero:

```text
Lambda_4111=(bk-2c)(bk+c/2)(bk+2c)^3,
Lambda_3211=(bk-2c)(bk+2c/3)(bk+c)(bk+2c)^2.             (C16)
```

Every zero of `(C16)` is an actual finite critical point, not an elimination
artifact.

### 7.3 The only post-clutch boundary escape

Assume `b Lambda_i !=0`.  At any root of `T`, put

```text
V=v t^m+...,             A=a t+...,             a=2/(2-m),  m>=3. (C17)
```

The coefficient of order `3m-1` in `K` is

```text
16m(m-1)/(m-2) b^4 v^3 (bk-ca).                         (C18)
```

It is nonzero by `(C16)`.  Thus `H_i` meets no root of `T` after the finite
clutch passes.

Let `s` be the root of `S`, put `t=x-s`, and expand

```text
V=v_1t+v_2t^2+v_3t^3+v_4t^4+v_5t^5+...,
v_1!=0.                                                   (C19)
```

The response identity determines the jet of `A`; in particular

```text
A=2t+(2v_2/(3v_1))t^2
   +4(3v_1v_3-v_2^2)/(15v_1^2)t^3+....                   (C20)
```

Direct substitution in `(C6)` gives

```text
[t^3]K=-(8/3)b^2v_1^2(bk-2c)
        (4b^2v_2+3(bk+2c)v_1^2).                        (C21)
```

Define

```text
eta_i=-2v_2/(3v_1^2)=-V''(s)/(3V'(s)^2),
tau_i=bk+2c-2eta_i b^2.                                  (C22)
```

In the two accessory fields,

```text
eta_4111=-34937(100u^2+244u+137)/131072,
eta_3211=49(1031225u^2+6942913u-5694514)/3359232.        (C23)
```

Equations `(C18)--(C23)` give the exact boundary resultant, up to a nonzero
accessory-field constant `gamma_i`:

```text
Res_x(g,H_i)=gamma_i b^18 Lambda_i tau_i.                (C24)
```

Thus `tau_i=0` is the only new boundary wall after the finite clutch passes.
It represents a root escaping through the `S` chart with `z` unbounded, not a
finite critical point over `S`.

### 7.4 The escape multiplicity is at most three

Assume `b Lambda_i!=0` and `tau_i=0`.  Put

```text
r=c/b^2,                     h=b^3,                       (C25)
```

so `k=b(2eta_i-2r)`.  Let

```text
L_0=3rv_1^2+v_2,
A_2(r)=30r^2v_1^4+10rv_1^2v_2+18v_1v_3-16v_2^2.         (C26)
```

The `S`-clutch condition `bk-2c!=0` is exactly `L_0!=0`.  The next coefficient
is

```text
[t^4]K=64b^3/(45v_1) L_0 (hA_2(r)-15v_1^3).             (C27)
```

If `(C27)` vanishes, then `A_2(r)!=0` and

```text
h=15v_1^3/A_2(r).                                        (C28)
```

After this substitution, the next two coefficients are

```text
[t^5]K=-(640/7)v_1^4 L_0 J_3(r)/A_2(r)^2,
[t^6]K=-(128/63)v_1^3 J_5(r)/A_2(r)^2,                  (C29)
```

where

```text
J_3=315r^3v_1^6+105r^2v_1^4v_2
    +r(315v_1^3v_3-280v_1^2v_2^2)
    -90v_1^2v_4+219v_1v_2v_3-128v_2^3,                 (C30)

J_5=113400r^5v_1^10+354375r^4v_1^8v_2
    +r^3(59535v_1^7v_3+148680v_1^6v_2^2)
    +r^2(12150v_1^6v_4+256770v_1^5v_2v_3-204165v_1^4v_2^3)
    +r(-15750v_1^5v_5-29700v_1^4v_2v_4+11529v_1^4v_3^2
       +202554v_1^3v_2^2v_3-157474v_1^2v_2^4)
    -5250v_1^3v_2v_5-7290v_1^3v_3v_4-4770v_1^2v_2^2v_4
    +13077v_1^2v_2v_3^2+30177v_1v_2^3v_3-25712v_2^5.   (C31)
```

Exact Euclidean division in each accessory field gives

```text
(deg J_5,deg J_3)=(5,3)->(3,2)->(2,1)->(1,0),
gcd(J_3,J_5)=1.                                          (C32)
```

Therefore `(C27)--(C29)` cannot all vanish.  Since `B_i` contains exactly
`S^3`,

```text
ord_S(H_i)<=3.                                           (C33)
```

This is the decisive hostile boundary: a constant clutch can move at most
three units of resultant multiplicity to infinity through `S`.

### 7.5 Uniform critical-point obstruction

There are now three exhaustive cases.

1. `b=0`: equation `(C14)` gives a finite `g`-fibre critical point.
2. `b!=0` and `Lambda_i=0`: equation `(C15)` gives a finite `g`-fibre critical
   point.
3. `b Lambda_i!=0`: `H_i` meets no root of `T`; it either misses `S`, or has
   `S`-multiplicity at most three by `(C33)`.  After full saturation it has
   degree at least `52-3=49` when `k!=0` and at least `28-3=25` when `k=0`.
   Hence it has a root with `V!=0`.  The leading `y`-coefficient of `R_1` is
   the constant `4`, so that root supports an affine common zero of `(C4)` and
   therefore a critical point of `P`.

This proves the theorem for both accessory fields and all four covers.
At any critical point of `P`, `Jac(P,Q)` vanishes for every polynomial `Q`,
so `(C2)` is impossible.

### 7.6 Generic Morse counts and hostile controls

Let

```text
Delta_i(b,c,k)=Disc_x(H_i).                              (C34)
```

It is not the zero polynomial: the control `(b,c,k)=(1,0,1)` is the
squarefree degree-52 control already audited in THM-3212.  Consequently

```text
k b Lambda_i tau_i Delta_i !=0                            (C35)
```

is a nonempty Zariski-open set with exactly `52` distinct reduced Morse
critical points away from `g`.  For constants common to all four covers, take
the intersection of these opens, equivalently the nonvanishing of the
accessory-field norms of `tau_i Delta_i` together with the rational factors in
`(C35)`.

On `k=0`, the exact control `(1,1,0)` gives a squarefree degree-28 `H_i`
coprime to `g` in both accessory fields.  Thus `Disc_x(H_i|_(k=0))` is also
nonzero, and a nonempty open subset of `k=0` has exactly `28` reduced Morse
points.

Two deliberately hostile escape controls show why boundary tuning cannot
repair the source:

```text
(b,c,k)=(1,eta_i,0):       H_i squarefree degree 28,
                           gcd(H_i,g)=S, hence 27 Morse points away;

(b,c,k)=(1,0,2eta_i):      H_i squarefree degree 52,
                           gcd(H_i,g)=S, hence 51 Morse points away.       (C36)
```

Good reductions are `(p,u,eta)=(113,85,101)` and `(101,64,43)`.  They preserve
all displayed degrees and give derivative gcd zero in all four controls.
The discriminant locus `Delta_i=0` merely collides or degenerates critical
points; it does not remove the existence conclusion in Section 7.5.

### 7.7 Lawful Faber and differential-owner handoff

For this constant family the formal Faber data remain

```text
T_F=F,
d_F=c-b^2/(4V),
s_F=bG/2-kx,
T_Fd_F=cF-(b^2/4)G^2.                                   (C37)
```

THM-3133 separately excludes the response in every normalized nonsplit
polynomial exact-square-prefix `R>=7` Faber bank at the common simple `S`-zero
of `V` and `M`.  The present theorem is a stronger source-level pre-owner
screen for constant `b,c,k`; it does not
construct a mate or prove `JC(2)`.

Because every source `(C1)` has a critical point, none supplies the complete
marked finite separable inverse pair required by THM-3172.  Feeding any of
them into the `B_1` differential-owner filtration would remain a type error.
The lawful next search must leave the constant-clutch family—for example by a
nonconstant `B,C_0,E_0` deformation that changes `(C4)--(C6)`—while retaining
all Faber fluxes.  A point on `Delta_i=0` or `tau_i=0` is not a lawful owner
candidate merely because the generic Morse count drops.


## Exact reproduction

Run

```text
python3 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
python3 -O 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
python3 04-computation/jc_heptic_offcenter_global_critical_audit_thm3212.py
python3 -O 04-computation/jc_heptic_offcenter_global_critical_audit_thm3212.py
python3 04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py
python3 -O 04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py
```

All three mode pairs must reproduce their declared outputs byte for byte.  The
companions use exact rational, polynomial, resultant, gcd, and accessory-
algebra arithmetic only.  They import no repository executable and use
explicit `require` gates, so optimized execution retains every check.  QED.
