---
id: THM-2872
title: "Four-slot shared multipole quartic norm and response-secant reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In adjacent-block coordinates, every candidate
  for the last four-slot factorial branch is exactly the common zero of two
  cubic remainder septics and one nonnegative quartic remainder norm.  The
  norm has a chart-free Hermitian normalization and an Euler-transverse
  normalization whose denominator contains one explicit gap conic.  On a
  cone-cutting cubic line above the bottom boundary, the quartic endpoint
  condition is equality of two weighted response secants; after that equality
  only one midpoint/Jensen defect remains.  The TP3/secant sign,
  cone-avoiding chamber, bottom boundary, and unbounded Euler avoidance remain
  open.
source: root/shared-multipole-quartic-norm-2026-07-28
audit: >
  audit-shared-secant-2026-07-28 (independent rederivation of both
  remainders, resultant/norm factors, Hermitian and Euler normalizations,
  all response-secant factors, the midpoint iff, and every scope boundary:
  ACCEPT)
depends_on:
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
  - THM-2860-euler-tangent-mobius-chord-and-factorial-cubic-avoidance-box
  - THM-2866-positive-factorial-difference-semiring-and-cubic-pascal-response-ladder
related:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
script: 04-computation/gmc_shared_multipole_quartic_norm_secant_thm2872.py
output: 05-knowledge/results/gmc_shared_multipole_quartic_norm_secant_thm2872.out
script_sha256: 3c4cf51eaa4c820a1d271e71c2f1c5db080d9b6c080209338fc700c5faecea89
output_sha256: 3fc436c21f5f62ef2626086594247de4b79dc1a7a4129d76dd41ab192fe4da13
hash_basis: LF-normalized bytes
---

# THM-2872 -- shared multipole quartic norm and response secants

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2866 removes the `F^o=0` alternative from THM-2848.  The sole
remaining first-window four-slot branch is therefore a real plane which is
simultaneously a cubic and quartic Maxwell multipole line.  This theorem
puts that branch into three exact forms:

```text
algebraic:       two cubic remainders and one positive quartic norm;
Euler:           one explicit gap conic separates tangent/transverse charts;
positive-cone:   one endpoint secant holonomy and one midpoint defect.
```

This is a strict reduction, not an exclusion of the branch.

## 1. Signed adjacent coordinates

Fix

```text
0<=a_0<a_1<a_2<a_3,             f_a=s^a/a!,
L(s^m)=m!,

B_i=f_(a_i)-f_(a_(i-1)),        i=1,2,3.               (1)
```

Each `B_i` is a nonzero positive adjacent-difference cone, and every
mean-zero polynomial on the four slots is uniquely

```text
Z=z_1B_1+z_2B_2+z_3B_3.                               (2)
```

If

```text
T_m(i_1,...,i_m)=L(B_(i_1)...B_(i_m)),
```

THM-2841 gives

```text
T_m(i_1,...,i_m)>0,                    m>=2.           (3)
```

The remaining branch is nonetheless signed: it asks for a nonzero complex
`z` satisfying

```text
P_m(z)=sum T_m(i_1,...,i_m)z_(i_1)...z_(i_m)=0,
                                             m=2,3,4. (4)
```

The quadratic `P_2` is positive definite over `R`.  Hence a nonzero
`P_2`-null vector is not a complex multiple of a real vector.  Its real
plane has projective normal

```text
n=Re(z) cross Im(z)
 =(
   Im(conj(z_2)z_3),
   Im(conj(z_3)z_1),
   Im(conj(z_1)z_2)
  ).                                                   (5)
```

Multiplication of `z` by a nonzero complex scalar multiplies `n` by its
squared modulus; conjugation reverses `n`.  Thus `[n]` is the intrinsic
unoriented real plane coordinate.

## 2. Two cubic equations and one quartic norm

Let

```text
E_n={x:n_1x_1+n_2x_2+n_3x_3=0}.                       (6)
```

On the chart `n_3!=0`, take the real basis

```text
U=n_3B_1-n_1B_3,               V=n_3B_2-n_2B_3.       (7)
```

It is independent because the first two coordinates are respectively
`n_3` and `n_3`.  Put

```text
q(z)=g_0+2g_1z+g_2z^2,
c(z)=t_0+3t_1z+3t_2z^2+t_3z^3,
f(z)=A_0+4A_1z+6A_2z^2+4A_3z^3+A_4z^4,               (8)
```

where these are the quadratic, cubic, and quartic moments of `U+zV`.
Real factorial positivity gives

```text
g_0>0,             g_2>0,
D_g=g_0g_2-g_1^2>0.                                   (9)
```

Exact Euclidean division by `q` gives

```text
c(z) mod q=(C_0+C_1z)/g_2^2,                          (10)

C_0=2g_0g_1t_3-3g_0g_2t_2+g_2^2t_0,

C_1=-g_0g_2t_3+4g_1^2t_3-6g_1g_2t_2+3g_2^2t_1,
```

and

```text
f(z) mod q=(R_0+R_1z)/g_2^3,                          (11)

R_0=A_0g_2^3-6A_2g_0g_2^2+8A_3g_0g_1g_2
    +A_4g_0^2g_2-4A_4g_0g_1^2,

R_1=4[
 A_1g_2^3-3A_2g_1g_2^2-A_3g_0g_2^2+4A_3g_1^2g_2
 +A_4g_0g_1g_2-2A_4g_1^3].
```

Define

```text
N_4=g_2R_0^2-2g_1R_0R_1+g_0R_1^2.                    (12)
```

Its exact square decomposition is

```text
N_4=[(g_2R_0-g_1R_1)^2+D_gR_1^2]/g_2.                 (13)
```

Equations `(9)` and `(13)` prove

```text
N_4>=0,
N_4=0 iff R_0=R_1=0 iff q divides f.                  (14)
```

Consequently the last four-slot branch is exactly

```text
C_0=C_1=N_4=0                                         (15)
```

in any valid cyclic chart.  If `n_k!=0`, the other two indices `i,j`
give the basis

```text
U=n_kB_i-n_iB_k,                  V=n_kB_j-n_jB_k.     (16)
```

The three charts cover `P^2(R)`.  Their zero predicates agree on overlaps
because each is precisely the intrinsic divisibility predicate.
The ordinary binary resultant records the same norm with the exact
normalization

```text
Res_z(q,f)=N_4/g_2^3.                                  (16a)
```

## 3. Intrinsic and Euler normalizations

Let `zeta` be the upper-half-plane root of `q`.  From

```text
zeta+conj(zeta)=-2g_1/g_2,
zeta conj(zeta)=g_0/g_2,
```

equations `(11)--(12)` give

```text
N_4=g_2^7 |L((U+zeta V)^4)|^2,                        (17)

L((U+zeta V)(U+conj(zeta)V))=2D_g/g_2.                (18)
```

Thus

```text
Omega(E)
 =|L((U+zeta V)^4)|^2
   /L((U+zeta V)(U+conj(zeta)V))^4
 =N_4/(16g_2^3D_g^4).                                 (19)
```

This is unchanged by chart, real basis, complex scale, orientation, or
conjugation.  Under `n -> lambda n`,

```text
g:lambda^2,  t:lambda^3,  A:lambda^4,
C_i:lambda^7, R_i:lambda^10, N_4:lambda^22,
```

and the denominator of `(19)` also has degree `22`.

Now put

```text
p=a_1-a_0,       qgap=a_2-a_1,       r=a_3-a_2.       (20)
```

If `(b_0,b_1,b_2,b_3)` is the plane normal in THM-2860's four-slot
coordinates, then

```text
n_i=b_i-b_(i-1).                                      (21)
```

The Euler-tangent Möbius conic becomes

```text
T_A(n)
 =p r n_2(n_1+n_2+n_3)
  -qgap(p+qgap+r)n_1n_3.                              (22)
```

Direct four-by-four expansion in the basis `(7)` gives

```text
det[U,V,AU,AV]=-n_3^2T_A(n).                          (23)
```

For `H=U+zeta V`, THM-2860's signed determinant is therefore

```text
Delta_A(E;H)=4D_g n_3^2T_A(n)/g_2^2.                  (24)
```

On the Euler-transverse locus,

```text
|L(H^4)|^2/Delta_A(E;H)^2
 =N_4/[16g_2^3D_g^2n_3^4T_A(n)^2].                   (25)
```

The denominator in `(25)` is not globally available: for consecutive
gaps and `n=(1,1,1)`, `T_A(n)=0`.  THM-2860 excludes this tangent conic
for cubic multipole lines only in the finite box `a_3<=30`.

## 4. The cone-cutting response-secant form

Assume `[n]` has mixed signs, so `E_n` cuts the nonnegative adjacent cone
in a two-ray wedge.  Choose its nonzero positive extreme rays `U,V`.
Assume also

```text
U(0)=V(0)=0.                                          (26)
```

Condition `(26)` holds uniformly when `a_0>=1`.  The `a_0=0` boundary
requires the correction terms already exposed in THM-2866 and is not
included below.

Here and below `U'=dU/ds` and `V'=dV/ds` are ordinary polynomial
derivatives.  Because `a_0>=1`, differentiation lowers every occupied
adjacent ray `d_j` to `d_(j-1)`.  Thus `U'` and `V'` remain nonzero
positive adjacent cones.

For a nonzero positive cone `W` and another positive cone `Y`, define

```text
R_W(Y)=L(YW^2)/L(YW),
S_W(Y)=L(YW^3)/L(YW).                                  (27)
```

Put

```text
alpha=t_0/g_0=R_U(U),              beta=t_3/g_2=R_V(V). (28)
```

Integration by parts, `(26)`, and `C_0=C_1=0` give

```text
delta_(2,U):=R_U(V)-R_U(U')
            =t_3g_0/(3g_1g_2)>0,

delta_(2,V):=R_V(U)-R_V(V')
            =t_0g_2/(3g_1g_0)>0.                     (29)
```

Every denominator and numerator in `(29)` is positive by THM-2841.
Thus cubic divisibility itself forces two strict quadratic-response
transport reversals.

Define

```text
delta_(3,U)=S_U(V)-S_U(U'),
delta_(3,V)=S_V(U)-S_V(V'),

kappa_U=delta_(3,U)/delta_(2,U),
kappa_V=delta_(3,V)/delta_(2,V).                      (30)
```

Write a possible quartic quotient as

```text
f=q(r_0+2r_1z+r_2z^2).                                (31)
```

Its left and right endpoint determinations are

```text
r_1^(L)=(2A_1g_0-A_0g_1)/g_0^2,
r_1^(R)=(2A_3g_2-A_4g_1)/g_2^2.                      (32)
```

Exact substitution of `(29)--(30)` gives

```text
beta kappa_U =(3/2)r_1^(L),
alpha kappa_V=(3/2)r_1^(R).                           (33)
```

Hence the two endpoint determinations agree exactly when

```text
beta kappa_U=alpha kappa_V.                           (34)
```

This is not a new invariant parallel to THM-2866's top-face
determinant.  It is exactly that determinant in response coordinates:

```text
J:=(2A_1g_0-A_0g_1)g_2^2-(2A_3g_2-A_4g_1)g_0^2
  =g_0^2g_2^2(r_1^(L)-r_1^(R))
  =(2/3)g_0^2g_2^2(beta kappa_U-alpha kappa_V).       (34a)
```

When `(34)` holds,

```text
r_1=(2/3)beta kappa_U=(2/3)alpha kappa_V.             (35)
```

The only coefficient not fixed by the two endpoints is the middle one.
Evaluating `(31)` at `z=1` proves that this coefficient agrees exactly
when

```text
S_(U+V)(U+V)-S_U(U)-S_V(V)
 -(4/3)beta kappa_U=0.                                (36)
```

Therefore, under `(26)`, a cone-cutting cubic line is shared with a
quartic line exactly when both `(34)` and `(36)` hold.

THM-2830 and THM-2866 prove that the atomic response coordinates

```text
n -> R_W(d_n),                 n -> S_W(d_n)
```

are separately strictly increasing.  What `(34)` needs is a stronger
comparison of their secant slopes.  A TP3 response-curve theorem is the
natural missing input.  It cannot be replaced by the top-face stochastic
order used in THM-2866: that order would make
`R_V(V')>=R_V(U)`, while cubicity forces the strict reverse inequality in
`(29)`.

## 5. Sharp chart and positivity boundaries

Four controls prevent overextension.

1. **Chart boundary.**  At `n=(1,1,0)`, formula `(7)` gives
   `U=V=-B_3` and is invalid.  The cyclic `n_1` chart instead gives the
   independent basis `B_2-B_1,B_3`.
2. **Complex isotropy.**  Real `V` in a valid chart has `g_2>0`, but this
   does not survive careless complexification.  On `{0,1,2,3}`,

   ```text
   V=B_1+(-1+i)B_2/2 !=0,             L(V^2)=0.       (37)
   ```

3. **Cone boundary.**  If all `n_i` have the same strict sign, the plane
   misses the nonnegative cone except at zero; equivalently the four-slot
   normal values are strictly monotone.  Zero coordinates give boundary
   faces.  Neither case receives the positive extreme-ray proof in
   Section 4.
4. **Euler boundary.**  `T_A=0` leaves `(12)` valid but makes `(25)`
   undefined.

Three algebraic degeneracies explain why those hypotheses cannot be
suppressed.  If `D_g=0`, take `q=(1+z)^2` and the nonzero remainder
`1+z`: then `N_4=0` although divisibility fails.  If `g_1=0`, as for
`q=1+z^2`, the response secants in `(29)--(30)` are undefined; positive
extreme rays exclude this because THM-2841 gives `g_1=L(UV)>0`.  Finally,
for `U=d_0`,

```text
L(U'U)=0,
```

so the bottom-boundary derivative response is genuinely unavailable.

There is also no stability or Lorentzian shortcut.  On `{0,1,2,3}`, put

```text
D(u,v,w)=primitive part of
 det[L(B_i(uB_1+vB_2+wB_3)^m)]_(m=1,2,3;i=1,2,3).     (38)
```

All `28` coefficients of `D` are positive, but

```text
D(1,v,1)
 =217v^6+4802v^5+43835v^4+210800v^3
  +560815v^2+779518v+441213                           (39)
```

has discriminant

```text
-2^21 3^6 5^10 7^2 2281 7508739397<0.                (40)
```

Thus `(39)` is not real-rooted and `D` is not stable.  The Hessian of
`partial_u partial_w^3 D` has determinant

```text
-43366523385600000
```

and positive trace `38634480`, so it has two positive eigenvalues.  This
also rules out the required Lorentzian signature.

## 6. Exact verification and scope

Run

```text
python3 04-computation/gmc_shared_multipole_quartic_norm_secant_thm2872.py
python3 -O 04-computation/gmc_shared_multipole_quartic_norm_secant_thm2872.py
```

Both modes byte-match the stored transcript.  The companion verifies:

1. the two Euclidean remainders `(10)--(11)`;
2. the positive norm, resultant, root norm, intrinsic normalization, and
   gauge degree;
3. the Euler gap conic and both determinant normalizations;
4. both transport reversals, endpoint secants, and midpoint defect;
5. the chart, complex-isotropy, degenerate-Gram, zero-cross, bottom-`d_0`,
   and Euler hostiles; and
6. the exact `28`-term stability/Lorentzian hostile.

The theorem leaves open:

```text
the general TP3/mixed-secant inequality
  (THM-2879 closes only the THM-2855 shifted high-tooth family);
the sign of the midpoint defect if endpoint holonomy vanishes;
the cone-avoiding and cone-boundary chambers;
the a_0=0 integration-by-parts correction;
unbounded Euler-tangent avoidance.
```

It does not close the shared-line branch, four-slot SFC, GMC(2), the
Gaussian Moment Conjecture, or the planar Jacobian conjecture.

## 7. Independent hostile audit

The independent audit rederived both Euclidean remainders and the exact
resultant factor, completed the square in `(13)`, and recomputed the
Hermitian root norm and `Omega`.  It independently expanded the Euler
determinant in all three pivot charts and checked the factors in
`(23)--(25)`.  The companion additionally checks exact coefficient-matrix
rank in the repaired cyclic chart, rather than merely comparing its two
basis polynomials.

For the positive chamber it rederived every integration-by-parts response,
the two strict reversals, both endpoint determinations, and the midpoint
iff.  It specifically required `D_g>0`, `g_1>0`, `(26)`, a valid pivot
chart, and a cone-cutting plane, and supplied the degenerate-Gram,
zero-cross, `d_0`, chart-boundary, and cone-avoiding hostiles above.  It
found no remaining factor, sign, or scope defect.

**QED.**
