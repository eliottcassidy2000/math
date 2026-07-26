---
id: THM-2376
title: "Degree-eighteen H2 coprime cube-locus closure"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On the
  degree-ten H_2 S_4^2 locus of THM-2357, THM-2371 makes THM-2360's
  coprime Laurent factorization unconditional. Returning that
  factorization to the moving-root coordinate gives a complete
  norm/trace coefficient system. Its y^4 equation has two factors.
  The norm scale excludes the first; the second reduces to a quadratic
  scale equation and one linear pivot in ratio coordinates. The pivot
  chart is empty explicitly, while the main chart ends in coprime
  univariate obstructions of degrees 62 and 53, certified by a
  degree-preserving nonzero resultant modulo 17. Consequently the
  degree-eighteen H_2 S_5^2 stratum is empty. The H_4 stratum, degree
  eighteen as a whole, JC(2), and DC(2) remain open.
source: codex-2026-07-25-degree-eighteen-h2-cube-closure
depends_on:
  - THM-2357-degree-eighteen-h2-moving-root-reduction
  - THM-2360-degree-eighteen-quadratic-ring-cube-descent
  - THM-2371-degree-eighteen-h2-common-root-elimination
related:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
  - THM-2359-degree-eighteen-perfect-quartic-wall-closure
script: 04-computation/jc2_degree18_h2_cube_locus_closure_thm2376.py
output: 05-knowledge/results/jc2_degree18_h2_cube_locus_closure_thm2376.out
script_sha256: 623a6eaa41da742e1de45d95965d6707c44b8545f353b5c79cd1d300f4cc82fa
output_sha256: 18c3328bbd4602d87c6ec28a3df24d2883c64eac6434ccd6421302b6266a90f4
hash_basis: working-tree bytes (LF)
---

# THM-2376 -- close the coprime H2 cube locus

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2357 reduces the mixed degree-eighteen branch to

```text
R_10
 =4(y-1)p_3^3+49q_5^2
 =H_2 S_4^2,                                      (1)

deg(H_2,p_3,q_5,S_4)=(2,3,5,4),

H_2 squarefree.
```

THM-2371 proves throughout this locus that

```text
q_5(1)!=0,                   Res_y(p_3,q_5)!=0.    (2)
```

Thus THM-2360's coprime linear-times-cube factorization applies
unconditionally. This theorem proves that its coefficient locus is
empty. Together with THM-2371, it therefore closes the complete
`H_2 S_5^2` stratum.

## 1. Return the Laurent factorization to the moving-root coordinate

The explicit THM-2357 covariants are

```text
p_3
 =35(y+1)(54B+7y^2+7),                            (3)

q_5
 =7y[
    1620By^2+1620By+2430B+26244C(y+1)
    +77(y^4+y^3+y^2+y)+182
   ].                                              (4)
```

In particular,

```text
lc(p_3)=245,

lc(q_5)=[y^4]q_5=539,

[y^5](7q_5)=[y^4](7q_5)=3773.                    (5)
```

Because `S_4` is monic,

```text
lc(H_2)=L=73060029
             =4*245^3+3773^2.                    (6)
```

Write

```text
H_2=Lh,

h=y^2+ay+b=z^2-Delta,             z=y+a/2.        (7)
```

Squarefreeness of `H_2` says `Delta!=0`. Choose `d` with `d^2=Delta`
and put

```text
s=(z+t)/d,                    t^2=h.              (8)
```

Then `s` has norm one. The selected Laurent root over `y=1` can be
written

```text
R=z(1)+t(1)=d r,

Delta=R(a+2-R).                                   (9)
```

Both `R` and `a+2-R` are nonzero because their product is `Delta`.
Moreover, (2) and

```text
R_10(1)=49q_5(1)^2
```

show that `h(1)!=0`; the two lifts over the moving root are distinct.

In the ring

```text
C[y,t]/(t^2-h)
```

set

```text
ell
 =(z+t)(z+t-R)
 =[2z^2-Delta-Rz]+[2z-R]t.                       (10)
```

This is the polynomial-coordinate form of `Delta*s(s-r)`, and

```text
N(ell)=-2R Delta (y-1).                           (11)
```

THM-2360 gives

```text
7q_5+sqrt(L)S_4 t
 =c_0 s^(-5)(s-r)U_3(s)^3,                       (12)
```

where `deg(U_3)=3` and `U_3(0)!=0`. Indeed, the latter follows from the
nonzero constant of THM-2360's degree-ten polynomial. Multiplying
`s^(-2)U_3` by a nonzero constant, and conjugating if necessary, gives
the lossless normal form

```text
u
 =z^2+mz+n+(-z+w)t.                               (13)
```

The constant adjustment is absorbed into `c`. Equations (10)--(13)
therefore give

```text
7q_5+sqrt(L)S_4 t=c ell u^3.                     (14)
```

Conjugation merely exchanges the two reciprocal orientations; it does
not create another chart.

## 2. Norm coefficients and the trace factor split

Taking norms in (14), using (1) and (11), and cancelling `y-1` in
`C[y]` gives

```text
c^2(-2R Delta)N(u)^3=-4p_3^3.                    (14a)
```

Unique factorization in `C[y]` now forces

```text
N(u)=lambda p_3,                    lambda!=0.     (14b)
```

Indeed, every irreducible valuation on the two sides of (14a) agrees
after division by three, so the two polynomials are associates. This is
the exact structural source of the norm equation; it is not imposed as
an extra ansatz.

Put

```text
k=m+w.                                             (15)
```

The norm of `u` has exact degree three and leading coefficient `2k`.
comparison of `y^3` gives

```text
lambda=2k/245,                         k!=0.       (16)
```

The `y^2` coefficient solves

```text
n
 =[R^2-Ra-2R-3ak+k^2-2km+2k]/2.                  (17)
```

The `y^1` coefficient then solves

```text
B=-7/(216k) [
   -4R^2k+2R^2m+4Rak-2Ram+8Rk-4Rm
   +3a^2k+6akm-4ak-2k^2m+4km^2-4km+4k
 ].                                                (18)
```

The two exact pivots are respectively `2` and `-108k/7`. Thus (17)
divides only by `2`, while (18) divides only by the already recorded
nonzero parameter `k`.

Only the norm constant remains.

Let `Sc` denote the scalar coordinate in the basis `{1,t}`. After
(15)--(17), direct multiplication gives

```text
[y^5]Sc(ell u^3)=2(k^3-R Delta),                  (19)

[y^4]Sc(ell u^3)-[y^5]Sc(ell u^3)
 =-(k^3+R Delta)(R+4a-3k+6m-4).                  (20)
```

By (5), equation (14) first gives

```text
c=3773/[2(k^3-R Delta)].                          (21)
```

Thus

```text
k^3-R Delta!=0.                                   (22)
```

The equality of the two target coefficients in (5), together with
(20), gives the exact split

```text
(k^3+R Delta)(R+4a-3k+6m-4)=0.                   (23)
```

The scalar compatibility in (14a) gives

```text
c^2(-2R Delta)lambda^3=-4.
```

Using (16) and (21), this is

```text
245^3(k^3-R Delta)^2
 =3773^2 R Delta k^3.                             (24)
```

On the first branch of (23), `k^3=-R Delta`. Equation (24) becomes

```text
(4*245^3+3773^2)(R Delta)^2
 =73060029(R Delta)^2=0,
```

contrary to `R Delta!=0`. Hence only

```text
m=(4-R-4a+3k)/6                                  (25)
```

can occur.

## 3. Exhaust the remaining coefficients

After (16)--(18), (21), and (25), retain these four primitive residuals:

```text
N_0:  constant coefficient of N(u)-lambda p_3;

T_3:  y^3 coefficient of c Sc(ell u^3)-7q_5;

T_0:  constant coefficient of c Sc(ell u^3)-7q_5;

T_1:  y coefficient after the y^2 equation solves C.                (26)
```

The coefficient of `C` in the `y^2` equation is the fixed nonzero
constant

```text
-49*26244=-1285956.                               (27)
```

The raw `T_1` numerator contains one factor `-(k^3-R Delta)`, which is
removed using (22). No other factor is discarded. The exact signatures

```text
(total degree, deg_a, deg_R, deg_k, terms)

N_0: (4,3,4,2,24)
T_3: (6,3,6,5,39)
T_0: (8,5,8,3,84)
T_1: (8,5,7,5,78)                                (28)
```

are reconstructed directly from (3)--(14) by the companion.

These equations exhaust all coefficients: the four norm coefficients
give (16)--(18) and `N_0`; the six scalar-trace coefficients give
(21), (23), `T_3`, the `C` pivot, `T_1`, and `T_0`.
Only necessity is used below: every Laurent factorization supplies this
system, so proving the system empty proves the factorization locus empty.
No converse reconstruction is assumed.

## 4. Ratio coordinates expose a linear pivot

Since `R,k!=0`, put

```text
x=k/R,

delta=(a+2-R)/R.                                  (29)
```

Thus

```text
Delta=R^2 delta,

R Delta=R^3 delta.
```

After cancelling only the inherited powers of `R`, equation (24)
becomes

```text
S(delta,x)
 =125delta^2-371delta x^3+125x^6=0.              (30)
```

The four residuals in (26) acquire signatures

```text
(total degree, deg_delta, deg_R, deg_x, terms)

N_0: ( 5,3,1,2,17)
T_3: ( 8,3,2,5,24)
T_0: (13,5,5,3,70)
T_1: (12,5,4,5,56).                               (31)
```

Most importantly, `N_0` is linear in `R`:

```text
R A(delta,x)+B_0(delta,x)=0,                      (32)
```

where

```text
A
 =(delta^2+5delta x+5delta+x)
  (delta x+2delta+2x+1),                          (33)

B_0
 =-2(delta+3x+2)(2delta x+3delta+x).             (34)
```

### 4.1. The pivot chart is empty

If `A=0`, equation (32) also forces `B_0=0`. The primitive exact
resultant is

```text
Res_delta(A,B_0)=x(x+1)^8.                        (35)
```

Since `x!=0`, necessarily `x=-1`. At this value,

```text
A=(delta-1)^2(delta+1),

B_0=-2(delta-1)^2.
```

Thus `delta=1`, but

```text
S(1,-1)=621!=0.                                   (36)
```

So the pivot chart has no point.

### 4.2. The main chart is empty

On `A!=0`, solve

```text
R=-B_0/A.                                         (37)
```

Substitute (37) into `T_3,T_0,T_1`, take primitive numerators, and
reduce in `Q(x)[delta]` modulo (30). Before reduction their signatures
are

```text
(total degree, deg_delta, deg_x, terms)

(12, 7, 9, 52)
(23,16,13,192)
(20,13,13,146).                                   (38)
```

The three primitive remainders have signatures

```text
(24,1,24,38)
(53,1,53,92)
(44,1,44,74),                                     (39)
```

and factor respectively as

```text
x^3 L_21(delta,x),

x^5 L_48(delta,x),

x^5 L_39(delta,x),                                (40)
```

where every `L_j` is linear in `delta`. Because `x!=0`, a common zero
must make the two coefficient determinants

```text
det_delta(x^3L_21,x^5L_48),

det_delta(x^3L_21,x^5L_39)                        (41)
```

vanish. Exact factor extraction gives

```text
det_delta(x^3L_21,x^5L_48)=x^12 D_62(x),

det_delta(x^3L_21,x^5L_39)=x^12 D_53(x),          (42)
```

where the subscripts are the exact degrees. The companion computes

```text
gcd_Q[x](D_62,D_53)=1.                            (43)
```

There is also a compact independent finite-field certificate. Modulo
`17`, both degrees are preserved and

```text
Res_x(D_62,D_53)=11 mod 17.                       (44)
```

A custom `115 x 115` Sylvester determinant and SymPy's independent
resultant implementation agree on (44). Therefore (43) holds, and no
`x!=0` can satisfy (42). The main chart is empty.

## 5. Saturation, hostile controls, and consequence

Every division used above is accounted for:

```text
Delta!=0          squarefree H_2;

R!=0              Delta=R(a+2-R);

k!=0              exact cubic norm degree;

k^3-R Delta!=0    nonzero y^5 trace target;

A!=0              main chart only, with A=0 treated separately.     (45)
```

The remaining denominators are nonzero rational constants. The
coefficient comparison does not divide by an unrecorded discriminant,
root separation, or coefficient variable.

Two hostile controls guard the delicate steps:

1. `A=B_0=0` really has the raw point `(delta,x)=(1,-1)`; it is the norm
   scale, not an accidental resultant cancellation, that removes it by
   the value `621`.
2. The prime `11` is deliberately rejected: the two final obstruction
   degrees drop to `61,52` and acquire a spurious gcd of degree `11`.
   The proof uses the degree-preserving prime `17`.

The modular-resultant routine also checks both a known common-factor
pair and a known coprime pair before certifying (44).

The two branches of (23), the pivot split of (32), and the main
coefficient system are all empty. Therefore the coprime cube locus of
THM-2360 is empty. THM-2371 says every degree-eighteen `H_2 S_5^2`
point lies on that locus. Hence

```text
the degree-eighteen H_2 S_5^2 stratum is empty.    (46)
```

This does **not** eliminate the `H_4 S_4^2` stratum. It does not close
degree eighteen as a whole and proves no instance of `JC(2)` or
`DC(2)`.

## 6. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_h2_cube_locus_closure_thm2376.py
python3 -O 04-computation/jc2_degree18_h2_cube_locus_closure_thm2376.py
```

Both transcripts are byte-identical to

```text
05-knowledge/results/jc2_degree18_h2_cube_locus_closure_thm2376.out
```

The companion reconstructs (3)--(44), checks all polynomial signatures
and factor orders, audits the pivot hostile and bad-prime control, and
compares custom and library modular resultants. No executable check uses
Python `assert`.

Independent audit is pending. QED.
