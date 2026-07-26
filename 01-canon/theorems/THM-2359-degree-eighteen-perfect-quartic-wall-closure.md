---
id: THM-2359
title: "Degree-eighteen perfect-quartic wall closure"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On the
  degree-eighteen wall 504D=115B^2, the structured quartic is the
  perfect square 5(27B+7y^2)^2. Every point whose Mordell polynomial
  has square-class degree at most four would have gcd(F,F') of degree
  at least four. In the B=1 chart, exact terminal subresultants show
  that the C=0 edge is impossible and that the C!=0 multiplicity locus
  is the single orbit C^2=-175/8748, W/C=-682/735. At that orbit F has
  one sixth-order root but its residual sextic is squarefree and
  coprime to that root, so the squarefree part has degree six. Hence the
  complete wall contains no Keller trajectory. Other parameter walls,
  the residual off-wall H_2/H_4 strata, and JC(2) remain open.
source: codex-2026-07-25-perfect-quartic-wall
depends_on:
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2347-degree-eighteen-double-zero-wall-saturation
  - THM-2357-degree-eighteen-h2-moving-root-reduction
script: 04-computation/jc2_degree18_perfect_quartic_wall_thm2359.py
output: 05-knowledge/results/jc2_degree18_perfect_quartic_wall_thm2359.out
script_sha256: 49f2bdce283d08f6c38b654445684ca047abd754a93ffd18f1cddb7affeb4153
output_sha256: 3d20b5fbdcd3088c16c0ea4e66c11f03b5283d14d291353e7542b900caeca575
hash_basis: working-tree bytes (LF)
---

# THM-2359 -- the perfect-quartic wall has no low square class

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2332 reduces every residual degree-eighteen Keller trajectory to

```text
F=4P^3+49Q^2=H S^2,                  deg(H) in {0,2,4}.       (1)
```

THM-2345 and THM-2347 close the two coordinate walls detected by the
values of `P` and `Q` at the origin.  The next intrinsic divisor is
detected by the discriminant of `P` as a quadratic in `y^2`.  On this
wall `P` itself becomes a perfect square.  The resulting factorization of
`F` exposes a single spectacular sixth-order collision, but that collision
is entirely an index square: the residual squarefree carrier still has
degree six.  This distinction closes the wall.

## 1. The perfect-quartic wall

Retain THM-2332's exact covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                            (2)
```

The discriminant of `P` as a quadratic in `y^2` is

```text
238140(115B^2-504D).                                (3)
```

Thus on

```text
504D=115B^2                                         (4)
```

one has the exact identity

```text
P=5(27B+7y^2)^2.                                    (5)
```

If `B=0`, equation (4) also gives `D=0`.  The complete plane `B=D=0`
is already empty by THM-2297.  Hence assume `B!=0` and use weighted
scaling to normalize

```text
B=1,                         D=115/504.              (6)
```

Write

```text
L=27+7y^2.
```

Then

```text
P=5L^2,

Q
 =539y^6+11340y^4+183708Cy^3-10935y^2
  +(2361960C+2480058W)y.                             (7)
```

Over `Q(sqrt(-5))`, equation (1)'s Mordell polynomial has the structural
factorization

```text
F
 =(7Q+10sqrt(-5)L^3)(7Q-10sqrt(-5)L^3).             (8)
```

This factorization is a useful view, not by itself a square-class
conclusion.  The common-root sidecar is

```text
Q mod L
 =(39366/7)(y(294C+441W)+32),

Res_y(L,Q)
 =3720786376356
  (333396C^2+1000188CW+750141W^2+1024).             (9)
```

The final quadratic in (9) is equivalently

```text
(27/7)(294C+441W)^2+1024.                           (10)
```

## 2. Low square class forces a four-dimensional gcd

The leading coefficient of `F` is the nonzero constant

```text
73060029,

deg_y(F)=12.                                        (11)
```

If (1) holds, then

```text
deg(S)=(12-deg(H))/2>=4,

S divides gcd(F,F').                                (12)
```

Consequently

```text
deg gcd(F,F')>=4.                                   (13)
```

This necessary condition remains valid at every higher collision.  Since
the leading coefficients of `F` and `F'` are constant and nonzero, the
subresultant specialization theorem has no hidden degree-drop chart.  The
complete generic subresultant degrees are

```text
12,11,10,9,8,7,6,5,4,3,2,1,0.                      (14)
```

Write

```text
Sres_3=a_3y^3+a_2y^2+a_1y+a_0,

Sres_2=b_2y^2+b_1y+b_0.                              (15)
```

Condition (13) forces every coefficient in (15) to vanish.  Only

```text
a_3=a_0=b_2=0                                       (16)
```

will be needed.

## 3. The C=0 edge is empty

First put `C=0`.  The primitive specializations of `a_3` and `a_2` are
polynomials in `W` with signatures

```text
(degree,number of terms)=(14,8), (13,7).             (17)
```

Their exact polynomial gcd over `Q[W]` is `1`.  Hence they cannot vanish
simultaneously over `C`, contradicting (16).  Thus every point satisfying
(13) has

```text
C!=0.                                               (18)
```

## 4. The ratio chart has one nilpotent orbit

On (18), introduce the involution-invariant ratio coordinates

```text
X=C^2,                      Z=W/C,                   (19)

U=8748X+175,                V=735Z+682.              (20)
```

Substitute `W=ZC` into the three coefficients in (16).  Their exact
`C`-orders are respectively

```text
0,1,0.                                               (21)
```

Divide the middle coefficient by `C`, replace every `C^(2j)` by `X^j`,
make the shift (20), and discard only nonzero rational contents.  The
three resulting generators in `Q[U,V]` have signatures

```text
(total degree, terms, U-degree, V-degree)

=(21,151,11,14), (18,109,9,12), (27,219,12,18).      (22)
```

Their exact grevlex basis in variable order `(U,V)` is

```text
G_0
 =-212914386392781422592U^2
  -776250367057015603200UV
  +1873417064878541875V^4
  +19823498739475200000V^3
  -707519865807175680000V^2,

G_1
 =845652816485908992U^3
  +25755182917786828800U^2
  +93899104387764480000UV
  +2724675804767984375V^3
  +85585121186764500000V^2,

G_2
 =15660237342331648U^2V
  -8718368619626496U^2
  -31785718925721600UV
  -51227634911866875V^3
  -28971358395840000V^2,

G_3
 =-8477454513299521536U^2
  +3083109226771543200UV^2
  -30907386246404505600UV
  +6409548769616444375V^3
  -28170794755837440000V^2.                         (23)
```

Exact division by (23) gives

```text
U^5=0 mod (G_0,G_1,G_2,G_3),

V^5=0 mod (G_0,G_1,G_2,G_3).                        (24)
```

Thus every complex common zero of the original three generators has

```text
U=V=0.                                              (25)
```

Equations (19)--(20) identify the complete multiplicity locus:

```text
C^2=-175/8748,               W/C=-682/735.           (26)
```

The two signs in (26) are exchanged by the weighted scaling `lambda=-1`,
so they form one weighted-projective orbit.

## 5. The unique orbit has squarefree degree six

Put `r^2=-21` and choose the representative

```text
C=5r/162,                       W=-341r/11907.        (27)
```

It satisfies (26).  Exact factorization in `Q(r)[y]` gives

```text
F=73060029(y+3r/7)^6 R_6(y),                         (28)
```

where

```text
R_6
 =y^6-(18r/7)y^5-(8721/161)y^4+(1356r/49)y^3
  +(203472/1127)y^2-(243000r/7889)y-364500/7889.
                                                            (29)
```

The residual checks are

```text
Disc(R_6)
 =1242517308060047668936704000000000000
  /8551083049093133387426609
 !=0,

R_6(-3r/7)=-1093500/343!=0.                         (30)
```

Therefore `R_6` is squarefree and coprime to the displayed repeated
root.  In particular,

```text
gcd(F,F')=(y+3r/7)^5

up to a nonzero scalar,                                  (31)
```

but the sixth power in (28) is itself a square.  The squarefree part of
`F` is exactly the degree-six polynomial `R_6`, not a polynomial of degree
at most four.  Thus the sole point surviving the coarse multiplicity test
(13) fails the actual square-class condition (1).

## 6. Close the complete wall

Let a degree-eighteen Keller trajectory satisfy (4).  The `B=0` chart is
closed by THM-2297.  On `B!=0`, normalize as in (6).  THM-2332 gives
(1), hence the necessary multiplicity (13).  Sections 3--4 force the
unique orbit (26), while Section 5 proves that orbit has squarefree degree
six, contradicting (1).

Consequently

```text
no degree-eighteen Keller trajectory satisfies 504D=115B^2.           (32)
```

This closes a third codimension-one divisor in the higher-support
degree-eighteen cone.  It does not close the residual off-wall `H_2` or
`H_4` strata, any other parameter wall, `JC(2)`, or `DC(2)`.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_perfect_quartic_wall_thm2359.py
python3 -O 04-computation/jc2_degree18_perfect_quartic_wall_thm2359.py
```

The ordinary transcript is byte-identical to the stored output.  The
optimized path contains no Python `assert` and executes the same checks,
but its final replay was deliberately stopped when the shared disk fell
below the session safety margin; byte comparison under `-O` remains an
independent-audit item.  The companion checks (3)--(11), the complete
subresultant profile (14), the `C=0` edge gcd, every chart division in
(19)--(22), every basis element in (23), both nilpotence reductions in
(24), and the exact orbit factorization, discriminant, and coprimality
checks (27)--(31).  Two off-orbit hostile controls have gcd degree zero.

Independent audit is pending. QED.
