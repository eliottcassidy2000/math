---
id: THM-2891
title: "Zero-normal-coordinate consecutive cone-boundary cubic emptiness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  integer n>=0, every real two-plane
  in span{d_n,d_(n+1),d_(n+2)} whose normal has a zero coordinate and
  whose intersection with the nonnegative adjacent-difference cone is
  two-dimensional is cubic-empty: its binary quadratic factorial moment
  does not divide its binary cubic factorial moment.  Two boundary
  families have a coefficientwise-negative division invariant.  The sole
  interlaced family has a nonzero exact resultant for every integer
  depth; its only sign-indefinite factor is negative at n=0 and becomes
  coefficientwise positive after n=m+1.  Together with THM-2879/2890,
  this closes the full closed consecutive cone-cutting plane atlas for
  n>=1 against shared cubic--quartic lines.  Cone-avoiding planes and the
  general mixed four-slot branch remain open.
source: codex/gmc-boundary-newton-2026-07-29
audit: root/gmc-boundary-audit-2026-07-29
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
  - THM-2890-discrete-newton-closure-of-strict-consecutive-gmc-wedges
related:
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2844-sharp-signed-adjacent-ray-orientation-boundary
  - THM-2869-gamma-positive-cone-y-zero-exit-portal
script: 04-computation/gmc_consecutive_zero_normal_boundary_thm2891.py
output: 05-knowledge/results/gmc_consecutive_zero_normal_boundary_thm2891.out
script_sha256: 42d247e0aa6f3e83a8907087a3e6d6c462b45c5a537af0500eb3d702a254b872
output_sha256: b852df3e6faac0218eb913f64227331f0e16e3ecc7f3779ffc1855b51a094285
hash_basis: LF-normalized bytes
---

# THM-2891 -- zero-normal-coordinate consecutive cone-boundary cubic emptiness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Put

```text
L(s^q)=q!,                  f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

Fix an integer `n>=0` and write

```text
e_1=d_n,              e_2=d_(n+1),          e_3=d_(n+2).  (2)
```

Let `P` be a real two-plane through the origin in the coefficient space
of `(2)`.  Suppose:

1. a normal vector of `P` has a zero coordinate; and
2. `P` meets the nonnegative `(e_1,e_2,e_3)` cone in a
   two-dimensional cone.

For any real basis `U,V` of `P`, define

```text
Q(alpha,beta)=L((alpha U+beta V)^2),
C(alpha,beta)=L((alpha U+beta V)^3).                   (3)
```

Then

```text
Q does not divide C.                                   (4)
```

Equivalently, every zero-normal-coordinate boundary of the consecutive
positive-cone atlas is cubic-empty.

## 1. The complete boundary atlas

Let the plane equation be

```text
r_1 c_1+r_2 c_2+r_3 c_3=0.                            (5)
```

If exactly one `r_i` is zero, the other two coefficients must have
opposite signs for the intersection in the hypothesis to be
two-dimensional.  Permuting the zero coordinate and rescaling `(5)`
gives exactly the three families

```text
A_t=span{e_1,             e_2+t e_3},
B_t=span{e_1+t e_3,       e_2},
C_t=span{t e_1+e_2,       e_3},             t>0.      (6)
```

If two normal coordinates vanish, `P` is one of the three coordinate
faces

```text
span{e_1,e_2},       span{e_1,e_3},       span{e_2,e_3}. (7)
```

Thus `(6)--(7)` are the whole projective boundary atlas.  In particular,
the apparent six chart boundaries of the high, middle, and low affine
charts pair up into three geometric families.

The `A` and `C` families, and the coordinate faces, also follow from the
transport-ordered positive-cone theorem THM-2830.  The genuinely
interlaced family is `B_t`.  It is the `Y=0` specialization

```text
U=e_1+X e_3,                    V=e_2                 (8)
```

of the high chart in THM-2879.  THM-2812/2824 do not cover this general
four-slot plane, and THM-2844 treats a different fixed-bottom signed ray.

## 2. Exact divisibility invariants

For every nonzero real polynomial `p`,

```text
L(p^2)=integral_0^infinity p(s)^2 e^(-s) ds>0.        (8)
```

Since `U,V` are independent real polynomials, `(8)` makes `Q` positive
definite over `R`.  The real binary quadratic `Q` is therefore
irreducible.  As in THM-2824, if `Q` and the real binary cubic `C` share
one complex projective root, conjugation supplies the other root of `Q`,
so `Q` divides `C`.  The converse is immediate.  Thus

```text
Q divides C.                                           (9)
```

Write

```text
g_0=L(U^2),       g_1=L(UV),       g_2=L(V^2),

t_0=L(U^3),       t_1=L(U^2V),
t_2=L(UV^2),      t_3=L(V^3).                         (10)
```

The two division-free remainders are

```text
I_1=3t_1g_0g_2-t_3g_0^2-2t_0g_1g_2,
I_2=3t_2g_0g_2-2t_3g_1g_0-t_0g_2^2.                  (11)
```

Thus `(9)` is equivalent to

```text
I_1=I_2=0.                                            (12)
```

The normalized adjacent tensors of THM-2879 divide every order-`k`
tensor by the same positive factor.  All three terms of each invariant
in `(11)` acquire the same positive scale.  Hence the normalization
preserves `(12)` for every `n>=0`.

## 3. Two signed families and all coordinate faces

For the `A_t` family, exact expansion gives

```text
I_1(A_t)
 =-(n+1)^2(n+2)^2 A(n,t)/D(n),                        (13)

D(n)=72(2n+1)^2(2n+3)^2
        (3n+1)(3n+2)(3n+4)(3n+5)>0,                  (14)
```

where the coefficient rows of `A`, in descending powers
`n^4,...,n^0` and columns `t^3,...,t^0`, are

```text
(1700,1380, 348, 28)
(6997,5893,1607,143)
(9769,8473,2475,243)
(5282,4678,1434,152)
( 912, 816, 256, 28).                                 (15)
```

Every entry is positive.  Therefore

```text
I_1(A_t)<0                         for n>=0,t>=0.      (16)
```

The endpoint `t=0` is the face `span{e_1,e_2}`.

For the `C_t` family,

```text
I_1(C_t)=-(n+2)^2 C(n,t)/D(n),                        (17)
```

where the coefficient rows of `C`, in descending powers
`n^6,...,n^0` and columns `t^4,...,t^0`, are

```text
( 1700, 11680, 29088, 31104, 12096)
(10397, 65016,150456,156384, 61776)
(25463,142420,303060,304128,121488)
(31817,156456,304512,293304,118068)
(21245, 90628,162000,149232, 60396)
( 7106, 26208, 43476, 38280, 15552)
(  912,  2952,  4608,  3888,  1584).                 (18)
```

Hence

```text
I_1(C_t)<0                         for n>=0,t>=0.      (19)
```

The endpoint `t=0` is the face `span{e_2,e_3}`.
The remaining coordinate face has

```text
I_1(e_1,e_3)
 =-(n+1)^2(n+2)^2
   (1700n^4+6997n^3+9769n^2+5282n+912)/D(n)<0.       (20)
```

Equations `(16),(19),(20)` close `A_t`, `C_t`, and every face already at
one cubic remainder.

## 4. The interlaced boundary resultant

For `B_t`, clear the positive denominators in `(11)` and call the two
numerators

```text
F_4(n,t),                       G_3(n,t),              (21)
```

where the subscripts are their `t`-degrees.  Their exact profiles are

```text
                         deg_t       terms in (n,t)
F_4                         4                45
G_3                         3                32.       (22)
```

Exact elimination gives

```text
Res_t(F_4,G_3)
 =1728 (n+1)^2(n+2)^17(2n+1)^4
        (3n+1)^3(3n+2)^2 P_10(n)P_11(n),             (23)
```

where

```text
P_10(n)=
 250000n^10+7215000n^9+54955825n^8+203222633n^7
 +436912597n^6+587900003n^5+508206858n^4
 +280435412n^3+94789320n^2+17740224n+1398528         (24)
```

is coefficientwise positive, and

```text
P_11(n)=
 1536n^11+111488n^10+2387360n^9+13553352n^8
 +33678576n^7+40665048n^6+18881565n^5
 -8228003n^4-14821853n^3-7575441n^2
 -1768896n-159876.                                    (25)
```

The last factor records a genuine continuous-depth sign change.  Exact
Sturm counting gives one root in

```text
(7/10,71/100).                                        (26)
```

Thus ordinary coefficient positivity in `n` is unavailable.  Integer
depth has a sharp two-piece certificate:

```text
P_11(0)=-159876!=0,                                   (27)
```

while, on writing `n=m+1`, the descending coefficient list is

```text
(1536,
 128384,
 3586720,
 40309992,
 241935792,
 880569288,
 2058721629,
 3148367062,
 3106290857,
 1875172704,
 610895628,
 76724856).                                           (28)
```

Every entry in `(28)` is positive.  Therefore

```text
P_11(n)>0                         for every integer n>=1. (29)
```

All other factors in `(23)` are nonzero at `n=0` and positive for
`n>=1`.  Hence the resultant never vanishes at an integer `n>=0`.
The two invariants cannot have a common complex `t` root.  In particular,

```text
I_1(B_t),I_2(B_t) are not both zero
                                    for n>=0,t>0.      (30)
```

This is stronger than positivity-domain emptiness: at each integer depth,
the two specialized division invariants have no common complex `t` root.

## 5. Completion of the closed consecutive cone atlas

Sections 1, 3, and 4 prove `(4)` on every zero-normal-coordinate
boundary.  Since a shared cubic--quartic multipole line first requires
`Q|C`, every such boundary misses the shared-line locus already at the
cubic stage.

For `n>=1`, THM-2879 and THM-2890 treat the three strict mixed-normal
charts:

1. the strict shared-high chart has a cubic branch but exits at quartic
   endpoint holonomy;
2. the strict shared-middle chart is cubic-empty; and
3. the strict shared-low chart is cubic-empty by discrete Newton
   positivity.

The present theorem supplies their common projective boundary.  Thus

```text
every closed consecutive cone-cutting plane for n>=1
misses the shared cubic--quartic multipole-line locus. (31)
```

At `n=0`, the present boundary theorem remains valid, but `(31)` is not
claimed because THM-2890's strict low chart starts at `n=1`.

## 6. Mechanism, failure boundary, and scope

The closest inherited mechanisms separate cleanly:

- THM-2830 already closes the transport-ordered `A/C` portions;
- THM-2879 supplies the exact divisibility coordinates;
- THM-2890 shows why the integer depth, rather than a continuous
  Maclaurin parameter, is the natural positivity variable.

The Gamma boundary portal of THM-2869 is the hostile control.  The same
interlaced `Y=0` shape acquires a genuine positive cubic-null point when
the Gamma shape is allowed to vary continuously to approximately
`23.2446`.  Therefore boundary emptiness is a factorial/integer-depth
fact, not a formal cone or total-positivity tautology.

This theorem does **not** treat:

- cone-avoiding planes, whose normal has one strict sign;
- nonconsecutive adjacent-difference triples;
- arbitrary four-slot supports or the general mixed shared-line branch;
- quartic midpoint exit outside the closed consecutive atlas; or
- any new instance of GMC(2), DvdK removal, `SFC(4)`, `JC(2)`, or
  `DC(2)`.

## 7. Exact verification

Run

```text
python3 04-computation/gmc_consecutive_zero_normal_boundary_thm2891.py
python3 -O 04-computation/gmc_consecutive_zero_normal_boundary_thm2891.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_consecutive_zero_normal_boundary_thm2891.out.
```

The companion pins the exact THM-2879 dependency hash and checks:

1. all eight positive clearing denominators;
2. every coefficient and the canonical digests of `(13)--(20)`;
3. the complete resultant factorization `(23)` and both factor digests;
4. the exact integer shift `(27)--(29)` and the continuous hostile
   interval `(26)`;
5. reconstruction of the resultant as the actual consequence object;
   and
6. independent literal-factorial resultants at `n=0,1,8`.

The coefficient-dictionary digests are frozen in the matching output.

## 8. Independent hostile audit

An independent root audit rederived the projective boundary
classification, the real-irreducible quadratic/conjugate-root step, both
division invariants, and the positive normalization scale.  It checked
the exact five linear resultant factors and their multiplicities, the
`P_10/P_11` split, the `n=0` exception, the `n=m+1` positivity
certificate, and the continuous hostile root.  It independently replayed
normal and optimized modes and added literal-factorial controls at
`n=2,3,5,13`.  No mathematical or scope defect remained after the script
was hardened to pin the exact linear factors.

**QED.**
