---
id: THM-3062
title: "Four-slot physical transverse augmentation and affine tail holotopy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For two fixed
  lower offsets, one affine moving lower offset, and a fixed terminal gap,
  the rank-two all-large physical face cancels but its first transverse
  augmentation is explicit and positive.  The actual normalized intrinsic
  resultant is asymptotic to
  G2^12(4^h-3^h)^24 U3,C^8 U4,C^6, hence to a positive constant times
  12^(24C)C^(-17).  After the exact four-slot width flag is restored, every
  fixed finite generalized-Hankel bank is eventually positive uniformly
  along the straight carrier-to-physical homotopy and every affine clock.
  This is not a selected-chart persistence or one all-order Stieltjes tail.
source: kind-pasteur-2026-08-01-four-slot-physical-augmentation
audit: >
  Three independent immutable-file audits rederived the high-layer bases,
  length-twelve localization with multiplicity, Poisson normalization and
  orientation, full Gamma asymptotic, exact carrier and corrected width flag,
  affine alternant, and uniform finite-bank holotopy.  All ordinary and
  optimized replays matched each other, the stored companion, and the
  declared LF hashes.  The final audit also supplied the all-width
  nonnegative norm boundary and caught the repaired THM-2843 slug and
  THM-2942 evidence dependency.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
  - THM-3060-three-slot-physical-terminal-face-and-affine-tail-holotopy
related:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-3058-k4-hafnian-initial-face-augmentation-and-unbounded-cancellation-jet
script: 04-computation/gmc_four_slot_physical_transverse_tail_thm3062.py
output: 05-knowledge/results/gmc_four_slot_physical_transverse_tail_thm3062.out
script_sha256: 5b22fc39fa132127800c39718aec9a67b538f31ba30e87bdf493281dd9d5b836
output_sha256: 27a28c55c2f61f91286acdc626815dd47192dbf9ae5e9cb83d83b1976177535b
hash_basis: LF-normalized bytes
---

# THM-3062 -- the four-slot physical face survives one transverse level later

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3060 finds the actual terminal face in three slots and proves that its
two powered lines have nonzero resultant.  It also identifies the exact
four-slot failure: the first all-large forms depend on only two normal
coordinates and share a transverse projective point.  The correct next move
is not to discard the physical specialization.  Keep the cubic and quartic
normal forms, and evaluate the fixed-low quadratic on that transverse point.
The resulting augmentation is positive and dominates every nominally larger
term after the equal-base cancellations are contracted.

## 1. Physical four-slot forms

Let `L(s^j)=j!` and `f_j=s^j/j!`.  Fix integers

```text
n>=1, 0<=a<b<C, h>=1, M=C+h.                          (1)
```

For homogeneous coordinates `(x_0,x_1,z)`, put

```text
Y_C=x_0(f_(n+a)-f_(n+M))
    +x_1(f_(n+b)-f_(n+M))
    +z   (f_(n+C)-f_(n+M)),                           (2)

F_(r,C)=L(Y_C^r)/L(f_n^r),             r=2,3,4,
R_C=Res(F_(2,C),F_(3,C),F_(4,C)).                     (3)
```

These are exactly THM-2925's normalized physical inclusion forms at the
integer width `M`; the resultant in `(3)` is the intrinsic standard
multivariate resultant.  No `2^M=3^M=4^M=0` corner has been taken.

Write

```text
w=x_0+x_1,                    P=[1:-1:0].              (4)
```

Split `F_(r,C)` according to the number `q` of factors drawn from the high
pair `{f_(n+C),f_(n+M)}`:

```text
F_(r,C)=sum_(q=0)^r F_(r,q,C).                        (5)
```

Every high factor is a linear form in `(w,z)`, hence

```text
F_(r,q,C) belongs to (w,z)^q.                         (6)
```

At the transverse point all high layers vanish and the terminal factor
cancels exactly.  Therefore

```text
F_(r,C)(P)=G_r
 :=L((f_(n+a)-f_(n+b))^r)/L(f_n^r),                  (7)

G_2>0.                                                (8)
```

The strict sign is the positive-definite Laguerre Gram norm of a nonzero
difference.

## 2. The two normal lines

For `r=3,4`, define the all-high scale

```text
U_(r,C)=L(f_(n+C)^r)/L(f_n^r)
       =(rn+1)_(rC)/(n+1)_C^r.                        (9)
```

The all-high blocks have complete coefficientwise inverse-`C` expansions,
whose leading forms are

```text
F_(r,r,C)/U_(r,C) -> L_r(w,z)^r,
L_r(w,z)=-r^h w+(1-r^h)z.                            (10)
```

The two lines are independent:

```text
det(L_3,L_4)=4^h-3^h>0.                              (11)
```

Put

```text
g_(3,C)=F_(3,C)/U_(3,C),
g_(4,C)=F_(4,C)/U_(4,C).                              (12)
```

The lower-layer exponential ratios, up to fixed powers of `C`, are

```text
r=3: q=2 -> (4/27)^C,      q=1,0 -> (1/27)^C;
r=4: q=3 -> (27/256)^C,    q=2 -> (1/64)^C,
     q=1,0 -> (1/256)^C.                              (13)
```

## 3. All twelve intersections concentrate at P

The limiting pair in `(10)` has only the common projective point `P`.
Compactness first forces all `3*4=12` points of

```text
V(g_(3,C),g_(4,C))                                   (14)
```

into one fixed affine chart around `P`, counted with intersection
multiplicity, for every sufficiently large `C`.

There is also an effective rate.  On that chart put
`r_P=||(w,z)||`.  Independence in `(11)` makes the all-high pair uniformly
directionally coercive.  Combining that coercivity with `(6)` and `(13)`
gives, at every common zero, one of the two inequalities

```text
c r_P^3 <=poly(C)[(4/27)^C r_P^2+27^(-C)r_P+27^(-C)],
c r_P^4 <=poly(C)[(27/256)^C r_P^3+64^(-C)r_P^2
                   +256^(-C)r_P+256^(-C)].           (15)
```

Their elementary root bounds imply uniformly

```text
r_P<=poly(C)3^(-C).                                   (16)
```

The quadratic now loses every high contribution.  Its fixed layer is smooth
at `P`, its one-high layer is `poly(C)r_P`, and its two-high layer is
`poly(C)4^C r_P^2`.  Hence at each point of `(14)`, with multiplicity,

```text
F_(2,C)=G_2+O(poly(C)(4/9)^C).                        (17)
```

This is the transverse mechanism: the large quadratic normal face vanishes
fast enough on the cubic--quartic intersection that the fixed positive Gram
value reappears.

## 4. Poisson product and the first nonzero augmentation

Let

```text
E_C=Res_bin(g_(3,C)|_(x_0=0),g_(4,C)|_(x_0=0)).      (18)
```

By `(10)--(11)`,

```text
E_C -> Res_bin(L_3^3,L_4^4)=(4^h-3^h)^12.            (19)
```

For large `C`, `E_C!=0`, so the homogeneous Poisson formula on the
`x_0!=0` chart is exact:

```text
R_C=U_(3,C)^8 U_(4,C)^6 E_C^2
    product_(p in V(g_3,g_4)) F_(2,C)(p).             (20)
```

Here `8=2*4` and `6=2*3` are the resultant multidegrees.  Equations
`(17),(19),(20)` prove

```text
R_C=G_2^12(4^h-3^h)^24 U_(3,C)^8 U_(4,C)^6
    (1+O(C^-1)).                                      (21)
```

Every object in `(18)--(20)` is a finite Gamma-ratio/resultant expression,
so `(21)` has a complete inverse-power expansion.  The root-product error
in `(17)` is exponentially smaller.  In particular

```text
R_C>0 for every sufficiently large C.                 (22)
```

There is also a sharp all-width boundary.  The quadratic is positive
definite on the real mean-zero plane, so the six points of its intersection
with the cubic occur in conjugate pairs.  Poisson evaluation of the quartic
therefore gives

```text
R_C>=0 for every admissible C,
R_C=0 iff F_(2,C),F_(3,C),F_(4,C) have a common point. (22a)
```

Thus eventual positivity in `(22)` is stronger than an orientation choice,
while strict positivity at every width is exactly the still-open arbitrary
four-slot first-window nullity problem.  Nothing here assumes that closure.

Equation `(20)` also proves the tropical optimality which termwise base
sorting misses.  Once all equal-base contributions are contracted, the
normalized resultant has the finite nonzero limit in `(21)`.  Therefore no
larger exponential base survives.  The first physical augmentation is
exactly the fixed quadratic plus the two highest powered lines; it is not an
arbitrary choice among the many equal-base Macaulay terms.

Using

```text
U_(r,C)=r^(rC) product_(j=1)^(r-1)
         [(n+j/r)_C/(n+1)_C],                         (23)
```

and THM-3060's constant `kappa_r`, `(21)` becomes

```text
R_C=A_(n,a,b,h) 12^(24C) C^(-17)(1+O(C^-1)),         (24)

A_(n,a,b,h)=G_2^12(4^h-3^h)^24 kappa_3^8 kappa_4^6
             >0.                                     (25)
```

The power is `8*(-1)+6*(-3/2)=-17`.

## 5. Exact carrier and affine tail holotopy

Define the positive carrier

```text
B_C=G_2^12(4^h-3^h)^24 U_(3,C)^8 U_(4,C)^6.          (26)
```

Formula `(23)` makes `(B_C)` a scaled product of Beta moments: eight copies
of each `Beta(n+j/3,1-j/3)`, `j=1,2`, and six copies of each
`Beta(n+j/4,1-j/4)`, `j=1,2,3`.

The exact four-slot width flag is

```text
W_M=(n)_M^26(n+1)_M^20/n^(46M).                      (27)
```

Thus

```text
Q_C=W_(C+h)R_C,                 G_C=W_(C+h)B_C       (28)
```

have the same positive leading symbol, while `(G_C)` is an exact strict
Stieltjes product with full support.

Fix integers `s>=1,b_0>=0`, put `C_c=sc+b_0`, and retain only tail indices
for which `a<b<C_c`.  For

```text
H_c^(theta)=(1-theta)G_(C_c)+theta Q_(C_c),
0<=theta<=1,                                           (29)
```

and every fixed order bound `d` and offset bound `B_0`, there is an effective
`c_0` such that every generalized minor

```text
det[H_(c+u_i+v_j)^(theta)]_(i,j=1)^r>0                (30)
```

whenever `c>=c_0`, `1<=r<=d`, and `u_i,v_j` are strictly increasing in
`[0,B_0]`.  Indeed `(21)` gives a full symbol for the actual/carrier ratio,
uniformly along `(29)`.  The common width flag contributes net Gamma order
`46s`, so the first confluent alternant is

```text
(46s)^N_r V(u)V(v)/product_(ell=0)^(r-1)ell!>0.       (31)
```

The finite bank makes the remainder uniform.  This is a literal finite-order
holotopy from an exact Beta--Gamma law to the actual physical resultant.

## 6. Exact evidence and scope

The exact companion performs:

- `48` normal-ideal layer and `12` transverse restriction checks;
- eight powered-line and `64` Gauss/Beta factorization checks;
- `33` exact physical intrinsic resultants from one selected `36x36`
  Macaulay chart after division by its proved universal flag;
- exact ratio intervals at `C=40,60,80`, descending toward the limit one;
- `39` carrier and `195` straight-homotopy Hankel controls.

Run

```text
python 04-computation/gmc_four_slot_physical_transverse_tail_thm3062.py
python -O 04-computation/gmc_four_slot_physical_transverse_tail_thm3062.py
```

Both modes equal the stored ten-line transcript after LF normalization.
An independent replay obtained script hash
`5b22fc39fa132127800c39718aec9a67b538f31ba30e87bdf493281dd9d5b836`
from `12,716` bytes in `392` lines and output hash
`27a28c55c2f61f91286acdc626815dd47192dbf9ae5e9cb83d83b1976177535b`
from `548` bytes in `10` lines.  Ordinary, optimized, and stored outputs were
byte-identical.  THM-2942 supplies the universal flag divided out by the
selected-chart checks; those checks do not assert chart persistence.

This theorem proves an intrinsic normalized physical-resultant asymptotic for
two fixed lower offsets, an affine moving offset, and a fixed terminal gap.
It does not prove that the displayed selected Macaulay chart persists at all
widths, that a wall-stripped norm core has the same sign, or that one tail is
Stieltjes at every order.  It does not handle moving fixed lows, varying gap,
or nonlinear clocks, and proves no new NC2, GMC, SFC, LRC, or Jacobian
statement.  THM-3058's DVR matching contraction is only a related filtered-
face analogy; it is not a dependency or a physical map used here.

**QED.**
