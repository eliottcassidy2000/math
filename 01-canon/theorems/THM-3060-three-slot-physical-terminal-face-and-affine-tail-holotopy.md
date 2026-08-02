---
id: THM-3060
title: "Three-slot physical terminal face and affine tail holotopy"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENTLY HOSTILE-AUDITED.  For the actual
  normalized three-slot inclusion system with
  support (0,C,C+h), the physical all-large face is a pair of distinct
  powered lines.  Three-slot factorial detection makes its intrinsic
  resultant positive at every width; it is asymptotic to a positive
  constant times 46656^C C^(-7/2).  After multiplication by the exact width
  flag, every fixed finite generalized-Hankel bank is eventually positive,
  uniformly along an explicit Beta-Gamma-carrier homotopy and every affine
  clock.  From four slots onward this first physical face has rank at most two
  and its multivariate resultant vanishes, isolating the next cancellation
  jet rather than claiming a higher-slot result.
source: kind-pasteur-2026-08-01-three-slot-physical-tail
audit: >
  Two independent audits reconstructed the physical inclusion pencil,
  Stirling face, powered-line resultant, all-width THM-2824 sign, corrected
  THM-3047 width flag, Beta-Gamma carrier, affine alternant, and higher-slot
  rank collapse.  Both replayed normal and optimized execution against the
  stored transcript; the first candidate's one-step flag shift was repaired
  before promotion and both final delta audits accepted the canonical blobs.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
related:
  - THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary
script: 04-computation/gmc_three_slot_physical_terminal_tail_thm3060.py
output: 05-knowledge/results/gmc_three_slot_physical_terminal_tail_thm3060.out
script_sha256: 4021e3cc59c40bc338ba7640108aaedba607b20b7e46b3ab9707755a3b961b8f
output_sha256: 14b5e527271e683c5bb54ffe350f73d3280b6de4b3aeb72755e1d8f5dddd7aa7
hash_basis: LF-normalized bytes
---

# THM-3060 -- the first actual terminal face survives exactly through three slots

**PROVED + VERIFIED-EXACT + TWO INDEPENDENTLY HOSTILE-AUDITED.**

THM-3054 controls the coefficientwise formal terminal corner.  The physical
specialization does not set the terminal exponentials to zero, so its
dominant face is different.  In three slots that face is nevertheless
nonzero and positive.  In four or more slots it acquires a common projective
zero and cancels.  This gives both a positive physical theorem and an exact
explanation for the first higher-slot obstruction.

## 1. The actual normalized inclusion forms

Fix integers

```text
depth n>=1, moving lower offset C>=1, fixed gap h>=1,
physical terminal width M=C+h.                         (1)
```

For an order `m` tuple `e=(e_1,...,e_m)`, put

```text
T_m(e;n)=(mn+1)_(sum e_i)/product_i (n+1)_(e_i).       (2)
```

This is THM-2925's normalized factorial tensor before denominator clearing.
For `m=2,3` and `0<=j<=m`, start with `j` copies of `C` and `m-j`
copies of `0`, replace any selected entries by `M`, and take signed
inclusion--exclusion:

```text
A_(m,j)(C)=binom(m,j) sum_(J subset [m])
             (-1)^|J| T_m(e(J);n).                    (3)
```

Define the binary forms and their standard Sylvester resultant by

```text
P_(m,C)(x,z)=sum_(j=0)^m A_(m,j)(C)x^(m-j)z^j,
R_C=Res(P_(2,C),P_(3,C)).                              (4)
```

These are the actual normalized physical inclusion forms at width `C+h`.
No formal substitution `2^M=3^M=0` has been made.

There is already an exact all-width sign theorem.  In the normalized
factorial basis `f_j=s^j/j!`, set

```text
H=x(f_n-f_(n+M))+z(f_(n+C)-f_(n+M)).                 (4a)
```

Then `L(H)=0`, and, up to the positive factors `L(f_n^m)`, the forms in
`(4)` are exactly `L(H^2)` and `L(H^3)`.  THM-2824 says that a nonzero
polynomial on any three distinct factorial slots cannot have its first
three moments all zero.  Thus `P_(2,C)` and `P_(3,C)` have no common point
in `C P^1`.

Moreover `P_(2,C)` is a positive-definite real binary quadratic.  If
`alpha,conj(alpha)` are its dehomogenized roots and `a>0` its leading
coefficient, the standard Sylvester orientation gives

```text
R_C=a^3 P_(3,C)(alpha)P_(3,C)(conj(alpha))
   =a^3 |P_(3,C)(alpha)|^2>0.                         (4b)
```

Consequently

```text
R_C>0 for every n>=1, h>=1, C>=1.                    (4c)
```

## 2. The physical all-large face

In a term of `(3)`, suppose exactly `t` offsets are of the form `C+delta`
and the rest remain bounded.  Stirling's formula gives exponential base
`t^(tC)`.  Hence the unique maximal base in an order-`m` coefficient comes
from replacing every bounded offset by `M`; the original `C` entries may be
kept or replaced.

For fixed shifts `delta_i`, the full Gamma-ratio expansion is

```text
T_m(C+delta_1,...,C+delta_m;n)
 =kappa_m m^(mC+sum delta_i) C^((1-m)/2)
  (1+c_1/C+c_2/C^2+...),                              (5)

kappa_m=(2pi)^((1-m)/2) m^(mn+1/2)
        Gamma(n+1)^m/Gamma(mn+1)>0.                   (6)
```

The expansion exists to every fixed inverse power, with a remainder of the
next order.  The all-large finite difference in the coefficient of
`x^(m-j)z^j` is exactly

```text
sum_(r=0)^j (-1)^(m-j+r) binom(j,r)m^((m-j+r)h)
 =(-m^h)^(m-j)(1-m^h)^j.                              (7)
```

Every `t<m` stratum is exponentially smaller.  Coefficientwise,

```text
P_(m,C)=kappa_m m^(mC) C^((1-m)/2)
 [L_m(x,z)^m+O(C^-1)+O(C^D rho_m^C)],                 (8)

L_m(x,z)=-m^h x+(1-m^h)z,
rho_m=(m-1)^(m-1)/m^m<1.                              (9)
```

The `O(C^-1)` term in `(8)` has a complete Poincare expansion, not merely a
size estimate.

The two leading lines are distinct:

```text
det[[-2^h,1-2^h],[-3^h,1-3^h]]=3^h-2^h>0.           (10)
```

Since `Res(L_2^2,L_3^3)=det(L_2,L_3)^6`, resultant
multihomogeneity and `(8)` yield

```text
R_C=K_(n,h) 46656^C C^(-7/2)(1+O(C^-1)),              (11)

K_(n,h)=kappa_2^3 kappa_3^2(3^h-2^h)^6>0.            (12)
```

Again `(11)` has a full inverse-power expansion.  Its leading sign agrees
with, and independently explains the tail of, the all-width theorem `(4c)`.

This is a genuinely physical terminal-face result.  The exponent `-7/2`
comes from three quadratic coefficient blocks of weight `-1/2` and two cubic
blocks of weight `-1`.

## 3. An exact positive carrier

For three slots, THM-3047's exact width flag is

```text
F_M=(n)_M^5(n+1)_M^2/n^(7M).                          (14)
```

Set

```text
Q_C=F_(C+h)R_C.                                       (15)
```

The asymptotic power in `(11)` is itself carried by a positive law.  Let

```text
B_C=46656^C (1)_C/(9/2)_C,                            (16)
G_C=F_(C+h)B_C.                                       (17)
```

`(B_C)` is the moment sequence of `46656*Beta(1,7/2)`.  The shifted width
flag in `(14)` is a product-Gamma moment sequence, so `(G_C)` is strict
Stieltjes with full positive support.  Equations `(11)` and `(16)` have the
same exponential base and power; their positive leading constants need not
be equal.

## 4. Physical affine tail holotopy

Fix integers `s>=1,b>=0`, write `C_c=sc+b`, and define

```text
H_c^(theta)=(1-theta)G_(C_c)+theta Q_(C_c),
0<=theta<=1.                                          (18)
```

For every fixed order bound `d` and offset bound `B_0`, there is an effective
threshold

```text
c_0=c_0(n,h,s,b,d,B_0)                                (19)
```

such that for all `c>=c_0`, `theta in [0,1]`, and

```text
1<=r<=d,
0<=u_1<...<u_r<=B_0,
0<=v_1<...<v_r<=B_0,                                  (20)
```

one has

```text
det[H_(c+u_i+v_j)^(theta)]_(i,j=1)^r>0.               (21)
```

Here is the point of retaining the full expansions in `(5)` and `(11)`.
Both summands in `(18)` have the same positive factorial scale and admit
complete inverse-`c` expansions.  Their leading constants are positive, so
the leading constant of the convex combination is bounded away from zero
uniformly in `theta`.  For a fixed shift `t`, after the exponential scale is
removed,

```text
H_(c+t)^(theta)/H_c^(theta)
 =Lambda^t c^(7s t)(1+a_1(theta,t)/c+...),             (22)
```

uniformly on `theta in [0,1]`.  Confluent alternant expansion gives, with
`N_r=r(r-1)/2`, the first nonzero coefficient

```text
(7s)^N_r V(u)V(v)/product_(ell=0)^(r-1)ell!>0.         (23)
```

Only finitely many minors occur in `(20)`, so their asymptotic remainders are
uniform and one threshold proves `(21)`.  Thus the normalized physical
three-slot resultant eventually lies in every prescribed finite strict
total-positivity chamber, and the explicit path `(18)` stays in that chamber.

The quantifier order is still

```text
for every fixed (d,B_0), there exists c_0(d,B_0).      (24)
```

No single tail is asserted to be Stieltjes at all orders.

## 5. Why the same face stops at four slots

The calculation also identifies the exact higher-slot boundary.  With
`k-2` fixed lower variables `x_1,...,x_(k-2)`, one moving variable `z`, and
orders `r=2,...,k`, the all-large face is

```text
L_r^r,
L_r=-r^h(x_1+...+x_(k-2))+(1-r^h)z.                  (25)
```

For `k>=4`, choose `z=0` and a nonzero vector with
`x_1+...+x_(k-2)=0`.  It is a common projective zero of every `L_r`, so

```text
Res(L_2^2,...,L_k^k)=0.                               (26)
```

Equivalently, the coefficient matrix of the leading lines has rank at most
two because all fixed-lower columns coincide.  This is the first failed
implication in any attempt to extend `(11)` verbatim: the physical dominant
coefficient strata survive, but their resultant cancels.  A higher-slot
theorem must compute the first transverse augmentation of this rank-two
face.  Equation `(26)` is a stopping theorem, not evidence that the next
face vanishes.

## 6. Exact evidence and scope

The dependency-free exact companion checks:

- `56` all-large finite-difference coefficients and eight powered-line
  resultants;
- `375` exact physical resultants for `1<=n,h<=5`, `1<=C<=15`, all positive
  as independently required by THM-2824;
- `16` exact ratio-slope controls around the predicted `-7/2` exponent;
- `1,008` strict carrier minors and `9,072` rational homotopy minors;
- the rank collapse for every `4<=k<=9`.

Run

```text
python 04-computation/gmc_three_slot_physical_terminal_tail_thm3060.py
python -O 04-computation/gmc_three_slot_physical_terminal_tail_thm3060.py
```

Both modes equal the stored eleven-line transcript after LF normalization.

This theorem concerns the normalized intrinsic inclusion resultant at fixed
depth.  It does not prove persistence of a selected denominator-cleared
Macaulay chart, positivity of a wall-stripped core, a global jet polynomial,
or an arbitrary nonlinear moving clock.  It does not prove a new NC2, GMC,
SFC, LRC, or Jacobian statement.  All-width scalar positivity `(4c)` does
not upgrade the sequence to one all-order Stieltjes tail.

Two independent hostile audits rederived the asymptotic and holotopy by
separate paths.  A final delta audit also checked the THM-2824 all-width
strengthening and the standard Sylvester orientation.  No defect remains.

**QED.**
