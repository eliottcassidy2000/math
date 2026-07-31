---
id: THM-2890
title: "Discrete-Newton closure of the strict consecutive GMC wedges"
status: >
  PROVED + VERIFIED-EXACT.  For every integer n>=1, the strict low
  consecutive wedge c1=a*c2+b*c3 with a,b>0 is cubic-empty: its binary
  quadratic factorial moment does not divide its binary cubic factorial
  moment.  Exact low-chart elimination leaves a degree-fifteen factor
  R(n,a).  Its ordinary monomial coefficients are not sign-definite, but
  all twenty-three Gregory--Newton coefficient polynomials are positive
  on a>0, except the terminal one which is an explicit nonnegative square.
  The sole quadratic resultant factor is saturated exactly: it is positive
  for n>=2 and at n=1 supports only the excluded b=0 boundary.  Combined
  with THM-2879, all three strict consecutive cone-cutting wedge charts
  miss a shared cubic--quartic multipole line.  Cone boundaries and the
  general mixed four-slot branch remain open.
source: codex/gmc-low-sector-newton-2026-07-29
depends_on:
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
related:
  - THM-2855-shifted-positive-cone-transverse-asymptotic-family
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
  - THM-2873-two-ray-factorial-response-tp3-curvature
script: 04-computation/gmc_strict_consecutive_low_wedge_newton_thm2890.py
output: 05-knowledge/results/gmc_strict_consecutive_low_wedge_newton_thm2890.out
script_sha256: cca87f78f8c250a6ade5c5dbb6bf3c7a390db03fe79dff1a5710768b0f14fd78
output_sha256: ead7e19aeb032571b4db8f1e386c668a3edcb78f7eaa73251cf9abf29166a1f7
hash_basis: LF-normalized bytes
---

# THM-2890 -- discrete-Newton closure of the strict consecutive GMC wedges

**PROVED + VERIFIED-EXACT.**

Put

```text
L(s^q)=q!,                  f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

Fix an integer `n>=1` and write

```text
e_1=d_n,              e_2=d_(n+1),          e_3=d_(n+2).  (2)
```

For every `a,b>0`, consider the strict low-chart plane

```text
c_1=a c_2+b c_3.                                      (3)
```

Its two positive extreme rays are

```text
a e_1+e_2,                         b e_1+e_3.          (4)
```

For any basis `U,V` of this plane, define

```text
Q(alpha,beta)=L((alpha U+beta V)^2),
C(alpha,beta)=L((alpha U+beta V)^3).                   (5)
```

Then

```text
Q does not divide C                                    (6)
```

for every integer `n>=1` and every `a,b>0`.

Consequently the sole strict cone-cutting orientation left by THM-2879 is
empty already at the cubic stage.  Together with that theorem's strict
quartic exit in the high chart and cubic emptiness in the middle chart,
all three strict consecutive wedge charts miss a shared cubic--quartic
multipole line.

## 1. The low rechart and its exact consequence object

In the high coordinates of THM-2879, a plane has equation

```text
c_3=X c_1+Y c_2.                                      (7)
```

Solving `(3)` for `c_3` gives

```text
X=1/b>0,                       Y=-a/b<0.               (8)
```

Let `I_1(n,X,Y),I_2(n,X,Y)` be the two exact cubic remainders used in
THM-2855/2879.  Their simultaneous vanishing is equivalent to `Q|C`.
Substitute `(8)`, cancel the common rational factors, and retain the
numerators

```text
L_1(n,a,b),                       L_2(n,a,b).          (9)
```

Every cancelled denominator is positive for `n>=1,a,b>0`.  The exact
profiles are

```text
             deg_b   deg_a   deg_n   terms
L_1             5       3       8      153
L_2             5       4       7      136.           (10)
```

Therefore a low cubic-divisible plane would force

```text
Res_b(L_1,L_2)=0.                                    (11)
```

The resultant itself, rather than an intermediate selector statistic, is
the consequence object checked below.

## 2. Exact factorization and the two easy factors

Over `Z[n,a]`, exact factorization gives

```text
Res_b(L_1,L_2)
 =432 P_+(n) a^2 Q_2(n,a) Q_4(n,a) R_15(n,a),        (12)
```

where:

- `P_+` is the product of six `n`-only factors, each coefficientwise
  positive;
- `Q_4` has bidegree `(4,4)` and all `25` coefficients positive;
- `Q_2` has bidegree `(5,2)` and `18` terms; and
- `R_15` has bidegree `(22,15)` and `368` terms.

The complete factor profile, including multiplicities, is

```text
(deg_a,deg_n,multiplicity,terms)

(0,1,2,2), (0,1,3,2), (0,1,3,2), (0,1,18,2),
(0,5,1,6), (0,10,1,11),
(1,0,2,1), (2,5,1,18), (4,4,1,25), (15,22,1,368).   (13)
```

Thus `P_+`, `a^2`, and `Q_4` cannot vanish in the strict domain.

As a cross-control on the rechart, the degree-two high resultant factor is

```text
G_n(Y)
 =(n+2)+2(2n+3)Y+2(2n+3)Y^2,                         (14)
```

and

```text
disc_Y G_n=-4(2n+3)<0.                               (15)
```

The exact linear-selector content is a positive `n`-factor times `G_n`.
Hence the content removed in THM-2879 is nonzero on the whole real line,
including the negative-`Y` sector considered here.

## 3. Saturating the quadratic boundary factor

The factor `Q_2` is not discarded by a sign guess.  It is the sole
boundary factor requiring saturation.

### 3.1. Integer tail `n>=2`

Put `n=m+2`.  All `18` coefficients of

```text
Q_2(m+2,a)                                             (16)
```

as a polynomial in `(m,a)` are strictly positive.  Since integer `n>=2`
means integer `m>=0`, `(16)` is positive for every `a>0`.

### 3.2. Exceptional depth `n=1`

At `n=1`, orient the primitive factor as

```text
q(a)=a^2-(435/109)a-595/109.                          (17)
```

It has one positive root.  Directly,

```text
gcd(L_1(1,a,0),L_2(1,a,0))
 =nonzero_rational_constant * a * q(a).              (18)
```

Thus the apparent positive root of `q` already supplies the projective
boundary root `b=0`; this boundary is not in `(3)`.  To prove that it
does not hide another finite root, work in the quadratic field

```text
K=Q[a]/(q).                                           (19)
```

In `K[b]`, equation `(18)` gives

```text
L_i(1,a,b)=b F_i(a,b),               deg_b F_i=4.     (20)
```

Exact reduction computes `Res_b(F_1,F_2)` modulo `q`.  The residue is
linear in `a`, and its field norm is nonzero.  Therefore it is a unit in
`K`; the residual quartics have no common root.  The only common root on
`q=0` is `b=0`, so no strict `b>0` point is lost.

The canonical digest of the linear residue is

```text
efb3a831ee7a5d2429e1e48919a2397702bdd85c524df1db6fe255ab27823e76.
                                                               (21)
```

This is the load-bearing saturation.  Merely observing that `Q_2` comes
from a boundary would not suffice.

## 4. The discrete Newton certificate for `R_15`

Ordinary coefficient positivity fails:

```text
R_15(n,a) in the (n,a) monomial basis:       326 positive,
                                               42 negative;

R_15(m+1,a) in the (m,a) monomial basis:     364 positive,
                                                4 negative.    (22)
```

Thus a Maclaurin-style continuous-parameter certificate misses the
integer structure.  Use the Gregory--Newton basis instead.  Define

```text
D_k(a)=Delta_n^k R_15(1,a),                 0<=k<=22. (23)
```

Since `deg_n R_15=22`, finite-difference interpolation is the exact
identity

```text
R_15(n,a)
 =sum_(k=0)^22 binom(n-1,k) D_k(a)                    (24)
```

for every `n`, not an asymptotic expansion.

For each `0<=k<=21`, the exact Sturm remainder sequence gives

```text
deg D_k=15,
D_k(0)>0,
LC_a(D_k)>0,
#{a>0:D_k(a)=0}=0.                                   (25)
```

Hence

```text
D_k(a)>0                         for a>0, 0<=k<=21.   (26)
```

The terminal coefficient factors exactly as

```text
D_22(a)
 =1024*22!*(a+2)^6*(a+3)^3
   *(a^2-12a-36)^2*(7a^2+24a+18)>=0                 (27)
```

for `a>0`.  If `n>=1` is an integer, every
`binom(n-1,k)` in `(24)` is nonnegative, while the `k=0` term is strictly
positive.  Equations `(24)--(27)` prove

```text
R_15(n,a)>0                  for integer n>=1, a>0.   (28)
```

This is the mechanism: the natural positivity basis belongs to the
discrete depth variable.  The continuous monomial basis records genuine
sign changes that are irrelevant at integer depth.

The canonical exact digests are

```text
R_15 coefficient dictionary:
c9c35f414c1411800080d561e8df2b7ebedab505414275021bf15ed858b3e81d,

twenty-three Newton coefficient dictionaries:
b9f87fb39354603d3976a538e2f82510433decd38e7ece8230cb472885c8db63.
                                                               (29)
```

## 5. Completion and scope

Every factor in `(12)` is now nonzero for integer `n>=1,a,b>0`.
Therefore `(11)` is impossible, proving `(6)`.

THM-2879 supplies the other two strict charts:

1. the strict shared-high chart has one cubic-divisible branch, but its
   quartic endpoint holonomy is strictly negative; and
2. the strict shared-middle chart is cubic-empty after an invertible
   rechart.

The present theorem proves:

3. the strict shared-low chart is cubic-empty.

Therefore no strict consecutive cone-cutting plane is a shared
cubic--quartic multipole line.

The strict inequalities are load-bearing.  This theorem does **not** treat
zero-coordinate cone boundaries, cone-avoiding planes, nonconsecutive
four-slot planes, or the general mixed shared-line branch.  It proves no
new instance of SFC beyond this chart, does not replace the unrestricted
GMC(2) proof, removes no DvdK dependency, and implies nothing about the
planar Jacobian conjecture.

## 6. Exact verification

Run

```text
python3 04-computation/gmc_strict_consecutive_low_wedge_newton_thm2890.py
python3 -O 04-computation/gmc_strict_consecutive_low_wedge_newton_thm2890.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_strict_consecutive_low_wedge_newton_thm2890.out.
```

The companion pins the exact THM-2879 dependency hash and checks:

1. the high `G_n^2 P_15` resultant and globally nonvanishing selector
   content;
2. every low clearing denominator and both numerator profiles;
3. the complete resultant factorization `(12)--(13)`;
4. `Q_2` tail positivity and the exact `n=1` quotient-field saturation;
5. all twenty-three exact Sturm certificates and the terminal square;
6. reconstruction of `(24)` as a polynomial identity; and
7. fixed-depth `R_15` factor-level positive-root controls at
   `n=1,2,8,64`.

**QED.**
