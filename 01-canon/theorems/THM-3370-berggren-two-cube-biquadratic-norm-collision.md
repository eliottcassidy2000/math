---
id: THM-3370
title: "Berggren U-spine and two-cube biquadratic norm collision"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.  A positive two-distinct-
  cube representation of the Berggren U-spine scalar Q_r=(2r+1)^2+2 is
  equivalent to a positive-distinct good divisor in THM-463.  After removing
  the common cube gcd, its Eisenstein cofactor is supported on 3 (to exponent
  at most one) and primes congruent to 1 or 19 modulo 24.  The depth lies in
  15 explicit classes modulo 63.  The first such scalar is exactly
  13,712,211=107^3+232^3=3703^2+2=Q_1851.  Five positive
  intersections occur in the exact box 1<=x<y<=5000.  Separately, the
  norm-one unit 23+4sqrt(33) gives an infinite fixed-sum Pell family with
  distinct signed integer cubes, starting at Q_31=16^3+(-5)^3.  Positive
  infinitude is not proved here and is subsequently closed by THM-3375;
  density, LRC, FC and JC remain open.
source: kps-s177-berggren-two-cube-2026-08-14
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
related:
  - THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy
  - THM-3368-weighted-berggren-horn-defect-tariff-and-complement-clock-separation
  - THM-3375-berggren-positive-two-cube-pell-ray
script: 04-computation/berggren_two_cube_norm_collision_kps_s177.py
output: 05-knowledge/results/berggren_two_cube_norm_collision_kps_s177.out
script_sha256: e0321576b71764fb65aa15cfd843869a9846b7bea1310d7cd6de8faed0cfc665
output_sha256: 8d61cd5dd4255ad2bd58f4bbc7172b7fb98cbe97c2943e09be31774aea11d60a
hash_basis: LF-normalized bytes
---

# THM-3370 -- the Berggren/two-cube intersection is a mod-24 norm collision

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**

This theorem joins two proved repository objects without identifying their
carriers.  THM-3334 supplies the consecutive-parameter Berggren `U`-spine;
THM-463 supplies the complete good-divisor criterion and Eisenstein cofactor
law for sums of two cubes.  Their intersection has a new biquadratic split
filter, a sharp first point, and an equally important coprime-carrier no-go.

## 1. The U-spine scalar

For `r>=1`, put

```text
P_r=(A_r,B_r,C_r)
   =(2r+1, 2r(r+1), 2r(r+1)+1),                         (1)
Q_r=2C_r+1=2B_r+3=(2r+1)^2+2=4r(r+1)+3.                (2)
```

These are the positive integral points on `C-B=1` in the Pythagorean cone:

```text
A_r^2+B_r^2=C_r^2.                                      (3)
```

For the standard Berggren matrix

```text
U=((1,-2,2),(2,-1,2),(2,-2,3)),                         (4)
```

direct multiplication gives `U P_r=P_(r+1)`.  Thus `(2)` is exactly the
scalar sequence `11,27,51,83,...` attached to the `U`-spine in THM-3334.
It displays two adjacent quadratic norms:

```text
C_r=(r+1)^2+r^2=N_(Z[i])((r+1)+ri),                     (5)
Q_r=(2r+1)^2+2=N_(Z[sqrt(-2)])((2r+1)+sqrt(-2)).        (6)
```

They are affine neighbours but arithmetically disjoint:

```text
gcd(C_r,Q_r)=gcd(C_r,2C_r+1)=1.                         (7)
```

## 2. Positive distinct cubes are positive distinct good divisors

By THM-463, unordered integer representations `x^3+y^3=n` correspond
bijectively to positive divisors `d|n` for which

```text
s=(d^3-n)/(3d) is integral,
Delta=(4n-d^3)/(3d)=e^2 is a square.                    (8)
```

Under the bijection,

```text
d=x+y,             s=xy,             e=|x-y|.           (9)
```

Consequently `Q_r` is a sum of two **distinct positive** cubes if and only if
it has a good divisor satisfying

```text
0<e<d.                                                   (10)
```

Indeed, the inverse roots are `(d+e)/2,(d-e)/2`; THM-463 proves their parity
is automatic.  Condition `(10)` is exactly positivity and distinctness.
This also fixes the convention here: zero and negative summands are excluded.

## 3. The mod-24 cofactor theorem

Suppose

```text
Q_r=x^3+y^3,                    x,y distinct integers.    (11)
```

Let `g=gcd(x,y)`, `X=x/g`, `Y=y/g`, and

```text
d=x+y,
q=x^2-xy+y^2=g^2 q_0,
q_0=X^2-XY+Y^2=N_(Z[omega])(X+Y omega),                 (12)
```

where `omega^2+omega+1=0`.  Then

```text
Q_r=dq=g^3(X+Y)q_0.                                     (13)
```

THM-463's primitive split lemma applies to `(X,Y)`: `q_0` is odd,
`v_3(q_0)<=1`, and every prime `p!=3` dividing `q_0` satisfies

```text
p=1 (mod 3).                                            (14)
```

But `q_0|Q_r`.  For such an odd `p`, equation `(2)` gives

```text
(2r+1)^2=-2 (mod p).                                    (15)
```

The left side cannot vanish modulo `p`, since that would make `p|2`.
The standard quadratic-character calculation

```text
(-2/p)=(-1)^((p-1)/2) (-1)^((p^2-1)/8)                 (16)
```

therefore shows

```text
p=1 or 3 (mod 8).                                       (17)
```

Combining `(14)` and `(17)` by the Chinese remainder theorem leaves exactly

```text
p=1 or 19 (mod 24).                                     (18)
```

Hence, for arbitrary common cube gcd `g`, the **primitive Eisenstein
cofactor** has the sharp form

```text
q_0=3^epsilon product_i p_i^(a_i),
epsilon in {0,1},             p_i=1 or 19 (mod 24).      (19)
```

For `g=1`, this is the cofactor `q` itself.  Primes in `g` need not satisfy
the Eisenstein condition `(14)`, so omitting the primitive reduction would
make `(19)` false in scope.  The sharp stored witness is

```text
Q_591=120^3+(-69)^3=1183^2+2,
g=3,             q=27,441=3^2*3,049,
q_0=3,049=1 (mod 24).                                   (20)
```

Thus `v_3(q)=2`, while the normalized `q_0` obeys `(19)`.  Equation `(7)`
supplies another boundary: the
Gaussian prime-toggle fibre of the native hypotenuse `C_r` cannot literally
be the Eisenstein factor fibre of `Q_r`.  The two norms meet through the
affine map `Q=2C+1`, not through shared prime factors.

## 4. A cheap local depth sieve

Cubes modulo `7` and modulo `9` are `0,+1,-1`; their pairwise sums are
`0,+1,-1,+2,-2`.  Substituting `(2)` gives the exact necessary conditions

```text
r mod 7 in {2,3,4},
r mod 9 in {1,2,4,6,7}.                                 (21)
```

The CRT intersection is the following `15` classes:

```text
r mod 63 in
{2,4,10,11,16,24,25,31,37,38,46,51,52,58,60}.           (22)
```

This removes `16/21` of depths, but it is only a local necessary sieve.  It
does not encode the square condition in `(8)` and cannot establish either
finiteness or infinitude of the intersection.

## 5. First point and exact box census

The companion scans every one of the

```text
binom(5000,2)=12,497,500
```

pairs `1<=x<y<=5000` using integer cube arithmetic and `isqrt`.  Ordinary and
optimized runs agree.  Exactly five pairs satisfy `(11)`:

| `x` | `y` | `r` | `Q_r=x^3+y^3` | `g` | primitive cofactor `q_0` |
|---:|---:|---:|---:|---:|---:|
| 107 | 232 | 1,851 | 13,712,211 | 1 | `40,449=3*97*139` |
| 360 | 3,003 | 82,352 | 27,127,737,027 | 3 | `896,281` |
| 1,907 | 3,472 | 110,441 | 48,789,299,691 | 1 | `9,070,329=3*3,023,443` |
| 107 | 4,960 | 174,660 | 122,025,161,043 | 1 | `24,082,329=3*19*409*1,033` |
| 4,403 | 4,576 | 212,825 | 181,178,773,803 | 1 | `20,178,057=3*19*354,001` |

Every displayed non-`3` cofactor prime is `1` or `19 mod 24`, as `(19)`
requires.  The first row factors more fully as

```text
13,712,211=3^2*97*113*139,
d=x+y=339=3*113,
q=q_0=40,449=3*97*139.                                  (23)
```

It is also the **global first positive-distinct intersection**.  If a smaller
`n=x^3+y^3` existed with `x<y`, then `y^3<n<13,712,211`, hence `y<=239`.
The companion exhausts those `28,441` pairs exactly and finds none.  Therefore

```text
min {Q_r=x^3+y^3 : r>=1, 1<=x<y}
 =13,712,211
 =107^3+232^3
 =3703^2+2
 =Q_1851.                                                (24)
```

The other four rows are a coordinate-box census only.  No height-complete
claim beyond the minimum is made.

## 6. An infinite signed Pell fibre

The positivity restriction is genuinely load-bearing.  For arbitrary
integer summands, put

```text
d=x+y,                 e=x-y,                 a=2r+1.
```

Then

```text
x^3+y^3=d(d^2+3e^2)/4,
Q_r=x^3+y^3 iff (2a)^2-3d e^2=d^3-8.                    (25)
```

Thus every fixed `d` is a Pell-type conic.  The fibre `d=11` has seed

```text
(a_0,e_0)=(63,21),
(x_0,y_0)=(16,-5),
Q_31=3,971=16^3+(-5)^3=63^2+2.                          (26)
```

The unit

```text
23+4 sqrt(33),                    23^2-33*4^2=1,         (27)
```

gives the integral recurrence

```text
a_(m+1)=23a_m+66e_m,
e_(m+1)=8a_m+23e_m.                                     (28)
```

Indeed, `(28)` is multiplication of `2a_m+e_m sqrt(33)` by `(27)`, so it
preserves `(25)`.  It also preserves odd parity, and both coordinates grow
strictly.  Hence

```text
x_m=(11+e_m)/2,              y_m=(11-e_m)/2,
r_m=(a_m-1)/2                                             (29)
```

are integers for every `m`, with `x_m^3+y_m^3=Q_(r_m)`.  The first iterate is

```text
(a_1,e_1,x_1,y_1,r_1,Q_(r_1))
=(2,835,987,499,-488,1,417,8,037,227).
```

This proves infinitely many **natural scalars** `Q_r` that are sums of two
distinct signed integer cubes.  Since already `e_0=21>d=11`, every point in
this forward orbit has one negative summand.  Distinct positive pairs are
exactly the Pell-window points `0<|e|<d`; their infinitude is not proved.

## 7. Cross-frontier interpretation and loss ledger

The actual junction is

```text
r
 -> (C_r,Q_r)                         Gaussian / sqrt(-2) norms
 -> (d,q_0,g) when Q_r=x^3+y^3        Eisenstein factorization. (30)
```

It preserves the depth, scalar and a witnessed cube pair.  Passing to `Q_r`
destroys THM-3346's Gaussian factor-choice labels because of `(7)`; passing
to the local residue set `(22)` destroys the good-divisor square condition;
passing to the unordered pair destroys its orientation.  Thus the necessary
sidecar for any later use is `(r,x,y,g,d,q_0)`, not merely a congruence class
or a CM-field label.

Homogenizing `(11)` gives the smooth cubic surface

```text
X^3+Y^3-A^2 Z-2Z^3=0.                                  (31)
```

Its four partial derivatives cannot vanish simultaneously in projective
space, so it is smooth.  Section 6 gives one infinite integral curve on it;
the honest next question is to classify the **positive chamber**, rather
than search for another modulus after `63`.  This also explains why the lane
resembles FC(3)'s diagonal cubics and
the Hessian/Jacobian work.  The observer is different, however: `(31)` asks
for integer **values**, while factorial moments ask for coefficient-weighted
cancellation.  No implication between them is supplied here.

Similarly, `(22)` is a periodic modular filter, but the exact positive support
need not be periodic; THM-3359's value-support warning applies.  No
density or harmonic pole follows from `15/63`.

## 8. Scope

Equations `(1)--(31)`, the positive minimum `(24)`, and the infinite signed
family `(26)--(29)` are proved.  The five-row positive box census is
FINITE-EXACT and the algebra/minimum were independently audited.  The theorem
does not classify the positive integral points of `(31)`, prove infinitely
many positive intersections, or show that the five displayed rows are the
only ones below a scalar bound larger than `(24)`.  The first of these
residuals is subsequently closed by THM-3375's moving-sum Pell ray; it is not
a consequence of the fixed-sum orbit proved here.  This theorem supplies
no physical runner row, owner, phase, factorial-moment tower, Keller mate, or
AMM integer flow.  LRC(14), FC(3), the planar Jacobian conjecture and AMM
12592 remain open.

## Reproduction

```text
python 04-computation/berggren_two_cube_norm_collision_kps_s177.py
python -O 04-computation/berggren_two_cube_norm_collision_kps_s177.py
```

Runtime checks remain active under optimized Python.

**End of proof.**
