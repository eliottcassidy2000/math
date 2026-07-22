---
id: THM-2138
title: "All-depth unit-annulus extremality closes the six-plus-one scalar tail"
status: >
  PROVED AND HOSTILE-REFEREED. For every m>=3, the largest
  intersection of one terminal unit mask with the 13-primitive guard-safe
  annulus modulo 13^m is attained by exactly the two sign classes a=6 and
  12a=-1. Its size is (10*13^m+130*(-1)^m)/91. For m>=5 every other mask
  loses at least twelve points. Consequently six such masks never cover the
  unit annulus. Together with THM-2135's exact mod-169 certificate, this
  eliminates the complete six-unit/one-deep scalar tail from THM-2133.
source: codex-2026-07-22-LRC-unit-annulus-all-depth
depends_on:
  - THM-2133
  - THM-2135
related:
  - THM-2080
  - THM-2128
script:
  - 04-computation/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.py
output:
  - 05-knowledge/results/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.out
---

# THM-2138 -- the all-depth unit-annulus law

Put `N=13^m`, and let

```text
U_N={z mod N:13 does not divide z and 7||z||_N>N},
S_a={z in U_N:14||az||_N<N},                         (1)
```

where `a` is a unit modulo `N` and `||x||_N=min(x mod N,-x mod N)`.
The sign classes `a` and `-a` define the same mask.

For every `m>=3`, set `epsilon=(-1)^m`. Then

```text
max_a |S_a| = M_m=(10N+130epsilon)/91.               (2)
```

The only maximizing sign classes are

```text
a=6,                    12a=-1 mod N.                (3)
```

For `m>=5`, every other sign class satisfies

```text
|S_a|<=M_m-12.                                        (4)
```

The finite exact census in THM-2135 proves (2)--(3) at `m=3,4` (and also
at `m=5`).  The proof below is uniform for `m>=5`.

## 1. Invert the mask and retain the primitive lattice layer

Let `b=a^(-1) mod N` and define the index-`2N` lattice

```text
Lambda_N(b)={(x,2y) in Z x 2Z:x=by mod N}.            (5)
```

Let `Q_N` be the open square

```text
Q_N={(x,Y):|x|<N/7 and |Y|<N/7},                     (6)
```

and write `G_N(b)=#(Lambda_N(b) intersect Q_N)`.  Multiplication by `a`
changes a terminal point `z` into the unique representative
`y=az` with `|y|<N/14`; its preimage is guard-unsafe exactly when the least
representative `x=by` has `|x|<N/7`.  Therefore the number removed from the
terminal mask by the guard is the number of points of (5)--(6) whose two
coordinates are not both divisible by thirteen.

The congruence in (5), and the fact that `b` is a thirteen-unit, give

```text
13|x iff 13|y iff 13|2y.                              (7)
```

Division by thirteen identifies the points of `Lambda_N(b) intersect Q_N`
with both coordinates divisible by thirteen with
`Lambda_(N/13)(b mod N/13) intersect Q_(N/13)`. Hence the exact primitive
subtraction is

```text
P_m(b)=G_N(b)-G_(N/13)(b mod N/13).                   (8)
```

If

```text
R_m=floor((N-1)/14),
T_m=(N+13epsilon)/91,                                 (9)
```

then the full terminal band contains

```text
B_m=2(R_m-R_(m-1))=12T_m                              (10)
```

thirteen-primitive points. Consequently

```text
|S_a|=B_m-P_m(b).                                     (11)
```

It remains to prove `P_m(b)>=2T_m`, with equality only in (3).

## 2. The one-line branch

Choose a nonzero vector `u=(r,2s)` of `Lambda_N(b)` with minimum sup norm

```text
d=max(|r|,2|s|).                                      (12)
```

Suppose first that `d<=12`.  Then `u` is nonzero modulo thirteen, and the
thirteen-primitive multiples `ku` lying in `Q_N` already number

```text
2C_d(N),
C_d(N)=floor((N-1)/(7d))-floor((N-1)/(91d)).          (13)
```

The function in (13) counts the positive integers up to `(N-1)/(7d)` that
are not divisible by thirteen, so it is nonincreasing in `d`.  Moreover,
because `13^m` is `1,13` modulo `84` according as `m` is even, odd, while
modulo `1092=13*84` the second floor is determined by
`13^(m-1)=13,1 mod 84` in the two parities,

```text
C_12(N)
 =floor((N-1)/84)-floor((N-1)/1092)
 =(N+13epsilon)/91=T_m.                               (14)
```

Thus the central line alone proves `P_m(b)>=2T_m`. If `d<=11`, the floor
bounds give the quantitative strictness

```text
C_11(N)-C_12(N)>=(N-1)/1001-2,                       (15)
```

which is greater than ten for `N>=13^5`.

It remains to classify equality when `d=12`. Since every integer divisor of
`r,s` is then prime to `N`, minimality implies `gcd(r,s)=1`. Up to signs, the
possibilities are

```text
(|r|,|s|)=(12,1),(12,5),(1,6),(5,6),(7,6),(11,6).   (16)
```

For `(12,1)` and `(1,6)`, every point `z` of the lattice in `Q_N` is
collinear with `u`. Indeed, `det(u,z)` is a multiple of `2N`, whereas

```text
|det(u,z)|<N(|r|+2|s|)/7<=2N.                        (17)
```

Thus these two directions have exactly the primitive points in (13), and
give equality.

For each of the other four directions in (16), `|r|+2|s|>14`. Extend the
shortest lattice vector to a basis. The first noncentral parallel lattice
line has an open chord in `Q_N` of parameter length at least `N/140`.
For `N>=13^5` this is greater than fourteen. Its lattice points are
consecutive integral parameters; modulo thirteen at most one parameter
class gives a point divisible by thirteen. Any fourteen consecutive
parameters therefore contribute at least twelve new primitive points. Hence

```text
P_m(b)>=2T_m+12                                      (18)
```

for every nonextremal `d=12` direction.

Finally, `(r,s)=(12,1)` says `b=12` up to sign, while `(r,s)=(1,6)` says
`6b=1` up to sign. Since `a=b^(-1)`, these are exactly the two sign classes
in (3).

## 3. The two-dimensional slicing branch

We now suppose `d>=13`. A shortest lattice vector is primitive as an element
of the free abelian group `Lambda_N(b)`, so extend it to a basis of determinant
`2N`. Reflect and interchange the square coordinates so that the absolute
coordinates of `u` are `d>=e>=1`.

Here `e` cannot vanish.  An axis vector of (5) has sup norm at least `N`,
whereas the vector obtained from `y=1` and the least representative of `b`
has sup norm less than `N`.

The parallel lattice lines are indexed by `j in Z` and have determinant
`2Nj` from `u`. Their open chord lengths, measured in the lattice parameter
along `u`, are

```text
ell_j = 2N/(7d)                              if 14|j|<=d-e,
ell_j = N(d+e-14|j|)/(7de)                  if d-e<14|j|<d+e,
ell_j = 0                                   otherwise.          (19)
```

This follows by writing the line as `dy-ex=2Nj` inside
`|x|,|y|<N/7`. An open interval of length `ell` contains at least `ell-1`
integers. Therefore, with

```text
Phi(d,e)=sum_j ell_j/N,
L(d,e)=#{j:ell_j>0}<=1+2d/7,                         (20)
```

we have the rigorous lattice-point lower bound

```text
G_N(b)>=N Phi(d,e)-L(d,e).                            (21)
```

There is a uniform arithmetic bound

```text
Phi(d,e)>=6/175.                                      (22)
```

Here are the details. The continuous extension of `ell_j/N` is even,
nonnegative and nonincreasing on the positive half-line. Its integral is
the area of `Q_N`, divided by the lattice covolume and by `N`, namely
`2/49`, and its value at zero is `2/(7d)`. Comparing the integral with the
left endpoint sums on the two half-lines gives

```text
Phi(d,e)>=2/49-2/(7d).                                (23)
```

This proves (22) for `d>=44`. For the remaining finite interval, put

```text
p=floor((d-e)/14),               J=floor((d+e-1)/14).
```

Summing (19) gives the exact expression

```text
Phi=2(1+2p)/(7d)
 +2[(J-p)(d+e)-7(J-p)(J+p+1)]/(7de).                 (24)
```

One coordinate of `(d,e)` is even, and (7) says that thirteen divides both
coordinates or neither. Substitution in (24) gives the following minima;
`d=13` has no admissible `e`.

```text
d       14    15     16      17      18      19      20    21
minPhi  2/49  5/126  13/336  23/595  4/105   26/665  11/280 2/49

d       22    23     24    25     26     27    28    29
minPhi  3/77  6/161  1/28  6/175  7/169  1/27  2/49  7/174

d       30     31       32      33       34      35    36
minPhi  5/126  43/1085  11/280  46/1155  19/476  2/49  5/126

d       37      38     39        40    41      42    43
minPhi  10/259  5/133  146/3549  1/28  11/287  2/49  73/1806. (25)
```

The minimum in (25) is `6/175`, proving (22).

Minkowski's first theorem applied to squares of halfwidth tending down to
`sqrt(2N)` gives

```text
d<=sqrt(2N),                                          (26)
```

because the lattice covolume is `2N`. Equations (20)--(22) now give

```text
G_N(b)>=6N/175-1-(2/7)sqrt(2N).                       (27)
```

At the lower level we need no discrepancy estimate: projection to the
second coordinate gives the trivial sharp upper bound

```text
G_(N/13)(b mod N/13)
 <=1+2 floor((N/13-1)/14).                            (28)
```

The right side of (28), plus the desired `2T_m`, is

```text
3N/91-4epsilon/7.                                    (29)
```

In the worse parity, subtracting (29) from (27) leaves at least

```text
3N/2275-11/7-(2/7)sqrt(2N).                           (30)
```

This expression is increasing for `N>=13^5`. Since
`sqrt(2*13^5)<862`, its value there is greater than

```text
3*13^5/2275-1735/7=6044/25>10.                       (31)
```

Thus in the two-dimensional branch `P_m(b)>2T_m+10`. The difference is
even, so it is at least twelve. Together with Sections 1--2, this proves
(2)--(4) for every `m>=5`, and the exact THM-2135 rows finish `m=3,4`.

## 4. Six masks cannot cover

The exact annulus size from THM-2135 is

```text
|U_N|=(60N-130epsilon)/91.                            (32)
```

If `m` is odd, equations (2) and (32) give

```text
6M_m=|U_N|-10,                                        (33)
```

so total capacity already forbids a six-mask cover. If `m` is even, then
`6M_m=|U_N|+10`. Any nonmaximal member loses at least twelve by (4), so its
packet has total capacity at most `|U_N|-2`. If all six are maximal, (3)
says that only two distinct masks occur, and their union has size at most
`2M_m<|U_N|`.  At the one even depth not covered by (4), namely `m=4`,
THM-2135's exact runner-up loss is `3140-3084=56>10`. Thus no six unit
masks cover `U_N` for any `m>=3`.

In the six-plus-one scalar tail of THM-2133, write the deep coefficient as
`v=13^s w` and take `m=s+2`. Normalization by the unit guard turns a proposed
cover into exactly such a six-mask cover of `U_N`. THM-2135 separately
excludes `m=2`. Hence the entire six-unit/one-deep scalar tail is empty.

## 5. Scope and preserved information

The proof uses the congruence lattice rather than a tournament. Its source
object is a unit terminal phase, the map is inversion followed by
`y -> (by,2y)`, and the preserved predicate is simultaneous terminal danger
and guard ineligibility. Quotienting by thirteen destroys precisely the
primitive layer, so (8) retains it by subtraction. The least-used sidecar is
the shortest lattice direction together with the noncentral chord lengths;
the extremizers are the two one-line directions and not merely the two
largest floor sums.

This closes case (I) of THM-2135. It does not address the five-unit/two-deep
tail, the fivefold nonblocker-guard pencil, or higher relation ranks, so
LRC(14) remains open.

## Hostile audit

An independent referee checked the lattice bijection and primitive
subtraction, the complete `d<=12` equality classification, every entry of
the `d=13..43` rational table, the Minkowski and error constants, and the
even-depth union argument.  The audit found no counterexample or substantive
gap and requested only the three local clarifications incorporated above.
QED.
