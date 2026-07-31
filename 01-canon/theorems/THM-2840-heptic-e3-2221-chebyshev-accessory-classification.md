---
id: THM-2840
title: "Heptic e=3 (2,2,2,1) Chebyshev accessory classification"
status: >
  PROVED + VERIFIED-EXACT.  In the first heptic maximal-pole response
  layer N=7,e=3,s=1,h=4, the normalized ordered pole passport
  (2,2,2,1) has a radical saturated accessory algebra of length six.
  Its six marked charts form one unmarked class, exactly the affine
  normalizations of (T_7+1)/(T_7-1).  The neighboring pole multisets
  (4,1,1,1) and (3,2,1,1), Keller-chart entry, JC(2), and DC(2) remain
  open.
source: root/heptic-e3-2221-chebyshev-accessory-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2816-maximal-pole-clean-nielsen-ribbon-tree-prufer-vandermonde-count
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
  - THM-2821-sextic-e3-power-wreath-chebyshev-dihedral-block-lattice
script: 04-computation/jc_heptic_e3_2221_chebyshev_accessory_thm2840.py
output: 05-knowledge/results/jc_heptic_e3_2221_chebyshev_accessory_thm2840.out
script_sha256: e1577aa95f61bddfff836eadf5d9b3f09b434e21b927e6f4de93d8696ec98879
output_sha256: 8d55733ff6897817b2be81657911e8e8f946e066e6ab0662fccd727a27e9191c
hash_basis: LF-normalized bytes
---

# THM-2840 -- the heptic `(2,2,2,1)` accessory stratum is Chebyshev

**PROVED + VERIFIED-EXACT.**

The first maximal-pole layer beyond the sextic atlas already contains a
classical rigid point.  The most symmetric heptic pole multiset has six
marked normalizations but only one underlying cover: the degree-seven
Chebyshev map with its two finite branch fibres exchanged by a target
Mobius transformation.

## 1. Exact normalized universe

Work over an algebraically closed field of characteristic zero.  In the
balanced response notation of THM-2796 fix

```text
N=7,             e=3,             s=1,             h=4.       (1)
```

Normalize the totally ramified third point to infinity and the ordered
pole passport to

```text
(beta_1,beta_2,beta_3,beta_4)=(0,1,lambda,mu),
(p_1,p_2,p_3,p_4)=(2,2,2,1).                         (2)
```

Put

```text
D=x^2(x-1)^2(x-lambda)^2(x-mu),
Theta=x(x-1)(x-lambda)(x-mu),                        (3)

K=
 2(x-1)(x-lambda)(x-mu)
+2x(x-lambda)(x-mu)
+2x(x-1)(x-mu)
+ x(x-1)(x-lambda),

E=K/7.
```

Then `E` is a monic cubic and

```text
D'=7(D/Theta)E.                                      (4)
```

Divide by the prospective double-zero square:

```text
D=S E^2+R,                 deg S=1,   S monic,

S=x-alpha,
alpha=(4lambda-5mu+4)/7.                            (5)
```

Identity `(4)` makes `R'` divisible by `E`.  Its quotient is linear, so
`R` is constant exactly when its two coefficients vanish.  After removing
nonzero rational scalars these equations are

```text
f_1=
 2lambda^2+2lambda mu-3lambda-3mu^2+2mu+2,

f_0=
 4lambda^2 mu+3lambda^2-5lambda mu^2+6lambda mu
 +3lambda-5mu^2+4mu.                                (6)
```

Thus `(6)` is not merely a necessary critical-value condition: it is
equivalent to the exact factorization `D+A=S E^2` for a constant `A`.

## 2. Radical saturated algebra of length six

Let `I=<f_1,f_0>` in `Q[lambda,mu]`.  Its lexicographic basis is

```text
mu+2lambda^5-5lambda^4-6lambda^3
   +14lambda^2+6lambda-6,

h(lambda)=p(lambda)q(lambda),                        (7)

p=lambda^3-2lambda^2-lambda+1,
q=lambda^3-lambda^2-2lambda+1.
```

More transparently,

```text
I
= <mu-lambda^2+lambda, p>
  intersect
  <mu+lambda^2-lambda-1, q>.                         (8)
```

The exact scalar invariants

```text
Disc(p)=Disc(q)=49,        Res(p,q)=-1,
Disc(h)=2401                                             (9)
```

show that `(7)` is a reduced algebra of length six.

For the actual accessory chamber define

```text
A=-R
 =-4lambda^2 mu^2(4lambda-5mu+4)/343,                (10)

Gamma=
 lambda mu(lambda-1)(mu-1)(lambda-mu)
 Disc_x(E) A Res_x(S,E).                              (11)
```

The factors in `(11)` respectively exclude pole collisions, confluent
double zeros, zero critical value, and collision of the simple zero with a
double zero.  Reduction by `(7)` gives

```text
pole collision:
 -lambda(lambda-2)(lambda-1)(lambda+1)(2lambda-1),

Disc_x(E):
 lambda^2(lambda-1)^2,

A:
 4(2lambda-1)
  (4lambda^4-8lambda^3-2lambda^2+6lambda-3)/343,

Res_x(S,E):
 (lambda-2)(lambda+1)(2lambda-1).                    (12)
```

Every polynomial in `(12)` is coprime to `h`.  The Jacobian of `(6)`
reduces to

```text
-42(2lambda-1)(lambda^2-lambda-1)
   (2lambda^2-2lambda-5),                            (13)
```

which is also coprime to `h`.  Therefore

```text
I:Gamma^infinity=I,                                  (14)
```

and all six admissible accessory points are simple.  This proves the exact
candidate count, rather than counting roots of an eliminant that might lie
on a discarded boundary.

## 3. Full response reconstruction

On `(7)`, equations `(5)--(6)` give exactly

```text
B=D+A=S E^2.                                         (15)
```

Set

```text
C=-7A,
F=B/D,
G=C E/(2D Theta),
V=4 S D Theta^2/C^2.                                (16)
```

Using `(4)` and `(15)` gives coefficientwise

```text
F'=2G,                  F=VG^2,
2VG'+V'G=2.                                           (17)
```

The units in `(12)` ensure that the two finite fibres have the exact
partitions

```text
zero fibre:          (2,2,2,1),
pole fibre:          (2,2,2,1).                      (18)
```

Moreover

```text
F-1=A/D,                                                (19)
```

so `F(infinity)=1` and the third fibre is `(7)`.  The simple factor `S`
also makes the potential genuinely nonsplit in THM-2796's squareclass
sense.

## 4. Exact Chebyshev identification

The additional affine denominator

```text
6mu-2lambda-2                                        (20)
```

reduces modulo `(7)` to

```text
-2(2lambda-1)
  (3lambda^4-6lambda^3-12lambda^2+15lambda+17),
```

which is coprime to `h`.  Hence the affine coordinate

```text
Y=1+7(x-mu)/(6mu-2lambda-2)                          (21)
```

is defined at every accessory point.  Exact reduction in the length-six
algebra gives

```text
A T_7(Y)=2D+A,                                       (22)
Y(mu)=1,                       Y(alpha)=-1.
```

Consequently

```text
F=(D+A)/D
 =(T_7(Y)+1)/(T_7(Y)-1).                             (23)
```

The classical factorizations

```text
T_7(Y)-1=(Y-1)(8Y^3+4Y^2-4Y-1)^2,
T_7(Y)+1=(Y+1)(8Y^3-4Y^2-4Y+1)^2                   (24)
```

explain both copies of `(2,2,2,1)` intrinsically.  The two cubic components
in `(8)` are the two orientation choices among the six affine
normalizations, not two different unmarked maps.

## 5. Marked and unmarked Nielsen count

There is also a self-contained permutation count.  Fix the total
seven-cycle for the fibre `(7)`.  A branch permutation of type
`(2,2,2,1)` is a three-edge matching on seven vertices, so there are `105`
raw matchings.  Exactly `35` are genus zero.  Their pole-passport
histogram is

```text
(4,1,1,1):   7,
(3,2,1,1):  21,
(2,2,2,1):   7.                                     (25)
```

For the last row, labelling the three double-pole cycles contributes `3!`
choices and quotienting by the seven rotations gives

```text
7*3!/7=6                                             (26)
```

marked Nielsen classes.  Forgetting those labels gives one class.  Thus
the six algebraic points in `(7)` exactly saturate the dessin count, and
all are the single unmarked Chebyshev cover `(23)`.

The same result for another placement of the unique simple pole follows
by relabelling the four pole points and affinely renormalizing two of them
to `0,1`.

## 6. Exact control and scope

The companion independently:

1. derives `(4)--(6)` by exact polynomial division;
2. verifies `(7)--(14)`, radicality, every boundary unit, and the
   Jacobian unit;
3. reconstructs `(15)--(19)` coefficientwise;
4. proves the affine Chebyshev identity `(20)--(24)`; and
5. enumerates all `105` matchings, the planar histogram `(25)`, and the
   marked/unmarked counts `(26)`.

It contains no Python `assert` node.  Run

```text
python 04-computation/jc_heptic_e3_2221_chebyshev_accessory_thm2840.py
python -O 04-computation/jc_heptic_e3_2221_chebyshev_accessory_thm2840.py
```

The normal, optimized, and stored transcripts agree byte for byte.

This theorem closes only the highest-symmetry pole multiset in the
`N=7,e=3,h=4` response layer.  The multisets `(4,1,1,1)` and
`(3,2,1,1)` remain open, as do higher layers.  It proves an abstract
one-variable response classification.  It does not supply a second
coordinate, satisfy the Faber/Keller chart-entry equations, or prove
`JC(2)` or `DC(2)`.
