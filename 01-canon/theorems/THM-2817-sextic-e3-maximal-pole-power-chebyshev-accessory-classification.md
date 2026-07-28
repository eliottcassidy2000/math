---
id: THM-2817
title: "Sextic e=3 maximal-pole accessory atlas: cubic power and Chebyshev carriers"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The minimal e=3 maximal-pole response chamber N=6,h=4 is
  completely explicit.  For every ordered positive four-part pole
  partition, the saturated two-variable equal-critical-value ideal is
  radical of length two; all ten ideals are listed.  Their twenty marked
  charts agree with an independent noncrossing-matching census and reduce,
  after forgetting equal-part labels and affine source normalization, to
  exactly two maps: a cubic-power carrier for (3,1,1,1) and a Chebyshev T_3
  carrier for (2,2,1,1).  This closes one response chamber, not
  Keller-chart entry, Faber-flux compatibility, JC(2), or DC(2).
source: root/sextic-e3-maximal-pole-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2808-three-pole-e2-maxwell-polynomial-and-finite-accessory-classification
  - THM-2816-maximal-pole-clean-nielsen-ribbon-tree-prufer-vandermonde-count
script: 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
output: 05-knowledge/results/jc_sextic_e3_maximal_pole_accessory_thm2817.out
script_sha256: 4d2c26143ef720013e183453dc8ae4234af5c2604bb899bc1f2276bba8ed8d6b
output_sha256: 52037f5fdd25ec270f678d73edbe180e4509245955e095f49ba8f18828ea137a
hash_basis: LF-normalized bytes
---

# THM-2817 -- sextic `e=3` maximal-pole accessory atlas

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The first chamber beyond the now-complete `e=2` layer has two accessory
coordinates.  In its smallest degree those coordinates do not produce a
new positive-dimensional family.  They form ten radical schemes of length
two, and the apparent diversity is only the marked form of two highly
symmetric covers.

## 1. The normalized accessory problem

Work over an algebraically closed field of characteristic zero.  Fix

```text
N=6,                    e=3,                    h=4.    (1)
```

Let

```text
p=(a,b,c,d),            a,b,c,d>=1,
a+b+c+d=6.                                        (2)
```

Order the four pole points and use the remaining affine source freedom to
put the first two at `0,1`.  Write the other two as `lambda,mu` and put

```text
D=x^a(x-1)^b(x-lambda)^c(x-mu)^d,
T=x(x-1)(x-lambda)(x-mu),
P=D/T,

E=D'/(6P).                                           (3)
```

The last expression is a monic cubic.  Divide

```text
D=H E+r_2 x^2+r_1 x+v.                              (4)
```

The three roots of `E` are the non-pole critical points of `D`.  When they
are distinct, they have one common critical value exactly when

```text
r_2=r_1=0.                                          (5)
```

The ordered normalization is rigid: an affine source map fixing the two
ordered points `0,1` is the identity.  This clause is necessary before a
finite accessory count can be stated.

## 2. Saturation and the ten exact ideals

The raw equations `(5)` include collisions.  Remove precisely the declared
boundary with

```text
Omega
 =lambda mu(lambda-1)(mu-1)(lambda-mu)
  Disc_x(E) v.                                       (6)
```

Define

```text
J_p=<r_2,r_1>:Omega^infinity
       in Q[lambda,mu].                              (7)
```

The following table gives a Groebner basis of every `J_p`.  Each row means
the displayed linear equation together with the displayed quadratic.

| ordered `p` | linear equation | quadratic in `lambda` |
|---|---|---|
| `(1,1,1,3)` | `3mu-lambda-1` | `lambda^2-lambda+1` |
| `(1,1,2,2)` | `mu+lambda-1` | `(4lambda-3)(4lambda-1)` |
| `(1,1,3,1)` | `mu-3lambda+1` | `3lambda^2-3lambda+1` |
| `(1,2,1,2)` | `mu-lambda+1` | `(lambda-4)(3lambda-4)` |
| `(1,2,2,1)` | `mu-lambda-1` | `(lambda-3)(3lambda-1)` |
| `(1,3,1,1)` | `mu+lambda-3` | `lambda^2-3lambda+3` |
| `(2,1,1,2)` | `mu-lambda-1` | `(lambda+3)(3lambda+1)` |
| `(2,1,2,1)` | `mu-lambda+1` | `(lambda+2)(3lambda-2)` |
| `(2,2,1,1)` | `mu+lambda-1` | `(2lambda-3)(2lambda+1)` |
| `(3,1,1,1)` | `mu+lambda+1` | `lambda^2+lambda+1` |

Every quadratic is separable, and `Omega` is a unit in every displayed
quotient.  Consequently

```text
J_p is radical of length two for every ordered p.     (8)
```

There are exactly two admissible normalized charts per passport and twenty
over the ten ordered passports.

### Exact finite proof of the table

Equation `(3)` constructs `E` without roots or resultants.  For each of the
ten compositions `(2)`, introduce an inverse variable `u` and compute the
lexicographic elimination ideal

```text
<r_2,r_1,1-u Omega> intersect Q[lambda,mu].            (9)
```

Exact polynomial reduction gives precisely the corresponding two
generators in the table.  Both containments are checked: every elimination
generator reduces to zero by the table row, and both table generators
reduce to zero by `(9)`.  Coprimality of `Omega` with the displayed
quadratic proves that no admissible point was discarded.

This is a finite exact symbolic proof, not a floating-point root census.
The companion reproduces all ninety ideal, radical, boundary, derivative,
and reconstruction gates.  Independent hostile audit is still required
before promotion.

## 3. Exact response reconstruction

At a root of `J_p`, equation `(4)` becomes

```text
D=H E+v.                                             (10)
```

Every root of `E` is a critical point of `D`, so it is a double root of
`D-v`.  The boundary `(6)` makes the three roots distinct.  Hence
`E^2` divides `D-v`; both are monic of degree six, and therefore

```text
D-v=E^2,                    v!=0.                    (11)
```

Conversely `(11)` implies the equal-critical-value equations `(5)`.  Thus
the table is already the complete ordered accessory classification.

Define

```text
F=E^2/D,
G=3v E/(D T),
V=D T^2/(9v^2).                                      (12)
```

From `(3)` and `(11)`,

```text
D'=6DE/T,
F=1-v/D,
F'=6vE/(DT)=2G,
F=VG^2,
2VG'+V'G=2.                                          (13)
```

The zero fibre has type `2^3`, the pole fibre has the ordered partition
`p`, and the third fibre is totally ramified at infinity.  All squarefree,
disjointness, nonzero-value, and pole-collision gates are exactly the
factors retained in `(6)`.

## 4. Only two unmarked carriers remain

For pole multiset `(3,1,1,1)`, choose a cubic root of unity
`omega!=1`.  The last row of the table gives the symmetric representative

```text
D=x^3(x^3-1),
E=x^3-1/2,
v=-1/4,
E^2=D-v.                                             (14)
```

The two accessory roots exchange `omega` and `omega^2`; they are the two
orientations of the same equal-part-labelled configuration.

For pole multiset `(2,2,1,1)`, use the centered coordinate `y`.  The
symmetric representative is

```text
D=(y^2-1/4)^2(y^2-1),
E=y^3-(3/4)y=T_3(y)/4,
v=-1/16,
E^2=D-v.                                             (15)
```

For example, normalizing its two double poles to `0,1` puts the two simple
poles at `-1/2,3/2`, exactly the two roots in the `(2,2,1,1)` row.

Affine renormalization and ordering of the four poles sends `(14)` or
`(15)` to every row of the corresponding pole multiset.  Since `(8)` gives
only two marked points per row, these constructions exhaust them.  Thus:

```text
twenty ordered charts,
two unmarked affine/target maps,
one map for each unordered pole multiset.             (16)
```

The `x^3` and `T_3` forms explain why this first ternary chamber is so
small.  They are two different symmetry mechanisms for forcing three
critical values to coincide; the conclusion is not that power and
Chebyshev carriers are generally interchangeable.

## 5. Independent Nielsen agreement

The order-two inertia in degree six is a perfect matching on a hexagon.
Exactly five of the fifteen matchings are noncrossing.  Compose each with
the fixed oriented six-cycle, label its four complementary cycles, and
quotient the labelled state by the six rotations.

The exact census has

```text
20 labelled rotation classes,
10 ordered pole profiles,
2 classes for every profile.                          (17)
```

This agrees with every algebraic length in `(8)` and with the specialization

```text
(e-1)! binom(N-e-1,e-1)=2
```

of the ribbon-tree candidate THM-2816.  The permutation census is included
as a separate control; this theorem's algebraic completeness does not rely
on THM-2816 being promoted first.

## 6. What changes at the next degree

The minimal equality `N=2e` forces the simple factor `S` in
`B=S E^2` to be constant, which is why the accessory problem becomes the
perfect-square identity `(11)`.  At `N=7,e=3,h=4`, `S` is linear and the
same Nielsen count rises to six per ordered passport.  The new coordinate
is not another pole: it is the simple zero of the numerator.  Retaining
that zero while deriving the multivariate Maxwell ideal is the next
algebraic test.

This identifies a precise recursion:

```text
minimal chamber:       D-v=E^2;
next chamber:          D-v=S E^2,       deg S=1;
general maximal pole:  D-v=S E^2,       deg S=N-2e.    (18)
```

The ribbon tree predicts the number of solutions; `(18)` is the accessory
coordinate problem that must realize them.

## 7. Exact controls

Run

```bash
python 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
python -O 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
```

Both executions byte-match

```text
05-knowledge/results/jc_sextic_e3_maximal_pole_accessory_thm2817.out.
```

The companion uses exact rational polynomial and permutation arithmetic,
contains no Python `assert` node, verifies both saturation-ideal
containments, checks every boundary unit and response reconstruction, and
independently enumerates the twenty marked Nielsen charts.  The finite
universe is exactly the ten positive ordered compositions of six into four
parts.

## 8. Scope and failure boundaries

This theorem closes the first accessory chamber with `e>=3`.  It does not
classify `N>=7,e=3`, any nonmaximal pole count, or the full `e>=3` response
layer.  It does not prove that a response enters a Keller chart or satisfies
the inherited Faber flux equations, and it proves neither `JC(2)` nor
`DC(2)`.

The saturation is load-bearing.  Keeping a pole collision, a multiple
critical point, or `v=0` creates false solutions.  Forgetting the ordered
source normalization also turns each isolated point into an affine orbit.
Finally, the two marked solutions in a row must not be divided by a
stabilizer order unless the actual equal-part label action is computed.
