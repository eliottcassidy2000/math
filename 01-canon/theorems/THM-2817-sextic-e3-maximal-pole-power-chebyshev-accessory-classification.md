---
id: THM-2817
title: "Sextic e=3 maximal-pole accessory atlas: cubic power and Chebyshev carriers"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The minimal
  N=6,e=3,h=4 maximal-pole response chamber has ten
  ordered pole passports and exactly two simple accessory points for each.
  Saturating the two equal-critical-value remainder equations gives an
  explicit linear-plus-quadratic ideal in every passport.  Up to unmarked
  pole permutation, the only carriers are the cubic power map and the
  cubic Chebyshev map, both followed by z^2/(z^2-1).  This is an abstract
  response classification, not Keller-chart entry, JC(2), or DC(2).
source: root/sextic-e3-maximal-pole-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2808-three-pole-e2-maxwell-polynomial-and-finite-accessory-classification
  - THM-2816-maximal-pole-clean-nielsen-ribbon-tree-prufer-vandermonde-count
script: 04-computation/jc_sextic_e3_power_chebyshev_accessory_thm2817.py
output: 05-knowledge/results/jc_sextic_e3_power_chebyshev_accessory_thm2817.out
script_sha256: 3745f150524fdd2bd068c486de57e9b66469eb3fd5eda4ffba0ed582097d414d
output_sha256: bca8f6825e92bb98913e6f330dd462a15850e0fb24b0fbfffbe08a9cd8c82343
secondary_script: 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
secondary_output: 05-knowledge/results/jc_sextic_e3_maximal_pole_accessory_thm2817.out
secondary_script_sha256: 4d2c26143ef720013e183453dc8ae4234af5c2604bb899bc1f2276bba8ed8d6b
secondary_output_sha256: 52037f5fdd25ec270f678d73edbe180e4509245955e095f49ba8f18828ea137a
independent_script: 04-computation/jc_sextic_e3_maximal_pole_accessory_independent_audit_thm2817.py
independent_output: 05-knowledge/results/jc_sextic_e3_maximal_pole_accessory_independent_audit_thm2817.out
independent_script_sha256: 8462778c730d03513506d106b54966cdfb186522a1d90546061fa57d1c014f42
independent_output_sha256: c3e9ce10c7ea93af1ed1f6201d12bc378e719aaf8712625fe3d4177706f480db
hash_basis: LF-normalized bytes
---

# THM-2817 -- sextic e=3 maximal-pole accessory atlas: cubic power and Chebyshev carriers

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The first genuinely ternary maximal-pole layer is already rigid enough to
solve completely.  Its two unordered answers are not accidental isolated
maps: they are the two classical cubic geometries, power and Chebyshev,
whose critical values are identified by the same even degree-two outer
response.

## 1. The two-variable equal-critical-value ideal

Fix

```text
N=6,                 e=3,                 h=4,
(a,b,c,d)>0,         a+b+c+d=6.                         (1)
```

Normalize the ordered poles to `0,1,lambda,mu`, the totally ramified third
point to infinity, and `F(infinity)=1`.  Put

```text
D=x^a(x-1)^b(x-lambda)^c(x-mu)^d,
T=x(x-1)(x-lambda)(x-mu).                            (2)
```

The non-pole critical polynomial is the cubic

```text
K=
 a(x-1)(x-lambda)(x-mu)
 +b x(x-lambda)(x-mu)
 +c x(x-1)(x-mu)
 +d x(x-1)(x-lambda),                                (3)

D'=D K/T,                  E=K/6.
```

Write

```text
D = K H+R_2 x^2+R_1 x+r_0.                           (4)
```

The three roots of `K` have one common critical value exactly when

```text
R_2=R_1=0.                                           (5)
```

Raw equations `(5)` contain colliding poles, confluent critical points,
and zero critical value.  Let

```text
Gamma=
 lambda mu(lambda-1)(mu-1)(lambda-mu)
 Disc_x(K)(-r_0).                                    (6)
```

The genuine accessory ideal is

```text
I_(a,b,c,d)
 =<R_2,R_1>:Gamma^infinity.                           (7)
```

Equivalently, adjoin `u Gamma-1` and eliminate `u`.

## 2. Complete saturated atlas

For every ordered positive four-part partition of six, `(7)` is radical
and has the following lexicographic basis.  Each displayed quadratic is
monic up to the harmless scalar shown.

```text
(a,b,c,d)   linear relation                 quadratic in lambda

(1,1,1,3)   3mu-lambda-1                    lambda^2-lambda+1
(1,1,2,2)   lambda+mu-1                     (4lambda-3)(4lambda-1)/16
(1,1,3,1)   mu-3lambda+1                    (3lambda^2-3lambda+1)/3
(1,2,1,2)   mu-lambda+1                     (lambda-4)(3lambda-4)/3
(1,2,2,1)   mu-lambda-1                     (lambda-3)(3lambda-1)/3
(1,3,1,1)   lambda+mu-3                     lambda^2-3lambda+3
(2,1,1,2)   mu-lambda-1                     (lambda+3)(3lambda+1)/3
(2,1,2,1)   mu-lambda+1                     (lambda+2)(3lambda-2)/3
(2,2,1,1)   lambda+mu-1                     (2lambda-3)(2lambda+1)/4
(3,1,1,1)   lambda+mu+1                     lambda^2+lambda+1.
                                                               (8)
```

Every quadratic in `(8)` is squarefree, and `Gamma` is a unit in its
quadratic residue algebra.  Thus each ordered passport has exactly two
simple admissible accessory points.

At either point put

```text
A=-r_0,                B=D+A,                C=-6A.   (9)
```

Equations `(4)--(8)` give

```text
B=E^2.                                                (10)
```

There is no simple-zero factor because `N=2e`.  With

```text
F=B/D,
G=C E/(2DT),
V=4 D T^2/C^2,                                       (11)
```

one obtains exactly

```text
F'=2G,                    F=VG^2,
2VG'+V'G=2.                                         (12)
```

Therefore `(8)` is a complete algebraic response atlas, not merely a
necessary eliminant.

## 3. Independent Nielsen saturation

The zero inertia is a noncrossing matching of three chords on the fixed
six-cycle.  There are five raw Catalan matchings.  Labelling the four pole
cycles and quotienting by rotation gives:

```text
ten ordered pole passports,
two Nielsen classes for every ordered passport.       (13)
```

The companion reconstructs `(13)` directly from permutations, independently
of the Groebner calculation.  Thus the two algebraic points in each row of
`(8)` exactly saturate the dessin count.  No accessory branch is missing.

This is the `e=3,N=6` specialization of the candidate general ribbon-tree
formula in THM-2816, but the present proof and exact permutation census are
self-contained and do not use that candidate as a dependency.

## 4. The two unmarked cubic carriers

There are only two unordered pole multisets.

### 4.1 Pole multiset `(3,1,1,1)`: cubic power

Use the last row of `(8)`.  Then

```text
lambda^2+lambda+1=0,          mu=-1-lambda,

(x-1)(x-lambda)(x-mu)=x^3-1.                         (14)
```

Consequently

```text
D=x^3(x^3-1),
E=x^3-1/2,
A=1/4.                                               (15)
```

With

```text
P_3=2x^3-1,
```

the response map is

```text
F=P_3^2/(P_3^2-1).                                   (16)
```

The two ordered accessory points exchange the two labelled simple poles
and give one unmarked map.

### 4.2 Pole multiset `(2,2,1,1)`: cubic Chebyshev

Use `(lambda,mu)=(3/2,-1/2)` in the ninth row of `(8)` and center
`y=x-1/2`.  Then

```text
D=(y^2-1/4)^2(y^2-1),
E=T_3(y)/4,
A=1/16,                                               (17)

T_3(y)=4y^3-3y.
```

Thus

```text
F=T_3(y)^2/[T_3(y)^2-1].                              (18)
```

Again the other accessory point swaps the two labelled simple poles, so
there is one unmarked map.

All other rows of `(8)` are the pole-relabelled affine normalizations of
`(16)` or `(18)`.

## 5. Why two and three are special here

Both sextic responses have the same exact decomposition

```text
degree 6 = degree 3 followed by degree 2,

F=phi(P_3),                  phi(z)=z^2/(z^2-1),       (19)
```

but their cubic inner maps are the power carrier and the Chebyshev carrier.
This is a rigorous co-occurrence of binary and ternary operations on one
map.  It is **composition**, not the free product
`C_2*C_3=PSL_2(Z)`: the intermediate cubic coordinate is a retained
sidecar, and forgetting it loses the distinction between the two carriers.

Both maps are nonsplit in THM-2796's sense because their pole partitions
contain odd parts.  They remain decomposable as rational covers; those are
different predicates.

## 6. Exact controls

The primary and secondary companions use exact rational Groebner,
polynomial, and permutation arithmetic to:

1. build `(3)--(7)` for all ten ordered passports;
2. verify the complete saturated bases `(8)`, radicality, and every gate;
3. check `B=E^2` and the response derivative `(12)`;
4. reconstruct all ten passports and two Nielsen classes per passport;
5. prove the power and Chebyshev identities `(14)--(18)`; and
6. verify normal, optimized, and stored transcript identity.

The hostile audit is methodologically separate: it factors the Sylvester
projection resultant, inspects every irreducible branch in its exact
quotient field, proves that every branch outside `(8)` is killed by
`Gamma`, and proves that `Gamma` and the Jacobian are units on the surviving
branches.  It then reconstructs the response coefficientwise, normalizes
every one of the ten rows directly to the power or Chebyshev carrier, and
independently enumerates the fifteen perfect matchings and twenty marked
Nielsen charts.  It imports neither other companion.

The companions have no Python `assert` node.  Run

```text
python 04-computation/jc_sextic_e3_power_chebyshev_accessory_thm2817.py
python -O 04-computation/jc_sextic_e3_power_chebyshev_accessory_thm2817.py
python 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
python -O 04-computation/jc_sextic_e3_maximal_pole_accessory_thm2817.py
python 04-computation/jc_sextic_e3_maximal_pole_accessory_independent_audit_thm2817.py
python -O 04-computation/jc_sextic_e3_maximal_pole_accessory_independent_audit_thm2817.py
```

The finite computation is exhaustive in the stated sextic universe.

## 7. Scope and boundaries

This is an abstract response classification.  It does not prove that either
map enters a Keller chart, satisfies the inherited Faber flux equations, or
comes from a Weyl-algebra endomorphism.  It proves neither `JC(2)` nor
`DC(2)`.

Omitting the saturation in `(7)` admits pole collisions, confluent critical
points, and zero critical value.  Forgetting pole labels before normalizing
collapses the two ordered points in each row.  Finally, the power/Chebyshev
classification is special to the minimal `N=2e=6` layer; no such two-map
claim is made for higher `N` or `e`.
