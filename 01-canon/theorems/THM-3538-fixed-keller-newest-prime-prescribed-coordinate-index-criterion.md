---
id: THM-3538
title: "Fixed Keller newest-prime prescribed-coordinate index criterion"
status: >
  PROVED + VERIFIED-EXACT.  At every depth, after the preceding newest-prime
  cover is split, the local indices of the prescribed integral coordinates
  y_n,z_n are exactly the valuations of an explicit product of one-step
  square factors and all cross-block resultants.  The same statement holds
  for u_n=1/x_n when its reciprocal chart is integral.  Hence local
  maximality is equivalent to a squarefree/pairwise-coprime special-fibre
  test, not to primitivity.  The carrier is exactly a unit for n=2,3,4;
  all-level unitness for these literal coordinates remains OPEN.  Literal
  x_n is nonintegral at the newest prime: its monic discriminant exponent is
  3-2*3^n+2i, while its primitive raw-cleared reversal and the reciprocal
  order both have exponent 1+2i.
source: codex/primitive-coordinate-audit/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
related:
  - THM-3532-fixed-keller-conjugacy-covariance-and-two-sided-one-step-boundary
  - THM-3537-fixed-keller-level-two-old-L-inertia-and-x-index
  - THM-3539-fixed-keller-newest-prime-decomposition-centralizer-and-lca-packet-floor
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_prescribed_coordinate_newest_prime_probe_20260816.py
output: 05-knowledge/results/keller_prescribed_coordinate_newest_prime_probe_20260816.out
script_sha256: 7e03534ef3ff1bd05d5108c327b3c40eae929f0a45baa28aa2980a9a9d632892
output_sha256: 5498fcdd2a35f86abb1bf80edea7a39ad5597520d4da5949b56f5e29e3187f93
semantic_sha256: 54e24e2b22edbc0882be21fc1c64f8b451a6c8767af779be52a3e16adcb6db34
hash_basis: LF-normalized bytes
---

# THM-3538 -- prescribed coordinates are maximal exactly when their split blocks do not collide

**PROVED + VERIFIED-EXACT, WITH ALL-LEVEL EQUALITY OPEN.**

Retain the fixed sporadic Keller map `F`, the degree-`3^n` inverse tower
`K_n/K_0`, the normalization `B_n` of

```text
A=Q[a,b,c]
```

in `K_n`, and THM-3530's raw image primes `P_j`, with `P_0=L`.  Fix `n>=1`
and write

```text
p=(P_(n-1)),       R=A_p,       N=3^n,       m=3^(n-1).       (1)
```

THM-3533 proves that the normalization discriminant at `p` has exponent one.
THM-3535 proves that each of `x_n,y_n,z_n`, and hence `u_n=1/x_n`, is a
primitive field generator.  The question left open by those theorems is
whether these particular primitive observations generate the maximal local
integral order.

## 1. The split newest-prime blocks

Let `R^sh` be a strict henselization of `R`.  The preceding normalization is
finite etale at `p`, so

```text
B_(n-1) tensor_R R^sh = product_(i=0)^(m-1) R^sh.        (2)
```

Write `q_i=(a_i,b_i,c_i)` for the resulting sections.  The degree-one
ancestry `V(L)->V(P_(n-1))` selects a unique section `q_0`, and THM-3533 gives

```text
v(L(q_0))=1,                 v(L(q_i))=0 for i>0.        (3)
```

The last inverse step over each section is one cubic block.  For the two
globally integral coordinates use the monic cubics

```text
f_(y,i)(Y) = (1/2) P_(r,q_i)(b_i-Y),
f_(z,i)(Z) = (1/8) Q_(z,q_i)(Z),                         (4)
```

where `P_r,Q_z` are the exact constant-leading cubics of THM-2546.  For the
reciprocal coordinate put

```text
T_i=4-3b_i c_i,
S_i=27a_i c_i^2-9b_i c_i+8,
f_(u,i)(U)=-(L(q_i)+T_i U^2-2c_i U^3)/(2c_i).           (5)
```

Equation `(5)` is a monic integral block exactly on the reciprocal-chart
gate

```text
G_(u,n):       every c_i is a unit in R^sh.              (6)
```

Equivalently, the norm of the preceding `c`-coordinate is a `p`-unit.  This
gate cannot be dropped: if an unramified block has `c_i=0` in the residue
field, its x-cubic has a zero root and `u=1/x` acquires a pole.

## 2. Exact one-block discriminants

The exact coordinate discriminants factor as

```text
Disc(f_(theta,i)) = -L(q_i) h_theta(q_i)^2,              (7)
```

where

```text
h_y(a,b,c) = 27a/2,
h_z(a,b,c) = 27M(a,b,c)/32,
h_u(a,b,c) = (27ac^2-9bc+8)/(2c^2).                     (8)
```

Here the explicit z-square factor is

```text
M = 729a^3b^3c^4 + 4374a^3b^2c^3 + 8748a^3bc^2 + 5832a^3c
    -729a^2b^4c^3 - 2754a^2b^3c^2 - 2268a^2b^2c + 648a^2b
    +27ab^6c^3 + 297ab^5c^2 + 360ab^4c - 268ab^3
    -9b^7c^2 - 28b^6c + 28b^5.                         (9)
```

For `y`, `(7)` is THM-2546's
`Disc(P_r)=-4(27a)^2L` divided by the fourth power of the leading
coefficient `2`.  For `z`, its raw discriminant
`-4(27M)^2L` is divided by `8^4`.  For `u`, THM-3533's raw reciprocal
identity `Disc(L+TU^2-2cU^3)=-4LS^2` is divided by `(-2c)^4`.  Thus the
three square factors in `(8)` include every possible internal block defect;
they are not decorative units.

## 3. The all-level carrier and index formula

For `theta=y,z`, and for `theta=u` under `(6)`, define after choosing the
labelling in `(2)`

```text
C_(theta,n)
 = product_i h_theta(q_i)
   product_(i<j) Res(f_(theta,i),f_(theta,j)).            (10)
```

Changing the split labelling changes `(10)` at most by a unit/sign, so its
valuation is intrinsic.  Since `theta_n` is primitive, its monic minimal
polynomial becomes

```text
m_(theta,n) = product_i f_(theta,i)                      (11)
```

over `R^sh`.  The product-discriminant identity gives

```text
Disc(product_i f_i)
 = product_i Disc(f_i) product_(i<j) Res(f_i,f_j)^2.     (12)
```

Equations `(3)`, `(7)`, and `(12)` therefore imply the exact all-level
formula

```text
v_p(Disc(m_(theta,n))) = 1+2v_p(C_(theta,n)).            (13)
```

Let

```text
i_(theta,n)=length_R(B_(n,p)/R[theta_n]).                (14)
```

THM-3533 and the order-index identity compare `(13)` with the normalization
exponent one and give

```text
i_(theta,n)=v_p(C_(theta,n)).                            (15)
```

This is the promised exact reduction.  It retains both sources of a
coordinate defect:

1. `h_theta(q_i)` detects an internal collision within one cubic block;
2. `Res(f_(theta,i),f_(theta,j))` detects two different predecessor blocks
   with a common residue value.

Full primitivity supplies the degree in `(11)` but says nothing about either
factor in `(10)`.

## 4. Equivalent local-maximality tests

Over the separably closed residue field of `R^sh`, the following are
equivalent:

```text
(a) i_(theta,n)=0;
(b) every h_theta(q_i) and every cross-block resultant in (10) is a unit;
(c) every i>0 block is squarefree modulo p, the q_0 block has only its
    forced double shadow, and all distinct block reductions are coprime;
(d) gcd(mbar_(theta,n), d mbar_(theta,n)/dT) is exactly the one linear
    shadow factor forced by q_0.                                  (16)
```

For `u`, all four statements include the chart gate `(6)`.  The baseline
linear gcd in `(d)` is load-bearing.  If a further block collides with the
same shadow, its multiplicity rises and the gcd degree rises as well; hence
the degree-one test does not accidentally hide a collision at the already
repeated root.

Thus `y_n,z_n` generate locally maximal integral orders exactly when `(16)`
passes.  The theorem proves this equivalence for every `n`; it does not
assert that the test passes for every `n`.

## 5. The literal x-coordinate: reciprocal order and raw clearing

Assume `(6)` and put `I=i_(u,n)`.  Viete's formula on the split x-blocks
gives

```text
Norm(x_n)=product_i (2c_i/L(q_i)),       v_p(Norm(x_n))=-1.       (17)
```

Hence `u_n` is integral but `x_n` is not.  Let `m_u(U)` and `m_x(X)` be the
monic minimal polynomials of `u_n` and `x_n`, respectively, and define the
primitive raw-cleared reversal

```text
E_(x,n)(X)=X^N m_u(X^(-1)).                              (18)
```

Its constant coefficient is one, its leading coefficient has valuation one,
and

```text
E_(x,n)=m_u(0)m_x,                  v_p(m_u(0))=1.       (19)
```

Reciprocal root differences and scalar discriminant covariance now give

```text
v_p(Disc(m_u))       = 1+2I,
v_p(Disc(m_x))       = 3-2N+2I,
v_p(Disc(E_(x,n)))   = 1+2I.                            (20)
```

Therefore “newest-prime exponent one for x” is correct only for the integral
reciprocal order or for the primitive nonmonic reversal `(18)`.  The literal
monic x-polynomial has a negative exponent already when `I=0`, and there is
no finite integral `R[x_n]` order whose index could be read from it.

## 6. Fixed levels where the carrier is exactly a unit

The criterion `(16)` is closed here only at newest levels two, three, and
four.  No finite row is promoted to an induction.

### Level two: characteristic-zero rank nine

Use the smooth ramified-block point

```text
q_0=(2/27,1,1) in V(L),
(h_y,h_z,h_u)=(1,49/32,1/2),                            (21)
```

and put

```text
eta=F(q_0)=(113651/19683,14399/6561,2584/19683).        (22)
```

The preceding x-fibre factors exactly as

```text
(27X-2)(118818549X^2+8801374X+5651208)/43046721.        (23)
```

Multiplying the three last-step coordinate blocks gives degree-nine
special-fibre polynomials whose derivative gcds are exactly

```text
y_2: Y-1/3,                  z_2: Z,                  u_2: U.     (24)
```

The other degree-six factor is squarefree and coprime to the ramified block
in all three rows.  Hence the rank-nine product orders have index zero.

### Level three: exact good reduction at 53

Over `F_53`, the smooth point and its image are

```text
q_0=(2,23,4),                 F^2(q_0)=(4,38,38).       (25)
```

The preceding fibre splits into nine distinct blocks, exactly one on `L=0`.
The total degree-27 polynomials have derivative gcds

```text
Y-7,                         Z,                        U,          (26)
```

and all distinct blocks are pairwise coprime.  The quotient-algebra route
agrees with the complete split product.

### Level four: exact good reduction at 41

Over `F_41`, use

```text
q_0=(7,26,3),                 F^3(q_0)=(15,11,7).       (27)
```

Nine rational parent blocks followed by exact cubic quotient norms cover all
`27` algebraic predecessor blocks.  Their coordinate products all have
degree `81`; the derivative gcds are exactly

```text
Y-18,                        Z,                        U.          (28)
```

Every inverse denominator is a unit, all nine parent norms have degree nine,
and exactly one algebraic predecessor block lies on `L=0`.

The witnesses `(25)`--`(28)` are exact nonidentity certificates, not
floating-point evidence.  Each `q_0` is a smooth point of the integral
`L`-model.  Hensel lifting therefore gives a characteristic-zero p-adic point
on `L`; `det J_F=-2`, the squarefree inverse equations, and the explicit
unit-denominator gates lift the split or quotient finite-etale tree with it.
The displayed nonzero resultants remain units after lifting.  Thus the
characteristic-zero carrier is nonzero on the relevant image divisor and is
a unit at its generic point.  An independent split-quotient crosscheck at
level two over `F_31` gives the same verdict.

Consequently, for `2<=n<=4`,

```text
i_(y,n)=i_(z,n)=i_(u,n)=0,
v_p Disc(m_y)=v_p Disc(m_z)=v_p Disc(m_u)=v_p Disc(E_x)=1.        (29)
```

The monic literal-x exponents in these three rows are, respectively,

```text
-15,                              -51,                              -159.
                                                                    (30)
```

## 7. The old-L hostile and exact boundary

THM-3537 is the canonical hostile to replacing `(16)` by primitivity.  On
the transverse old-`L` depth-two DVR it proves

```text
normalization exponent 4,       inertia (4)(2)(1)^3,
literal x_2 order exponent 8,   index length 2.                    (31)
```

The quartic and quadratic Newton packets both reduce to `x_2=0`; their
cross-packet resultant contributes the extra four in the order
discriminant.  This does not conflict with `(29)`.  At an old prime the
preceding cover is already ramified/nonproper, so the finite-etale split
hypothesis `(2)` and the unique last-step `L=0` model `(3)` are unavailable.
The old-`L` example instead demonstrates exactly why cross-block collisions
must be audited.

The following remain **OPEN**:

1. whether `C_(y,n),C_(z,n),C_(u,n)` are units for every `n>=5`;
2. an induction or recurrence controlling those cross-block resultants;
3. analogous prescribed-coordinate indices at old primes beyond `(31)`;
4. transport of the literal standard-coordinate statement under nonlinear
   conjugacy--THM-3532 transports the observation with its chart, not the
   standard linear family unchanged.

THM-3539 supplies the exact group-theoretic sidecar to item 2.  At the newest
prime it proves only `I<=D<=C_(W_n)(I)`.  If the image of `D` on predecessor
blocks has the full marked-leaf point-and-pair orbits, the factors in `(10)`
collapse to exactly `n^2` LCA valuation packets.  That decomposition-group
saturation is not implied by full global wreath monodromy and remains open;
without it the raw exponential factor census can persist.

The numerical old-`L` sequence `1,4,14,46,142` is discovery data only.  In
particular, the formula fitting its first four rows predicts `146` at depth
five and is refuted.  Nothing here proves an arbitrary-Keller-map
classification, a Jacobian-conjecture counterexample family theorem,
`JC(2)`, `DC(2)`, LRC(14), or a tournament/tree identification.

## 8. Exact companion

Reproduce with

```text
python -B 04-computation/keller_prescribed_coordinate_newest_prime_probe_20260816.py
python -B -O 04-computation/keller_prescribed_coordinate_newest_prime_probe_20260816.py
```

The ordinary and optimized transcripts match the stored output after the
declared LF normalization.  The
companion uses explicit failures rather than executable assertions, checks
the rational rank-nine product directly, independently crosschecks split and
quotient-algebra products, and records every finite chart, degree, unique
`L`-block, and derivative-gcd gate.

**QED.**
