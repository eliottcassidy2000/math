---
id: THM-3913
title: "Moving-triple-root one-place elliptic decic normal nonmonogenic cubic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  coefficient-depth-three binary cubic
  (AU+CV)^3+C(AU+CV)U^2+A^2V^3 defines a normal, globally nonmonogenic S3
  cubic order.  Its discriminant is an absolutely irreducible decic with one
  smooth projective infinity place, and its quadratic-resolvent affine class
  group has nonzero 3-torsion.  The projective normalization is the elliptic
  curve y^2=w(w^2-1), however, so the affine branch is not polynomially
  uniruled and the cubic etale surface admits no plane atlas.  This pays the
  order, confluence, and three-class invoices at the first moving-triple-root
  depth but is not a Jacobian counterexample; rational depth-three
  deformations and JC(2) remain open.
source: root / first post-THM-3908 coefficient-depth-three construction, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS after two proof-completeness clarifications.
  The C=0 irreducibility specialization now records that the U^3 leading
  coefficient A(A^2+C) is not divisible by C, so neither primitive factor
  can lose degree on specialization.  The no-atlas step now explicitly uses
  the deleted-ramification valuation for containment in the nonproperness set
  and purity plus irreducibility for the full-component conclusion.  The
  audit independently rechecked Delone--Faddeev normality, absolute
  irreducibility at the unique smooth infinity point, the elliptic inverse
  and pole ledger, boundary units, the nontrivial Kummer class, and the
  polynomial-uniruledness contradiction.  Normal and optimized runs
  byte-match the frozen output; raw hashes and documentation pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3908-quadratic-depth-common-zero-one-point-sextic-two-place-obstruction
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3912-even-one-place-split-boundary-a2-three-torsion-design-sieve
script: 04-computation/jc2_moving_triple_root_one_place_elliptic_decic_thm3913.py
output: 05-knowledge/results/jc2_moving_triple_root_one_place_elliptic_decic_thm3913.out
script_sha256: 49503a5400ad6ac08ea824e0a38548185d64ef7fe60b1ee5caddded735530efa
output_sha256: e677adc21fddb894030a526d876fd097b58ea1ba880b74fe6455ea91b3a074d9
semantic_sha256: 8a591b8a768f89173c9f426e28fc5467fb5e5d6ff87525fca3b10035c99d1805
hash_basis: raw LF bytes
---

# THM-3913 -- the first moving triple root pays every invoice except rationality

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Put `R=k[A,C]`
and consider the binary cubic

```text
ell=AU+CV,
Phi=ell^3+C ell U^2+A^2V^3.                                (1)
```

Let `T/R` be its Delone--Faddeev cubic algebra.  Then:

1. `T` is a normal finite-flat domain of rank three, is globally
   nonmonogenic, and has generic Galois group `S3`;
2. its discriminant curve is absolutely irreducible of degree ten and has
   exactly one projective infinity point, with one smooth normalization
   place;
3. for the affine quadratic resolvent

   ```text
   Q=Spec R[W]/(W^2-Delta),                                  (2)
   ```

   one has `Q*=k*` and `Cl(Q)[3]!=0`;
4. the normalization of the projective discriminant curve is elliptic, and
   the corresponding maximal etale cubic surface admits no polynomial-plane
   atlas.

Thus coefficient depth three crosses the complete quadratic-depth barrier of
THM-3908 and realizes the one-place/nonmonogenic/normal-`S3`/three-class
packet simultaneously.  It fails at polynomial uniruledness, not at the
quadratic-resolvent or cubic-order layers.

## 1. A sparse normal nonmonogenic cubic order

Expanding `(1)` gives coefficients

```text
a=A^3+AC,
b=3A^2C+C^2,
c=3AC^2,
d=C^3+A^2.                                                (3)
```

They have gcd one, so `Phi` is primitive, but they all belong to the proper
ideal `(A,C)`.  Every value of the binary index form therefore belongs to
`(A,C)` and cannot be a scalar unit.  By the intrinsic index criterion of
THM-3801, `T` is globally nonmonogenic.

The generic binary cubic is irreducible.  Indeed, after dehomogenizing with
`V=1`, a factorization over `k(A,C)` would, by Gauss's lemma, give primitive
positive-degree factors over `k[A,C]`.  Their leading coefficients multiply
to

```text
A(A^2+C),                                                   (3a)
```

which is not divisible by `C`.  Hence neither factor loses degree after
setting `C=0` (and neither specialization is zero).  But

```text
Phi|_(C=0)=A^2(AU^3+V^3),                                  (4)
```

and `AU^3+V^3` is irreducible over `k(A)`: a root would make `-1/A`
a cube, while its `A`-valuation is `-1`.  This contradiction proves generic
irreducibility.

The exact discriminant is

```text
Delta=-27A^10-54A^8C-27A^6C^2
      -36A^4C^5-4A^2C^6-4C^9.                             (5)
```

Section 2 proves that `(5)` is absolutely irreducible.  Hence it has
height-one valuation one.  The Delone--Faddeev algebra is finite free and
therefore `S2`; away from `(Delta)` it is etale, and at `(Delta)` the
square-index discriminant formula forces index zero in the maximal DVR
order.  Thus it is `R1`, hence normal.  Its generic fibre is a field and
`Delta` is nonsquare, so the Galois group is `S3`.

## 2. The moving triple root creates a smooth one-place decic

The degree-three leading coefficient row is the moving triple root `ell^3`.
The quadratic perturbation evaluates at its repeated-root address
`[U:V]=[-C:A]` as

```text
(C ell U^2+A^2V^3)|_(-C,A)=A^5.                           (6)
```

Accordingly the degree-ten discriminant form is `-27A^10`.  This is the
first possible moving-triple-root degree after the degree-two exclusion in
THM-3908.

Let `Delta_h` be the degree-ten homogenization.  At infinity,

```text
Delta_h(A,C,0)=-27A^10.                                   (7)
```

There is only the point `P_infinity=[0:1:0]`, and the surviving term
`-4C^9Z` gives

```text
partial_Z Delta_h(P_infinity)=-4.                         (8)
```

Thus the projective curve is smooth at its sole infinity point.  This also
proves absolute irreducibility without relying on a rational factorizer.  If
`Delta_h` had two positive-degree geometric factors, each projective
component would meet the line at infinity; uniqueness in `(7)` would force
both through `P_infinity`, making the product singular there, contrary to
`(8)`.

## 3. The hidden normalization is elliptic

Choose constants `k_0,l_0 in k*` satisfying

```text
6k_0^2=1,                         3l_0^2+k_0=0.            (9)
```

On the smooth elliptic curve

```text
E: y^2=w(w^2-1)                                             (10)
```

define

```text
A=l_0 w(w^2+2)y,
C=k_0 w(w^2-1)(w^2+2).                                    (11)
```

Direct substitution gives `Delta(A,C)=0`.  The map is birational, not an
uncontrolled cover.  In the scaled coordinates

```text
A_0=A/l_0,                  C_0=C/k_0,
R_0=A_0^2/C_0,                                              (12)
```

put

```text
zeta=(C_0^2/R_0+2R_0-2)/(R_0+1).                          (13)
```

On `(11)`, one has

```text
R_0=w^4+2w^2,                 zeta=w^2,
w=C_0/((zeta-1)(zeta+2)),
y=A_0/(w(w^2+2)).                                            (14)
```

These are rational inverse formulas.  At the elliptic origin `O`, the pole
orders are

```text
ord_O(C)=-10,                       ord_O(A)=-9.            (15)
```

Consequently `(11)` extends to the projective normalization, maps only `O`
to infinity, has image degree ten, and has local parameter `A/C` there.  The
affine normalization is exactly

```text
E minus {O}.                                                 (16)
```

## 4. The quadratic resolvent really has a three-class

On a smooth resolution of the projective double plane `W^2=Delta_h`, the
preimage of the line at infinity splits as rational curves `B_+,B_-`.
Their local equations differ first in order five, and the standard
double-plane adjunction gives

```text
B_+^2=B_-^2=-4,                       B_+B_-=5.             (17)
```

Thus the split-boundary Gram matrix and Smith data are

```text
[-4  5],                    SNF=(1,9),                     (18)
[ 5 -4]
```

in exact agreement with the first three-primary boundary degree allowed by
THM-3912.  The matrix is nonsingular.  The affine double plane `Q` is normal:
its reduced irreducible plane branch has only finitely many singular points.
As in THM-3911, the divisor of a unit on `Q` pulls back to a combination of
`B_+,B_-`; intersecting with the two
boundary curves forces both coefficients to vanish.  Therefore `Q*=k*`.

The normal `S3` cubic has a connected cyclic cubic Galois-closure layer over
the quadratic resolvent.  Generic transposition inertia is absorbed by the
quadratic base change, so this layer is etale on `Q_reg`.  Normality gives

```text
Gamma(Q_reg,O)^*=k*,                    Pic(Q_reg)=Cl(Q).   (19)
```

The Kummer sequence over the algebraically closed field identifies the
nontrivial cyclic cover with a nonzero class in `Cl(Q)[3]`.  Thus the
three-primary permission in `(18)` is genuinely used here; it is not merely
a lattice possibility.

## 5. Why this is still not a plane atlas

Suppose an affine plane were a dominant etale degree-three atlas for this
cubic function field.  THM-3801 identifies it with an open of the finite
normalization.  The ramified divisor over the sole branch `(5)` is absent.
The deleted-divisor valuation argument of THM-3841 therefore puts
`V(Delta)` inside the nonproperness set of the resulting generically finite
polynomial plane map.  That set is pure of dimension one when nonempty;
since `V(Delta)` is irreducible, it is an entire irreducible component, not
merely a subset of one.

The polynomial-uniruledness gate used in THM-3841 would then give a dominant
polynomial curve `A1 -> V(Delta)`.  Normality of `A1` lifts it through `(16)`
to a nonconstant morphism

```text
A1 -> E minus {O}.                                          (20)
```

Properness extends `(20)` to `P1 -> E`, but a nonconstant map from `P1` to
an elliptic curve is impossible by Riemann--Hurwitz.  Hence no plane atlas
exists.

The failure coordinate is now unusually sharp.  Quadratic coefficient depth
could not confluence the infinity places at all.  This cubic-depth example
confluences them, remains nonmonogenic and normal, and globalizes the
resolvent three-class; it pays for those gains with genus one.  The next
counterexample search should deform `(1)` so that the elliptic normalization
degenerates rationally while preserving the smooth one-place infinity and
the nonmonogenic `S3` order.  `JC(2)` remains **OPEN**.

## 6. Exact replay

Run

```bash
python3 04-computation/jc2_moving_triple_root_one_place_elliptic_decic_thm3913.py
python3 -O 04-computation/jc2_moving_triple_root_one_place_elliptic_decic_thm3913.py
```

Both streams must byte-match the frozen output.  The companion checks the
compact binary cubic, primitive/common-zero coefficient packet, exact
discriminant, moving-root quintic, unique smooth infinity, Kummer
irreducibility specialization, elliptic map and inverse, pole orders, and
split-boundary Smith data in 24 active gates.  The Delone--Faddeev normality,
absolute-irreducibility, Kummer-class, and Jelonek bridges have been
independently hostile-audited.  **QED.**
