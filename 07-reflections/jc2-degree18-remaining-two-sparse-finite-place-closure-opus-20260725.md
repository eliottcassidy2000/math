---
source: codex-2026-07-25-degree18-remaining-two-sparse
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Within the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, one good finite place for each irreducible
  B--W, C--D, and D--W ratio packet forces ten simple finite branch
  values, a geometrically irreducible trigonal curve, and normalization
  genus at least three.  The genus/deck argument empties all sixteen
  ratios on those planes.  Together with the B--C and B--D packets and
  THM-2297's C--W closure, this empties the complete 31-ratio
  exactly-two-sparse locus.  Higher-support strata,
  split/even descent, JC(2), and DC(2) remain open.
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-bc-algebraic-finite-place-genus-floor-opus-20260725
  - jc2-degree18-bc-rational-nodal-genus3-closure-opus-20260725
  - jc2-degree18-bd-4075-smooth-cubic-closure-opus-20260725
  - jc2-degree18-bd-quartic-nodal-atlas-opus-20260725
script: 04-computation/jc2_degree18_remaining_two_sparse_modular_floor.py
output: 05-knowledge/results/jc2_degree18_remaining_two_sparse_modular_floor.out
script_sha256: 3eb2da72c29a4e6ab4d0127952e2910fca7bcf1ac8de911e710c3458f900c236
output_sha256: 13025ec869e816b81f25bed6dc993a2daf7f3f1249ba238bd28626f8c2636b6e
supporting_script: 04-computation/jc2_degree18_bc_algebraic_modular_floor.py
supporting_script_sha256: 262d99129916bb6904507322b824830e3c6ab140c6c996f76a0d0d31f12b2d4b
hash_basis: LF-normalized working-tree bytes
---

# The remaining two-sparse planes have a uniform genus floor

## 1. Three rationalizing coordinates

The finite-place mechanism from the `B`--`C` packet extends because the
other weighted ratios also admit a constant algebraic coordinate in
which the spectral cubic is defined over `Q(t)`.

### The B--W plane

Normalize `B=1`, write `W^2=t`, put

```text
x=Wy,
```

and multiply by `t^3`.  The curve is

```text
H_BW
 =-26040609t^3u^3
  +(49601160t^3+1607445t^2x^2)u^2
  +(-20995200t^3-2857680t^2x^2-138915tx^4)u
  -5878656t^3x+777600t^2x^2+78120tx^4+1127x^6.
                                                            (1)
```

### The C--D plane

Normalize `C=1`, write `D^3=t`, then put

```text
y=Dx,                      u=D^2v.
```

After division by the nonzero constant `t`, the curve is

```text
H_CD
 =-26040609tv^3+1607445tx^2v^2
  +(-52907904-138915tx^4)v
  +1959552x^2-435456x^3+1127tx^6.                  (2)
```

### The D--W plane

Normalize `D=1`, write `W^4=t`, and put

```text
y=Wx,                      u=W^2v.
```

After division by `W^2`, the rational curve is

```text
H_DW
 =-26040609tv^3+1607445tx^2v^2
  +(-52907904-138915tx^4)v
  -5878656x+1959552x^2+1127tx^6.                  (3)
```

Every transformation is an invertible constant change over the
algebraic closure and preserves constancy of a Keller trajectory.

## 2. Exact good-place atlas

For each factor of THM-2311's three ratio banks, the script searches
good primes and proves the parameter reduction irreducible by exhaustive
monic-factor division through half its degree.  The accepted atlas is:

```text
plane packet     deg  prime  monic parameter polynomial       x_*

BW    quadratic   2    19    5+4t+t^2                          6
BW    sextic      6    13    4+12t+11t^2+12t^3+12t^4+t^5+t^6 1

CD    linear      1    11    5+t                              3
CD    cubic       3     5    3+4t^2+t^3                       1

DW    linear      1    11    1+t                              8
DW    cubic       3    13    12+5t+9t^2+t^3                   3.       (4)
```

Here `x_*` is a residue-field specialization at which the cubic in the
cover coordinate has no root.  The test is the exact gcd with
`v^Q-v`, where `Q` is the residue-field size.

At every row of (4), independently,

```text
the leading cubic coefficient is a unit;

the parameter polynomial is irreducible;

deg_x Delta_t=12;

deg_x gcd(Delta_t,Delta_t')=1;

the leading cubic at infinity is separable;

H_t(cover_coordinate,x_*) is irreducible.           (5)
```

The residue-field orders, respectively, are

```text
19^2, 13^6, 11, 5^3, 11, 13^3.                    (6)
```

## 3. Characteristic-zero consequences

Each irreducible parameter factor in THM-2311 divides the exact
repeated-branch resultant.  The constant algebraic changes in
(1)--(3) only rescale the branch polynomial and its base coordinate by
nonzero constants, so repeated-root existence is unchanged.  Therefore
its characteristic-zero
discriminant has a nonconstant gcd with its derivative.  At a good
place the gcd degree can only increase.  Equations (5) give

```text
deg gcd(Delta_t,Delta_t')=1                         (7)
```

over every parameter field in (4).

Since `deg Delta_t=12`, every conjugate curve has one double
discriminant value and ten simple values.  The ten simple values are
ten smooth index-two branch points.  The double value may be nodal,
colliding, or more ramified; no classification is needed because it
cannot subtract ramification from the ten simple points.

The irreducible specialization in (5) proves the trigonal polynomial
irreducible over the parameter function field: after division by the
unit leading coefficient, a factorization over the local fraction field
would have monic integral factors and would reduce without losing
degree.

Geometric irreducibility follows from a rational smooth point.  In all
three normal forms,

```text
H(0,0)=0,                     H_cover(0,0)!=0.       (8)
```

Indeed the linear cover-coordinate coefficient is a nonzero constant
multiple of `t` or of `1`.  A rational smooth point cannot lie on
several conjugate geometric components of an irreducible curve.

## 4. Genus floor and closure

The normalization is a connected degree-three cover of `P^1`.  Its ten
simple finite branch points give ramification degree at least ten, so

```text
2g-2>=-6+10=4,

g>=3.                                               (9)
```

The infinity separability in (5) confirms that no cover degree is lost
there; it is not needed for the lower bound.

A rational Keller trajectory into the normalization is constant.
Equations (1)--(3) then make the original `u,y` constant after restoring
the fixed algebraic scaling constants.  For `y!=0`, THM-2262's first
flux makes `T^2,T,q` constant and contradicts the genuine nonsplit
deck.  For `y=0`, THM-2262 Section 4 closes the exceptional center.

Thus the packet closes

```text
BW: 2+6=8 ratios;

CD: 1+3=4 ratios;

DW: 1+3=4 ratios;                  total 16.         (10)
```

Together with the nine `B`--`C` ratios and six `B`--`D` ratios, all
`31` off-axis ratios in THM-2311 are empty.  THM-2297 already removes
the raw `C`--`W` plane.  The finite-place packets are independently
hostile-audited.  The genuine nonsplit polynomial exact-square-prefix
branch therefore has no exactly-two-sparse degree-eighteen survivor.

## 5. Information ledger

```text
source:
  the 16 B--W/C--D/D--W candidate ratios;

map:
  weighted normalization, rationalizing source coordinate, then one
  inert finite place per irreducible parameter packet;

preserved:
  full conjugacy packets, cover degree, the branch gcd inequality,
  and a rational smooth point;

destroyed:
  exact locations and local types of the double branch values;

sidecar:
  THM-2311 supplies gcd-degree>=1 in characteristic zero, while the
  good reduction supplies gcd-degree<=1 and an irreducible fibre;

hostile boundary:
  an irreducible parameter reduction alone does not prove the spectral
  cubic irreducible; the separate no-root specialization is essential;

next target:
  move from exactly two nonzero translation invariants to the first
  three-sparse face, using finite-place simple-branch floors before
  attempting full normalization.                                  (11)
```

Run

```text
python3 04-computation/jc2_degree18_remaining_two_sparse_modular_floor.py
python3 -O 04-computation/jc2_degree18_remaining_two_sparse_modular_floor.py
```

Both executions must match

```text
05-knowledge/results/jc2_degree18_remaining_two_sparse_modular_floor.out
```

after LF normalization.
