---
id: THM-4431
title: "Colored lattice basis and three-direction LRC network closure"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4386/4414 + FINITE-EXACT +
  INDEPENDENTLY AUDITED. Every eligible complete live carrier body with
  exactly three primitive unoriented directions has all three degree-zero
  raw network projections strictly below 6/77. Two live primitive vectors
  form a basis of the owner lattice, and the third direction has only the
  circuit normal forms z=u+v or z=-2u+v. This is not chart entry,
  synchronization, or a proof of LRC(14).
source: projection_inequality + root / next-sharp continuation, 2026-09-06
depends_on:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4428-lrc14-two-direction-network-closure-and-sharp-one-direction-gap
script: 04-computation/lrc14_colored_basis_three_ray_h99_independent_referee.py
output: 05-knowledge/results/lrc14_colored_basis_three_ray_h99_independent_referee.out
script_sha256: 227576fef415844edb01dd2ec66873caa4de77102059322548f6f69ad46fedbb
output_sha256: 7dc3cfe799cebbf3b4337a6114c008d60fb0fa566380dc3e16b4d87f6e5160eb
semantic_sha256: 6ecbb5a2e079338af842a2f5d6c474ccdb8eb359941d9a8b192490d677115729
hash_basis: raw LF bytes
audit: >
  PASS. An import-free middle-coordinate enumerator checks the complete
  5,409-row eligible universe below height 99, including all 1,791 exactly
  three-direction rows and both circuit normal forms. It reconstructs every
  carrier and projection, checks the determinant and reciprocal inequalities,
  and exercises wide controls. Its 65,789 explicit gates survive optimization;
  normal and optimized outputs byte-match the frozen LF artifact. The
  all-height colored-basis descent and circuit classification are proved
  symbolically below and were separately audited line by line. The independent
  overnight hexagon proof supplies an additional no-three-in-line cap theorem.
---

# THM-4431 -- Colored lattice basis and three-direction LRC network closure

**PROVED ELEMENTARY RELATIVE TO THM-4386/4414 + FINITE-EXACT +
INDEPENDENTLY AUDITED.** This theorem closes the exact three-direction local
network. It supplies no physical chart entry or synchronization;
`LRC(14)` remains **OPEN**.

## 1. Statement and native colored lattice

Let `w=(a,b,c)` be primitive, sorted, distinct, positive, odd, and nonzero
modulo three. Let `Lambda(w)` be the complete THM-4414 live carrier set and
`E_1,E_2,E_3` its raw network projections. If `Lambda(w)` has exactly three
primitive unoriented directions, then

```text
E_i(w)<6/77 for i=1,2,3.                               (1)
```

Put

```text
L={C in Z^3:C dot w=0},
Gamma={C in L:aC_1=bC_2=cC_3 (mod 3)},
H=3L.                                                   (2)
```

The common residue defines `tau:Gamma->F_3` with kernel `H`. Modulo `3L`,
`L` is the two-dimensional kernel of `w`, while `Gamma/3L` is its
one-dimensional equal-product line. Hence

```text
[L:Gamma]=3,       Gamma/H is cyclic of order 3.        (3)
```

The live set is exactly `(Gamma\H) intersect K`, where `K` is the open,
bounded, centrally symmetric convex roof body in the real kernel plane.
Projection to the first two coordinates identifies `L` with

```text
{(x,y) in Z^2:ax+by=0 (mod c)},
```

of determinant `c`; thus `Gamma` has determinant `3c` in these coordinates.
If a live vector is `C=g d` with `d` primitive, then `3` does not divide `g`.
Division preserves its nonzero owner color and central convexity preserves
membership in `K`. Therefore every live ray has a primitive live
representative.

## 2. Colored-basis descent

**Lemma.** If the complete live body has two nonparallel directions, two live
primitive representatives form a basis of `Gamma`.

Choose a nonparallel live primitive pair `u,v` of least positive index `D`
in `Gamma`. If `D>1`, choose a nonzero representative

```text
x=t u+s v,       |t|,|s|<=1/2                          (4)
```

of a nontrivial coset of `Zu+Zv` in a centered fundamental parallelogram.
Both coefficients are nonzero: an integral point on the line of a primitive
integer vector is an integral multiple of it. Since `|t|+|s|<=1`, central
convexity puts `x` in `K`.

If `x` is live, primitive reduction gives a live pair of determinant
strictly between zero and `D`. If `x` is dead, then `x in H`. When
`|s|<=|t|`, set

```text
y=u-sign(t)x=(1-|t|)u-sign(t)s v.                      (5)
```

The absolute coefficient sum is `1-|t|+|s|<=1`, so `y in K`; also
`y=u (mod H)`, so `y` is live, and

```text
0<|det(u,y)|=|s|D<D.                                   (6)
```

For `|t|<|s|`, the symmetric repair `y=v-sign(s)x` gives determinant
`|t|D`. Primitive reduction again preserves color and decreases determinant.
Both alternatives contradict minimality. Thus `D=1` in `Gamma`, or

```text
|det(u,v)|=3c                                          (7)
```

in the original first-two coordinates. The subgroup sidecar is essential:
subtracting a dead representative preserves the live coset.

## 3. Complete three-direction classification

Assume there are exactly three unoriented primitive live directions. Choose
the basis above and orient `u,v` so `tau(u)=tau(v)=1`. Write the third as

```text
z=m u+n v,       gcd(m,n)=1,       m+n !=0 (mod 3).     (8)
```

At least one of `|m|,|n|` equals one. If instead both are at least two and
have the same sign, orient `z` so `m,n>=2`; then

```text
u+v=[z+(n-1)u+(m-1)v]/(m+n-1)                          (9)
```

is in the convex body, has nonzero owner color, and after primitive reduction
gives a fourth ray. If the signs differ, write `z=m u-r v`, `m,r>=2`.
Coprimality gives `m!=r`. For `m>=r+1`,

```text
2u-v=[2z+2(r-1)u+(m-r-1)(-v)]/(m+r-1);                (10)
```

for `r>=m+1`,

```text
u-2v=[2z+(r-m-1)u+2(m-1)(-v)]/(m+r-1).                (11)
```

Each is a convex combination, is live, and yields a fourth direction.

Swap the basis and reverse `z` if needed to write `z=m u+v`. The segment
from `v` to `z` contains every integral `(k,1)` between `0` and `m` in basis
coordinates. Its owner color is `k+1`, so exactly `k=2 (mod 3)` is dead.
Distinctness excludes `m=0`; `m=-1` and `m=2` are dead; `m>=3` forces the
extra live direction `(1,1)`; and `m<=-3` forces `(-2,1)`. Consequently the
only normal forms are

```text
z=u+v                 or                 z=-2u+v.       (12)
```

Equivalently, the absolute primitive circuit coefficients are `(1,1,1)` or
`(1,1,2)`. This is a classification of complete colored convex bodies, not
of arbitrary selected triples of rays.

## 4. Reciprocal tail

Every direction in a multi-direction body has `l1` norm at least sixteen by
THM-4386; otherwise the zero-defect rigidity would make the whole body
collinear. Hence its maximum coefficient `M_i` is at least seven. For any two
distinct directions `u,v`,

```text
u cross v=k w,       3 divides k,       k!=0.           (13)
```

The factor three follows because all live owner vectors are parallel modulo
three. Taking the third coordinate gives, with `P=3c/2`,

```text
M_i M_j>=P for every pair.                              (14)
```

Let `m=min M_i` and `R=sum_i 1/M_i`. If `P<=m^2`, then
`R<=3/sqrt(P)`. If `P>m^2`, the other two maxima are at least `P/m`, so
`R<=1/m+2m/P`. When `P>=14m`, this is at most `1/7+14/P`. When
`m^2<P<14m`, put `rho=P/m^2`; then `1<rho<=2` and
`(rho+2)^2<=9rho`, giving `R^2<=9/P`. Thus

```text
R<=1/7+14/P       or       R^2<=9/P=6/c.               (15)
```

The exact residue-deleted multiplier lists on the three rays give

```text
E_i<(12/49)R+12/(7c).                                  (16)
```

In the linear branch this is at most `12/343+4/c`; in the radical branch it
is at most `(12/49)sqrt(6/c)+12/(7c)`. Both decrease with `c`. At `c=99`,

```text
12/343+4/99=2560/33957<6/77,
96/26411<4/1089                                        (17)
```

where the second exact inequality is the positive squared form of the
radical comparison. Therefore `(1)` holds for every `c>=99`.

## 5. Complete finite head and audit

The import-free referee exhausts

```text
1<=a<b<c<99, all odd and nonzero mod 3, gcd(a,b,c)=1.  (18)
```

It enumerates `5,409` eligible rows by solving for the middle carrier
coordinate, independently reconstructs all projections, and finds `3,500`
multi-direction rows and `1,791` exactly-three-direction rows. The normal
forms split as

```text
z=-2u+v: 684 rows;       z=u+v: 1,107 rows.             (19)
```

All three projections are strict. The head leader is

```text
w=(5,37,43),
E=(240/11137,18/301,2822/55685).                       (20)
```

The wide row `(5,191,199)` has sixteen carriers on exactly three directions,
so the theorem is not a six-point statement. The `(7,611,613)` control has
thirteen directions and lies outside the hypothesis. Normal and optimized
runs pass `65,789` explicit checks and byte-match the frozen raw-LF output.

## 6. No-three-in-line cap

For a complete colored carrier set with no three distinct collinear points,
parity zero is absent: otherwise `-C,C/2,C` are three live collinear points.
Each of the two live owner colors therefore uses at most one point in each of
the three other parity classes, so the total size is at most six, sharply at
`(23,29,37)`. The same proof gives `2(2^d-1)` in rank `d` with an index-three
dead subgroup; no higher-rank sharpness is claimed. This is a cap theorem,
not a bound of six on arbitrary dictionaries or on several points of one ray.
The independent proof and raw-carrier controls are recorded in
`05-knowledge/results/lrc14_colored_basis_three_ray_overnight_hexagon_sep05.md`.

## 7. Reproduction and remaining scope

```powershell
python -B 04-computation/lrc14_colored_basis_three_ray_h99_independent_referee.py
python -B -O 04-computation/lrc14_colored_basis_three_ray_h99_independent_referee.py
```

Four or more directions are not covered by this theorem's structural
classification. THM-4434 subsequently closes their local projections, but
does not supply chart entry, synchronization, or `LRC(14)`.
