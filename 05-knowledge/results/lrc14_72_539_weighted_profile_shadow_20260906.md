# The two `72/539` occurrences: a scalar shadow, not a conjugacy

**Verdict.** There is a **PROVED exact scalar bridge**, stronger than a bare
numerical coincidence:

```text
THM-4437 auxiliary coefficient-envelope slope at (4,7,11)
 =2/49+4/77+2/49
 =THM-4449 pair energy of (1,7,11)
 =THM-4449 physical owner-cut mass of (1,7,11)
 =72/539.                                                (1)
```

The last equality uses a special zero-overlap owner mechanism. However,
there is **no affine label equivalence, no THM-4449 common-dilation ray
equivalence, and no coefficient-as-speed owner transport** between the two
triples. The faithful typing is therefore:

> **PROVED unlabelled weighted-profile seventh-wall shadow; REFUTED affine,
> common-dilation, and coefficient-as-speed owner identifications.**

This yields no new THM-4449 equality shape and no LRC(14) closure.

## 1. What the earlier occurrence actually was

In the THM-4437 coefficient box, `(4,7,11)` is a **relation-coefficient
magnitude pattern**, not a speed triple. Its maximizing signed sector and
normalized speed-polygon vertex are

```text
v=(-11,4,7),       w=(1,1,1).                            (2)
```

Thus `v dot w=0` is just the coefficient identity `-11+4+7=0`. The point
`w=(1,1,1)` is an equal-speed boundary vertex of the closed normalized
polygon, not a physical distinct-speed row.

The exact-token repository search found the pre-THM-4449 value in the
THM-4437 coefficient outputs and their independent audit:

- `05-knowledge/results/lrc14_coefficient_box_empty_core_three_ray_sep06.out`;
- `05-knowledge/results/lrc14_all_parity_coefficient_box_thm4437.out`;
- `05-knowledge/results/lrc14_all_parity_coefficient_box_thm4437_independent.out`.

The self-contained replay confirms that `(4,7,11)` is the unique pattern
with slope `72/539` among all `750` patterns in THM-4437's complete
max-coefficient-18 parity-free box. Apparent earlier July hits returned by a
substring search were pieces of larger denominators such as `...72/5399`
and `...72/5392695`; they are not occurrences of the rational number.

## 2. Analytic reconstruction of the coefficient value

Retain THM-4434/4437 notation with error-cube radius `R=3/14`. Because all
coordinates of `v=(-11,4,7)` are units modulo three, the allowed defects are

```text
delta=-3,0,3,
```

and the owner-residue weight is `rho=2/3`. At `w=(1,1,1)`, use the scalar
coordinate

```text
s=(e_0-e_2)/4
```

on the plane slices `v dot e=delta` inside `[-3/14,3/14]^3`. Direct vertex
evaluation gives

```text
width s(v dot e=-3)=3/49,
width s(v dot e= 0)=6/77,
width s(v dot e= 3)=3/49.                               (3)
```

For example, at `delta=-3` the scalar extrema are `9/196` and `3/28`,
whose difference is `3/49`; at `delta=0` they are `-3/77` and `3/77`.
The positive defect is the central reflection of the negative one. Hence

```text
F^(3)_v(w)
 =(2/3)(3/49+6/77+3/49)
 =2/49+4/77+2/49
 =72/539.                                                (4)
```

Convexity of each section width reduces the coefficient maximum to polygon
vertices as in THM-4437. The new audit reconstructs all vertices with exact
rational clipping, obtains (2)--(4), and independently rescans the entire
750-pattern box.

## 3. Analytic reconstruction of the dyadic value

In THM-4449, `(1,7,11)` is instead a genuine primitive odd-3-unit **speed
triple**. Its three quotient pair cross-combs have the exact components

```text
Sigma_(1,7):  (6/49,1/7),       (6/7,43/49);
Sigma_(1,11): (6/77,8/77),      (69/77,71/77);
Sigma_(7,11): (13/49,2/7),      (5/7,36/49).             (5)
```

Their respective masses are

```text
2/49, 4/77, 2/49.                                       (6)
```

The six open intervals in (5) are pairwise disjoint. Thus the mixed
three-owner correction `Omega_T` in THM-4449 is empty up to endpoints, and

```text
mu(F_(1,7,11))
 =sum_{{a,b}} mu(Sigma_(a,b))
 =2/49+4/77+2/49
 =72/539.                                                (7)
```

This proves the exact bridge (1). It also identifies its source: both
calculations resolve the same `1/14` wall through seventh-denominator interval
widths. The central coefficient defect has weighted width `4/77`, matching
the `1:11` pair, while the two reflected outer defects each have weighted
width `2/49`, matching the other two pairs.

This is not an owner bijection. The signs `delta=+3,-3` do not canonically
distinguish the equal-mass `1:7` and `7:11` owners, and coefficient-slice
coordinates retain neither their six circle locations nor their tail labels.

## 4. Three exact obstructions to a stronger identification

### No relevant projective equivalence

THM-4449's projective equivalence is permutation plus common positive
dilation. Both triples are primitive, and

```text
{1,7,11} != {4,7,11},
```

so they lie on different projective rays. An abstract `PGL_2(Q)` map can
send any three rational points to any other three, but it does not preserve
the multiplication combs and is not an LRC symmetry.

### No affine label map

An affine bijection of label sets would scale every pairwise difference by
one constant. The sorted positive differences are

```text
{4,6,10} for {1,7,11},       {3,4,7} for {4,7,11}.
```

Their largest-to-smallest ratios are `5/2` and `7/3`, so no affine map can
exist. Exhausting all six target permutations gives zero affine maps.
Likewise a circle covering `x -> kx` changes every speed by the same factor,
and a translation changes phases but not frequency magnitudes; neither can
bridge the two supports.

### Owner-conjugacy hostile

If the coefficient pattern is incorrectly reinterpreted as the speed triple
`(4,7,11)`, then

```text
D_4+1/2=D_4.
```

The even speed can own both half-lifts by itself, whereas every odd tail in
THM-4449 can own at most one. Exact wall integration gives

```text
mu(F_(4,7,11))=9/49 !=72/539.                            (8)
```

This is a direct owner-invariant refutation, not merely a mismatch of labels.

The tempting transformation

```text
(1,p,q) -> (q-p,p,q)
```

also fails immediately nearby:

```text
coefficient slope at (6,7,13) =233/1911,
physical mass at (1,7,13)     =80/637=240/1911.          (9)
```

Thus (1) is a special scalar alignment, not a transport theorem for a family.

## 5. Exact scope and replay

The connection map is

```text
source: THM-4437 defect-slice widths at v=(-11,4,7), w=(1,1,1)
target: THM-4449 pair energies for T=(1,7,11)
map: (-3,0,+3) weighted widths -> (2/49,4/77,2/49)
preserved: total scalar mass and the central-versus-outer split
destroyed: circle address, tail labels, sheet owners, scale, distinctness
sidecar making target physical: pair-comb disjointness, Omega_T=0
hostile: coefficient tuple as speeds has mass 9/49
next use: search only for termwise slice/energy identities, never infer a
          phase or coefficient-as-speed owner map from equal totals.   (10)
```

Run:

```text
python -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906.py
python -O -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906.py
python -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
python -O -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
```

The script is self-contained, uses active `need` gates rather than asserts,
and ends with

```text
status=PASS;checks=261633;typed=SCALAR_BRIDGE_NOT_AFFINE_PROJECTIVE_OWNER_CONJUGACY
```

```text
04-computation/lrc14_72_539_weighted_profile_shadow_20260906.py
  SHA256: 143CCA9032AE9690FF33069D93A1011AD246027F1C9AF2677C1078F00A0BE91F
05-knowledge/results/lrc14_72_539_weighted_profile_shadow_20260906.out
  SHA256: 439843AF0B3C9C7496C2419C50ED7576673524136024074C70760EFD508F8ACB
04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
  SHA256: 0063EBF7FB63743323D849FC0FB24A8CBBE7991E451F30CEB3126A6E7E5D4485
05-knowledge/results/lrc14_72_539_weighted_profile_shadow_20260906_independent.out
  SHA256: 64A3D57DE8C2DFE54915E3599E9198793048159C073107D8488D53B8E124AD6D
```

The clean-room report is
`05-knowledge/results/lrc14_72_539_weighted_profile_shadow_20260906_independent_audit.md`.

Canonical sources read:

- `01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md`;
- `01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md`;
- `01-canon/theorems/THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry.md`.
