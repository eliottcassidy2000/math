---
id: THM-2804
title: "Fourteen-rail V4 configuration completion and eleven-carry anchor defect"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT HOSTILE AUDIT.
  On each of the first fourteen THM-2749 rails, the three maximal
  source-twelve configurations of THM-2797 and one canonical adjacent
  configuration form the even-parity V4 in (sector,edge,kappa), with
  h=6(1-sector).  Their carry capacities are exactly 12,12,12,11.  The
  fourth vertex is the carry-support meet, restores private root 12, and
  omits target labels one and two.  Its source-twelve facet is positive on
  every rail and shares a marked positive clock, yet all 23 nonempty marked
  comparisons die at the present anchor, with zero open overlap and zero
  closure seams.  Abstractly AGL(2,2)=S4 and its translation-V4 quotient is
  S3, but capacity singles out the fourth vertex, so no S3-equivariant
  physical or semantic transport follows.
source: lrc-a12-carry-bridge/v4-configuration-completion-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2797-fourteen-rail-source-twelve-configuration-switch-semantic-base-no-go
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2795-tree-star-circuit-dwell-time-v4-diamonds-and-ternary-move-boundary
script: 04-computation/lrc14_v4_configuration_completion_anchor_defect_thm2804.py
output: 05-knowledge/results/lrc14_v4_configuration_completion_anchor_defect_thm2804.out
script_sha256: 0143f83e0762946caee23043275cf28a7aa863f92a0a033fa325c3908062a3f4
output_sha256: c22b12b909f939363494eff78590ec4bb8f00945c47dbaf3a592e5438dd12d85
hash_basis: LF-normalized bytes
---

# THM-2804 -- the V4 completion restores root 12 but loses a carry

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT HOSTILE
AUDIT.**

THM-2797 finds the same three maximal source-twelve configurations on every
one of the first fourteen semantic rails.  Their bit pattern is not
accidental.  They are three vertices of the even-parity `V4` inside
`F_2^3`.  This theorem inserts the fourth vertex and determines exactly what
is gained and lost.

The result gives a literal `2`/`3` symmetry boundary.  The four configurations
have the abstract affine symmetry

```text
AGL(2,2)=S4,                 S4/V4=GL(2,2)=S3.                     (1)
```

But the physical carry-capacity function is not equivariant: it assigns
`12,12,12,11`.  Thus the abstract ternary quotient cannot be used as a
semantic transport law.

## 1. The four configuration vertices

Write a configuration bit as

```text
x=(sector,edge,kappa) in F_2^3.                                   (2)
```

On every rail `j=0,...,13`, use

```text
vertex   x       h
  A     000      6
  B     011      6
  C     110      0
  D     101      0.                                                (3)
```

These are exactly the even-parity triples

```text
sector+edge+kappa=0 in F_2,                                       (4)
```

and their height is the affine sector law

```text
h=6(1-sector).                                                     (5)
```

Coordinatewise xor therefore makes the four points a literal `V4`.

The companion independently enumerates all invertible affine maps of
`F_2^2`, using `(sector,edge)` as coordinates and
`kappa=sector+edge`.  It obtains all `24` permutations of the four points.
The four translations form a normal subgroup, and the induced action of the
six linear maps on the three nonzero directions gives every permutation of
those directions.  Hence `(1)` is exact.

This is an abstract symmetry of the four configuration labels.  No claim is
made that all affine permutations are lawful physical carrier operations.

## 2. Exact carry-support diamond

Let `U_X` be the set of carries on which vertex `X` has a normalized
THM-2640 unit.  On every one of the fourteen rails:

```text
U_A=U_B=F_13\{6},
U_C=F_13\{0},
U_D=F_13\{0,6}.                                                     (6)
```

Thus

```text
|U_A|,|U_B|,|U_C|,|U_D| = 12,12,12,11,                             (7)
U_D=U_A intersect U_C=U_B intersect U_C,
U_A union U_C=F_13.                                                 (8)
```

The first three vertices are precisely the three maximal signatures in
THM-2797.  The fourth is their exact carry-support meet.

At source carry `c0=12`, the private roots are

```text
(r_A,r_B,r_C,r_D)=(12,12,11,12).                                   (9)
```

Thus `D` repairs the root mismatch of `C`: it has the desired private root
`12`.  Under the source-twelve label rebase

```text
delta(c)=2(c-12) mod 13,                                           (10)
```

the omitted target labels are

```text
A:{1},        B:{1},        C:{2},        D:{1,2}.                 (11)
```

The price of the root repair is exactly one additional missing carry/label.

Equations `(7)` and `(11)` distinguish `D` from the other three vertices.
Consequently neither the translation `V4` nor the quotient `S3` preserves the
physical capacity data.  The abstract relation `(1)` cannot yield an
`S3`-equivariant marked-current or semantic-attachment theorem.

## 3. The fourth vertex is honestly positive

For `D`, fix

```text
(sector,edge,kappa,h)=(1,0,1,0),
source carry=12,
active labels=F_13\{1,2}.                                         (12)
```

On every rail, intersect its eleven translated rail supports, eleven
clock-matched present supports, and root-`12` private half, then apply the
exact sector-one delayed prefix.  All fourteen facets are positive.  Their
positive clock supports are

```text
(1), (6), (1,2,3), (0,5), (2,3), (0,4,6), (0,2,3),
(0,4,5,6), (0,1,2,3), (0,4,5), (0,1,3), (4,5), (0,2), (4,5,6).    (13)
```

Each intersects the positive-clock support of its matching THM-2749 marked
rail.  The total exact facet numerator is

```text
103771702146779222198643000.                                      (14)
```

Thus `D` is not a formal or zero-mass boundary object.

## 4. Root repair does not repair the semantic anchor

Let `R_D,j` be the elevenfold common rail of vertex `D` and
`M_j,ell` the matching THM-2749 marked source base.  Insert the label-zero
anchor

```text
F_(ell,0),                                                         (15)
```

because `h=0`.  Exactly as in the height-zero branch of THM-2797, every
nonempty marked comparison satisfies

```text
R_D,j intersect M_j,ell is nonempty,
(R_D,j intersect M_j,ell) intersect F_(ell,0) is empty.            (16)
```

The first-failure census over `14*7=98` clock comparisons is

```text
anchor-present:23,                         marked-empty:75.         (17)
```

After the anchor, later translated presents, the restored root-`12` private
half, and the delayed prefix are only further intersections.  Therefore the
full `D` base has

```text
strict open overlap with the marked base =0.                       (18)
```

Unlike the middle maximal vertex of THM-2797, `D` has no boundary rescue:

```text
facet-right/marked-left seams=0,
marked-right/facet-left seams=0                                    (19)
```

on every rail.  Restoring the correct root has moved the packet away from
the marked base rather than opening a seam.

## 5. Consequence and modular-group boundary

The four configurations realize exactly the local algebra suggested by the
recurring dyadic/ternary grammar:

```text
four affine points = a V4 torsor,
three directions   = the S3 quotient action.                       (20)
```

What fails is equally exact.  The carrier invariant

```text
X -> |U_X|                                                          (21)
```

is not constant on the torsor.  It singles out `D`, and the support labels in
`(6)` distinguish the other vertices further.  Therefore:

- the abstract `S4/V4=S3` quotient is real;
- it organizes the configuration moves;
- it does not preserve unit-carry capacity, root/label data, or semantic
  attachment;
- it does not produce a physical `PSL_2(Z)=C_2*C_3` action.

The next successful move must carry a sidecar that restores the information
destroyed by the quotient—at minimum the carry-support label and anchor
polarity.  A bare `V4`/`S3` symmetry argument cannot close the attachment.

## 6. Relation to current LRC work

THM-2795 finds a different exact `V4`/ternary boundary in tree star moves.
The common lesson is now two-thread evidence:

```text
V4/S3 can organize local moves,
but a distinguished capacity/connected-component invariant can break the
putative ternary transport.                                        (22)
```

Here the breaker is carry capacity and the present anchor; in THM-2795 it is
fixed-tree move-graph connectivity.  This is an explanatory cross-link, not
a dependency or a new Graceful Tree consequence.

## 7. Scope

The theorem does **not** prove:

- that the four configurations carry a lawful physical `S4` or `S3` action;
- a modular-group action on LRC carriers;
- an attachment after mixing different vertices labelwise;
- a no-go for a different fifth configuration or different rail bank;
- a twelve-label result for `D`, which has only eleven active labels;
- a row exclusion or LRC(14).

The theorem proves the exact local symmetry and the exact obstruction to
using it.

## 8. Exact companion

Run

```bash
python 04-computation/lrc14_v4_configuration_completion_anchor_defect_thm2804.py
python -O 04-computation/lrc14_v4_configuration_completion_anchor_defect_thm2804.py
```

Both modes byte-match the stored transcript.  The companion uses exact
integer/group/interval arithmetic, explicit exception gates, no floating
point, and no truth-bearing Python assertions.  It pins the proved THM-2672
and THM-2749 scripts; rebuilds all four unit-support rows on fourteen rails;
enumerates the affine group and its normal translation subgroup; verifies
the carry diamond, root and missing-label laws; reconstructs every `D`
facet and marked rail; and checks the first failed layer, exact mass, zero
open overlap, and zero seams.

The theorem remains a candidate until an immutable independent hostile audit
rebuilds the `V4/S3` typing, carry and root laws, `D` positivity, anchor
failure, replay, and LF hashes.

QED, conditional only on candidate status promotion after independent audit.
