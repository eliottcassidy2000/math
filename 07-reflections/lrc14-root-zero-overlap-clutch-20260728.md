# The zero-root clutch is an open chart overlap

> **STATUS: FINITE-EXACT + VERIFIED on one canonical semantic pair, with
> uniform adjacent-chart geometry and a fourteen-rail coefficient audit.**
> The relative-present scalar lift has no edge-preserving nonzero-root edge,
> but it does admit a partial physical chart isomorphism through the open
> overlap of adjacent half-tooth charts.  On rail `8` the overlap coefficient
> is recomputed identically on both sides and is a primitive unit at source
> root `12` and target root `1`.  This is not a global action, endpoint
> current, row exclusion, or proof of LRC(14).

## 1. The old no-go was correctly typed but not terminal

The relative-present computation leaves two semantic residue classes:

```text
residue 7: carry 12, right edge, root 12;
residue 8: carry  6, left  edge, root  1.
```

A scalar lift from residue `7` to residue `8` has root increment `+1`.
If its right edge is retained, it sends root `12` to root `0`; the reverse
left-edge lift sends root `1` to root `0`.  Therefore the edge-preserving
private-root relation is indeed empty.

The missing object is the overlap between two charts of the same physical
point.  In the conventions of THM-2640, write

```text
H^L_r=(14r-13,14r)/182,       epsilon=1,
H^R_r=(14r,14r+13)/182,       epsilon=0.                 (1)
```

The script field called `edge` uses the opposite numeric convention:
`edge=0` is left and `edge=1` is right.  Keeping those two conventions
separate prevents a cosmetic edge switch from being mistaken for a physical
identity.

For every `r mod13`, the adjacent charts overlap by

```text
H^R_r intersect H^L_(r+1)
  =(14r+1,14r+13)/182,                                 (2)
```

an open strip of width `12/182=6/91`.  Each half-tooth has only a one-unit
outer flank which does not admit the opposite-edge chart.

## 2. The canonical lift lands far inside the overlap

Put

```text
R=13^6=4826809,                 tau=7/R.
```

Since `c3=2*13^5`, translation by `tau` moves the deep coordinate
`z=182{c3 x}` by exactly `+14 mod182`.  On the relative-present fibre,

```text
y={Rx}=11195237/20792408,
eta=2y-1=799033/10396204.
```

The source and target coordinates are

```text
z =168+14 eta=125553481/742586
  =169+56447/742586,

z'=     14 eta=799033/742586
  =  1+56447/742586.                                  (3)
```

Thus

```text
H^R_12 intersect T_tau^(-1)H^L_1=(169,181)/182,

T_tau((169,181)/182)=(1,13)/182
 =H^R_0 intersect H^L_1.                              (4)
```

The inverse translation gives the reverse overlap law.  The physical point
is fixed by the rechart at the target, while its chart label changes from
`(right,root0)` to `(left,root1)`.

The lower overlap margin in the original circle coordinate is

```text
56447/100360982066072
 =56447*(1/100360982066072).                           (5)
```

The factor in parentheses is the inherited open-cylinder radius.  Hence the
entire common cylinder, not just its midpoint, lies in the chart overlap.
Equivalently, the switch threshold is `eta>1/14`; here

```text
eta-1/14=56447/10396204>0.                             (6)
```

This also explains why THM-2672's old maximal-face cap is not contradicted.
Its canonical missing-label strip begins strictly beyond deep coordinate
`181`, on the one-edge flank outside `(169,181)`.  The new E3 fibre lies on
the other side of exactly that boundary.

## 3. A strict semantic and full-target witness

The exact adjacent address pair is

```text
n=6715:
q =47850889647341/100360982066072,

n=6716:
q'=47851035194197/100360982066072=q+7/R.              (7)
```

Both points have semantic record

```text
E3 -> D^6 -> Q_(3,{1,2}).                              (8)
```

They lie on source-one rail `8`, whose metadata is `(1,4,12)`, in the same
weighted segment

```text
(141869054887560,142120818960120),
weight 27581135604.                                    (9)
```

Their reverse-clock pair is `(shallow,owner)=(1,4)`, their delayed coordinate
is identical, and the active delayed clocks are exactly `1,...,6`.
Both lie in the clock-`1`, label-`7` relative present complement.

The two full target-label banks are close but not identical:

```text
source s=(0,1,2,3,8,9,10,11,12),
source t=(3,4,5,6,7,8,9,10,11),

target s=(0,1,2,3,8,9,10,11,12),
target t=(3,4,5,6,7,8,9,10,11,12).                    (10)
```

Their common bank contains all source labels, in particular the lawful
choice `(s,t)=(0,3)`.  Each endpoint's own label bank is constant on its
whole inherited open cylinder.  The extra target label `t=12` is retained as
a sidecar difference rather than silently identified with the source bank.

## 4. The coefficient is recomputed before the determinant

Let `rho_j` be the weighted profile of rail `j`, let `P_ell^c` be the
complement of the old one-target present packet, and let `T_tau` be
translation by `7/R`.  In the source coordinate, restrict to

```text
support(rho_j)
 intersect T_tau^(-1)support(rho_j)
 intersect P_ell^c
 intersect T_tau^(-1)P_ell^c
 intersect H^R_12
 intersect T_tau^(-1)H^L_1.                            (11)
```

Source and pulled-target rail weights are retained separately.  Apply the
exact THM-2640 delayed-carry functional to `(11)` with source carry `12`.
Translate the restricted target profile forward by `tau` and independently
apply the same functional with target carry `6`.  Because

```text
{R(x+7/R)}={Rx},                                       (12)
```

the delayed phase is unchanged; equality is nevertheless checked on the two
interval banks before any determinant or unit test.

On rail `8`, every overlapping source/target step has the same weight, so no
extra weight-equality cut is needed.  The two exact seven-clock vectors are

```text
A^-=A^+=(0,a,a,a,a,a,a),

a=5359949020444386606638400.                           (13)
```

Every entry is divisible by the inherited canonical content `26`, and

```text
(a/26) mod13=11.
```

After the THM-2640 root normalization, the reduced `Phi_7` coordinates are

```text
source root 12: (11,0,0,0,0,0),
target root  1: ( 2,0,0,0,0,0).                        (14)
```

Both are nonzero constants in `F_13[z]/(Phi_7)`, so both multiplication
determinants are nonzero.  The overlap vector is therefore a private unit at
both nonzero endpoint roots.

There is also a transparent characteristic-zero clock statement.  If

```text
P(X)=a(X+X^2+...+X^6),
```

then

```text
P(X)=-a mod Phi_7(X),
P(zeta_7^k)=-a for k=1,...,6.                           (15)
```

Thus every nontrivial seven-clock character survives on the overlap carrier.
This is a clock character, not the deepest `t`-character of THM-2742 and not
a physical deck character.

## 5. Exact extent of the construction

The result has three different uniformity levels which must not be blended.

1. The adjacent-chart identity `(2)` holds uniformly for all thirteen roots.
2. On the fourteen source-one rails, the uncut overlap vectors agree on
   exactly

   ```text
   (0,4,5,6,7,8,9,10,11,12),                           (16)
   ```

   ten rails in all.
3. If each rail is further restricted to its exact equal-weight locus, all
   fourteen rails give a positive vector of shape `(0,a_j,...,a_j)` and a
   unit at roots `12` and `1`.  Rail `8` is stronger: its whole overlap has
   equal weights and contains the explicit semantic pair `(7)`.

The third statement is a lawful subcarrier statement, but the equal-weight
cut on four rails is derived data.  It is not promoted to a canonical global
action.

## 6. Hostile controls and the sharp surviving no-go

Three controls prevent overstatement.

- Without the overlap restriction, the canonical rail-`8` relative rows
  have support sizes `7` and `6` and are unequal.  Full-half transport still
  does not intertwine.
- The strips `H^R_0\H^L_1=(0,1)/182` and
  `H^L_1\H^R_0=(13,14)/182` admit no opposite-edge rechart.  The map is only
  partial.
- The overlap map is the identity on the physical point but not on chart
  labels: it changes right to left and root `0` to root `1`.  It therefore
  does not refute the edge-preserving zero-root no-go.

The strongest proved conclusion is

```text
semantic scalar lift
 + relative-present endpoint support
 + open adjacent-tooth overlap
 -> partial physical chart isomorphism
 -> equal recomputed overlap vector
 -> private units at roots 12 and 1.                    (17)
```

No endpoint transition amplitude, global `C_13` action, target/physical
diagonal, row exclusion, or LRC(14) conclusion follows.

## 7. Cheapest next refinement

The next computation should make `(17)` semantic inside the integral rather
than only at the explicit endpoint.  Retain one common lawful full-target
label, for example `(s,t)=(0,3)`, and intersect both overlap carriers with

```text
E3 intersect F_(ell,0,3) intersect D^(-6)Q_(3,{1,2})    (18)
```

before applying the delayed-carry functional.  Then recompute source/target
equality and both content-`26` unit tests.  The current companion verifies
that all factors in `(18)` coexist on a common open cylinder, but it does not
yet prove that the globally integrated vectors after this extra cut remain
equal or units.

## 8. Exact reproduction

Run normally and with optimization:

```bash
python3 04-computation/lrc14_root_zero_overlap_clutch_20260728.py
python3 -O 04-computation/lrc14_root_zero_overlap_clutch_20260728.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_root_zero_overlap_clutch_20260728.out.
```

SHA-256:

```text
script  e27981478cd30c8e3cceada128049b145b254410c8d0b6d525a8a1830545d55f
output  ba9d0a67dfede0b64cf97ff55af7e86c9bb46c962c669d33015df7e574e8e91e
```

The companion pins the six direct carrier dependencies, checks all thirteen
adjacent overlaps, both translation directions, the exact semantic witness
and cylinder margins, all fourteen source-one rails, equality before the
determinant, both unit tests, and the three hostile controls.  It uses no
truth-bearing Python `assert`.
