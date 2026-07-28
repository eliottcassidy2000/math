---
id: THM-2749
title: "Two-sided fully marked root-zero clutch and cyclotomic-unit target window"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  rail-8 root-zero overlap, intersect the full E3/lawful-target mask with its
  translated-target pullback and insert the actual D^6 Q_(3,{1,2}) terminal
  fork plus both relative-present factors before delayed integration.  The
  source and target vectors are then exactly equal and primitive units for
  precisely t=3,...,11, with root-normalized gain -1.  Their nine-label
  target window is a norm-one cyclotomic unit with the positive three-term
  inverse z^2+z^6+z^10 modulo Phi_13, so all twelve primitive target
  characters survive.  This is a two-sided common-section coefficient
  theorem, not naturality of the single target sheet, a global transporter,
  endpoint current, row exclusion, or LRC(14).
source: root/fully-marked-root-zero-target-profile-2026-07-28
audit: root-zero-two-sided-clutch-hostile-audit-2026-07-28 (independent interval/witness, direct 239-piece forward-coordinate coefficient reconstruction, single-sheet shear boundary, unit/gain/window/Bezout/rail-hostile audit, and normal/-O/hash replay: ACCEPT)
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2742-full-two-target-present-sheet-deepest-source-semantic-current
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
related:
  - THM-830-b3-deletion-deck-mirror-current-calculus
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
script: 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
output: 05-knowledge/results/lrc14_fully_marked_root_zero_target_profile_thm2749.out
script_sha256: 12f150dc8e0fc543cc36fafaa2b84dd57a2dde6e40ce3cbadd8d057817bce3dc
output_sha256: 72fed42be733aca63fc0ccd0a907eadcb02d224c8832d2cd5c42208b34a18048
hash_basis: LF-normalized bytes
---

# THM-2749 -- the fully marked clutch exists on a two-sided common section

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2744 proves an open physical rechart through the root-zero landing, but
its last section leaves one decisive test: put the semantic source, lawful
target sheet, and actual terminal word inside the overlap integral before
testing the coefficient.  The answer depends on which categorical object is
used.  One natural target sheet has unequal source and target coefficients.
The two-sided common-section fibre product instead has an exact signed clutch,
and its full target-label profile is an integral cyclotomic unit.

## 1. The minimal physical seam

Use the canonical row

```text
(H,q1,...,q5,c1,c2,c3)
  =(1,14,27,40,53,66,13,13^3,2*13^5),
R=13^6=4826809,
T=297836897838480.                                      (1)
```

For THM-2640's open half-tooth charts

```text
H^L_r=(14r-13,14r)/182,
H^R_r=(14r,14r+13)/182,                                 (2)
```

put

```text
Sigma_-=H^R_12 intersect H^L_0=(169,181)/182,
Sigma_+=H^R_0  intersect H^L_1=(  1, 13)/182.           (3)
```

These are intervals in the scaled coordinate `182{c3 x}`.  Translation by

```text
tau=7/R                                                  (4)
```

adds `14` modulo `182`, because `c3 tau=14/13`, and maps `Sigma_-`
bijectionally to `Sigma_+`.  It fixes the delayed phase because

```text
{R(x+tau)}={Rx}.                                         (5)
```

The predecessor carry changes by `7`, hence `12 -> 6`; the edge-preserving
deep root changes `12 -> 0`; and the overlap then recharts that same target
point from right root `0` to lawful left root `1`.  Thus the cospan is

```text
(carry12,R12) --T_tau--> (carry6,R0) --rechart--> (carry6,L1). (6)
```

The corrected strict address pair

```text
q_-=47850889647341/100360982066072,
q_+=47851035194197/100360982066072=q_-+7/13^6           (7)
```

has scaled deep coordinates

```text
125553481/742586=169+56447/742586,
   799033/742586=  1+56447/742586.                       (8)
```

It lies strictly inside both seams.  Both points have semantic record
`E3 -> D^6 -> Q_(3,{1,2})`, source-one rail `8` with metadata `(1,4,12)`,
the same delayed phase, common `s`-labels
`(0,1,2,3,8,9,10,11,12)`, and common `t`-labels `3,...,11`.

## 2. The two-sided fully marked carrier

Let `F_(1,0,t)` be THM-2742's lawful full present section and put

```text
A_t=E3 intersect F_(1,0,t),
M_t=A_t intersect T_tau^(-1)A_t.                        (9)
```

The second factor in `M_t` is load-bearing: `(9)` is a two-sided
source/target fibre product, not the natural coefficient of the single sheet
`A_t`.

For delayed clock `j`, let `P^c_(j,7)` be THM-2744's relative complement and
let `rho_8` be the weighted rail-`8` profile.  In the source coordinate,
restrict to

```text
support(rho_8) intersect T_tau^(-1)support(rho_8)
 intersect P^c_(j,7) intersect T_tau^(-1)P^c_(j,7)
 intersect M_t intersect Sigma_-.                       (10)
```

Inside the delayed prefix, before applying the carry functional, intersect
with the actual terminal fork `D^6 Q_(3,{1,2})`, delayed sector zero, `h=6`,
and `kappa=1`.  Thus every semantic, target, terminal, relative-present,
rail, and seam factor occurs inside the integral.

Apply the THM-2640 delayed-carry functional with source carry `12`.  Retain
the pulled target weight separately, and independently recompute the target
after pushing `(10)` to `Sigma_+` and selecting carry `6`.  Equation `(5)`
identifies the delayed prefix.  Exact interval integration verifies the
source-coordinate and direct target-coordinate answers before any content
division or determinant.

## 3. Exact nine-label profile and mirror gain

For fixed `(ell,s)=(1,0)`, the exact result is

```text
V_t^-=V_t^+=0                         for t=0,1,2,12,

V_t^-=V_t^+=(0,C,C,C,C,C,C)          for t=3,...,11,    (11)

C=339633525654239542165440.
```

Every supported label has common grid mass `6320326320`.  The amplitude has

```text
v_13(C)=1,                    (C/26) mod13=9.            (12)
```

After division by the inherited content `26`, root normalization, and
reduction in `F_13[z]/(Phi_7)`, the two classes are

```text
source root12: (9,0,0,0,0,0),       determinant 1,
target root1:  (4,0,0,0,0,0),       determinant 1.     (13)
```

The target constant is the negative of the source constant, and

```text
4/9=12=-1 mod13.                                        (14)
```

Thus the fully marked two-sided cospan preserves the raw amplitude exactly,
while the standard root normalization contributes the mirror sign.  This is
a THM-830-type signed seam current and the central `J=-1` shadow of the
THM-2716 `C4` transporter.  It is not the whole transporter groupoid: the
objects here are half-tooth charts and only one common section has been
selected.

## 4. The target window is an integral cyclotomic unit

Let `u` denote the target-label variable and write the scalar profile of
`(11)` as `C W(u)`, where

```text
W(u)=u^3+u^4+...+u^11.                                  (15)
```

Since `W` is nonzero of degree eleven, `W(zeta_13^k)` is nonzero for every
`k=1,...,12`: otherwise the degree-twelve polynomial `Phi_13` would divide
`W`.  More strongly, put

```text
V(u)=u^2+u^6+u^10.                                      (16)
```

Direct multiplication gives the integral Bezout identity

```text
W V-1=(u^9+u^5+u-1) Phi_13(u).                          (17)
```

Hence `W` is a unit in `Z[u]/(Phi_13)` with inverse `V`, and

```text
Res(Phi_13,W)=Norm_{Q(zeta_13)/Q}(W(zeta_13))=1.        (18)
```

Equivalently, in the circular group algebra `Z[C_13]`,

```text
W*V=(3,2,2,...,2)=delta_0+2(1+u+...+u^12).             (19)
```

Thus three positive **coefficient sections** indexed by `2,6,10`
demodulate the nine-window to one spike after quotienting the uniform
target-null sector.  This is a coefficient-derived recombination of lawful
`t`-tables.  It is not three physical translates of one Boolean packet:
the fixed `E3`, seam, rail, relative, and terminal marks have not been proved
whole-packet target-covariant, and `(19)` has multiplicities two and three.

There is also an exact characteristic-`13` coefficient-algebra corollary.
Since

```text
Phi_13(u)=(u-1)^12 mod13,       W(1)=9,       V(1)=3,   (20)
```

multiplication by `W` and by `V` on
`F_13[u]/(Phi_13)` has determinant `1`.  Tensoring with either constant
root class in `(13)` gives a formal `72`-dimensional bigraded unit, again of
determinant `1`.  This is a new coefficient-algebra unit.  It is not the
inherited THM-2640 unit, a semisimple physical target current, or a canonical
endpoint decoder.

## 5. Sharp naturality and global-transporter boundaries

The two-sided factor in `(9)` cannot be dropped.  The FINITE-EXACT single-
sheet computation in
`07-reflections/lrc14-semantic-root-zero-clutch-refinement-20260728.md`
uses the natural `A_3` carrier without intersecting its translated pullback.
It gives unequal raw vectors

```text
(0,a,a,a,a,a,a),       a=1812281403506324508838080,
(0,b,b,b,b,b,b),       b=1826551436254490256030720,     (21)
```

with root-normalized profiles `5` and `9`.  Both are units, but the natural
single sheet does not intertwine.  Its normalized target/source ratio is
`9/5=7=6/12`, the carry ratio.  That multiplicative scalar cannot by itself
encode the additive odometer class: `Hom(C_13,F_13^*)=0`.  THM-2749 therefore
proves a common-section clutch, not single-sheet naturality.

The inherited THM-2744 hostile is also unchanged.  On the unmarked common
seam, target rails `1` and `3` are nonunits; eleven unit-unit rails have gain
`-1`; rail `13` has the exceptional gain

```text
8=5^3 in <5>={1,5,12,8}.                                (22)
```

On the unrestricted relative-present rows, the thirteen unit rail pairs
have thirteen distinct ratios in `F_13[z]/(Phi_7)`, none a scalar clock shift
or scalar dihedral clock transport.  The outer flanks `(0,1)/182` and
`(13,14)/182` have no opposite-edge rechart.  Finally THM-2657 still forbids
a physical order-thirteen section of the odometer extension.

The exact connection contract is:

| item | retained or lost content |
|---|---|
| source | rail `8`, carry `12`, right root `12`, `E3`, `F_(1,0,t)`, relative present, and terminal fork |
| target | translated pullback of the same marks, carry `6`, recharted left root `1` |
| map | `x -> x+7/13^6`, followed by the identity rechart on `H^R_0 intersect H^L_1` |
| preserved | physical measure, delayed phase, all factors in `(10)`, raw vector, and primitive clock unit |
| changed | carry `12 -> 6`, chart `(R,12)->(R,0)->(L,1)`, normalized gain `-1` |
| not supplied | single-sheet naturality, all-rail functor, canonical arm, global `C_13` action, owner/root endpoint current, or row decrement |
| needed sidecar | a coherent common-selector/gain section into the canonical endpoint current, or a lawful physical realization of the coefficient decoder `(19)` |

No global transporter, target/physical diagonal, endpoint current, scalar-
ledger decrement, row exclusion, or proof of LRC(14) follows.

## 6. Exact reproduction

Run

```bash
python 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
python -O 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_fully_marked_root_zero_target_profile_thm2749.out.
```

The companion checks seven direct dependencies by LF-normalized SHA-256,
then reconstructs the two-sided semantic carrier, inserts the terminal fork
inside every delayed prefix, recomputes the pulled and direct target vectors,
scans all thirteen target labels, checks both determinants and the gain, and
verifies `(17)`--`(20)` by exact integer arithmetic.  The factored target-
label scan is independently compared with the slower full Boolean
construction at `t=3`.  No truth-bearing Python `assert` is used.

An independent hostile audit separately reconstructed the target coefficient
from `239` exact forward-coordinate pieces, checked the strict overlap witness,
all marked factors, content and valuation, both private-root determinants,
the nine-label support, all primitive target characters, the Bezout decoder,
the unrestricted-rail hostile, and the distinction from the natural
single-sheet shear.  Its normal and optimized runs byte-match the stored
`22`-line transcript and the declared LF hashes.

QED.
