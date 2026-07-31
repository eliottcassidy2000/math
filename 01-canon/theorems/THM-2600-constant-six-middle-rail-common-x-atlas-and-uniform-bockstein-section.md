---
id: THM-2600
title: "Constant-six middle-rail common-x atlas and uniform primitive Bockstein section"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, including the
  projective owner-cycle coefficient-support addendum.  On the
  canonical typed row, the two arrival-root-six middle edges of THM-2584's
  depth-five toothpick have disjoint three-cell zero sets and together cover
  all 84 nonzero-displacement/owner-clock cells.  Pulling both edges against
  all thirteen THM-2585 target sections on THM-2592's literal common-x
  carrier gives one globally primitive 162x13 bank.  Every base cell admits
  a positive unit-Bockstein attachment.  More strongly, the single target
  section q=0 works in all 84 cells: take the theta-one edge t=0 exactly for
  s in {6,11} and (s,ell)=(8,2), and the theta-zero edge t=12 otherwise.
  The selected 84-slice bank has the same global content as the full bank and
  all 84 septimal factors are units.  This keeps arrival=future digit six on
  one physical fibre product, but sacrifices the deep diagonal on 15 cells,
  does not identify the two clock labels, and supplies neither a semantic old
  head nor a transition cochain.  No row is excluded and LRC(14) remains open.
source: codex-2026-07-28-constant-six-middle-rail
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2592-fallback-rail-digit-diagonal-pullback-and-primitive-bockstein
  - THM-2603-hurwitz-projective-root-owner-atlas-and-nonabelian-seven-edge-trivialization
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2586-depth-five-arrival-to-future-root-diagonal
  - THM-2590-boolean-bockstein-and-theta-selector-incidence-spectrum
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
script: 04-computation/lrc14_constant_six_middle_rail_pullback_thm2600.py
output: 05-knowledge/results/lrc14_constant_six_middle_rail_pullback_thm2600.out
script_sha256: c85ecd26df053ad68362970d3d4056ad09d3668ba717fba5356a6968dff2e95e
output_sha256: b89d10055b53b6cb43de216c11fce3bd198f9a4af08de6fa86a2e1af1cd61e41
hash_basis: working-tree bytes (LF)
---

# THM-2600 -- the middle vertex attaches every depth-five cell

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2592's priority theta-zero selector attaches only three fallback cells.
Its other `81` cells use arrival digit zero and are orthogonal to the entire
delayed word, not merely to its digit-zero piece.  That hostile is sharp for
those physical rail cells, but it does not exhaust THM-2584's four-edge
toothpick.  The shared middle arrival vertex has digit six and two incident
deep edges.  Keeping both edges changes the physical rail before the delayed
word is inserted:

```text
theta-zero middle edge:  (v,t)=(6,12),   w=7t=6=v;
theta-one middle edge:   (v,t)=(6, 0),   w=7t=0=v+7.       (1)
```

Their holes are disjoint.  Since the delayed word has positive digit-six
mass in every owner phase, the two-edge atlas supplies the uniform common
carrier that the one-edge priority rule misses.

## 1. The two middle edges cover all 84 cells

Use THM-2584's canonical typed row and depth-five packet

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),

sigma={b},                 d=13^5,                 R=13^6. (2)
```

Its fine positive tensor is

```text
K_ell(s,v,t)
 =13 integral U(x)V(x-s/13)
       1_(floor(13x)=v)
       1_(floor(26x)=t mod 13)
       1_(frac(169x) in win_ell) dx.                     (3)
```

Here `s in F_13^*` and `ell in F_7`.  Exact reconstruction of all `15,379`
fine cells gives

```text
Z_(6,0) ={(7,1),(7,2),(7,3)},

Z_(6,12)={(6,4),(6,5),(6,6)}.                            (4)
```

The sets in (4) are disjoint.  Consequently every one of the `12*7=84`
base cells supports at least one middle edge, and the complete positive
middle-edge bank has

```text
81 theta-zero edges + 81 theta-one edges = 162 edges.     (5)
```

This is not a resection of THM-2592's primary `v=0` rail.  It replaces that
rail by a different positive THM-2584 cell with `v=6`, before either tensor
is marginalized.

## 2. The literal common-x pullback

Retain THM-2592's present `sigma={a}` packet

```text
F_(ell5,s5)(x),             Delta_r(x),

Q^+_(ell5,h)(Rx)
 =Q^+_ell5(Rx) 1_(floor(13 frac(Rx))=h).                 (6)
```

Its exact digit support obeys, for every `ell5 in F_7`,

```text
Q^+_(ell5,0)=Q^+_(ell5,12)=empty,

measure(Q^+_(ell5,6))>0.                                 (7)
```

For every positive middle edge `(s,ell4,6,t)`, target position
`q in F_13`, present owner cell `ell5`, and deep probe `r`, define

```text
J_(s,ell4,t,q;ell5,r)
 =13 integral U(x)V(x-s/13)
       1_(ell4,6,t)(x)
       F_(ell5,-q)(x) Delta_r(x)
       Q^+_(ell5,6)(Rx) dx.                              (8)
```

All factors in (8) are evaluated at the same physical `x` before any sum or
Fourier transform.  Expanding the two Perron profiles gives a nonnegative
finite sum of literal Boolean intersections, exactly as in THM-2592.  The
two digit-six occurrences give the physical identity

```text
floor(13x)=6=floor(13 frac(13^6 x)).                     (9)
```

Thus positivity of one entry of (8) supplies an actual common-ancestry
Boolean sheet, not a coupling inferred from equal labels or marginals.  The
Bockstein below belongs to the globally primitive aggregate slice (8); it
does not canonically select that hidden positive Perron sheet.

## 3. One global primitive atlas

Put every `162*13*7*13` entry of (8), including the identically zero `r=0`
cells, over the common exact denominator.  Before the route-two factor `13`,
the nonzero integer numerators have global content

```text
g0=4,244,240,                 nu_13(g0)=1.                (10)
```

Including that factor, the physical raw content and primitive denominator
are

```text
g=55,175,120,

D_prim=607,037,547,933,467,614,874,742,741.              (11)
```

No edge, target section, or owner cell is reprimitivized separately.  The
full exact support census is

```text
positive edge/section pairs:            2,024 / 2,106,
positive fine (edge,q,ell5,r) entries: 61,248 / 176,904,

fine supports per pair:
  0:82, 12:154, 24:814, 36:902, 48:154.                  (12)
```

For each primitive slice define

```text
Y_ell5=sum_(r=1)^12 a_(ell5,r) r^(-1) mod 13,

Y(z)=sum_(ell5=0)^6 Y_ell5 z^ell5
       in R_7=F_13[z]/(Phi_7(z)).                         (13)
```

Exactly `1,740/2,106` slice polynomials in (13) are units.  There are `364`
zero-Bockstein slices and exactly two further nonzero nonunits.  Nevertheless
every one of the `84` base cells has at least one positive unit slice.  The
number of base cells admitting a unit at each fixed target position is

```text
q:        0  1  2  3  4  5  6  7  8  9 10 11 12
cells:   84 84 81 83 83 83 84 84 83 83 83 81 84.        (14)
```

Thus not every target section is uniform: `q=2` and `q=11` are sharp
three-cell hostiles.  But five fixed sections, including `q=0`, cover all
`84` cells.

## 4. A uniform q=0 unit selector

Fix the one target-shift position

```text
q=0.                                                       (15)
```

For every `(s,ell) in F_13^* x F_7`, select

```text
t_(s,ell)=0,   if s in {6,11} or (s,ell)=(8,2),

t_(s,ell)=12,  otherwise.                                 (16)
```

Every selected edge exists in (3), every selected pair (8) has positive
mass, and every selected polynomial (13) is a unit of `R_7`.  Hence its first
Bockstein

```text
beta=Omega Y(zeta_7^kappa)                               (17)
```

is nonzero for all six labelled owner colours
`kappa in F_7^*`.  In fact unitness is stronger: no nonzero downstream
septimal coefficient can annihilate a selected factor.

The `84` selected slices have global content exactly `g0`, the same content
as the entire `162*13` bank.  Thus (16) does not hide a second primitive
division.  Its primitive serialization digest is

```text
ab59b257dd3a9788d6108a46f1e01705caa2598fc9c47cee95f26165a0242110.
                                                                    (18)
```

The selector uses the theta-zero diagonal edge on `69` cells and the
theta-one edge on `15`.  On those `15` cells arrival and future digits still
equal six by (9), but the rescaled deepest root is `w=0`, not `w=v`.  This is
the exact price of uniform attachment.

### 4.1 Projective owner-cycle coefficient support

THM-2603 puts the order-seven owner map into the affine root chart

```text
A(q)=(7q+5)/(10q+11)
```

with projective cycles

```text
O0=(0,4,6,10,7,5,3),
O1=(1,8,infinity,2,9,12,11).                              (19)
```

There is no proved identification between that projective root coordinate
and the target-section coordinate `q` of this theorem.  Nevertheless, a
cheap hostile test is exact: treat the labels as a proposed interface, rotate
each cycle against the seven clock positions, and ask at every **finite**
vertex only whether the corresponding cell has a positive unit rail.  For
each physical displacement `s`, minimize the number of forced theta-one
choices over the seven phase rotations.

The first cycle has support at all seven vertices for every `s`.  The second
has support at all six finite vertices; `infinity` is not a zero but an absent
target section, so it is outside the test.  In the order `s=1,...,12`, the
minimum theta-one invoices are

```text
O0:        (0,0,0,0,0,5,0,1,0,0,7,0),
O1 finite: (0,0,0,0,0,4,0,1,0,0,6,0).                   (20)
```

Thus, for each cycle separately, nine displacement lanes admit theta-zero
unit support at every tested vertex.  The exceptional set is exactly

```text
{6,8,11},                                                 (21)
```

the same displacement set on which the uniform `q=0` selector (16) uses a
theta-one rail somewhere.  This is a structural positive signal: the
coefficient bank nearly follows both projective owner cycles, and its only
unrepresented projective vertex is the fourteenth boundary point already
isolated by THM-2603.

The quantifiers make the limitation equally sharp.  The phase is minimized
separately for each displacement and cycle; (20) is vertex support, not one
ordered seven-edge fibre product.  It supplies no `q -> q'` ancestry gluing,
no physical root/projective intertwiner, no value at `infinity`, and no
semantic endpoint.  It therefore does not define the transition kernels of
THM-2602 or change the LRC verdict.

An independent hostile audit rederived both projective cycles, checked every
phase rotation and both invoice vectors, treated `infinity` as absent rather
than zero, and independently byte-matched normal and optimized runs to the
stored transcript.

## 5. Sharp hostiles and loss ledger

The theorem closes the physical carrier gap between THM-2584's middle rail
and THM-2585's primitive target sections on this packet.  It does not close
the semantic or holonomy gaps.

| item | retained or lost content |
|---|---|
| source | one positive common-`x` Boolean fibre-product aggregate for each selected cell |
| selected arrival/future digit | the literal constant value `v=h=6` |
| target section | the fixed absolute shift position `q=0`, not a target character |
| deep sidecar | `t` is retained, but `w=v` fails on the 15 theta-one cells |
| clocks | `ell4` and `ell5` are both retained and are not identified |
| semantic loss | no selected empty old head, relation residue, or named endpoint |
| cohomological loss | (16) is vertex selection, not a mixed-square edge cochain |

The positive-mass witness sheet and the aggregate unit Bockstein are two
facts on the same pre-marginal carrier, but the latter does not identify the
former.  Reusing its coefficient as a sheetwise charge would repeat the
hidden-sheet error already excluded in THM-2592's scope ledger.

Four controls make the boundary sharp.

1. Digit `0` and digit `12` are empty in every delayed owner phase, while
   digit `6` is positive in all seven.
2. The theta-zero middle edge alone has two base cells with no target-section
   attachment; the theta-one middle edge alone also has two.  Their union is
   essential.
3. The complete bank contains zero Bocksteins and nonunit nonzero factors;
   arbitrary edge/section selection is not safe.
4. THM-2591 still applies to a root chosen chart by chart: it is a zero-cochain
   and cannot by itself pay the required seven-edge holonomy.  The fixed
   `q=0` atlas supplies a common physical input for a future mixed square, not
   that transition boundary.

The row (2) is a canonical typed non-cover control.  Nothing here transports
the construction to the other `164` rows, identifies a semantic old head,
produces a THM-2334 current, excludes a scalar row, or proves LRC(14).

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_constant_six_middle_rail_pullback_thm2600.py
python3 -O 04-computation/lrc14_constant_six_middle_rail_pullback_thm2600.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_constant_six_middle_rail_pullback_thm2600.out.
```

The companion independently rebuilds THM-2584's fine tensor, checks both
middle-edge zero sets, reconstructs the delayed word and digit gates, forms
every entry of (8), removes one global content, evaluates every first
Bockstein and multiplication determinant, checks (12)--(18), freezes the
positive and hostile controls, and records the primitive digests.  Its
weighted delayed sweep is cross-checked against the original grouped sweep
and its deep-comb restriction against an independent direct intersection.

An independent hostile audit rederived both disjoint middle-edge zero sets,
the `81+81` cover, the exact `15/69` selector split, the global content and
all `84` unit tests.  It separately checked that `q=0` is an absolute target
section, that `v=h=6` is physical on one common-`x` carrier, and that the
deep diagonal is lost on exactly the fifteen `t=0` cells.  Normal and
optimized executions byte-match the stored transcript after LF
normalization, with the declared hashes.  No mathematical or scope defect
remains.

No row is removed and LRC(14) remains open. **QED.**
