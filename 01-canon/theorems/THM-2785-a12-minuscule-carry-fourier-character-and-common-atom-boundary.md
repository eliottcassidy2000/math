---
id: THM-2785
title: "A12 minuscule carry Fourier character and common-atom boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT HOSTILE AUDIT.
  The thirteen THM-2672 missing-facet incidence vectors project to the
  regular minuscule omega_12 orbit in the A12 standard plane, but every one
  has the same class -gamma in P(A12)/Q(A12).  Thus the physical carry is
  not that quotient.  The selected positive component intervals nevertheless
  have a genuine ordinary spatial Fourier descent: whenever the base
  component coefficient is nonzero, an integer frequency N gives a C13
  character exactly when 13^5 divides N.  Frequencies
  N=a*13^5, 1<=a<=12, realize all twelve nontrivial dual characters; the
  physical deep speed c3=2*13^5 realizes -2gamma and multiplies by zeta^-1
  under the slope-seven step.  This works because the Fourier modes
  annihilate the THM-2657 1/13^5 extension kernel; it does not geometrically
  glue the components.  The THM-2749-typed rail-eight reconstruction has
  the same spectrum, but its selected source-root-12 component is exactly
  disjoint from E3 and S_(1,0,4), while both its actual translate and the
  selected carry-six component are disjoint from the fully marked target
  carrier.  Hence no common
  semantic atom, THM-2334 endpoint current, allocation K4, row exclusion,
  or LRC(14) conclusion follows.
source: lrc-a12-carry-bridge/minuscule-fourier-2026-07-28
depends_on:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects
related:
  - THM-2688-affine-facet-normalization-vertical-sphere-diagonal-simplex-and-lens-quotient
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
script: 04-computation/lrc14_a12_minuscule_carry_fourier_boundary_thm2785.py
output: 05-knowledge/results/lrc14_a12_minuscule_carry_fourier_boundary_thm2785.out
script_sha256: e1b81f994d194bce1f1c846b85e2cec6ecf96f5fc0568c35adcd943284175478
output_sha256: 8570a9a6a7a5e2f5c5004026b0e3cb156edc5fde7a9814f6b566e4276b7c1317
hash_basis: LF-normalized bytes
---

# THM-2785 -- the carry orbit is minuscule, while its physical character is dual

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT HOSTILE
AUDIT.**

THM-2774 makes `P(A12)/Q(A12)=Z/13` available as an exact lattice
coordinate.  THM-2672 independently has thirteen predecessor-carry strata.
The tempting identification is false in the most rigid possible way:

```text
the thirteen facet weights form the full minuscule A12 orbit,
but all thirteen lie in one and the same P/Q class.                    (1)
```

There is nevertheless a genuine physical bridge on the dual side.  The
selected positive component intervals in the thirteen carry strata are
equal translates modulo the THM-2657 extension kernel.  Ordinary spatial
Fourier modes divisible by `13^5` annihilate that kernel and recover all
twelve nontrivial `C13` characters.  At the actual deep blocker speed
`c3=2*13^5`, the character is `-2 gamma`.

The strongest natural attachment test then fails.  The rail-eight
configuration matching THM-2749's source/target roots has the same physical
character, yet its selected source component misses `E3` and the fully
marked semantic section exactly.  Fourier survival therefore does not
supply a common endpoint atom.

## 1. The thirteen physical facets

Use THM-2672's scales

```text
p=13,             R=13^6,             S=13^5,
c3=2S=742586.                                                        (2)
```

In the canonical fixed configuration

```text
(rail,sector,edge,kappa,h)=(0,0,0,0,6),                              (3)
```

let `F_c` be the positive twelve-label facet selected after choosing source
carry `c in F_13`.  Its omitted label is

```text
m(c)=2(6-c)=12-2c                         in F_13.                    (4)
```

Let `I_c` be the lexicographically first open physical component used by
THM-2672 in that facet.  Every `I_c` lies in delayed clock one and has length

```text
ell_0=3/(7*13^11)=3/12545122758259.                                  (5)
```

The words *facet*, *carry*, and *component* remain distinct below.
`F_c` is a twelve-label intersection, `c` is its predecessor-carry stratum,
and `I_c` is one selected open component of that intersection.  The coarse
virtual `S^11` of THM-2672 is not used.

## 2. The direct type-A12 map collapses

Let

```text
V={x in R^13:sum_i x_i=0},
Q={x in Z^13:sum_i x_i=0},
P={x in V:x_i-x_j in Z for all i,j}.                                  (6)
```

These are the standard plane, root lattice, and weight lattice of `A12`.
Put

```text
lambda_i=e_i-(1/13)1,                 gamma=[lambda_0] in P/Q.         (7)
```

The label-incidence vector of `F_c` is

```text
v_c=1-e_(m(c)) in Z^13,                sum_i (v_c)_i=12.               (8)
```

Projecting orthogonally to `V` gives

```text
w_c=v_c-(12/13)1=(1/13)1-e_(m(c))=-lambda_(m(c)).                     (9)
```

### 2.1 Minuscule regular simplex

As `c` ranges over `F_13`, `(4)` ranges over all labels, so the `w_c` are
the thirteen weights of the minuscule `omega_12` orbit.  They satisfy

```text
<w_c,w_d> = 12/13  if c=d,
           -1/13  if c!=d.                                           (10)
```

Thus they form a regular twelve-simplex centred at zero.  Every difference
is an `A12` root.  In particular, the slope-seven carry step obeys

```text
m(c+7)=m(c)-1,
w_(c+7)-w_c=e_(m(c))-e_(m(c)-1).                                    (11)
```

It is the Coxeter thirteen-cycle on this weight orbit.

For a primitive thirteenth root `zeta`, define the unnormalized
coefficient-space DFT

```text
W_k=sum_(c in F_13) zeta^(-kc) w_c.                                  (12)
```

The cyclic Gram sums from `(10)` are

```text
A_0=sum_c <w_c,w_c>=12,
A_d=sum_c <w_c,w_(c+d)>=-1                    (d!=0).                 (13)
```

Therefore, for every `k!=0`,

```text
||W_k||^2=sum_d zeta^(kd) A_d=12-sum_(d!=0)zeta^(kd)=13.             (14)
```

All twelve coefficient-space colours survive.  Relabelling the orbit by
`c -> c+7` multiplies `W_k` by `zeta^(7k)`.  This is an exact geometric
property of the regular simplex; it does not yet assert any physical
endpoint coefficient.

### 2.2 Constant central character

Equation `(9)` also gives

```text
[w_c]=-[lambda_(m(c))]=-gamma                    in P/Q               (15)
```

for every `c`.  The Weyl group `S_13` permutes the `lambda_i` and acts
trivially on `P/Q`, so no relabelling repairs this collapse.  Consequently

```text
c |-> [w_c]
```

is constant, not an isomorphism from the physical carry to
`P(A12)/Q(A12)`.  A formula such as `c |-> c gamma` is an external affine
marking: it requires a chosen carry origin and generator and is not induced
by the facet-incidence weights.

This is the direct-map no-go promised in `(1)`.

## 3. Exact interval section and physical Fourier descent

The selected components have more structure than their incidence vectors.
With carry zero as the displayed origin, exact interval arithmetic gives

```text
I_c=I_0+(c-13)/R,                       c=1,...,12.                    (16)
```

For an integer spatial frequency `N`, define the genuine component-indicator
Fourier coefficient

```text
C_c(N)=integral_(I_c) exp(-2 pi i N x) dx.                            (17)
```

Suppose `C_0(N)!=0`.  By `(16)`, the internal ratio is

```text
C_(c+1)(N)/C_c(N)=exp(-2 pi i N/R)                  (c=1,...,11).     (18)
```

The ratios in `(18)` extend consistently across `c=0` to a `C13` character
exactly when their thirteenth power is one:

```text
exp(-2 pi i 13N/R)=1
 iff S divides N.                                                       (19)
```

Conversely, if `N=aS`, then `(16)` gives the exact character law

```text
C_c(aS)=zeta^(-ac) C_0(aS),                 c in F_13.                (20)
```

For each `a=1,...,12`,

```text
0<aS ell_0=3a/(7*13^6)<1,                                      (21)
```

so the interval integral in `(17)` is nonzero.  Thus the physical component
bank realizes every nontrivial dual character:

```text
N=a*13^5             corresponds to the character class -a gamma.     (22)
```

Here `(22)` uses the standard dual pairing fixed by the displayed carry
origin and `zeta`.  Changing the spatial origin multiplies all thirteen
coefficients by one common phase; it does not change their character ratios.

At the actual deep frequency `N=c3=2S`,

```text
C_c(c3)=zeta^(-2c) C_0(c3),

|C_0(c3)|
 =sin(6 pi/(7*13^6))/(2 pi 13^5)>0.                                  (23)
```

Under the physical slope step `c -> c+7`, `(23)` is multiplied by

```text
zeta^(-14)=zeta^(-1).                                                  (24)
```

The sign in `(24)` is evaluation of the physical character along the carry
step.  It is the inverse convention to the relabelling eigenvalue of the
coefficient-space vector `(12)`; the two must not be silently identified.

## 4. Why the Fourier character exists despite the odometer cocycle

The interval section `(16)` is not geometrically equivariant under the
slope-seven lift.  Exact comparison gives

```text
I_c+7/R = I_(c+7)+epsilon_c/S,                                       (25)

(epsilon_0,...,epsilon_12)
 =(1,0,0,0,0,0,0,1,1,1,1,1,1).
```

Hence

```text
sum_c epsilon_c/S=7/S.                                                (26)
```

This is exactly THM-2657's thirteen-step odometer translation and cocycle
class seven.  It obstructs cyclic gluing of the actual components.

For `N=aS`, however,

```text
exp(-2 pi i N epsilon_c/S)=exp(-2 pi i a epsilon_c)=1.                (27)
```

Thus `(20)` is the Fourier quotient which annihilates the extension kernel.
This is controlled forgetting:

```text
preserved: component mass, carry character, slope eigenvalue;
forgotten: exact component origin modulo 1/S and semantic membership.  (28)
```

Equation `(27)` explains the positive result; it does not trivialize the
physical cocycle or construct transition maps between the components.

## 5. Rail eight carries the same spectrum

Now use the THM-2749-typed configuration

```text
(rail,sector,edge,kappa,h)=(8,0,1,1,6),
rail metadata=(1,4,12).                                                (29)
```

All thirteen carry-rebased facets again have positive selected components,
all in delayed clock one, and the root row is

```text
(r_0,...,r_12)=(1,3,5,7,9,11,0,2,4,6,8,10,12).                       (30)
```

Their common selected length is

```text
ell_8=1/(28*13^11)=1/50180491033036,                                  (31)
```

and they obey the same translation law `(16)`, minuscule typing
`(9)--(15)`, Fourier descent `(19)--(22)`, and odometer law `(25)--(27)`.
At `c3`,

```text
|C_0^(8)(c3)|
 =sin(pi/(14*13^6))/(2 pi 13^5)>0.                                   (32)
```

The source/root-12 and carry-six/root-zero selected components are

```text
I_12=(
 70792935919/148462991222,
 23928012340623/50180491033036),                                     (33)

I_6=(
 70792751371/148462991222,
 23927949963399/50180491033036).                                     (34)
```

The root zero in `(34)` is the private root in the pulled-back
fixed-configuration coordinate.  THM-2749's physical target translation
changes the local target-root label; `I_6` itself is not silently renamed
the translated source component.

Their positive full-facet numerators are respectively

```text
8028307708883516100046000,
8066356743663093598282560.                                            (35)
```

On all of `I_12`, both the source rail and the pulled-back target rail have
the same positive weight

```text
27581135604.                                                          (36)
```

Moreover

```text
I_12+7/R=I_6+1/S.                                                     (37)
```

The `1/S` discrepancy in `(37)` is invisible to `(32)` because `c3/S=2`.
Thus this is a particularly sharp positive Fourier control on the very rail
used by the semantic clutch.

## 6. The common semantic atom is absent

Let `E3` be THM-2749's exclusive deepest-source cell and

```text
S_(1,0,4)=E3 intersect F_(1,0,4)                                     (38)
```

its clock-one full two-target source section.  Let `C_source` and
`C_target` be the clock-one fully marked rail-eight carriers from
THM-2749 `(3)`, retaining both rail copies, both relative-present
complements, the two root halves, both semantic sections, and the delayed
fork.

THM-2749 supplies the positive control

```text
A^-_(8,1)(0,4)=A^+_(8,1)(0,4)
 =339633525654239542165440>0.                                        (39)
```

So the marked carrier itself is not empty.  Nevertheless direct exact
intersection gives

```text
mu(I_12 intersect E3)=0,
mu(I_12 intersect S_(1,0,4))=0,
mu(I_12 intersect C_source)=0,                                       (40)

mu(I_6 intersect C_target)=0,
mu((I_12+7/R) intersect C_target)=0.                                 (41)
```

In fact the first failure is already `I_12 intersect E3=empty` up to the
strict half-open boundaries.  Thus the rail, clock, source root, target root,
equal rail weight, carry character, and odometer phase can all match while
the physical semantic support does not.

This is the decisive hostile for reserved THM-2782's proposed rail-eight
common-ancestry attachment.  The cheapest transplant test is not another
Fourier or quotient calculation.  It is the literal intersection

```text
I_c intersect rail intersect T_tau^(-1)rail
    intersect P^c intersect T_tau^(-1)P^c
    intersect S intersect T_tau^(-1)S,                                (42)
```

before any coefficient integration.  For the canonical source/root-12
component, `(42)` fails already at `I_12 intersect S`.

The theorem does not claim that every component in the full rail-eight
configuration misses the semantic carrier.  It proves the exact failure of
the selected positive component section which carries `(20)` and `(32)`.

## 7. Comparison with the current right-wing and allocation theorems

The nearby `C13` objects are not interchangeable.

### THM-2771

THM-2771's target coordinate `t` is an endpoint-dipole label in a physical
`C7 x C13` coefficient table.  Its intrinsic Bockstein divides raw
coefficients by thirteen and decodes a chart-by-target square.  Here `c` is
the THM-2640 predecessor carry and `(17)` is ordinary spatial Fourier
integration of one component indicator.  No target convolution,
uniformizer division, root-deck decoder, or common ancestry is inherited.

### THM-2772

THM-2772 needs four bare/source/target/both allocation coefficients on one
common endpoint atom before Fourier projection.  The family `(17)` supplies
one positive component state only.  It has no endpoint pair `(L,R)`,
determinant sector, allocation bits, or Segre mixed face.  Equations
`(40)--(41)` show that the common-atom prerequisite fails in the most direct
rail-eight transplant.

### THM-2763 and THM-2774

THM-2763's address `(r,k,l)` retains carrier and endpoint harmonics but not
factor allocation.  Nothing in `(17)` selects such an exact relation
address.  THM-2774 supplies the lattice quotient `P(A12)/Q(A12)`, but
`(15)` proves that the direct facet map into it is constant.  Only the
Fourier-dual character after the kernel quotient survives.

## 8. Source, map, loss, sidecar, and stopping test

The complete connection contract is:

| item | exact content |
|---|---|
| source | THM-2672 selected positive carry components `I_c` |
| geometric map | facet incidence projection `v_c -> w_c=-lambda_(m(c))` |
| geometric target | one regular minuscule `A12` simplex |
| direct quotient | constant class `-gamma` in `P/Q`; no carry identification |
| physical map | ordinary spatial Fourier coefficient `C_c(N)` |
| physical target | the `C13` character `c -> zeta^(-ac)` for `N=aS` |
| preserved | positivity, component length, carry ratios, slope eigenvalue |
| destroyed | exact `1/S` component origin, semantic section, endpoint pair, allocation |
| needed sidecar | one literal common semantic endpoint atom before integration |
| decisive test | intersection `(42)`; it is zero for the selected rail-eight source |

This is stronger than a shared-prime analogy and weaker than an endpoint
current.  It identifies exactly where the lattice/Fourier bridge works and
exactly where the physical attachment fails.

## 9. Exact companion and audit boundary

Run

```bash
python 04-computation/lrc14_a12_minuscule_carry_fourier_boundary_thm2785.py
python -O 04-computation/lrc14_a12_minuscule_carry_fourier_boundary_thm2785.py
```

Both modes byte-match the stored transcript.  The companion uses explicit
exceptions, no floating-point truth checks, and no truth-bearing Python
assertions.  It:

1. pins the audited THM-2672 and THM-2749 helper hashes;
2. rebuilds all thirteen full positive facets in configurations `(3)` and
   `(29)`;
3. verifies their roots, missing labels, selected clocks, exact intervals,
   lengths, translation sections, and positive numerators;
4. verifies all `13^2-13` pairwise minuscule root differences, the regular
   simplex Gram matrix, constant `P/Q` class, slope roots, cyclic Gram
   profile, and all twelve exact cyclotomic DFT nonvanishings;
5. verifies the iff divisibility arithmetic in `(19)`, all twelve physical
   component characters, the `c3` character, and every odometer correction
   in `(25)--(27)`;
6. rebuilds THM-2749's positive rail-eight carrier and checks the exact rail
   weight, source/target interval transport, and all five zero intersections
   in `(40)--(41)`.

The status remains candidate until the new exact companion and the final
scope wording receive an immutable independent hostile audit.  The
underlying lattice signs, Fourier phases, component translation law, rail
typing, and selected common-atom failure were independently rederived before
this packaging.

```text
PROVED IN THE CANDIDATE:
  minuscule A12 orbit and regular-simplex DFT;
  constant direct P/Q class;
  all twelve physical component-indicator carry characters;
  c3 character -2gamma;
  Fourier annihilation of the odometer kernel;
  rail-eight positive spectral control;
  selected rail-eight semantic common-atom failure.

NOT PROVED:
  direct carry=P/Q identification;
  geometric cyclic gluing of the facets;
  full component-refined nerve homology;
  THM-2771 target/Bockstein identification;
  THM-2334 relation-address or endpoint current;
  THM-2772 allocation K4 or mixed face;
  failure of every rail-eight component;
  row exclusion or LRC(14).                                           (43)
```

QED, conditional only on candidate status promotion after independent audit.
