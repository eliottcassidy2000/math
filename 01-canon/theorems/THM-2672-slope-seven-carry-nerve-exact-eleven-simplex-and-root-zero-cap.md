---
id: THM-2672
title: "Fixed-configuration slope-seven carry facets, false coarse S11, and odometer twist"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; final promotion awaits the independent
  audit of the carry-rebase referee.  In every fixed THM-2640
  (rail,sector,edge,kappa,h) configuration, the slope-seven private-root law
  permits at most twelve unit carry charts.  Exactly 534 configurations attain
  twelve, and every one has a positive honest twelve-fold physical overlap.
  One canonical configuration has a positive twelve-face in every missing
  label after rebasing the source carry.  Forgetting that carry produces the
  coarse nerve boundary of a 12-simplex, with virtual reduced H_11 of rank one;
  retaining the predecessor carry instead gives a carry-refined label nerve
  of thirteen disjoint filled 11-simplices and H_11=0.  Full connected-arc
  gain refinement remains uncomputed.  Their cyclic rebase accumulates the
  nonzero THM-2657 odometer
  kernel translation 7/13^5.  Configuration switching can evade the fixed-edge
  root-zero cap in general.  The displayed component lies in a one-edge strip
  and cannot itself be extended by the missing label, so it is a maximal local
  simplex after the THM-2658 component refinement.  This does not make its
  twelve labels a maximal simplex in the unrestricted union nerve: another
  common component of those labels could still meet the thirteenth union
  label.  The existence of a thirteen-fold component elsewhere remains open.
  No row exclusion or LRC(14) conclusion follows.
source: root/physical-cech-higher-overlap-2026-07-28
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
related:
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2682-central-arrival-clock-trap-and-three-event-dilation-nilpotence
script: 04-computation/lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py
output: 05-knowledge/results/lrc14_slope7_fixed_configuration_carry_nerve_thm2672.out
script_sha256: 83ccf3a38660a92cc990bdf304fd4ea4475339731c3e7e92ad35383ef097f361
output_sha256: c2fee76c461263df43a7f12e972e6d12e8a6b15bc753694e2f04882a99247f26
secondary_scripts:
  - 04-computation/lrc14_slope7_twelve_chart_component_witness_thm2672.py
  - 04-computation/lrc14_slope7_rebase_facet_torsor_thm2672.py
secondary_outputs:
  - 05-knowledge/results/lrc14_slope7_twelve_chart_component_witness_thm2672.out
  - 05-knowledge/results/lrc14_slope7_rebase_facet_torsor_thm2672.out
secondary_script_sha256s:
  04-computation/lrc14_slope7_twelve_chart_component_witness_thm2672.py: f69d98386edae69234100130ad16c52bba862ab8072d4445414b1d9686a079a4
  04-computation/lrc14_slope7_rebase_facet_torsor_thm2672.py: 722c86b1e74df79cbc5d78da47be7c03317437c25d2783748e5789a9143c1347
secondary_output_sha256s:
  05-knowledge/results/lrc14_slope7_twelve_chart_component_witness_thm2672.out: d01a819562453e8ae77c9e6b9281c2e88a478c6650693f8f0a567f3929b43f5a
  05-knowledge/results/lrc14_slope7_rebase_facet_torsor_thm2672.out: 6b7e3687ed1e52f3b5d1e088088d281071636dd4fc9b85b1a75ca6fc0b7c6675
hash_basis: LF-normalized bytes
---

# THM-2672 -- twelve charts glue, while their coarse sphere does not

**PROVED CANDIDATE + VERIFIED-EXACT.**

THM-2640 makes every nonzero predecessor-carry root private, but its
slope-seven physical lift is not a group action.  The first genuine
higher-overlap calculation now gives both sides of that boundary:

```text
one fixed configuration has honest positive 12-fold components;
all thirteen missing-label facets occur after changing the base carry;
the thirteen facets do not live in one carry-refined label nerve.        (1)
```

Thus coarse label forgetting manufactures an `S^11`, while retaining the
physical carry dissolves it into thirteen contractible pieces.

## 1. Fixed-configuration slope-seven charts

Use THM-2640's notation

```text
p=13,                 R=13^6,                 S=13^5,
c3=2S.                                                               (2)
```

Fix one THM-2640 configuration

```text
e=(rail j, delayed sector sigma, half-edge epsilon, kappa, h),         (3)
```

including its base cell.  Here `epsilon=1` is the left deep half and
`epsilon=0` the right half.  Put

```text
b=floor((2h+kappa)/13),
r_e(c)=2c+b+epsilon                         in F_13.                    (4)
```

For a source carry `c0` and target label `delta`, take the canonical
integer-linear lift

```text
tau_delta=7delta/R,
c_delta=c0+7delta                           in F_13.                    (5)
```

Let `P_(e,c)` denote the full unintegrated THM-2640 packet in configuration
`e` and predecessor-carry stratum `c`, with its rail, present factors,
delayed word, future half-digit, private deep half-tooth, and global-content
unit test retained.  Its seven delayed-clock slices remain physical; they
are not replaced by a coefficient-only row.  When the row at `c_delta` is a
unit, define the pulled-back chart

```text
U^(c0)_delta={x:x+tau_delta belongs to P_(e,c_delta)}.                 (6)
```

The exact computations below intersect the displayed physical factors in
`(6)` before integrating.  In particular, they do not infer a common chart
from pairwise coefficient support.

## 2. The fixed-configuration cap is exactly twelve

Because `R tau_delta=7delta` is integral, the future coordinate `{Rx}`, and
hence `(h,kappa)`, is unchanged.  The source carry changes as in `(5)`, and
the private root obeys

```text
r_e(c_delta)=r_e(c0)+14delta=r_e(c0)+delta       in F_13.              (7)
```

Thus the thirteen fixed-edge roots traverse `F_13` exactly once.  The
THM-2640 chart bank retains only nonzero private roots, so one and only one
label in a fixed configuration is forbidden.  Consequently

```text
at most twelve labels U^(c0)_delta can be unit.                         (8)
```

This is an orbitwise statement.  It does not yet say that the surviving
twelve charts meet.

The exact full-carrier census rebuilds THM-2640's global primitive content
`26` and scans every canonical configuration `(3)`.  Exactly

```text
534                                                                    (9)
```

configurations have twelve unit carries.  Their unique missing-carry
histogram is

```text
missing 0: 133,              missing 6: 269,
missing 12:132.                                                       (10)
```

For every one of the `534` configurations, choose the least unit carry as
`c0`, pull all twelve charts back by `(5)`, and intersect:

```text
* every translated rail support;
* every translated present F_(ell5,-h) factor;
* one common delayed clock ell5;
* the aligned private deep half-tooth; and
* the predecessor-carry/future-half-digit delayed word.                 (11)
```

All `534/534` intersections have positive exact numerator.  Across this
census the raw overlap numerators have

```text
minimum  =442778772135550408718160,
maximum  =15728067643498917694609986,
distinct =168,                    gcd=9646.                            (12)
```

Equations `(8)` and `(9)`--`(12)` prove that the maximal
same-configuration physical nerve dimension is exactly eleven.

There is also an exact coefficient-level repair prefilter.  At the unique
root-zero carry of each of the `534` maximal configurations, keep the base
cell, `h`, and `kappa`, switch the deep half-edge, and scan every rail and
both sectors of that same cell for a nonzero THM-2640 unit.  The result is

```text
opposite-edge candidates =0/534.                            (12a)
```

Thus no same-cell edge/rail/sector switch repairs any maximal fixed
configuration.  Extending any one of these `534` fixed-configuration facets
therefore requires a base-cell change (hence a source/owner-chart change),
not merely a local half-edge switch.  This does not exclude an unrelated
same-cell thirteen-fold intersection assembled labelwise from several
nonmaximal configurations, with no fixed configuration supplying twelve of
its labels.

## 3. One explicit open twelve-chart component

The first canonical witness is

```text
base cell              (1,0),
rail                    0,
(sector,edge,kappa,h)   (0,left,0,6),
source carry            c0=0,
unit carries            F_13\{6},
target labels           F_13\{12}.                                  (13)
```

Both delayed clocks `ell5=1,3` contribute positively.  At `ell5=1` an exact
open component is

```text
( 4855397/10396204,
  23436073938185/50180491033036 ),                                    (14)
```

with

```text
length   =3/12545122758259,
midpoint =3348010562597/7168641576148.                                (15)
```

An independent rational replay at the midpoint checks the carry, future
digit, rail, present factor, delayed word, private half-tooth, and unit flag
for every label.  Labels `0,...,11` survive with roots `1,...,12`; label
`12` has carry `6` and root zero and fails.  Thus `(14)` is an actual open
common component, not only a positive integral.

This component is nonextendable inside the **full** union of configurations.
In
the pulled-back base deep coordinate

```text
q=182{c3 x},
```

its exact strip is

```text
13<q<13+12/371293.                                          (15a)
```

At the missing label `delta=12`, the deep translation adds `168`, so

```text
181<q_missing<181+12/371293.                                (15b)
```

On `(15b)`, the left edge has root zero and is never a THM-2640 unit.  The
only alternative edge is the right half of root `12`, whose half-open support
ends at `181`; it is absent throughout the open strip.  Carry, `h`, and
`kappa` are pointwise fixed, while the root law is independent of rail and
sector.  Therefore no rail, sector, edge, delayed clock, or other
configuration choice can extend `(14)` by the thirteenth label.  Thus `(14)`
is a maximal local eleven-simplex after splitting charts into the connected
physical components required by THM-2658.

This is deliberately not called a maximal face of the unrestricted
union-labelled nerve.  The same twelve labels may possess another disconnected
common component which does meet the thirteenth union label.  The present
strip calculation excludes extension of `(14)`, not extension of every common
component carrying that label set.

By THM-2658, each connected component in `(14)` determines a complete
triangle-balanced integer-gain section, and its mass is the least selected
pair-component overlap.  The present exact construction supplies the
component directly rather than reconstructing it from the pair graph.

## 4. Every missing label occurs after rebasing the carry

Keep the single physical configuration `(13)`, whose unique nonunit carry is
`c_*=6`, but let the source carry range over all `c0 in F_13`.  The omitted
label is

```text
m(c0)=2(c_*-c0)                         in F_13.                       (16)
```

For `c0!=c_*`, the active labels are the twelve deltas which transport `c0`
to a unit carry.  For `c0=c_*`, omit `delta=0`; all twelve nonzero deltas
land on unit carries, so the pulled-back base half has root zero but no
retained chart uses that root.  Exact interval intersection gives a positive
twelve-fold component in all thirteen cases.  The source/missing bijection is

```text
(0,12),(1,10),(2,8),(3,6),(4,4),(5,2),(6,0),
(7,11),(8,9),(9,7),(10,5),(11,3),(12,1).                              (17)
```

The lexicographically first open component in every row of `(17)` has the
same exact length

```text
3/12545122758259.                                                       (18)
```

This is a label-complete **rebase facet bank**.  It is not yet the boundary
of one carry-retained physical simplex.

## 5. Coarse forgetting creates a false sphere

The typing difference is visible before topology.  For fixed `c0`, every set
`U^(c0)_delta` in `(6)` lies in predecessor-carry stratum `c0`.  Different
source carries are disjoint away from the null half-open boundaries.  Hence
the **carry-refined positive-overlap label nerve** furnished by `(17)` is

```text
disjoint_union_(c0 in F_13) Delta^11.                                 (19)
```

Its `f`-vector is

```text
(156,858,2860,6435,10296,12012,
 10296,6435,2860,858,156,13),                                         (20)
```

it has thirteen connected components, reduced `H_0` of rank twelve, and

```text
H_11=0.                                                               (21)
```

Now deliberately forget the source carry and form the disconnected unions

```text
V_delta=union_(c0:m(c0)!=delta) U^(c0)_delta.                          (22)
```

Every proper label subset is contained in the positive facet with a missing
label outside that subset, while the full thirteen-label intersection is
empty in every carry stratum.  Therefore the coarse label nerve is exactly

```text
N((V_delta)_delta)=boundary Delta^12.                                  (23)
```

Its `f`-vector is

```text
(13,78,286,715,1287,1716,1716,1287,715,286,78,13),                    (24)
```

and its reduced `H_11` has rank one.  This is genuine homology of the coarse
nerve, but it is **not** a component-refined Cech equivalence.  Equation
`(22)` merged thirteen disjoint physical choices under each label.  THM-2658
requires a further split into proper connected arc components and full
integer gains.  The explicit components in Sections 3--4 supply at least one
balanced `Delta^11` inside every carry stratum, but the complete connected-arc
gain nerve has not been enumerated.  Already the coarser act of retaining
`c0` replaces `(23)` by `(19)` and kills this particular virtual sphere.

Thus `(23)` is an exact physical realization of the disconnected-label
hostile in THM-2658: high label-nerve homology can be manufactured solely by
forgetting even the predecessor-carry component.  No statement about all
homology of the fully component-refined gain nerve is made.

## 6. The facet atlas carries the odometer twist

One slope-seven step changes the source carry and omitted label by

```text
c0 -> c0+7,                       m -> m-1.                            (25)
```

The canonical lift of that step is `tau_1=7/R`.  Thirteen successive lifts
accumulate

```text
13 tau_1=91/R=7/13^5 !=0                         mod 1.                (26)
```

This is exactly THM-2657's extension-kernel translation and cocycle class
`7`.  The `f`-vectors in `(19)`--`(24)` follow independently from the
disjoint carry strata.  Equation `(26)` instead obstructs a tempting
**equivariant cyclic gluing** of those facets: repeated use of the fixed
generator lift rotates through all missing labels but ends at a nontrivial
translate, not a return to the starting component.  The speed-one guard has
trivial rotation stabilizer, so that residual translation cannot be erased
as a symmetry of the full packet.

Equation `(26)` is an obstruction to the naive cyclic gluing.  It does not
itself provide transition maps between the components in `(19)`; no such
maps are claimed.

## 7. Why the unrestricted thirteen-chart problem remains open

The cap `(8)` fixes the half-edge and the other configuration data.  It
cannot be applied after allowing a different branch at each target label.
Indeed the right half of root `r` and the left half of root `r+1` are

```text
[14r,14r+13)/182,             [14r+1,14r+14)/182,                      (27)
```

and overlap in positive length `12/182`.  A point at the fixed-edge zero
root may therefore geometrically switch edge/root at that label, so the
fixed-edge proof alone cannot give a global cap.  For each of the `534`
maximal fixed-configuration facets, the exact prefilter `(12a)` shows that no
such repair is a unit while its base cell is retained.  The
full union-labelled THM-2640 chart family may still have a positive same-cell
or inter-cell thirteen-fold component assembled from unrelated configurations;
the present theorem neither proves nor excludes it.  Equations
`(15a)`--`(15b)` prove that one explicit twelve-fold component cannot be so
extended even by changing the base cell.  They do not prove maximality of its
label set in the ordinary disconnected-union nerve.

This boundary is logically independent of THM-2680/2682.  Those theorems
show that the correctly oriented dilation handoff `D(x)={13x}` has positive
two-event fibre products but is support-nilpotent at three events on the
central-arrival carrier.  The slope-seven rebase atlas is a different
handoff geometry, not a refinement that repairs that nilpotent chain.

No lawful THM-2365 target action, common semantic owner, endpoint current,
configuration-switching thirteen-simplex, physical `S^11` transition,
holonomy trivialization, row exclusion, or LRC(14) conclusion follows.

## 8. Exact companions and audit boundary

Run

```bash
python 04-computation/lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py
python -O 04-computation/lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py
python 04-computation/lrc14_slope7_twelve_chart_component_witness_thm2672.py
python -O 04-computation/lrc14_slope7_twelve_chart_component_witness_thm2672.py
python 04-computation/lrc14_slope7_rebase_facet_torsor_thm2672.py
python -O 04-computation/lrc14_slope7_rebase_facet_torsor_thm2672.py
```

Each normal/optimized pair must byte-match its stored transcript in
`05-knowledge/results`.  The first referee rebuilds the full THM-2640 carrier
and proves the `534/534` census.  The second isolates `(13)`--`(15)` and
replays every physical factor at the rational midpoint.  The third proves
all thirteen rebased facets, `(17)`--`(26)`, and the exact coarse and
carry-refined label-nerve `f`-vectors, together with the one-edge local
nonextension strip `(15a)`--`(15b)`.

An independent hostile audit has accepted `(4)`--`(15)` with the mandatory
fixed-configuration scope and explicitly rejected the invalid full-union
root-zero inference.  Final promotion awaits the independent immutable audit
of the rebase/topology companion and all declared hashes.

QED (candidate pending final audit).
