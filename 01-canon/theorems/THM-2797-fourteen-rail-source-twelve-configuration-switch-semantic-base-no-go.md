---
id: THM-2797
title: "Fourteen-rail source-twelve configuration-switch semantic-base no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the first fourteen THM-2749 semantic rails there are exactly 55 maximal
  THM-2672 configurations, of which exactly 42 admit source carry 12.  Every
  one of those 42 source-twelve facets is positive and has at least one
  positive clock in common with its matching fully marked semantic rail.
  Nevertheless all 69 clock comparisons with nonempty marked bases become
  disjoint at the first anchor-present insertion; no later translated
  present, root-half, or delayed component can repair them.  Exactly fourteen
  middle configurations have closure contacts, totalling 7018, all oriented
  facet-right to marked-left; the other 28 have none.  This is a
  configurationwise no-go on the fourteen matching rails, not a
  labelwise-configuration union, row exclusion, or LRC(14).
source: lrc-a12-carry-bridge/fourteen-rail-semantic-switch-scan-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2785-a12-minuscule-carry-fourier-character-and-common-atom-boundary
  - THM-2794-rail-eight-fixed-configuration-all-component-semantic-atom-no-go
script: 04-computation/lrc14_fourteen_rail_source12_semantic_switch_no_go_thm2797.py
output: 05-knowledge/results/lrc14_fourteen_rail_source12_semantic_switch_no_go_thm2797.out
script_sha256: 389202f67d38cd7dbd81ce3d120d74039471251ca43b474a27e2bba989e91272
output_sha256: 5a25e91af51a6ec06894b42bb411bc9e447fec51cee4a2398be2e02c5afa2f61
hash_basis: LF-normalized bytes
---

# THM-2797 -- all fourteen fixed semantic-rail switches miss at the anchor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2794 shows that all `879,237` delayed components in one rail-eight
configuration miss THM-2749's marked semantic source.  Its independent audit
found that the real obstruction occurs earlier: the facet and marked bank
have opposite anchor-present polarity.

This theorem exhausts that mechanism over the entire fourteen-rail bank on
which THM-2749 constructs its uniform marked source/target profile.  Every
maximal fixed THM-2672 configuration that can use source carry `12` is
positive, but every one misses its matching semantic base before either the
private half or delayed prefix can matter.

“Configuration switch” here means choosing one of the admissible fixed
configurations.  It does not allow different target labels to choose
different configurations inside one intersection.

## 1. Exact universe and source-carry rebase

Use the THM-2672 scales

```text
p=13,                   R=13^6,                   T=297836897838480.  (1)
```

Restrict its rail bank to the first fourteen rails

```text
j=0,1,...,13,                                                   (2)
```

which are exactly the rails used in THM-2749's uniform marked profile.
For a configuration

```text
(j,sector,edge,kappa,h),                                        (3)
```

let `U` be its set of predecessor carries on which the normalized seven-clock
row is a THM-2640 unit.  A configuration is maximal when `|U|=12`.

There are exactly

```text
55 maximal configurations in the first fourteen rails,
42 of them with 12 in U.                                        (4)
```

For each of the latter, fix source carry `c0=12` and rebase the target labels
by

```text
delta(c)=2(c-12) mod 13,                c in U.                  (5)
```

This is a twelve-element set containing zero.  Intersect the twelve rail
translates by `7 delta/R`, the twelve clock-matched present translates, and
the exact private half determined by `(edge,kappa,h,c0)`.  Before inserting
the delayed prefix, call the resulting strict base `B_e,ell`.

Every rail has exactly the same three admissible signatures

```text
(sector,edge,kappa,h,missing carry)
   =(0,0,0,6,6),
     (0,1,1,6,6),
     (1,1,0,0,0).                                               (6)
```

Thus the 42 configurations consist of three per rail, with missing-carry and
height censuses

```text
missing: 0 -> 14, 6 -> 28;
height:  0 -> 14, 6 -> 28.                                     (7)
```

## 2. Matching semantic base and positive controls

For each rail `j`, let `M_j,ell` be the matching THM-2749 fully marked source
base before its delayed prefix.  It retains the equal-weight source/target
rail overlap, both relative-present complements, the exact deep overlap, and
the common semantic source/target section.

The comparison is literal:

```text
B_e,ell and M_j,ell live in the same source coordinate,
use the same rail j and clock ell,
and B_e,ell uses source carry 12.                                (8)
```

All `42` rebased facets have positive delayed mass.  Moreover every one has
at least one clock on which both its delayed facet mass and THM-2749's marked
rail coefficient are positive.  The common-positive-clock support histogram
is

```text
{1}:3,       {1,2,3}:3,   {1,3}:3,
{2}:3,       {2,3}:9,     {4,6}:3,
{5}:11,      {5,6}:4,     {6}:3.                                (9)
```

Their total exact source-twelve facet numerator is

```text
290130506081396373729700182.                                  (10)
```

Thus the no-go is not caused by an empty facet, an empty semantic rail, or
disjoint positive clock supports.

## 3. Universal first-failure law

Let `R_j` be the twelvefold common rail before present factors.  Insert only
the label-zero present anchor:

```text
A_e,ell = R_j intersect F_(ell,-h).                             (11)
```

For every clock on which the marked base is nonempty, exact interval
intersection gives

```text
R_j intersect M_j,ell is nonempty,
(R_j intersect M_j,ell) intersect F_(ell,-h) is empty.          (12)
```

There are exactly `69` such comparisons.  The remaining `225` of the
`42*7=294` clock comparisons have empty marked base.  Hence the exact
first-failure census is

```text
anchor-present: 69,                  marked-empty: 225.          (13)
```

For `h=6`, equation `(12)` contains THM-2794's literal polarity

```text
F_(ell,7) versus F_(ell,7)^c.                                  (14)
```

For `h=0`, the complement statement becomes true only after restriction to
the twelvefold common rail `R_j`; the theorem does not claim that
`M_j,ell` is globally contained in `F_(ell,0)^c`.

All remaining facet operations are intersections:

```text
B_e,ell subset A_e,ell.                                        (15)
```

Therefore `(12)` proves

```text
B_e,ell intersect M_j,ell=empty
for every one of the 42 configurations and every clock.        (16)
```

Every delayed component is a subset of `B_e,ell`, so no choice of delayed
component inside any of these configurations can repair the attachment.

## 4. Exact seam boundary

The companion independently sweeps the full private-half bases against the
marked bases.  It finds

```text
strict positive overlap clocks =0,
facet-right to marked-left closure seams =7018,
marked-right to facet-left closure seams =0.                    (17)
```

The first and third signatures in `(6)` have no seams on any rail.  Only

```text
(sector,edge,kappa,h,missing)=(0,1,1,6,6)                       (18)
```

has closure contacts.  Its fourteen railwise seam counts are

```text
278,51,504,394,790,477,886,659,1011,416,473,393,406,280,        (19)
```

whose sum is `7018`.  The rail-eight value `1011` is exactly the base-seam
count in THM-2794.

The aggregate strict-base piece counts are

```text
facet pieces=45449,                    marked pieces=30534.      (20)
```

Because all intervals are strict open, a closure seam carries no positive
mass.  Still, `(17)--(19)` identify a sharp boundary: only the middle
signature can be repaired infinitesimally, and only by changing or removing
the anchor polarity.  The other two signatures require a nonlocal change.

## 5. Consequence and stopping rule

Within the first fourteen THM-2749 rails, replacing the rail-eight
configuration in THM-2785/2794 by any other maximal fixed THM-2672
configuration admissible at source carry `12` cannot produce a common
semantic atom.  The obstruction is not the delayed word, private root,
component section, or lack of simultaneous positive clocks.  It is the first
present anchor after the twelvefold rail restriction.

The source-to-target contract is:

| item | exact content |
|---|---|
| source | one of 42 positive source-12 maximal fixed facets |
| target | matching THM-2749 marked source base on the same rail |
| map | literal clock-matched intersection |
| preserved | physical coordinate, rail, clock, source carry, twelve-label rebase |
| first failure | label-zero anchor present factor |
| sharp boundary | 14 middle configs, 7018 one-sided seams |
| needed sidecar | remove/change anchor polarity or mix configurations labelwise |

This turns the next search away from component enumeration.  The cheapest
unsettled operation is a lawful **labelwise** configuration mixture or a
factor-repair identity that replaces the anchor-present factor.  Merely
changing the fixed configuration is exhausted.

## 6. Scope

The theorem does **not** prove:

- a no-go for rails beyond the fourteen THM-2749 rails;
- a no-go for source carry other than `12`;
- a no-go for nonmaximal THM-2672 configurations;
- a no-go when different target labels choose different configurations;
- that a union of configuration bases cannot intersect the semantic base;
- that the target rather than source marked bank behaves identically;
- a row exclusion or LRC(14).

It is compatible with THM-2687: that theorem excludes positive
thirteen-label service even after arbitrary labelwise configuration
switching in a relaxed envelope, whereas the present theorem studies
twelve-label source-12 facets and a much narrower semantic attachment.

## 7. Exact companion

Run

```bash
python 04-computation/lrc14_fourteen_rail_source12_semantic_switch_no_go_thm2797.py
python -O 04-computation/lrc14_fourteen_rail_source12_semantic_switch_no_go_thm2797.py
```

Both modes byte-match the stored transcript.  The companion uses exact
integer interval arithmetic and explicit exception gates, with no
truth-bearing Python assertions or floating point.  It pins the proved
THM-2672 and THM-2749 scripts; rebuilds the fourteen-rail unit flags; checks
all 55 maximal and 42 source-admissible configurations; independently
recomputes every facet and marked vector; records the first failed layer;
and verifies every open intersection and oriented closure seam.

An immutable independent hostile audit rebuilt the source-carry rebase and
all three signatures.  It reproduced all `42` positive facets, their common
positive clocks, the `69/225` first-failure census, zero strict overlap,
all `7018` one-sided seams, normal/optimized/stored replay, dependency and
artifact hashes, and the documentation/uniqueness gates.  No repair was
required.

QED.
