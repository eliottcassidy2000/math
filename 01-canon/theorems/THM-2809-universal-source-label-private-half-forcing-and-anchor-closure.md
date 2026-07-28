---
id: THM-2809
title: "Universal source-label private-half forcing and anchor closure"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  On each of the first fourteen THM-2749 rails, scan the
  full 104-member THM-2640 configuration bank.  A source-label-zero packet
  at carry 12 can meet the marked source deep band through only two half
  geometries.  The exact unit gate leaves the unique configuration
  (sector,edge,kappa,h)=(0,1,1,6), whose label-zero present anchor is
  F_(ell,7), while the marked source retains F_(ell,7)^c.  Hence any
  arbitrary labelwise configuration mixture containing the source label is
  empty against the marked source.  The unique wrapped label-one
  configuration is retained as the sharp boundary for faces omitting label
  zero.  No target-endpoint, outside-rail, row, or LRC(14) claim is made.
source: lrc-a12-carry-bridge/universal-source-label-forcing-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2797-fourteen-rail-source-twelve-configuration-switch-semantic-base-no-go
  - THM-2804-fourteen-rail-v4-configuration-completion-and-eleven-carry-anchor-defect
script: 04-computation/lrc14_universal_source_label_private_half_closure_thm2809.py
output: 05-knowledge/results/lrc14_universal_source_label_private_half_closure_thm2809.out
script_sha256: d02cec64a45cb89a9d88ccdd68e4c2d03fae8490caf0a49cabdb78da8d4927b7
output_sha256: 964fd5edabddf80ccee9cba43ef9c3efcf360d0a5fef2263075639010afa1536
hash_basis: LF-normalized bytes
---

# THM-2809 -- the marked source label has one possible half and the wrong anchor

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

THM-2797 excludes all forty-two maximal fixed-configuration attempts to
attach a source-carry-twelve facet to THM-2749's marked semantic source, but
leaves two apparent freedoms:

1. use a nonmaximal configuration for a label; or
2. choose configurations independently from label to label.

Both freedoms disappear at the source label itself.  This theorem scans all
`104` configurations on every one of the first fourteen rails and proves
that only one source-label packet can reach the marked deep band.  Its
present anchor has exactly the forbidden polarity.

The result is stronger than the four-vertex `V4` corollary suggested by
THM-2804: no information about the other labels is needed.

## 1. Source-labelled attachment model

Use THM-2672's scales

```text
p=13,                 R=13^6,                 c0=12.                (1)
```

Fix one of the first fourteen THM-2749 rails `j` and one delayed clock
`ell`.  For a target label `delta`, use the slope-seven lift

```text
tau_delta=7 delta/R,
c_delta=12+7 delta                         in F_13.                  (2)
```

Allow that label to choose an arbitrary THM-2640 configuration

```text
e_delta=(sector,edge,kappa,h)
       in {0,1}^3 x F_13,                                         (3)
```

subject only to the exact THM-2640 unit flag at carry `c_delta`.  Pull its
packet back by `tau_delta` to the common source coordinate.

A **source-labelled attachment** contains the label `delta=0` packet.  This
is load-bearing terminology: a disconnected twelve-face which omits label
zero is not a source attachment and is outside the theorem.

Every other factor in a proposed common atom is an intersection.  It
therefore suffices to compare label zero's private half and present anchor
with the matching marked source base.

## 2. Only two half geometries can reach the marked band

For `e=(sector,edge,kappa,h)`, write

```text
r_e(c)=2c+floor((2h+kappa)/13)+1_(edge=0)       in F_13.           (4)
```

Since `14=1 mod 13`, equations `(2)` and `(4)` give

```text
r_e(c_delta)=r_e(12)+delta.                                      (5)
```

The physical translation changes the deep-root index by the same `delta`.
Thus the pulled-back private half is determined by the **source root**
`r_e(12)`, independently of the label.

For a nonzero root `r`, write strict half intervals in numerator coordinates
over denominator `182` as

```text
edge 0: (14r-13,14r),          edge 1: (14r,14r+13).              (6)
```

At root zero the corresponding circular intervals are

```text
edge 0: (169,182),             edge 1: (0,13).                   (7)
```

THM-2749's marked source base retains the strict deep overlap

```text
H_mark=(169,181)/182.                                             (8)
```

An exact endpoint comparison of all `2*13` half types gives

```text
H_e intersects H_mark is nonempty
iff (edge,r_e(12)) is (0,0) or (1,12).                            (9)
```

The two possible intersections are respectively all of `H_mark` and
`(169,181)/182`.  No null seam is being promoted: every failed half ends at
`168/182` or starts at `181/182` or later.

## 3. The full 104-configuration unit bank leaves one row

For each rail, the exact companion rebuilds THM-2640's primitive-content
unit table and scans

```text
2 sectors * 2 edges * 2 kappa values * 13 heights =104           (10)
```

configurations at source carry `12`.

The edge-zero/root-zero alternative in `(9)` is nonunit.  Among all
edge-one/root-twelve alternatives, exactly one configuration is a unit:

```text
B=(sector,edge,kappa,h)=(0,1,1,6).                               (11)
```

This classification is identical on all fourteen rails.  Hence every
source-label packet that can have positive open intersection with the
marked deep band is forced to use `B`.

The label-zero present factor of `B` is

```text
F_(ell,-h)=F_(ell,7).                                             (12)
```

THM-2749's marked source construction instead retains

```text
M_(j,ell) subset F_(ell,7)^c.                                    (13)
```

Therefore

```text
P_(B,12) intersect M_(j,ell)=empty                               (14)
```

for every one of the first fourteen rails and every clock.  All translated
labels, delayed prefixes, words, owners, and component selections are later
intersections and cannot repair `(14)`.

## 4. Universal labelwise consequence

Consider any family of label packets on one of the first fourteen rails.
The configurations may vary arbitrarily with the labels; they need not be
maximal and need not belong to THM-2804's `V4`.

If the family contains label zero and its common intersection with the
matching THM-2749 marked source were nonempty, then label zero's private
half would have to meet `(8)`.  Sections 2--3 force its configuration to be
`B`, after which `(12)--(14)` make the intersection empty.  Thus

```text
every arbitrary labelwise configuration mixture containing source label 0
has empty intersection with the matching marked source base.           (15)
```

This includes every twelve-of-thirteen mixture except the face which omits
label zero.  It also includes smaller source-labelled families.  The
four-vertex THM-2804 statement is an immediate corollary rather than a
dependency.

## 5. Sharp wrapped-label boundary

The excluded face omitting label zero has a genuine geometric distinction.
For label one,

```text
c_1=6.                                                            (16)
```

The full unit/half scan finds exactly one configuration whose pulled-back
half meets `(8)`:

```text
(sector,edge,kappa,h)=(1,0,1,12),
target root=1,
pulled-back source root=0,
private half=(169,182)/182.                                      (17)
```

Thus the wrapped root-zero half is not a fictitious endpoint artifact.  It
is an honest unit after translation because the **target** root is one.
Equation `(17)` does not assert that its rail, present factors, delayed
prefix, or the other eleven labels meet.  It records the cheapest boundary
which any theorem about faces omitting label zero must retain.

In particular, one may not silently call such a face a source attachment or
reuse the label-zero anchor `(12)`.

## 6. Information ledger and modular reframe

The connection from the configuration bank to the marked source is:

| item | exact content |
|---|---|
| source | one arbitrary THM-2640 label-zero packet at carry 12 |
| target | matching THM-2749 marked source base |
| map | literal source-coordinate intersection |
| preserved | rail, clock, source carry, root half, present anchor |
| first reduction | marked deep band leaves two edge/root types |
| second reduction | unit gate leaves the unique row `B` |
| first failure | `F_(ell,7)` versus `F_(ell,7)^c` |
| destroyed by omitting label zero | source anchor and source-attachment typing |

THM-2804's abstract `V4` and quotient `S3` organize four notable
configurations, but the physical marked-band observable collapses the entire
104-bank before those symmetries act.  The relevant sidecar is not another
abstract permutation: it is the distinguished source label together with
its strict half and anchor polarity.

## 7. Scope

The theorem does **not** prove:

- a no-go for a twelve-face omitting label zero;
- a no-go for the target marked endpoint;
- a no-go outside the first fourteen THM-2749 rails;
- a no-go for a different source carry or marked band;
- positivity or emptiness of the wrapped label-one row after all factors;
- a row exclusion or LRC(14).

It proves a complete, arbitrary-configuration no-go for attaching any
source-labelled THM-2640 family on the first fourteen rails to the
corresponding THM-2749 marked source.

## 8. Exact companion

Run

```bash
python 04-computation/lrc14_universal_source_label_private_half_closure_thm2809.py
python -O 04-computation/lrc14_universal_source_label_private_half_closure_thm2809.py
```

Both modes must byte-match the stored transcript.  The companion pins the
proved THM-2672 and THM-2749 scripts, rebuilds the first-fourteen-rail unit
flags, verifies the pulled-root law for all `52*13` edge/kappa/height/label
tuples, classifies all `26` strict half types, scans every one of the
`14*104` source configurations, and independently records the wrapped
label-one boundary.  It uses exact integer arithmetic, explicit exception
gates, and no truth-bearing Python assertions or floating point.

Promotion requires an immutable independent hostile audit of the
source-label typing, strict circular half endpoints, full unit-bank
classification, anchor polarity, wrapped-label boundary, normal/optimized
replay, dependency hashes, and documentation gates.

QED, pending independent hostile audit.
