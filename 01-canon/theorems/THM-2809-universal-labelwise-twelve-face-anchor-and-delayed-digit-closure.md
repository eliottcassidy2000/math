---
id: THM-2809
title: "Universal labelwise twelve-face anchor and delayed-digit closure"
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
  configuration has future digit 12, while the complete marked source has
  future digit 6; slope-seven pullback fixes that coordinate, so the two
  delayed prefixes are disjoint.  Consequently all thirteen arbitrary-
  configuration twelve-faces are empty against the full marked source.  An
  eleven-face omitting labels zero and one remains open.  No target-endpoint,
  outside-rail, row, or LRC(14) claim is made.
source: lrc-a12-carry-bridge/universal-labelwise-twelve-face-closure-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2797-fourteen-rail-source-twelve-configuration-switch-semantic-base-no-go
  - THM-2804-fourteen-rail-v4-configuration-completion-and-eleven-carry-anchor-defect
script: 04-computation/lrc14_universal_labelwise_twelve_face_closure_thm2809.py
output: 05-knowledge/results/lrc14_universal_labelwise_twelve_face_closure_thm2809.out
script_sha256: f6dd8ba6e1909cd24f0532b1822b5ac6697a18e9871703820ccce4bd43eccfa3
output_sha256: 3b2dda93c766c2186cd8b90871b70e7270dac528f8895c6a047232c6f9982408
hash_basis: LF-normalized bytes
---

# THM-2809 -- every twelve-face meets either the wrong anchor or the wrong digit

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

THM-2797 excludes all forty-two maximal fixed-configuration attempts to
attach a source-carry-twelve facet to THM-2749's marked semantic source, but
leaves two apparent freedoms:

1. use a nonmaximal configuration for a label; or
2. choose configurations independently from label to label.

The first freedom disappears at the source label itself.  This theorem scans
all `104` configurations on every one of the first fourteen rails and proves
that only one source-label packet can reach the marked deep band.  Its
present anchor has exactly the forbidden polarity.  The unique twelve-face
which omits that source label must contain label one; its sole
band-compatible configuration occupies a future half-digit disjoint from
the complete marked source.

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
zero is not a source attachment and is outside this first subargument.
Section 5 treats that unique face separately.

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
`(169,181)/182`.  No null seam is being promoted.  Both surviving
intersections equal `H_mark=(169,181)/182`.  The nearest failed half ends at
`168/182`; every other failed half lies farther from the open marked band.

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

## 4. Every family containing label zero dies

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

This includes twelve of the thirteen twelve-of-thirteen mixtures.  It also
includes smaller source-labelled families.  The four-vertex THM-2804
statement is an immediate corollary rather than a dependency.

## 5. The wrapped label-one row dies at the delayed digit

The sole twelve-face not covered by Section 4 omits label zero and therefore
contains label one.  That label has a genuine geometric distinction.
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

But this unique row has

```text
h=12,                    kappa=1,                                 (18)
```

so its ordinary delayed prefix is contained in the future half-digit

```text
25/26 < y={R x} < 1.                                             (19)
```

THM-2749's complete marked source prefix has `h=6,kappa=1`, hence is
contained in

```text
13/26 < y={R x} < 14/26.                                        (20)
```

The slope-seven pullback does not change this coordinate:

```text
R tau_1=7 in Z
implies {R(x+tau_1)}={R x}.                                      (21)
```

Therefore `(19)` and `(20)` are literal disjoint open intervals in the same
source `y`-coordinate.  Exact reconstruction of the ordinary and semantic
prefixes gives zero intersection on all seven clocks.  The label-one packet
and the **full delayed** marked source are disjoint.

This conclusion is intentionally not stated for the bare pre-prefix marked
base.  The delayed future digit is load-bearing.

## 6. Complete twelve-face classification

Let `L` be any twelve-element subset of `F_13`, and choose arbitrary
THM-2640 configurations independently for every label in `L`.

- If `0 in L`, Section 4 kills the intersection at label zero's anchor.
- If `0 notin L`, then `L=F_13\{0}` contains label one, and Section 5 kills
  the full intersection at the future half-digit.

Thus

```text
all thirteen arbitrary-configuration twelve-faces have empty intersection
with the matching full THM-2749 marked source.                    (22)
```

The exact failure census is

```text
12 faces: label-zero anchor F_(ell,7) versus its complement;
 1 face:  label-one digit (25,26)/26 versus (13,14)/26.           (23)
```

An eleven-face omitting both labels zero and one avoids both arguments.
It is the sharp residual; no positivity is asserted.

## 7. Information ledger and modular reframe

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
| repair after omitting label zero | label-one future digit remains |
| residual after both gates | an at-most-eleven face omitting labels zero and one |

THM-2804's abstract `V4` and quotient `S3` organize four notable
configurations, but the physical marked-band observable collapses the entire
104-bank before those symmetries act.  The relevant sidecar is not another
abstract permutation: it is the distinguished source label together with
its strict half and anchor polarity.

## 8. Scope

The theorem does **not** prove:

- a no-go for the target marked endpoint;
- a no-go outside the first fourteen THM-2749 rails;
- a no-go for a different source carry or marked band;
- a no-go for an eleven-face omitting labels zero and one;
- a row exclusion or LRC(14).

It proves a complete arbitrary-configuration no-go for every twelve-label
THM-2640 family on the first fourteen rails against the corresponding full
THM-2749 marked source.

## 9. Exact companion

Run

```bash
python 04-computation/lrc14_universal_labelwise_twelve_face_closure_thm2809.py
python -O 04-computation/lrc14_universal_labelwise_twelve_face_closure_thm2809.py
```

Both modes must byte-match the stored transcript.  The companion pins the
proved THM-2672 and THM-2749 scripts, rebuilds the first-fourteen-rail unit
flags, verifies the pulled-root law for all `52*13` edge/kappa/height/label
tuples, classifies all `26` strict half types, scans every one of the
`14*104` source configurations, independently records the wrapped label-one
row, reconstructs both delayed prefix banks, and proves their disjointness
on all seven clocks.  It uses exact integer arithmetic, explicit exception
gates, and no truth-bearing Python assertions or floating point.

Promotion requires an immutable independent hostile audit of the
label typing, strict circular half endpoints, full unit-bank classification,
anchor polarity, pulled-back future-coordinate identity, delayed-prefix
disjointness, twelve-face census, normal/optimized replay, dependency hashes,
and documentation gates.

QED, pending independent hostile audit.
