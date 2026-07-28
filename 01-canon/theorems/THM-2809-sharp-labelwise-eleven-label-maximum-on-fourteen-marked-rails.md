---
id: THM-2809
title: "Sharp labelwise eleven-label maximum on fourteen marked rails"
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
  configuration twelve-faces are empty against the full marked source.
  This bound is sharp: labels 2,...,12 form a positive full delayed
  eleven-face on every rail, and each of their 2^11 adjacent-edge assignments
  gives the same atom.  Thus maximum label cardinality is exactly 11
  (simplex dimension 10) in this marked-source model.  No target-endpoint,
  outside-rail, row, or LRC(14) claim is made.
source: lrc-a12-carry-bridge/sharp-labelwise-eleven-face-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2797-fourteen-rail-source-twelve-configuration-switch-semantic-base-no-go
  - THM-2804-fourteen-rail-v4-configuration-completion-and-eleven-carry-anchor-defect
script: 04-computation/lrc14_sharp_labelwise_eleven_label_maximum_thm2809.py
output: 05-knowledge/results/lrc14_sharp_labelwise_eleven_label_maximum_thm2809.out
script_sha256: d74ea1db38238dcd95a32598d1f728d2e37a08be7cc8c951c87bc3de15af6551
output_sha256: 0f4f9bc3747131f37c57534049c09543b5412fab2f95c654251d8debf8c9d844
hash_basis: LF-normalized bytes
---

# THM-2809 -- eleven labels survive exactly; twelve never do

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
THM-2804: no information about the other labels is needed for the upper
bound.  Better, the only residual label set is an honestly positive
eleven-face, so the upper bound is exact.

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
The next section proves that this residual is positive.

## 7. The residual eleven-face is positive on every rail

Set

```text
L_*={2,3,...,12}.                                                 (24)
```

For every `delta in L_*` and every one of the fourteen rails, the exact unit
bank contains precisely the following two marked-band-compatible rows with
the common data `sector=0,kappa=1,h=6`:

```text
edge 0: target root delta,   pulled root 0,  half (169,182)/182;
edge 1: target root delta-1, pulled root 12, half (168,181)/182.   (25)
```

Both target roots are nonzero because `2<=delta<=12`, so both rows are
honest THM-2640 units.  Both pulled halves contain all of `H_mark`.
Therefore the edge choice imposes no further cut once the marked source is
retained.

The remaining factors are common across all choices:

```text
sector=0,                  h=6,                  kappa=1.          (26)
```

In particular, THM-2749's semantic delayed prefix is its terminal-fork
refinement and hence a subset of every label's ordinary `(26)` prefix.
Also `R tau_delta=7delta` is integral, so each pullback uses the same
future coordinate and predecessor carry `12`.

For each rail and clock, form the exact physical intersection:

1. start with the full THM-2749 marked source base;
2. for every `delta=2,...,12`, intersect the rail pullback by
   `7delta/R`;
3. intersect the corresponding pullback of `F_(ell,7)`;
4. retain the common marked deep band and semantic section already present;
5. apply the common sector-zero, `h=6,kappa=1`, carry-`12` marked prefix.

The positive-clock supports and total exact delayed numerators are

```text
rail  0: (1)       399580256360672050023360
rail  1: (6)        74205644260590152069760
rail  2: (2,3)     724908063903933297548160
rail  3: (5)       565104521676801927300480
rail  4: (2,3)    1130171627188809393027840
rail  5: (4,6)     682117240653421081629120
rail  6: (2,3)    1267162127454119622485760
rail  7: (5,6)     941819893732588135224960
rail  8: (1,2,3)  1449825103908006680574720
rail  9: (5)       596479469905204957431360
rail 10: (1,3)     676409208657101856256320
rail 11: (5)       562231844838877400066880
rail 12: (2)       582228901121553500855040
rail 13: (5)       399555625773821502585600.                       (27)
```

All fourteen rows are strictly positive.  Since the edge affects only the
private half in `(25)`, every one of the

```text
2^11=2048                                                          (28)
```

labelled edge assignments realizes the same common atom on each rail.

Sections 4--6 exclude every label set of cardinality twelve, while `(24)`--
`(28)` exhibit cardinality eleven.  Hence the sharp result is

```text
maximum label cardinality=11,
maximum simplex dimension=10                                      (29)
```

for arbitrary labelwise THM-2640 configurations against the full marked
THM-2749 source on each of the first fourteen rails.

## 8. Information ledger and modular reframe

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
| residual after both gates | the positive eleven-face `L_*` |
| edge freedom on `L_*` | all `2^11` assignments give the same atom |

THM-2804's abstract `V4` and quotient `S3` organize four notable
configurations, but the physical marked-band observable collapses the entire
104-bank before those symmetries act.  The relevant sidecar is not another
abstract permutation: it is the distinguished source label together with
its strict half and anchor polarity.

## 9. Scope

The theorem does **not** prove:

- a no-go for the target marked endpoint;
- a no-go outside the first fourteen THM-2749 rails;
- a no-go for a different source carry or marked band;
- the full component/gain nerve or homology beyond the displayed atom;
- a row exclusion or LRC(14).

It proves the exact maximum label cardinality in the first-fourteen-rail,
full-marked-source attachment model: all twelve-label families fail and one
eleven-label family is positive.

## 10. Exact companion

Run

```bash
python 04-computation/lrc14_sharp_labelwise_eleven_label_maximum_thm2809.py
python -O 04-computation/lrc14_sharp_labelwise_eleven_label_maximum_thm2809.py
```

Both modes must byte-match the stored transcript.  The companion pins the
proved THM-2672 and THM-2749 scripts, rebuilds the first-fourteen-rail unit
flags, verifies the pulled-root law for all `52*13` edge/kappa/height/label
tuples, classifies all `26` strict half types, scans every one of the
`14*104` source configurations, independently records the wrapped label-one
row, reconstructs both delayed prefix banks, and proves their disjointness
on all seven clocks.  It then verifies both edge rows for all `14*11`
rail/label pairs, reconstructs the complete eleven-label physical
intersection, checks all fourteen support/mass rows in `(27)`, and proves
that the edge choices impose no extra cut.  It uses exact integer arithmetic,
explicit exception gates, and no truth-bearing Python assertions or floating
point.

Promotion requires an immutable independent hostile audit of the
label typing, strict circular half endpoints, full unit-bank classification,
anchor polarity, pulled-back future-coordinate identity, delayed-prefix
disjointness, twelve-face census, normal/optimized replay, dependency hashes,
the positive eleven-face construction, all `2^11` edge choices, sharp
cardinality/dimension wording, and documentation gates.

QED, pending independent hostile audit.
