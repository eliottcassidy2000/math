---
id: HYP-2651
title: LRC14 core-gap atlas - exact fixed-observer collar and the near-collar template target
status: OPEN; exact bounded atlas through max speed 19
source: codex-2026-06-19-S32
depends_on:
  - THM-523
  - HYP-2649
  - HYP-2648
  - HYP-2647
  - HYP-2644
related:
  - HYP-2638
  - HYP-2639
  - HYP-2646
  - THM-541
  - THM-542
  - THM-543
  - HYP-2653
  - HYP-2654
  - OPEN-Q-108
---

# HYP-2651 - LRC14 Core-Gap Atlas

## Claim

The fixed-observer core-gap crux in THM-523 should be attacked with exact
positive-core gap measure first, and only then routed through Freiman,
state-word, and far-element quotients.

For a 12-core `C` of positive speeds, define

```text
G_C = {t in [0,1): ||c t|| > 1/14 for all c in C}.
```

The exact bounded atlas through all primitive `C subset [1,19]` shows the
known drop-6 core

```text
C0 = (1,2,3,4,5,7,8,9,10,11,12,13)
```

is still the unique minimum, with

```text
meas(G_C0) = 7/858.
```

The next distinct safe measure in the atlas is

```text
426/35035
```

at the endpoint-drop core `(1,2,3,4,5,6,7,8,9,10,11,13)`, giving separation

```text
426/35035 - 7/858 = 841/210210.
```

## Evidence

Script:
`04-computation/lrc14_core_gap_atlas_codex_s32.py`

Stored output:
`05-knowledge/results/lrc14_core_gap_atlas_codex_s32.out`

The script scans every primitive positive 12-core in `[1,B]`, with exact
rational interval arithmetic at target `1/14`.  It deliberately does **not**
translate cores to contain zero: fixed-observer gap measure is scale-invariant,
but it is not freely translation-invariant.

For `B=19`, the exact scan covers `50,388` primitive positive cores.  The
cumulative minimum is:

```text
B=12: 6617/194040 at (1,2,...,12)
B=13..19: 7/858 at (1,2,3,4,5,7,8,9,10,11,12,13)
```

The single-hole AP-window ledger for `[1,13]\{e}` has a unique minimum at
`e=6`:

```text
drop 6  ->  7/858
drop 12 ->  426/35035
drop 4  ->  97/4004
drop 10 ->  1520/63063
```

The critical interval ledgers are small enough to become a hand-checkable local
lemma.  For the minimizer `[1,13]\{6}`:

```text
[29/182, 9/56]     length 1/728
[29/168, 27/154]   length 5/1848
[127/154, 139/168] length 5/1848
[47/56, 153/182]   length 1/728
```

The total is

```text
2*(1/728) + 2*(5/1848) = 7/858.
```

For the next competitor `[1,13]\{12}`:

```text
[15/182, 13/154] length 2/1001
[29/70, 41/98]   length 1/245
[57/98, 41/70]   length 1/245
[141/154,167/182] length 2/1001
```

The total is `426/35035`.  This gives the AP-window collar proof target in
concrete form: prove every other single-hole ledger contains at least this much
safe length or decomposes into the displayed larger intervals.

The near-minimum rows through `B=19` are AP-cluster templates with a small
number of holes, usually retaining the central `6` hole.  Higher sumset excess
does not by itself make rows dangerous; the near-collar rows are recognized
better by hole/state template than by scalar excess.

## Proof Reading

This suggests a three-step route for the OPEN-Q-108 uniform fattening lemma.

1. The AP-window single-hole collar lemma is now proved as THM-541:

```text
For C=[1,13]\{e}, meas(G_C) is minimized uniquely at e=6,
with value 7/858.
```

THM-541 gives the exact addressed wall certificate.  The drop-6 components are
the signed boundary gaps `R(13,2)->L(12,2)`, `R(12,2)->L(11,2)`,
`R(11,9)->L(12,10)`, and `R(12,10)->L(13,11)`, with determinant numerators
`[3,5,5,3]`.  This confirms that the collar is an addressed boundary-owner
phenomenon, not a scalar hole-position statistic.

2. Prove a near-collar state-template lemma, now refined by HYP-2654:

```text
If meas(G_C) < 426/35035, then C has the drop-6 mouth-retention template.
```

This would convert the observed separation into a structural dichotomy:
drop-6-like rows are handled by the explicit collar, and every other bounded
template has a uniform margin above it.

Important correction from HYP-2654: this template cannot mean the exact one-hole
AP row.  The AP-tail row `(1,2,3,4,5,7,8,9,11,12,13,20)` already has
`3859/420420 = 7/858 + 1/980 < 426/35035`.  It keeps the old drop-6 mouth
components undamaged and only adds new mouth mass.  The target is therefore
boundary-owner mouth retention, not exact-row rigidity.

THM-542 closes the first infinite subcase of this corrected target: for every
one-tail AP-drop-6 row `({1,...,13}\{6,h}) union {r}`, the only row below
`426/35035` is the `(h,r)=(10,20)` tail, and it keeps old mouth mass exactly
`7/858`.

THM-543 upgrades this from the drop-6 one-tail subcase to the entire
one-replacement AP-tail layer:

```text
({1,...,13}\{a,b}) union {r}, 1 <= a < b <= 13, r >= 14.
```

Across all `78` two-hole bases, exact periodic-comb cutoffs reduce the infinite
tail to `3277` finite rows, and the only row below `426/35035` is again
`(a,b,r)=(6,10,20)` with

```text
3859/420420 = 7/858 + 1/980.
```

Thus the live near-collar gap has moved past one-replacement AP tails.  Any
remaining below-second bounded AP-tail row must use multi-hole/multi-tail
state-word damage, or else the proof must route to HYP-2653/HYP-2644 style
far-element discrepancy.

3. Compose with the far-element plateau recursion:

Rows with genuinely large or dissociated elements should be routed through
HYP-2644's plateau estimate.  The exact atlas shows that merely appending a
tail through `19` does not beat the collar; the top tail rows still look like
AP clusters with a few holes, not independent strangers.

The concurrent KPS HYP-2653 strengthens this branch by proving a
sigma-dependent far-element decorrelation bound and reducing the uniform
far-element constant to a structured breakpoint discrepancy.  Thus the current
split is: HYP-2654 for bounded AP-tail mouth retention, HYP-2653/HYP-2644 for
the far-element plateau route.

The later exact HYP-2653 engine found a bounded resonant core where the proposed
small uniform constant is too optimistic.  This does not affect THM-543; it
clarifies the architecture.  Bounded resonant AP-tail layers should be certified
by exact rational cutoffs, while genuinely wide rows need a separate structured
breakpoint-discrepancy argument.

## Relationship To Freiman Small-Excess

HYP-2638 is valid for the `L_y`/sector-functional pocket where translation and
dilation normal forms are appropriate.  The THM-523 `G_C` crux is stricter:
translation changes a fixed observer's gap set.  Therefore HYP-2651 keeps the
positive core itself as the first vertex, then records sumset excess only as a
classifier.

In the `B=19` atlas the minimum has exact excess `2`, while the second row has
exact excess `1`; later near-minimum rows have excess up to `11`.  Scalar
excess is therefore a routing hint, not the collar order parameter.

## Tournament Analysis

Vertices are proof quotients:

```text
exact_positive_core_gap
AP_drop_profile
state_word_template
sumset_excess_classifier
far_element_plateau
raw_runner_bound
```

Pairwise observable: which quotient preserves the exact OPEN-Q-108 predicate
before routing to a coarser proof object.

Switch/gauge: keep fixed-observer positivity before using Freiman, state-word,
or far-element simplifications.

Hamiltonian path:

```text
exact_positive_core_gap
> AP_drop_profile
> state_word_template
> sumset_excess_classifier
> far_element_plateau
> raw_runner_bound
```

Fingerprint: transitive priority tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1}`, no directed 3-cycles.

## Assumption Challenge

Alternate vertex sets considered: runners, holes, safe components, exact
sumset excess, state words, far-element plateaus, residue/coset classes, and
proof obligations.

The challenged assumption is that the exact core-gap crux can inherit
translation-normalized Freiman tables directly.  It cannot without extra
argument, because the observer is fixed.  HYP-2651 therefore treats exact
positive-core gap measure as the primary object and uses the other structures
only after the predicate is preserved.
