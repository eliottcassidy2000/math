# HYP-2486 - A000568 hides in LRC as the source-deleted observer fiber

**Status:** OPEN proof program with exact source-cone lemma and exact small-clock audits.
**Source:** codex-2026-06-13.
**Extends:** HYP-1977, HYP-1982, HYP-2020, HYP-2120, HYP-2481, THM-381, THM-497.
**Artifacts:** `04-computation/lrc_a000568_source_fiber_codex.py`, `05-knowledge/results/lrc_a000568_source_fiber_codex.out`, `07-reflections/lrc-a000568-source-fiber.md`.

## Statement

Unlabelled tournament isomorphism classes, counted by A000568, do not decide
LRC directly.  They appear one layer lower:

```text
LRC witness state
  <=> observer is a source in the threshold source-lift tournament
  <=> rooted class is source_cone(deleted moving-runner class).
```

Thus the LRC target slice is a source-cone copy of A000568 on the moving
runners.  The arithmetic speed clock does not visit all of that slice; it walks
through a restricted source-deleted image.  This is why raw A000568 projections
mix safe and unsafe states, while the observer-source/deleted-source fiber is
the proof-relevant quotient.

## Exact Lemma

For every tournament class `C` on `m` vertices, add one distinguished source
vertex to obtain `source_cone(C)`.  Deleting that distinguished source recovers
`C`.  Therefore source-rooted tournament isomorphism classes on `m+1` vertices
are in canonical bijection with unlabelled tournament classes on `m` vertices.

The computation verifies the bijection by canonical enumeration through
`m=6`:

```text
m  A000568  enumerated  rooted_source_classes  collisions  deletion_failures
1        1           1                      1          0                  0
2        1           1                      1          0                  0
3        2           2                      2          0                  0
4        4           4                      4          0                  0
5       12          12                     12          0                  0
6       56          56                     56          0                  0
```

The proof is elementary, but the enumeration checks the implementation used by
the LRC audits.

## LRC Source Lift

Given a time state with stationary observer `0`, keep the usual half-turn
tournament among moving runners.  For observer-runner edges, orient

```text
0 -> i  iff  ||v_i t|| >= 1/N.
```

Then the observer is a source exactly when the LRC predicate holds at that
time.  At every good sampled state, the rooted observer class is exactly
`source_cone(deleted moving-runner class)`.

Small exact audits:

```text
N4 consecutive: states=25, good=2, runner_mixed=2, rooted_mixed=0, cone_exact=True
N5 prime-ish:   states=93, good=16, runner_mixed=2, rooted_mixed=0, cone_exact=True
N6 mixed:       states=153, good=24, runner_mixed=4, rooted_mixed=0, cone_exact=True
N7 mixed:       states=305, good=28, runner_mixed=5, rooted_mixed=0, cone_exact=True
```

Batch primitive scans also keep `rooted_mixed=0` and `cone_exact=True` through
the audited `N<=7` ranges.

## LRC14 Transfer

For `N=14`, full canonical enumeration of A000568(13) is intentionally too
large for this script, so the computation records fingerprints at first unit
witnesses instead.  Each hard row lands in the same formal source-cone slice:
observer outdegree is `13`, and deleting the observer leaves a 13-vertex
runner tournament fingerprint.

```text
AP13:                 witness 1/14, root_outdeg=13, H=2641713, c3=85
one-stranger-611:     witness 9/41, root_outdeg=13, H=3351471, c3=89
HYP2470-E1-q33:       witness 13/33, root_outdeg=13, H=1234135, c3=74
HYP2470-E2-q31:       witness 14/31, root_outdeg=13, H=2744463, c3=84
THM497-band2-refuter: witness 7/43,  root_outdeg=13, H=110233,  c3=43
```

This does not prove LRC14.  It identifies the structural coordinate that
survives the projection: a hard row must either hit the source-deleted fiber or
carry enough endpoint/band-cover debt to be forced into the Q31/7-ideal/13-clock
portal of HYP-2471/HYP-2480.

## Tournament Analysis

Primary vertices were quotient layers, not runners.  The pairwise observable
was the profile

```text
(source purity, A000568 proximity, n14 transfer, non-tautology,
 compression, low label cost).
```

Lexicographic dominance gives a transitive meta-tournament:

```text
source_deleted_fiber
> observer_source_lift
> gap_threshold_stack
> band_cover_unit_fiber
> blocking_height_dominance
> observer_marked_phase
> raw_runner_A000568.
```

This ranking is a method choice, not a theorem.  It says the best current
coordinate is the source-deleted A000568 fiber, while THM-497's band-cover
unit fiber and HYP-2481's blocking-height dominance are the LRC14 transfer
coordinates.

## Assumption Challenge

Alternate LRC tournament vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, matroid circuits, proof obligations.
```

Selected quotient:

```text
vertices = quotient layers and, inside the target, deleted moving-runner classes.
```

Preserves:

```text
observer-source predicate and deleted moving-runner unlabelled class.
```

Destroys:

```text
exact gap lengths, speed labels, owner intervals, endpoint pressure,
and multiplicative dilated-band cover loads.
```

Challenged assumption:

```text
A000568 should not be read as the raw full phase class deciding loneliness.
The useful A000568 object is the source-deleted base below a threshold-marked
observer source cone.
```

## Proof Route

The new proof target is:

```text
Any long LRC14 blocked-band walk either enters a source-cone deleted class,
or its failure to enter that slice forces balanced cover congruences, hence
the Q31/7-ideal/13-clock portal already isolated by HYP-2471/HYP-2480.
```

Concrete next steps:

1. Add a support-criticality scan to HYP-2481 rows, but record the
   source-deleted runner fingerprint at the first leak.
2. For below-eight-core rows, build a ledger whose state is
   `(deleted A000568 class, Q27 obligation class, divisor fiber, owner/Bprime debt)`.
3. Search for repeated source-deleted fingerprints among hard rows; if many
   hard rows share a fingerprint but diverge in witness shell, the missing
   coordinate is the multiplicative band-cover fiber, not the tournament base.
