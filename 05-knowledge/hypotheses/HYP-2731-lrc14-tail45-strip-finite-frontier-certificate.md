---
id: HYP-2731
title: LRC14 generated tail45 strip has a two-witness finite frontier certificate
status: OPEN; exact finite frontier certificate plus Lean cheap-side theorem
source: codex-2026-06-21
depends_on:
  - HYP-2728
  - HYP-2727
  - HYP-2722
  - HYP-2702
  - HYP-2698
related:
  - HYP-2726
  - HYP-2725
  - HYP-2721
  - THM-561
  - THM-558
  - OPEN-Q-108
---

# HYP-2731: Tail45 Strip Frontier Certificate

## Claim

The HYP-2728 generated-word separator can be sharpened to a finite frontier
certificate:

```text
182/2005 <= (q5 + 5 q6) / q0 <= 10910/21539
```

for every normalized generated atom move in the HYP-2722 sparse-tail
frontier.  The lower and upper endpoints are both attained by size-3
challengers:

```text
floor:   shape=(0,1,3),  context=[3+1],  tail45=182/2005
ceiling: shape=(0,1,13), context=[4],    tail45=10910/21539
```

The five HYP-2721 cheap abstract directions lie outside this strip:

```text
r=1:  1       above ceiling
r=2: -1       below floor
r=3:  3/2     above ceiling
r=4: -1       below floor
r=5: -1       below floor
```

Thus the generated-word exclusion target can be stated as a tail-strip theorem
before relation-code/Delsarte packet classification is invoked.

## Exact Evidence

Script:
`04-computation/lrc14_tail45_strip_certificate_codex_20260621.py`.

Stored output:
`05-knowledge/results/lrc14_tail45_strip_certificate_codex_20260621.out`.

The script recomputes the `318` HYP-2722/S71 generated frontier rows and checks
the strip with exact Fraction arithmetic:

```text
below_floor = 0
above_ceiling = 0
```

It also reports by-size and by-context extrema.  The size-3 layer owns both
global endpoints:

```text
size=3 min_tail=182/2005 at (0,1,3)/[3+1]
size=3 max_tail=10910/21539 at (0,1,13)/[4]
```

The remaining sizes are strictly internal:

```text
size=4: 21/188 <= tail45 <= 143/444
size=5: 8269/61661 <= tail45 <= 15749/37082
size=6: tail45 = 62/207
```

## Formalization State

`TournamentH7.LRCFactorialAtom` now contains the cheap-side strip constants
and cross-multiplied predicates:

```text
tailFloorNum = 182
tailFloorDen = 2005
tailCeilNum = 10910
tailCeilDen = 21539
belowTailFloor(q,scale)
aboveTailCeil(q,scale)
```

The theorem

```text
cheapScaled_outside_tailStrip_bool
```

proves by finite `native_decide` that the Boolean cross-multiplied strip
exclusion test is true for every scaled HYP-2721 cheap direction.  This does
not formalize the generated-side automaton inequality yet; it formalizes the
exact cheap-side handoff once that inequality is proved.

## Proof Target

The hand proof should now split into two finite-generated statements:

1. **Context-local tail floor:** among all generated sparse-tail rows, the
   minimum normalized `tail45` is `(0,1,3)` against context `[3+1]`.
2. **Context-local tail ceiling:** among all generated sparse-tail rows, the
   maximum normalized `tail45` is `(0,1,13)` against context `[4]`.

The context tournament in the scout is transitive under the observable
`larger floor, smaller ceiling, smaller width`, but the two endpoints live in
different context vertices.  This argues against a single monotone context
merging proof.  The likely route is a two-bound proof: a death-chain floor
lemma for `[3+1]`-type mixed context pressure, and a far-shape ceiling lemma
for the `[4]` context.

## Assumption Challenge

The vertices are not runners, arcs, or tournaments on shapes.  For this proof
obligation, the useful vertices are context partitions and sign profiles of
the normalized generated atom move.  This quotient preserves the tail-strip
predicate that excludes cheap abstract atom-cone moves.  It destroys relation
height and full row geometry, which should re-enter only after the generated
strip has done its filtering.
