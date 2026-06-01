---
id: THM-382
name: lrc-bounded-threshold-fiber-separation
status: PROVED (exact finite computation, bounded scope)
date: 2026-06-01
session: codex-2026-06-01-S513
depends_on:
  - THM-372
  - THM-381
---

# THM-382: Threshold-decorated fibers separate the bounded n=3 and n=4 LRC audits

## Statement

In the exact bounded audit

```text
04-computation/lrc_iso_class_constraint_s512.py
```

the sampled state space consists of all half-turn walls, all LRC endpoint
walls, and all midpoints between consecutive walls for the following primitive
speed systems:

```text
total n=3, moving speeds <= 16: 79 systems
total n=4, moving speeds <= 10: 109 systems
```

For each sampled state, the script computes LRC-good status and seven
tournament-class lifts:

```text
phase_half
phase_marked_observer
phase_safe_colored
gap_rank_marked_adjacent
gap_threshold_fiber
pair_deficit_origin_marked
pair_deficit_threshold_fiber
```

In this finite state space, the raw and rank-only lifts are mixed and certify no
sampled speed system, while the threshold-decorated gap and pair-deficit fibers
have no mixed classes and certify every sampled speed system.

For total `n=3`, the exact class/fiber table is:

```text
lift                         classes  good-only  bad-only  mixed  certifies
gap_rank_marked_adjacent           3          0         1      2        0/79
gap_threshold_fiber                8          2         6      0       79/79
pair_deficit_origin_marked         3          0         1      2        0/79
pair_deficit_threshold_fiber       9          2         7      0       79/79
phase_half                         2          0         0      2        0/79
phase_marked_observer              4          0         1      3        0/79
phase_safe_colored                14          3        11      0       79/79
```

For total `n=4`, the exact table is:

```text
lift                         classes  good-only  bad-only  mixed  certifies
gap_rank_marked_adjacent           6          0         2      4        0/109
gap_threshold_fiber               20          5        15      0      109/109
pair_deficit_origin_marked        17          0        13      4        0/109
pair_deficit_threshold_fiber      50          4        46      0      109/109
phase_half                         4          0         0      4        0/109
phase_marked_observer             11          0         1     10        0/109
phase_safe_colored                52         10        42      0      109/109
```

Here "certifies" means that the speed system visits at least one good-only class
for that lift.

## Scope

This is a finite theorem about the bounded S512 audit.  It does not claim that
the displayed bounds include all primitive speed systems at total `n=3` or
`n=4`, and it does not claim that the same class counts persist at larger
`n`.  Its content is the exact separation observed in the finite state space:
raw A000568-style phase classes and rank-only fibers are too coarse there,
while threshold-decorated fibers separate good from bad states there.

## Proof

The script enumerates the stated primitive speed systems and constructs the
exact rational event set for each one.  The event set includes:

1. half-turn tournament walls for every speed pair;
2. LRC endpoint walls `||v_i t|| = 1/n`;
3. midpoint samples between consecutive walls.

For each sampled time, all coordinates are rational `Fraction` values.  The
script computes LRC-good status by the condition

```text
forall i, ||v_i t|| >= 1/n.
```

It then canonicalizes each lift under color-preserving relabeling.  For every
lift and every class key, it counts whether the class occurs at good states,
bad states, or both.  Finally, for each speed system, it checks whether the
classes visited by that speed system intersect the good-only set.

The stored output reports exactly the two tables in the statement.  Therefore,
within the finite audit scope, the threshold-decorated fibers have zero mixed
classes and certify all sampled speed systems, while the raw/rank-only fibers
listed above have mixed classes and certify none.

## Verification Record

Stored output:

```text
05-knowledge/results/lrc_iso_class_constraint_s512.out
```

The relevant lines are the "TOTAL n=3" and "TOTAL n=4" lift-separation tables.
The audit also records the sampled state counts:

```text
total n=3: 10995 exact sampled states
total n=4: 13197 exact sampled states
```

## Significance

This theorem formalizes the finite positive/negative evidence behind
HYP-1982.  The A000568 base is real, but raw tournament class, observer mark
alone, and rank-only gap/pair data are too coarse in the bounded audits.  The
first separating objects are threshold-decorated fibers: class data plus the
actual `1/n` threshold colors.

For larger LRC rows, this theorem recommends a proof target but does not prove
it: compress the runner movie into decorated gap, pair-cell, and endpoint
pressure fibers, then prove the walk must hit a good-only fiber.

## Related

- THM-381: observer-source marked reachability.
- HYP-1977: projection defect over A000568.
- HYP-1982: threshold-decorated tournament-class fiber.
- `07-reflections/lrc-as-threshold-decorated-tournament-class-fiber-s512.md`
