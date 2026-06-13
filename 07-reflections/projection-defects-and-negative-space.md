# Projection Defects And Negative Space

**Session:** kind-pasteur-2026-05-29-S1
**Computation:** `04-computation/projection_defect_waggly_layers_s1.py`
**Results:** `05-knowledge/results/projection_defect_waggly_layers_s1.out`

## The Unlooked Place

The project has spent enormous attention on named structures: OCF, wiggly edges, blue/black lines, the merged tournament metagraph `G_n/Z_2`, and the even graph metagraph `E_n`. The negative space between them is smaller and stranger: the same tiling hypercube `Q_m` projects to both `G_n/Z_2` and `E_n`, but these two projections do not commute with waggly motion.

That noncommutation is not a bug. It is a new invariant.

For any pair of tilings at Hamming distance `d`, ask two yes/no questions:

- Does the merged tournament class change?
- Does the even graph class change?

This gives four outcomes:

| Outcome | Meaning |
|---------|---------|
| silent_both | both quotients absorb the move |
| tournament_only | tournament class changes, even graph class does not |
| even_only | even graph class changes, tournament class does not |
| joint | both quotients change |

The old bridge work mainly measured this at `d=1`. The new measurement checks every waggly layer `d=1,...,m` for `n=3..6`.

## What The Data Says

The headline is synchronization: joint changes dominate, and they dominate more strongly as `n` grows.

| n | all-layer joint | all-layer silent_both | tournament_only | even_only |
|---|-----------------|-----------------------|-----------------|-----------|
| 4 | 46.43% | 17.86% | 14.29% | 21.43% |
| 5 | 72.32% | 3.57% | 14.58% | 9.52% |
| 6 | 85.40% | 0.77% | 10.05% | 3.78% |

At `n=6`, the ordinary wiggly layer `d=1` has 80.57% joint changes, but the middle layers are even more synchronized:

| d | joint% at n=6 |
|---|---------------|
| 1 | 80.57% |
| 4 | 84.85% |
| 5 | 86.56% |
| 6 | 86.61% |
| 9 | 87.81% |
| 10 | 77.73% |

So `d=1` is not the whole story. The strongest tournament/even-graph agreement lives in the high-entropy middle of the hypercube, while the complement-tiling layer `d=m` is special again.

## Meta-Hypothesis

Many prior mistakes were caused by identifying two nearby projections:

- tiling complement vs tournament complement
- line vs metagraph edge
- blue/black line color vs class-level SC/NS type
- tournament class graph vs even graph class graph
- d=1 wiggly edges vs all waggly layers

The better principle is:

> Tournament theory is organized by projection defects: whenever two natural quotients of the same binary triangle almost agree, the places where they disagree mark the next theorem.

The theorem-looking part is the synchronization rate. The application-looking part is the defect signature. A class, dataset, or engineered codec can be fingerprinted not only by its invariants, but by how its perturbations split across several quotient lenses.

## Engineering Reading

For ranking TDA or tournament fingerprints, `G_n/Z_2` and `E_n` should be treated as two coupled feature channels:

- `joint` moves are robust structural mutations.
- `tournament_only` moves are orientation-sensitive but cycle-space-invisible; they look like score/cut-space changes.
- `even_only` moves are cycle-space changes absorbed by tournament symmetry; these are good candidates for compression redundancy or canonicalization shortcuts.
- `silent_both` moves are true neutral mutations across both lenses.

This suggests a practical feature extractor: for a tournament or local tiling neighborhood, sample perturbations at several Hamming radii and report the four-way defect profile. It is cheap, invariant-aware, and may separate ranking datasets that have the same `H`, score sequence, or Betti profile.

## Where To Look Next

1. Extend the all-layer defect table to `n=7` by sampling, then exact if a faster canonical cache is available.
2. Compute the same four-way table by structured move family: range flips, vertex-star flips, anti-diagonal flips, and random middle-layer flips.
3. Condition the defect profile on spine/ribs/sea position in `G_n/Z_2`.
4. Compare defect profiles for `H`-maximizers, Paley/interval circulants, and near-transitive tournaments.
5. Package the defect profile as a `tournament_tda.py` feature.

The illuminating flashlight here was not another invariant, but a commutator: take two maps we believe we understand and measure where they fail to tell the same story.
