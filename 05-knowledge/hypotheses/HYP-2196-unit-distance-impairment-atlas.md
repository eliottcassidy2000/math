# HYP-2196: Unit-Distance Progress Should Use Small Impairment Atlases To Find Load-Bearing Side Channels

**Status:** OPEN, supported by S623 small-carrier experiments.

**Claim.** Unit-distance construction and proof work should not only widen the best
beam or strengthen the graph-only sieve. It should deliberately impair small
carriers and record which damaged side channels change the optimum. The useful
side channels are the ones whose removal changes dense small examples: compact
frontier shape, high-gain extension packets, direction support, canonical orbit
budget, deletion-core resilience, and obstruction labels.

This complements HYP-2194's incoming impairment spectroscopy. HYP-2194 isolates
Moser direction-pair shadow prices and gain caps through `n=14`; HYP-2196
adds triangular controls, beam-width starvation, ranking-policy damage,
canonicalization damage, and a broader damage-response ledger.

This is the unit-distance analogue of several earlier repo lessons:

- LRC quotient compression only works after reattaching owner/carry/pinch side
  channels.
- Tournament scalar counts miss packet structure such as SCCs, cycles, and
  Hamiltonian-path fibers.
- Adversarial cauldron schedules reveal a forgotten resource by changing the
  turn word.
- Equidecomposability separates equal volume/count from predicate-preserving
  construction class.
- S617 frontier gains show that the right observable can be state-local even
  when the proof still needs geometry.

## Evidence

S623 adds `04-computation/unit_distance_impairment_atlas_s623.py` and stores
`05-knowledge/results/unit_distance_impairment_atlas_s623.out`, using the
incoming S622/HYP-2194 spectroscopy lab as the narrower direction-dropout
baseline.

### Width Impairment

On the Moser carrier at target `14`, the healthy policy is deliberately starved
of beam width:

| Width | Best true edges | Span | Frontier-gain histogram |
|-------|-----------------|------|-------------------------|
| `1` | `28` | `6` | `{1: 175, 2: 9, 3: 1}` |
| `3` | `29` | `6` | `{1: 174, 2: 10}` |
| `10` | `29` | `6` | `{1: 174, 2: 10}` |
| `30` | `33` | `4` | `{1: 111, 2: 34, 3: 1, 4: 1}` |
| `60` | `33` | `4` | `{1: 111, 2: 34, 3: 1, 4: 1}` |

The small carrier already shows a threshold: width is not just more samples, it
preserves future high-gain shape.

### Policy Impairment

At target `12`, the triangular carrier is robust under the tested ranking
impairments: all policies found `24` edges. The Moser carrier is not:

- healthy: `27` edges, span `4`, gain histogram `{1: 109, 2: 25, 3: 1}`
- edge-only: `27` edges, but span `10`, so compactness is lost even when the
  scalar edge count survives
- sprawl-bias: `24` edges, a loss of `3`
- balance-bias and future-gain-bias: `27` edges, span `4`

The impairment exposes compact frontier geometry as a retained invariant, not
mere aesthetics.

### Direction-Drop Jackknife

Dropping any one of the three triangular antipodal directions at target `12`
still reaches `24` edges. In the Moser carrier at target `10`, six of the nine
antipodal directions are edge-critical: dropping any of
`(-1,0,0,0)`, `(-1,1,0,0)`, `(0,-1,0,0)`, `(0,0,-1,0)`,
`(0,0,-1,1)`, or `(0,0,0,-1)` lowers the best from `20` to `19`, while the
other three lose `0`.

This suggests direction-support jackknife certificates as a construction/proof
tool: a dense candidate should carry a profile of which unit-vector directions
are load-bearing.

### Gain-Ceiling Impairment

On the Moser carrier at target `12`, forbidding high-gain moves produces:

| Gain ceiling | Best true edges |
|--------------|-----------------|
| `1` | `12` |
| `2` | `21` |
| `3` | `25` |
| `4` | `27` |

This makes high-gain extension packets a first-class proof obligation. For the
`n=22` frontier from S614/S617, a `61`-edge witness must be explained through
precisely this kind of gain-threshold extension structure over dense `21`-cores.

### Canonicalization Impairment

For the triangular carrier at target `10`, translation-only canonicalization
creates `1797` last-layer unique children, while D6 canonicalization creates
`1428`; both find `19` edges. The `1.26x` duplicate-work ratio is small here but
already measurable. Orbit budget is therefore a method impairment to track
before scaling.

## Technique Program

1. **Damage-response ledger.** For every small carrier and target size, record
   which impairments preserve the best edge count and which ones change it.
2. **Direction jackknife certificates.** Attach a per-direction loss vector to
   dense candidates and prove that high-edge constructions need a specific
   direction-support profile.
3. **Gain-threshold extension solver.** Replace raw widening with a solver that
   asks whether a core admits gain `g` or higher under all side-channel filters.
4. **Deletion-core resilience score.** For each dense core, record how many
   deletions preserve extendability, direction support, and obstruction status.
5. **Orbit-budget accounting.** Treat automorphism/canonicalization loss as a
   search resource, not just an implementation detail.

## Tournament Analysis

Vertices are impairment/repair lenses, not points or unit edges:

- extension gain ledger
- direction-drop jackknife
- canonical-orbit repair
- deletion-core impairment
- gain-ceiling adversary
- raw wider beam
- triangular-only baseline
- count-only quotient

The S623 lens tournament is transitive: score histogram
`{0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}`, zero directed
`3`-cycles, singleton SCCs, and one Hamiltonian path. The ranking places
extension gain ledger, direction jackknife, and canonical-orbit repair above
raw wider beam.

## Assumption Challenge

Potential tournament vertices include points, unit directions, frontier
candidates, gain packets, deletion cores, obstruction filters, construction
classes, and proof obligations. S623 uses impairment lenses because they
preserve the question "which method loss changes the dense small optimum?"

This quotient preserves small-size response to controlled loss. It destroys
exact planar embedding outside the finite carrier, full graph isomorphism, and
proof-grade totally-unfaithful obstruction certificates. The challenged
assumption is that the next unit-distance improvement is simply a wider beam.
The atlas says width matters, but only as one retained side channel among
direction support, gain packets, orbit budget, and deletion resilience.

## Next Test

Run the same impairment atlas around the S617 retained `21`-cores:

1. select each `56`/`57`-edge retained core,
2. compute gain-`4` and gain-`5` extension packets,
3. jackknife directions and orbit representatives,
4. attach totally-unfaithful obstruction labels, and
5. ask whether any side-channel-complete extension lane survives to `61`.
