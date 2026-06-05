---
id: HYP-2246
status: OPEN enumeration/invariant conjecture with exact evidence through n=8
source: codex-2026-06-05-S671
related:
  - HYP-2245
  - HYP-2236
  - HYP-2228
  - HYP-2187
  - HYP-2209
---

# HYP-2246: Endpoint Half-Filter Traces Enumerate Tournament Classes

## Claim

The HYP-2245 ultrafilter/descent picture gives a concrete A000568 speedup:
enumerate one-vertex endpoint extensions by automorphism orbits in the incident
Boolean cube, then classify candidates by a three-state deletion/filter trace.

For a tournament `T` on `n` vertices, define

```text
side_T(v) =
  L if 2*outdeg_T(v) < n-1
  M if 2*outdeg_T(v) = n-1
  U if 2*outdeg_T(v) > n-1
```

and the half-filter deletion trace

```text
Phi(T) = multiset_v ( iso_class(T-v), side_T(v) ).
```

S671 verifies that `Phi` is injective on tournament isomorphism classes through
`n=8`.

This is the user's old three-state automaton, now with a concrete job:

```text
left / middle / right = below-half / tied-half / above-half deleted-owner state.
```

## Endpoint-Orbit Enumeration

The script `04-computation/ultrafilter_endpoint_enumeration_s671.py` builds the
exact tournament-class tower through `n=8` by:

1. keeping one representative of each parent class on `n` vertices;
2. taking one representative from each `Aut(parent)`-orbit in the incident
   bit cube `{0,1}^n`;
3. extending by that incident pattern;
4. canonicalizing only for the audit.

It recovers the known A000568 values through `n=8`:

| child n | A000568 | endpoint-orbit candidates | duplicate factor | fixed-path tilings / candidates | labelled tournaments / candidates |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 2 | 2.000 | 0.500 | 1.000 |
| 3 | 2 | 4 | 2.000 | 0.500 | 2.000 |
| 4 | 4 | 12 | 3.000 | 0.667 | 5.333 |
| 5 | 12 | 48 | 4.000 | 1.333 | 21.333 |
| 6 | 56 | 296 | 5.286 | 3.459 | 110.703 |
| 7 | 456 | 3040 | 6.667 | 10.779 | 689.853 |
| 8 | 6880 | 54256 | 7.886 | 38.653 | 4947.572 |

The one-step `n=8 -> 9` estimate is:

```text
Aut-orbit endpoint candidates = 1,716,608
known A000568(9)              =   191,536
candidate/known ratio         = 8.962
raw endpoint patterns          = 1,761,280
```

So automorphism-orbit endpoint generation is already a compact replacement for
the fixed-path cube and a massive replacement for labelled enumeration.

## Trace Collision Evidence

On finished classes, the trace comparison through `n=8` is:

| n | classes | `score,c3,SCC` mixed | card multiset mixed | half-filter mixed | paired-score mixed |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 0 | 1 | 0 | 0 |
| 4 | 4 | 0 | 1 | 0 | 0 |
| 5 | 12 | 2 | 1 | 0 | 0 |
| 6 | 56 | 11 | 4 | 0 | 0 |
| 7 | 456 | 41 | 0 | 0 | 0 |
| 8 | 6880 | 136 | 2 | 0 | 0 |

The coarse scalar signature `(score histogram, c3, SCC profile)` collapses
badly by `n=8`: only `167` groups for `6880` classes, with a max bucket of
`577`.

The unpaired deletion-card multiset is nearly enough but not quite: at `n=8`
it has two collision buckets, four colliding classes.

The half-filter trace repairs all tested collisions, using only `L/M/U` rather
than exact deleted score.

The two `n=8` card-deck collisions are especially informative.  Both buckets
consist of regular strongly connected classes with identical coarse payloads:
score histogram `((3,4),(4,4))`, `c3=20`, and matching Hamiltonian-path counts
(`633` in one bucket, `645` in the other).  The card deck therefore loses a
pure orientation address, not a visible scalar property.  The half-filter
repairs the loss by detecting a symmetric transfer of deleted-card IDs between
`L` and `U`; for example one bucket differs only by exchanges such as
`(230,L) <-> (230,U)` and `(232,L) <-> (232,U)`.

This is a concrete local lemma target: paired deletion cards may fail only
when the visible deck forgets which side of the half-filter owns the card.
If every future card-deck collision has this form, `Phi` is not a lucky
refinement but a canonical orientation repair.

## Candidate Fallback Evidence

For each raw endpoint-orbit candidate, S671 computes traces before child
canonicalization and asks which trace buckets mix true child classes.  Through
the tested range, the half-filter trace has zero mixed buckets:

| transition | candidates | `score,c3,SCC` avoid | card multiset avoid | half-filter avoid |
|---:|---:|---:|---:|---:|
| 4 -> 5 | 48 | 52.1% | 83.3% | 100.0% |
| 5 -> 6 | 296 | 15.5% | 83.8% | 100.0% |
| 6 -> 7 | 3040 | 2.9% | 100.0% | 100.0% |
| 7 -> 8 | 54256 | 0.3% | 100.0% except 24 mixed candidates | 100.0% |

Thus a production enumerator can try:

```text
parent class
  -> Aut(parent)-orbit incident patterns
  -> raw child candidate
  -> Phi(raw child) using only n-vertex card canonicalization
  -> full (n+1)-canonicalization only if Phi bucket is mixed.
```

The audit used full child canonicalization to verify purity.  The proposed
production use computes `Phi` directly from the raw child, then calls the
expensive child canonicalizer only for mixed trace buckets.

## Property Payload

The same tower computes tournament properties while enumerating:

| n | distinct H | distinct c3 | SCC profiles | H range | c3 range |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 1..3 | 0..1 |
| 4 | 3 | 3 | 3 | 1..5 | 0..2 |
| 5 | 7 | 6 | 4 | 1..15 | 0..5 |
| 6 | 19 | 9 | 6 | 1..45 | 0..8 |
| 7 | 77 | 15 | 8 | 1..189 | 0..14 |
| 8 | 320 | 21 | 11 | 1..661 | 0..20 |

So the method is not only for A000568 counts; it can stream `H`, `c3`, SCC, and
later OCF/path-homology payloads.

## Tournament Analysis

Vertices are enumeration routes:

1. `paired_deletion_filter_trace`
2. `endpoint_aut_orbit_generation`
3. `modular_substitution_dp`
4. `half_filter_trace`
5. `THM410_interval_ledger`
6. `raw_fixed_path_cube`
7. `raw_labeled_enumeration`

Observable:

```text
(exactness,
 canonical-work avoidance,
 invariant retention,
 scaling,
 implementation readiness,
 LRC transfer)
```

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

Top route: `paired_deletion_filter_trace`.  The half-filter trace ranks behind
the exact paired trace but is the surprising minimal side channel.

## Challenged Assumption

Do not assume the enumeration vertices are arcs or tournaments.  In this
session the useful vertices were endpoint incident patterns, automorphism
orbits, deletion cards, deleted-owner side states, interval ledgers, module
blocks, and proof obligations.

Preserved predicate:

```text
tournament isomorphism class, and streamed H/c3/SCC properties.
```

Destroyed information:

```text
raw labels and exact deleted scores when only L/M/U side is retained.
```

The tested surprise is that exact deleted score was unnecessary through `n=8`.

## Next Tests

1. Push the half-filter trace to `n=9` with a faster canonical engine
   (`nauty/traces`, C++, or a cached Python/NumPy hybrid).
2. Extract and study the two `n=8` card-multiset collision pairs; prove why
   the `L/M/U` side bit separates them.
3. Combine the trace with modular decomposition: if a candidate has modules,
   classify blocks by half-filter traces before any full flattening.
4. Use the same `L/M/U` owner-state trace inside LRC14 owner/carry fibers.
