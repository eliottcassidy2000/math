---
id: HYP-1981
status: PARTLY CANONIZED
source: oracle-2026-06-01-S511; formalized by codex-2026-06-01-S513
related:
  - THM-381
  - HYP-1977
  - HYP-1978
  - HYP-1979
  - HYP-1982
  - HYP-1987
  - HYP-1988
  - THM-369
  - THM-385
---

# HYP-1981: LRC is source-reachability in the observer-marked A000568 quotient

## Statement

The exact LRC witness bit is a marked tournament source condition.

For speeds `{0,v_1,...,v_{n-1}}`, orient the observer incident edges by

```text
0 -> i  iff  ||v_i t|| >= 1/n.
```

Orient moving-runner pairs by a complete comparator such as the half-turn
runner phase clock.  THM-381 proves the canonized core:

```text
observer lonely at t  <=>  marked observer is a source.
```

It also proves that the source-marked target classes on `n` vertices are in
bijection with ordinary tournament classes on `n-1` vertices, hence counted by
`A000568(n-1)`.

The open part of the hypothesis is the reachability/classification layer:

```text
Which source-marked classes are reachable by arithmetic runner clocks?
Why must every primitive speed system reach one?
```

## Evidence

The S511 script

```text
04-computation/lrc_observer_source_tournament_s511.py
```

checked the direct equivalence on representative `n=4,5,6` families with zero
mismatches.  The computation is a sanity audit rather than the proof; THM-381
now supplies the proof.

The same session recorded the target count convention in

```text
05-knowledge/results/lrc_qr_paley_speed_tournament_s511.out
```

with

```text
n=4..10 source target sizes = A000568(3..9)
                           = 2,4,12,56,456,6880,191536.
```

S509/S510 supply the negative background: raw unmarked and even
observer-pointed half-turn classes can mix safe and unsafe states.  The
observer-source construction repairs that by using the LRC threshold on the
observer incident arcs.

S517 adds THM-385: the observer-source construction also records exact
distance-to-source.  The observer indegree is the number of runners currently
blocking the observer, so source, almost-source, and every observer-score layer
have direct LRC meanings.  The S517 audit found zero marked-class mixing by
blocker count in bounded `n=4,5,6,7` windows, while runner-subtournament
classes did mix blocker counts.

## Predictions

1. The reachable source classes are a small subset of `A000568(n-1)`, likely
   those whose runner sub-tournament is realizable by `n-1` points inside a
   safe arc of length `1 - 2/n`, with boundary strata included.
2. Initial segments and known tight examples should appear as boundary-only
   source touchers: the marked walk reaches source classes only at equality
   walls, not in open cells.
3. A genuine counterexample would be a closed marked tournament walk avoiding
   every source-marked class while also carrying the THM-380 endpoint-pressure
   obstruction data.
4. The operation-grid labels from HYP-1976/HYP-1980 should explain which
   source classes are reachable: addition drives wall crossings, while
   multiplication and odd-core/Burnside structure constrain the class menu.

## Next Tests

1. Enumerate reachable source classes for total `n<=6` and compare them with
   the full `A000568(n-1)` target.
2. Add source-class ids to the n=14/n=18 corridor scripts, including exact
   equality-wall hits.
3. Combine source reachability with the threshold-decorated fibers of HYP-1982:
   a source class is the witness target, while decorated bad/mixed fibers
   explain how an A000568 projection can hide or reveal that target.
4. Add observer-score, side-defect, and observer 2-king repair layers to the
   reachable-source classification, especially for tight and near-tight rows.

## Sources

- `01-canon/theorems/THM-381-lrc-observer-source-marked-reachability.md`
- `04-computation/lrc_observer_source_tournament_s511.py`
- `05-knowledge/results/lrc_observer_source_tournament_s511.out`
- `05-knowledge/results/lrc_qr_paley_speed_tournament_s511.out`
- `07-reflections/lrc-as-source-reachability-in-marked-A000568-s511.md`
- HYP-1977
- HYP-1979
- HYP-1982
- HYP-1987
- HYP-1988
