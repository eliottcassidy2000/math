---
id: THM-377
name: lrc-n14-n18-selected-pressure-acyclicity
status: PROVED (exact finite computation, selected-row scope)
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on:
  - THM-372
  - THM-375
  - THM-376
---

# THM-377: Selected n=14/n=18 pressure lifts have no strict SCC or 3-cycle

## Statement

In the selected-row audit of

```text
04-computation/lrc_n14_n18_tournament_pingpong_s492.py
```

the following eight rows were tested:

```text
initial n=14
initial n=18
n=14 lpd ladder, scale 7 skip 6
n=18 lpd ladder, scale 9 skip 8
n=14 gate ladder, scale 14 skip 6
n=18 gate ladder, scale 18 skip 8
n=14 single-gate repair: replace 6 by 14*16
n=18 single-gate repair: replace 8 by 18*18
```

At the selected exact endpoint/gap times, three incomplete pressure
tournaments were computed:

```text
k1 relief:      nearest-distance deletion relief
k2 relief:      two-nearest-distance-sum deletion relief
deficit relief: two-neighbor threshold-deficit deletion relief
```

For both `n=14` and `n=18`, none of the selected rows produced a directed
3-cycle or a strict strongly connected component of size greater than `1` in
any of the three pressure lifts:

```text
n=14:
  k1      cyclic-or-SCC rows 0/23
  k2      cyclic-or-SCC rows 0/23
  deficit cyclic-or-SCC rows 0/23

n=18:
  k1      cyclic-or-SCC rows 0/23
  k2      cyclic-or-SCC rows 0/23
  deficit cyclic-or-SCC rows 0/23
```

## Scope

This theorem is intentionally finite and selected-row scoped.  It does not
claim that all perturbations of the rows are pressure-acyclic, nor that every
time chamber in the full runner movie is pressure-acyclic.  It certifies the
exact selected audit used by S492.

## Verification Record

The stored output is:

```text
05-knowledge/results/lrc_n14_n18_tournament_pingpong_s492.out
```

The script was run through `run_and_save.sh`, and all arithmetic in the
underlying LRC endpoint/gap computations is exact over `Fraction`.

## Proof

For each selected time and each pressure lift, the script constructs the
directed graph by comparing exact rational relief values for every unordered
pair of runners.  Ties are left as missing arcs.

It then computes:

1. the number of strict directed 3-cycles;
2. the largest strict strongly connected component.

The stored output reports zero directed 3-cycles and largest strict SCC equal
to `1` for every selected row and every pressure lift.  This proves the finite
claim.

## Significance

This gives a formal negative certificate for the current n=14/n=18
near-counterexample rows: strengthening nearest-only pressure to k2 and
threshold-deficit pressure still does not produce a mobile cyclic core.  The
next disproof-like target is therefore sharply defined:

```text
find k2_largest_scc > 1 or deficit_largest_scc > 1
```

after private endpoint and pressure leaves are removed.

## Related

- HYP-1930
- HYP-1950
- THM-375
- THM-376
- `07-reflections/lrc-n14-n18-tournament-pingpong-s492.md`
