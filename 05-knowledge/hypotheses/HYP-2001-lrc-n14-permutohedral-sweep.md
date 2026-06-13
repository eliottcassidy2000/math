---
id: HYP-2001
status: OPEN
source: codex-2026-06-01-S526
related:
  - HYP-1986
  - HYP-1987
  - HYP-1990
  - HYP-1991
  - HYP-1997
  - HYP-1999
  - HYP-2000
  - HYP-2002
  - HYP-2003
  - THM-381
  - THM-384
  - THM-387
---

# HYP-2001: n=14 LRC is a labelled permutohedral sweep-hitting problem

## Statement

At total denominator `n=14`, the observer-marked LRC movie can be lifted from
scalar gap tracking to a path through the braid arrangement / permutohedral
normal fan of the `14` circle points.  A no-lonely counterexample would be a
closed one-parameter rational sweep whose chamber sequence keeps the two
observer-adjacent gaps out of the long-long source facet from THM-384.

The raw version of the hoped-for proof is false: blocker owners can hand off
inside the permutohedral fan without crossing the source facet.  The surviving
claim is labelled.  Every source-avoiding handoff should export endpoint debt
to quotient/depth labels, and after the private endpoint leaves are peeled no
closed blocker-owner circulation should remain realizable by a primitive n=14
speed system.

This refines:

- HYP-1997: the metagraph is faithful but non-reducing because free walks are
  too large;
- HYP-1999: the fixed-`n` source target is a tiny Ferrers interval menu;
- HYP-2000: arc/tile flips are dependent and should be read in recursive
  ranking/tiling coordinates.

It sits between the concurrent permutohedral summaries HYP-2002 and HYP-2003:
HYP-2002 identifies the closed geodesic and central-box formulation, while
HYP-2003 localizes the n=14 covering obstruction at the AP wall-only extreme.
HYP-2001 records the labelled handoff debt that survives inside that geometry.

## Evidence

`lrc_n14_permutohedron_s526.py` compares candidate families by:

1. chamber/facet sequence in the braid fan,
2. observer-adjacent gap sign state `(L,R)`,
3. CRT blocker owner and handoff order,
4. permutohedral descent or circulation certificates,
5. Tournament Analysis fingerprints of the row families.

Selected rows:

```text
initial 1..13: source_measure=0, wall_sources=6
seven-ladder:  source_measure=1/143, chambers=2034, raw_handoff_defects=142
S380 gate:     source_measure=1/143, chambers=4224, raw_handoff_defects=310
double gate:   source_measure=1/143, chambers=8506, raw_handoff_defects=600
```

The hard ladders have blocker-state max SCC `3`, so a pure acyclic handoff
proof is unavailable.  In contrast, all `560` primitive systems with
`max_speed<=16` still hit the source; only the initial segment is wall-only.

## Predictions

1. The hard ladder family preserves the open source measure `1/143` while
   scaling chamber count and raw handoff defects by the dyadic gate depth.
2. Raw blocker-owner SCCs are not proof obstructions by themselves; the useful
   obstruction is an SCC after endpoint-private leaves and quotient-depth debt
   are charged.
3. A successful n=14 proof should be expressible as a Hall/Farkas certificate
   on labelled permutohedral handoffs: every closed source-avoiding circulation
   violates a private endpoint row.

## Sources

- `04-computation/lrc_n14_permutohedron_s526.py`
- `05-knowledge/results/lrc_n14_permutohedron_s526.out`
- `07-reflections/lrc-n14-permutohedral-proof-attempt-s526.md`
- HYP-2002 / `07-reflections/lrc-permutohedron-geometry-and-n14-attempt-s521o.md`
- HYP-2003 / `05-knowledge/hypotheses/HYP-2003-lrc-n14-covering-wallonly-AP.md`
