---
id: HYP-1900
status: OPEN
source: codex-2026-05-31-S430
related:
  - THM-357
  - THM-358
  - THM-369
  - HYP-1866
  - HYP-1868
  - HYP-1880
  - HYP-1881
  - HYP-1890
  - HYP-1891
  - HYP-1892
  - HYP-1893
---

# HYP-1900: LRC and tournament recursion share a labelled incidence core

## Statement

Most useful formulations of the Lonely Runner Conjecture are shadows of one
finite incidence object:

```text
rows       = constraints that must be protected/covered
columns    = available protectors
witness    = uncovered row or positive gap
bad object = every row covered with no residual gap
```

For LRC, rows can be forbidden endpoints, coarse denominators, p-adic frontier
leaves, zonotope lattice cells, distance-graph color conflicts, or view
obstruction cells.  For tournaments, rows can be Hamiltonian-path cuts,
cycle-conflict rows, bucket transport rows, projection defects, or good-cut
constraints.

The conjectural bridge is:

```text
LRC counterexample = leafless labelled protection hypergraph
tournament SC/core = every cut protected by crossing arcs
```

The two theories differ in what a protector is allowed to be, but they share
the same obstruction grammar.

## Evidence

`lrc_lens_atlas_s430.py` records 20 lenses.  Eleven lie in a single
formulation/transport component, with lens status distinguishing strict
restatements from implication or necessary-constraint lenses:

```text
runner frame
Diophantine norm
forbidden intervals
endpoint hypergraph
coarse sieve
torus/subtorus
view obstruction
distance graph/circular coloring
flows/matroids
zonotope covering radius
lonely-runner spectrum
```

Ten proof technologies form another connected component around the same
obstruction:

```text
finite checking
prime-product sieve
coarse denominator sieve
integer program/set cover
endpoint hypergraph
Bruhat-Tits frontier
adelic gap-debt product
natural row/column modes
product-sum/Egyptian ledgers
distance graph coloring
```

The repo already has exact tournament-side analogues:

- THM-357 turns LRC into endpoint rows protected by speed intervals.
- THM-369 turns coarse rational times into denominator rows.
- S386 shows tournament good cuts are protected by backward arcs crossing a
  Hamiltonian-path cut.
- HYP-1890 turns LRC into an integer-program row/column system.
- HYP-1881/HYP-1891 split the row/column modes by the same `x+2` and `x*2`
  natural-number recursion that organizes tournament add-two and blowup modes.

## Predictions

1. The recent finite-checking prime-product sieve should have a dual
   row-weight formulation inside the S420 integer-program matrix.
2. LRC disproof searches should generate candidate labelled incidence cores
   before searching for speed sets; most candidate cores should fail
   realizability.
3. The Giri-Kravitz spectrum recursion should have a tournament analogue in
   vertex deletion, cut quotient, or root-spectrum accumulation.
4. View-obstruction and zonotope lenses should produce a tournament tiling
   polytope whose facets are good cuts and odd-cycle conflict rows.
5. A complete proof likely needs two commuting incidence invariants: row-mode
   p-adic/frontier mass and column-mode odd-payload transfer.

## Sources

- `04-computation/lrc_lens_atlas_s430.py`
- `05-knowledge/results/lrc_lens_atlas_s430.out`
- `07-reflections/lrc-lens-atlas-s430.md`
- THM-357, THM-358, THM-369
- HYP-1881, HYP-1890, HYP-1891, HYP-1892, HYP-1893
