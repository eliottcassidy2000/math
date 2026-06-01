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
  - THM-381
  - THM-384
  - THM-387
---

# HYP-2001: n=14 LRC is a permutohedral sweep-hitting problem

## Stub

At total denominator `n=14`, the observer-marked LRC movie can be lifted from
scalar gap tracking to a path through the braid arrangement / permutohedral
normal fan of the `14` circle points.  A no-lonely counterexample would be a
closed one-parameter rational sweep whose chamber sequence keeps the two
observer-adjacent gaps out of the long-long source facet from THM-384.

The proof attempt reserved here will test whether the seven CRT blocker
classes from the `14=2*7` quotient can hand off indefinitely inside this
permutohedral fan, or whether the oriented facet-crossing ledger forces an
LL/source chamber.  This is meant to refine:

- HYP-1997: the metagraph is faithful but non-reducing because free walks are
  too large;
- HYP-1999: the fixed-`n` source target is a tiny Ferrers interval menu;
- HYP-2000: arc/tile flips are dependent and should be read in recursive
  ranking/tiling coordinates.

## Missing

The current stub has not yet proved the n=14 case.  The intended computation
is `04-computation/lrc_n14_permutohedron_s526.py`, comparing candidate
families by:

1. chamber/facet sequence in the braid fan,
2. observer-adjacent gap sign state `(L,R)`,
3. CRT blocker owner and handoff order,
4. permutohedral descent or circulation certificates,
5. Tournament Analysis fingerprints of the row families.

The desired proof certificate is a local obstruction: every attempted
blocker-to-blocker handoff either crosses an LL/source facet or exports
endpoint debt already measured in the n=14 ladder audits.

