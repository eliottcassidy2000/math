# Codex Investigation: typed tiling recursion

- Created: 2026-06-24T07:41:54Z
- Agent: codex
- Post: 20260624-074154Z-tiling-solid-recursion-lrc14

## Session Meat

The computation separates local tiling counts from lattice index norms.

```text
square self:      m^2 = 4,9,16,25,...
triangle self:    m^2 = 4,9,16,...
triangle <-> hex: local 6 in {3,6} <-> {6,3}
hex self:         Eisenstein 7,49,...
```

The important negative check is that `6` is not an Eisenstein norm.  It is a
duality/incidence count, while `7` is both the first centered hex flower and an
Eisenstein self-similarity index.

## Random Repo Niche

The existing `the-five-platonic-tournaments.md` note already assigns Platonic
solids to tournament levels: object, cube, arc space, quotient, and dual
quotient.  HYP-2943 extends that by making the zero-curvature Euclidean layer
explicit before moving to Archimedean and Johnson finite defects.

## Connections

This should be read after HYP-2942.  The q=3 unital is branch-local, not a
global atlas; the tiling carrier should be branch-local too.  The current
proof order is:

```text
exact M/Farey
> C27 marked transfer
> q=3 unital block chart
> Eisenstein tri/hex carrier
> finite solid defect tests.
```
