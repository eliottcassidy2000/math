---
id: HYP-1824
status: EXPLORATORY
source: codex-2026-05-31-S369
related:
  - HYP-1823
  - HYP-1825
  - THM-363
---

# HYP-1824: The LRC 56 missed cells encode a six-vertex chirality boundary

## Statement

The `56` missed cells of the S367 fourteen-runner quotient extremal should be
compared to the six-vertex tournament layer through mirror/chirality structure,
not through a raw bijection of classes.

The proposed bridge is:

```text
LRC missed cells:
  56 = 7 odd shifts * 8 alpha stencils
     = 14 outer mirror-pair cells + 42 interior cells
     = 7 shifts * (2 outer caps + 6 interior positions)

six-vertex tournaments:
  56 = 12 self-converse classes + 44 chiral classes
     = 35 strong classes + 21 non-strong classes
     = 48 five-vertex perspectives + 8 perspective-gap classes.
```

The shared missing object is an eight-fold mirror/perspective residue that
appears exactly when moving from five-object structure to six-object structure.

## Evidence

`lrc_tournament_56_bridge_s369.py` exactly enumerates tournament isomorphism
classes on five and six vertices and reads the S367 missed-cell ledger.

Tournament side:

```text
T(5)=12
T(6)=56
self-converse/chiral split at n=6: 12 + 44
strong/non-strong split at n=6: 35 + 21
sum vertex orbits over T(5): 48
perspective gap: 56 - 48 = 8
```

LRC side:

```text
missed cells: 56
odd shifts: 7
alpha stencils: 8
mirror pairs: 4
outer mirror pair: 14 cells
interior stencils: 42 cells
```

Thus the user's decompositions are structurally meaningful:

- `56-12=44` matches the self-converse/chiral tournament split.
- `56-14=42` matches the LRC outer-pair/interior split.
- the unexplained `8` in the old five-to-six perspective failure matches the
  eight LRC stencils.

S370 extends this bridge to the Paley/Fano and base-42 layers.  Paley `T7`
has directed odd-cycle counts `14+42+24=80`, support counts `14+21+1=36`,
and support excess `44`.  Base 42 has `phi(42)=12`, with its coprime residue
classes split as `8` easy and `4` hard in the Erdos-Straus cover.  This
supports the broader HYP-1825 grammar: `12` is inherited symmetric core, `44`
is chiral/support residue, `42` is doubled boundary/interior, and `8` is
projection failure.

## Predictions

1. The eight LRC stencils should carry the same symmetry type as the eight
   extra six-vertex perspective classes not accounted for by the old
   five-vertex rooted-perspective heuristic.
2. The outer LRC mirror pair should correspond to a decomposable or boundary
   tournament layer, while the six interior stencils should correspond to the
   `42 = 2*21` doubled non-strong boundary.
3. A proof of HYP-1823 should split by mirror/chirality: first handle the
   outer fourteen-runner mirror pair, then handle the forty-two-cell interior.

## Sources

- `04-computation/lrc_tournament_56_bridge_s369.py`
- `05-knowledge/results/lrc_tournament_56_bridge_s369.out`
- `07-reflections/lrc-tournament-56-bridge-s369.md`
- `04-computation/chirality_perspective_atlas_s370.py`
- `05-knowledge/results/chirality_perspective_atlas_s370.out`
- HYP-1825
- HYP-1823
- THM-363
