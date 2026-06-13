---
id: HYP-1825
status: EXPLORATORY
source: codex-2026-05-31-S370
related:
  - HYP-1823
  - HYP-1824
  - THM-363
---

# HYP-1825: Chirality-perspective residues organize the 12/42/44/56 layer

## Statement

Several repo structures with counts near `56` are governed by the same
projection-residue grammar:

```text
12 = inherited symmetric core
44 = chiral or support-multiplicity residue
42 = doubled boundary/interior layer
8  = projection failure or stencil layer
```

In this grammar, a projection from a richer object to a symmetric/inherited
shadow loses a mirror-sensitive residue.  The residue is often measured as
`44`, while the boundary/interior part of the residue appears as `42`.

## Evidence

`chirality_perspective_atlas_s370.py` gives exact ledgers.

Six-vertex tournaments:

```text
56 = 12 self-converse + 44 chiral
42 = 2 * 21 non-strong classes
12 = T(5) = phi(42) = self-converse T(6)
8  = T(6) - sum(vertex orbits over T(5))
```

S367 Lonely Runner missed cells:

```text
56 = 7 odd shifts * 8 alpha stencils
14 = 7 shifts * 2 outer caps
42 = 7 shifts * 6 interior positions
```

Paley/Fano `T7`:

```text
directed odd cycles = 14 + 42 + 24 = 80
cycle supports      = 14 + 21 + 1  = 36
support excess      = 44
```

Base 42:

```text
phi(42)=12
units mod 42 split as 8 easy + 4 hard in the Erdos-Straus cover
```

## Predictions

1. The eight S367 LRC stencils should match an eight-object residue in the
   five-to-six tournament perspective failure.
2. The `42` interior LRC cells should behave like a doubled boundary layer,
   analogous to `2*21` non-strong six-vertex classes and Paley's two directed
   pentagons per five-subset.
3. HYP-1823 should become easier after quotienting by mirror/chirality: handle
   the outer cap pair first, then prove that every interior residue has a
   positive-margin cell.
4. Other repo occurrences of `44` should usually be support excess, chiral
   residue, or projection defect rather than a primary count.

## Sources

- `04-computation/chirality_perspective_atlas_s370.py`
- `05-knowledge/results/chirality_perspective_atlas_s370.out`
- `07-reflections/chirality-perspective-atlas-s370.md`
- HYP-1824
- HYP-1823
