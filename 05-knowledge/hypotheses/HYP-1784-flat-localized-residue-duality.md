---
id: HYP-1784
status: EXPLORATORY
source: codex-2026-05-30
related:
  - HYP-1779
  - HYP-1780
  - HYP-1781
  - HYP-1782
  - THM-212
  - THM-344
---

# HYP-1784: Flat-Versus-Localized Residue Duality

## Statement

Structured extremal tournaments separate into two residue regimes:

```text
Paley-like / flat:
  cycle mass, deletion loss, score, and code shadows are spread uniformly.

Core-like / localized:
  cycle mass is concentrated around one kill or near-kill vertex.
```

Both regimes can produce special Hamiltonian-path behavior, but for opposite
reasons.  Paley-like objects optimize by distributing residue evenly; core-like
objects optimize or unlock gaps by concentrating residue until the conflict
graph becomes almost complete.

## Evidence

Paley T7:

```text
scores = (3,3,3,3,3,3,3)
all vertices have deletion loss rho = 60/80 = 0.75
support_count = 36
support_excess = 44
even-graph projection = 14 edges, degree sequence 4^7
```

Interval T7 has the same score sequence and support count but a different
residue profile:

```text
cycles = 59
support_excess = 23
disjoint_pairs = 14
even-graph projection = 7 edges, degree sequence 2^7
```

The H=63 classes are the opposite extreme:

```text
one vertex has rho = 1
delete it and Omega vanishes
Omega = K31
```

THM-025 sits between these regimes:

```text
one vertex has rho = 92/94
delete it and a two-cycle residue remains
```

## Predictions

1. H-maximizer families should often be either flat-residue or localized-core
   families, with mixed regimes appearing near transitions.
2. Real-root failures should be enriched in near-localized residues with small
   but nontrivial deletion remainders.
3. Paley-to-interval crossover should be visible as a shift from flat
   multiplicity residue to better disjoint-cycle packing.

## Sources

- `07-reflections/residue-calculus-feedback-loop.md`
- `07-reflections/projection-defect-as-common-residue.md`
- `07-reflections/paley-gives-dual-codes.md`
- `04-computation/projection_defect_bridge_s12.py`
- `04-computation/paley_codes_s306.py`
