---
id: HYP-1798
status: EXPLORATORY
source: opus-2026-05-30-S5
related:
  - HYP-1795
  - HYP-1785
  - HYP-1796
  - THM-025
  - THM-344
---

# HYP-1798: Standard Contrast Lens Separability

## Statement

The standard tournament anomaly set is separable by three coarse lenses:

```text
localized residue kill,
broad phase/fiber contrast,
incidence/rank interface.
```

Concretely, the set

```text
transitive, Paley, interval, H=63 single-core, THM-025, n=6 H=37 trap
```

should not collapse to one explanation.  Exact kills, near-kills,
phase/fiber crossovers, and landscape traps occupy different regions of a
`residue_vector + phase_vector + incidence_vector` feature space.

## Evidence

`residue_phase_incidence_contrast_s5.py` gives the first small contrast table:

- the two `H=63` single-core classes have `rho=1` and deletion-residue
  `alpha=(0,0,0)`, hence exact residue kill;
- THM-025 has `rho=.979`, residue `alpha=(2,1,0)`, and residue rank `2`, hence
  a dangerous near-kill;
- Paley and interval `T7` have broad rank-2 residues, equal score variance, and
  equal `t3`, but different odd-cycle fiber/disjointness profiles;
- the `n=6`, `H=37` local trap has no exact kill and only rank-1 max-loss
  residue, so the interesting feature is landscape/phase rather than
  projection kill.

## Predictions

1. Adding lightweight phase features should separate Paley/interval and local
   traps better than adding more deletion-residue features.
2. Adding residue-rank features should separate THM-025-style failures better
   than phase features alone.
3. Endpoint-transfer failures will remain poorly explained until incidence
   features are included.

## Test Plan

1. Extend `tournament_tda.py` with explicit `residue_*`, `phase_*`, and
   operation-level `incidence_*` blocks.
2. Benchmark ablations on the standard contrast set.
3. Require every future broad synthesis to add or reuse a contrast table before
   promoting a single-lens explanation.

## Sources

- `04-computation/residue_phase_incidence_contrast_s5.py`
- `05-knowledge/results/residue_phase_incidence_contrast_s5.out`
- `07-reflections/applied-residue-phase-incidence-programs.md`
