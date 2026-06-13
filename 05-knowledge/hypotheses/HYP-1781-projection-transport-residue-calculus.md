---
id: HYP-1781
status: EXPLORATORY
source: codex-2026-05-30
related:
  - HYP-1779
  - HYP-1776
  - HYP-1777
  - HYP-1778
  - HYP-1780
  - HYP-1783
  - THM-350
  - THM-351
  - THM-353
---

# HYP-1781: Projection-Transport Residue Calculus

## Statement

Many tournament anomalies and extremal families are better predicted by a
combined residue vector than by any one projected invariant.

For a tournament or fixed-path tiling, define a feature block

```text
R(T) = (
  deletion-loss profile,
  support-multiplicity profile,
  Omega alpha-vector,
  even-graph projection where defined,
  quotient-transport leakage by selected mask families,
  endpoint-transfer parity row
)
```

The hypothesis is that `R(T)` separates the known exceptional regimes:

- complete-core H=63 classes;
- THM-025-style non-real-rooted independence polynomials;
- Paley/interval cycle-packing differences;
- good-cut bucket transport boundaries;
- path-homology rank-defect anomalies.

## Evidence

Existing computations already show the same defect shape in different
languages:

- H=63 at n=8 has an old-projection kill vertex: `rho=1`, `Omega=K31`.
- THM-025 at n=9 has a near-kill vertex: `rho=92/94`, with residue
  `alpha(T-v)=[1,2,1]`.
- Paley T7 and interval T7 have the same regular score sequence and the same
  36 odd-cycle supports, but different support multiplicity, disjoint-pair
  count, and even-graph projection.
- THM-350/351/353 prove exact quotient-transport row checksums, so escaping
  half-line mass is a formally audited projection residue.
- `residue_vector_codex_2026_05_30.py` now emits a compact seed table for the
  representative examples used in this synthesis.

## Test Plan

1. Consolidate the existing exact feature blocks from
   `projection_defect_bridge_s12.py`, `goodcut_transport_excess_s15.py`, and
   `tournament_tda.py`.
2. Sample n=9 and n=10 tournaments with non-real-rooted `I(Omega,x)` and test
   whether they have unusually high `max_v rho_v` or distinctive small
   deletion residues.
3. Compare Paley, interval, near-transitive, H-maximizer, and complete-core
   families on the two coarse axes:

```text
localization = max_v loss_v / |Omega|
mobility     = normalized crossHalf leakage
```

## Interpretation

The common object is not a single invariant.  It is a calculus for finite
projections:

```text
upstairs structure -> projected shadow + residue
```

This makes the residue, not the shadow, the first-class diagnostic object.

## Sources

- `07-reflections/residue-calculus-feedback-loop.md`
- `07-reflections/projection-defect-as-common-residue.md`
- `07-reflections/good-cut-interval-gas.md`
- `07-reflections/staircase-top-bucket-is-strong-connectivity.md`
- `04-computation/projection_defect_bridge_s12.py`
- `04-computation/residue_vector_codex_2026_05_30.py`
- `04-computation/goodcut_transport_excess_s15.py`
- `04-computation/tournament_tda.py`
- `05-knowledge/results/residue_vector_codex_2026_05_30.out`
