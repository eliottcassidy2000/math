---
id: HYP-1795
status: EXPLORATORY
source: opus-2026-05-30-S4
related:
  - HYP-1779
  - HYP-1780
  - HYP-1785
  - THM-025
  - THM-094
  - THM-136
---

# HYP-1795: Residue-Phase Bifurcation

## Statement

Tournament anomalies split into two main mechanisms:

```text
residue anomalies = a projection leaves a structured survivor;
phase anomalies   = an orthogonal/eigen/character channel survives or changes sign.
```

The hypothesis is that future feature extraction and proof searches should
separate these into two vectors:

```text
residue_vector(T) = deletion loss, SCC defect, support multiplicity,
                    Omega alpha residue, quotient gaps, fiber defects;

phase_vector(T)   = complement parity, Walsh/Krawtchouk energy by level,
                    trace signs, circulant eigenphase data, mod-p Taylor zeros.
```

Real-root failures and H-spectrum gaps should be enriched by residue-rank
features.  Paley/interval and circulant maximizer questions should be driven
primarily by phase-channel dominance.  Mixed examples require both.

## Evidence

- H=63 single-core classes are exact deletion kills: the residue is empty.
- THM-025 is a near-kill: a tiny deletion residue survives with rank 2.
- Paley and interval T7 share broad support shadows but differ in fiber
  multiplicity/disjointness, so they are residue-related but not near-kills.
- THM-094 is phase-like: tournament-dependent channels vanish mod 2, giving
  `F(T,x) = (1+x)^(n-1) mod 2`.
- THM-136 is phase-like: Paley and interval trace differences are controlled
  by eigenvalue phases around `pi/2` and a `k mod 4` sign pattern.
- Super-orthogonality is phase-like: complement symmetry, Walsh parity, OCF
  amplitudes, and homological balances reinforce each other through
  orthogonal decompositions.

## Predictions

1. Non-real-rooted `I(Omega(T),x)` examples should be easier to classify by
   residue-rank features than by trace or Walsh phase features alone.
2. Circulant H-maximizer transitions should be easier to classify by phase
   features than by deletion-residue features alone.
3. Examples that are hard for either single feature family should be mixed:
   a small structured residue plus a bad phase/root channel.
4. A combined `residue_vector + phase_vector` should separate the standard
   contrast set better than either vector alone:

```text
transitive, Paley, interval, H=63, THM-025, n=6 H=37 trap.
```

## Test Plan

1. Add phase-level summaries to `tournament_tda.py`: complement parity flags,
   low Krawtchouk/Walsh energy where feasible, trace sign profile, and
   circulant eigenphase summaries for circulants.
2. Join these with the existing residue-rank features.
3. Benchmark against real-root failures, H-maximizer families, local H traps,
   and Paley/interval circulant classes.
4. Use ablations: residue-only, phase-only, and combined classifiers.

## Sources

- `07-reflections/residue-phase-incidence-synthesis.md`
- `07-reflections/residue-calculus-theses.md`
- `07-reflections/super-orthogonality.md`
- THM-025, THM-094, THM-136, HYP-1785.
