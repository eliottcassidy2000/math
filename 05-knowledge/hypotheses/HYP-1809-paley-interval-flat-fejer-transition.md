---
id: HYP-1809
status: EXPLORATORY
source: codex-2026-05-30 S359
related:
  - HYP-1801
  - HYP-1808
---

# HYP-1809: Paley-Interval Crossover Is Flat-Versus-Fejer Phase Competition

## Statement

The Paley/interval circulant crossover is a competition between two structured
root-sign phase regimes:

```text
Paley    = multiplicative-character flatness
Interval = additive-chamber Fejer concentration
```

Fejer concentration is not itself a monotone predictor of `H` at small primes.
It becomes useful only after the Paley small-prime advantage stops dominating.

## Evidence

`fejer_root_sign_phase_probe_s359.py` computes `H` for all prime circulants
through `p=13` and compares it with phase features.

The sign of the correlation between `H` and spectral concentration changes:

```text
p=7:  corr(H,top_fraction)=-1.000000, corr(H,ipr)=-1.000000
p=11: corr(H,top_fraction)=-0.106581, corr(H,ipr)=-0.087860
p=13: corr(H,top_fraction)=+0.444395, corr(H,ipr)=+0.509255
```

At `p=7,11`, Paley wins `H` with a flat spectrum.  At `p=13`, the
interval-unit orbit wins with Fejer-level concentration.

Literal Fejer alignment is not orbit-invariant: unit rotations of an interval
can have the same `H` and same spectral concentration but much lower alignment
with the fixed interval representative.

## Predictions

1. A correct circulant phase profile must record orbit-invariant concentration
   features: sorted spectrum, top fraction, IPR, entropy, and peak orbit.
2. Fixed-coordinate features such as low-frequency pair mass or literal Fejer
   alignment will misclassify unit-rotated interval chambers.
3. A character-resolved OCF should show Paley winning early through total cycle
   mass, while interval wins later through disjoint packet packing.

## Sources

- `04-computation/fejer_root_sign_phase_probe_s359.py`
- `05-knowledge/results/fejer_root_sign_phase_probe_s359.out`
- `07-reflections/fejer-root-sign-phase-synthesis.md`
