---
id: HYP-1801
status: EXPLORATORY
source: opus-2026-05-30-S5
related:
  - HYP-1795
  - HYP-1798
  - T257
  - THM-136
---

# HYP-1801: Circulant Maximizer Transition Is Phase-Dominant

## Statement

The Paley/interval circulant H-maximizer transition is controlled primarily by
phase-channel data, not by localized deletion-residue rank.

For circulant tournaments `C_p(S)`, the decisive features should be Fourier
character phases, trace signs, additive-energy profiles, and cycle-packing
fiber multiplicities.  Deletion residue remains useful as a sanity check but
should not be the primary transition coordinate.

## Evidence

The contrast table shows Paley and interval `T7` are both regular, have equal
three-cycle count, and have broad rank-2 deletion residues.  Their separation
is instead:

```text
Paley:    alpha=(80,7,0),  support_excess=44
Interval: alpha=(59,14,0), support_excess=23
```

This matches earlier trace/eigenphase threads: Paley wins at small primes via
flat spectral distribution and high total cycle production; interval later
wins by different packing concentration.

## Predictions

1. A circulant `phase_profile` should predict the Paley/interval crossover
   better than max-loss deletion residue.
2. Trace-sign and additive-energy features should locate the transition before
   full Hamiltonian-path computation at larger primes.
3. The same phase profile should be useful for circulant path-homology symbol
   matrices.

## Test Plan

1. Implement `phase_profile(C_p(S))`: eigenvalue arguments/magnitudes, low
   trace signs, additive energy, and cycle-count deviations.
2. Compare against `H` for Paley, interval, and sampled connection sets.
3. Reuse the same data in the Paley/circulant path-homology program.

## Sources

- `04-computation/residue_phase_incidence_contrast_s5.py`
- `05-knowledge/results/residue_phase_incidence_contrast_s5.out`
- `07-reflections/applied-residue-phase-incidence-programs.md`
