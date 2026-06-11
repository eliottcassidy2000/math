# HYP-2435 - Naive Construction A is the wrong 72 bridge

**Status:** OPEN warning / route-selection hypothesis.
**Source:** codex-2026-06-11-P4.
**Companions:** HYP-2415, HYP-2425, HYP-2430, HYP-2434, OPEN-Q-066,
T784.

## Statement

The existence of an extremal 72-dimensional even-unimodular lattice should not
be read as direct evidence for the binary `[72,36,16]` code by naive binary
Construction A. That construction retains frame/root data and does not preserve
the extremal support obstruction in the way the code problem needs.

The promising bridge is instead a support/polarization question:

```text
Can the 72-dimensional lattice-side support be decomposed, framed, or polarized
so that the binary Type II support constraints reappear?
```

If not, the failure mode should itself be a proof-facing obstruction.

## Why This Matters

The scalar theta gate and scalar Gleason gate both work at level 72. The lattice
support object exists; the binary support object is open. The missing theorem is
therefore not a modular-form identity but a compatibility theorem between
support categories.

## Next Tests

1. Search for frame/polarization invariants in known 72-dimensional extremal
   lattice constructions that would force or forbid a binary Type II support.
2. Translate the `[72,36,16]` weight-16 `5-(72,16,78)` layer into a matroid/Tutte
   leakage test.
3. Build a route Tournament Analysis whose vertices are support bridges:
   binary matroids, skew-Hadamard gauges, lattice polarizations, Z4/code lifts,
   and design incidence systems.
