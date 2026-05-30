# HYP-1807: Exclusion diagram expansion

**Status:** EXPLORATORY proof-technology hypothesis.

## Statement

When a spectral, Fourier, transfer, or cover kernel counts walks or overlaps
too generously, the right proof object is a finite expansion:

```text
easy kernel count - illegal coincidence diagrams = actual simple object.
```

The obstruction is often not the kernel estimate but the organization of the
exclusion diagrams.

## Evidence

For Fejer/circulant tournaments:

- `J_3` has no repeated-vertex corrections and collapses to additive energy.
- `J_5` is the first place where repeated-vertex correction diagrams appear.

For other repo threads:

- path homology is all paths minus non-allowed face diagrams;
- endpoint transfer rank is support plus collision diagrams;
- Lonely Runner is forbidden interval mass minus endpoint-protection overlaps;
- Caccetta-Haggkvist is rooted walk growth minus return/collision diagrams.

## Prediction

Classifying the `J_5` correction diagrams will either produce a low-degree
autocorrelation formula or identify the first genuinely new invariant beyond
additive energy.  The same diagram grammar should then be portable to path
homology and endpoint-protection graphs.

## See

- `07-reflections/kernel-residue-trick-atlas-2026-05-30.md`
- HYP-1804 and HYP-1805

