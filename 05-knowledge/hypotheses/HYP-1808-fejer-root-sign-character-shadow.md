---
id: HYP-1808
status: EXPLORATORY
source: codex-2026-05-30 S359
related:
  - HYP-1795
  - HYP-1801
  - HYP-1809
---

# HYP-1808: Fejer Kernel Is A Root-Sign Character Shadow

## Statement

For prime circulant tournaments, the Fejer kernel is the Fourier character
shadow of the all-one cyclic root-sign chamber.

Let `p` be odd prime, `m=(p-1)/2`, and encode a circulant tournament on
`Z/pZ` by signs

```text
sigma_d in {-1,+1},  d=1,...,m,
```

where `sigma_d=+1` chooses the step `d` and `sigma_d=-1` chooses `-d`.
Then for each nonzero Fourier character `k`,

```text
lambda_k = -1/2 + i * sum_{d=1}^m sigma_d sin(2*pi*k*d/p),
```

and therefore

```text
|lambda_k|^2 = 1/4 + (sum_d sigma_d sin(2*pi*k*d/p))^2.
```

When `sigma=(1,...,1)`, this becomes

```text
|lambda_k|^2 = sin^2(pi*m*k/p) / sin^2(pi*k/p),
```

the Fejer kernel.

## Evidence

`04-computation/fejer_root_sign_phase_probe_s359.py` verifies the identity for
all circulant signs through `p=13` and samples signs for `p=17,19,23,29,31`.
The maximum numerical error in the root-sign eigenvalue formula is below
`4e-12`, and the interval/Fejer identity is below `5e-13` in the tested range.

The same probe shows that spectral concentration maximizers for
`p=7,11,13,17` are exactly interval-unit-orbit representatives in the tested
full samples.

## Predictions

1. The circulant `phase_profile` should be implemented as a character profile
   of the root-sign vector, not merely as raw eigenvalues.
2. Spectral-concentration maximizers should be the unit orbit of interval
   chambers.
3. A character-resolved OCF should separate interval-like additive chamber
   packing from Paley-like multiplicative flatness.

## Sources

- `04-computation/fejer_root_sign_phase_probe_s359.py`
- `05-knowledge/results/fejer_root_sign_phase_probe_s359.out`
- `07-reflections/fejer-root-sign-phase-synthesis.md`
