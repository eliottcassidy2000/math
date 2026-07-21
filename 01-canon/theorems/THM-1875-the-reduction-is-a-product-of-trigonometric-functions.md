---
id: THM-1875
title: "THE REDUCTION IS A PRODUCT OF TRIGONOMETRIC FUNCTIONS: a unifying trigonometric lens on the repo's reduction principles. (a) SPECTRAL REDUCTION: char_A is multiplicative under order-join (block-triangular), so a tournament's spectrum is the disjoint union of its strong components' spectra — the spectral companion of THM-1862/THM-1830. (b) The strong CIRCULANT atoms have eigenvalues that are explicit TRIGONOMETRIC sums (Gauss sums for Paley, Dirichlet/Chebyshev kernels for interval connection sets). (c) The signed tournament partition function Σ_T(−1)^{back}x^{score}=∏_{i<j}(x_j−x_i) evaluated on the unit circle x_k=e^{iθ_k} is a PRODUCT OF SINES ∏2i·sin((θ_j−θ_i)/2), so mac-mini's transitive-core involution IS a trigonometric factorization; at roots of unity it is the cyclotomic discriminant. This is the tournament mirror of the LRC covering sum ∏_j sinc(k_jδ) (Chebyshev-equioscillation / Fejér-certificate side). CLAIMED — boxeph-S194, IN PROGRESS."
status: >
  CLAIMED / IN PROGRESS (boxeph-2026-07-21-S194). Reserving the number. Synthesis + verification:
  char_A join-multiplicativity (block-triangular, extends THM-1862 to the spectrum), circulant-atom
  eigenvalues as trig closed forms, the sine-product evaluation of the signed partition function,
  and the bridge to the LRC product-of-sincs / Chebyshev-equioscillation side. To be filled with
  exhaustive small-n verification.
source: boxeph-2026-07-21-S194 (owner: another archeology session; pursue more reduction principles; think trigonometric)
depends_on: []
related:
  - THM-1862  # my order-join reduction principle (this is its spectral/trigonometric companion)
  - THM-1830  # reducible => char_A = product of char(SCC) (block-triangular); SCCs are the atoms
  - THM-456   # blow-up spectrum law (cycle-length spectra) — a different blow-up reduction
  - "07-reflections/the-covering-min-is-a-chebyshev-equioscillation-and-why-greedy-has-no-shortcut.md"
  - "07-reflections/the-cyclotomic-magic-function-is-the-fejer-kernel-kps.md"
  - "07-reflections/the-sign-reversing-tournament-involution-as-a-repo-wide-engine-macmini-S159.md"
script: 04-computation/trig_reduction_tournament_boxeph_S194.py (+ .out)
---

# THM-1875 — the reduction is a product of trigonometric functions (STUB, in progress)

Placeholder claimed by boxeph-S194. See script + reflection for the developing content.
