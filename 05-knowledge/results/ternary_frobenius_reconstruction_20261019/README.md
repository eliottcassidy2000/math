# Ternary Frobenius reconstruction and integral lift obstructions

Research packet dated 2026-10-19, with independent proof and coefficient review.

## Main results

- Complete classification of odd-degree cyclic two-variable-supported norm solutions over every field of characteristic p >= 5: dimension floor(log_p n).
- Complete supported kernel over Z/p^m for every m >= 2: it is p^(m-1) times the mod-p kernel. This gives an exact prime-to-six Smith-torsion description.
- Explicit all-interior norm-zero families of degrees 3p+3 for every p >= 5 and 2p+12 for every p >= 11. These give universal upper bounds on failure of coordinate-plane reconstruction.
- Independent complete minimum calculations for p = 2,3,5,7,11,13,17,19,23,29,31.
- A characteristic-29 exception: the minimum is degree 42, disproving equality with the two-family bound 70. Its 23-alternant certificate has no cyclic interior lift modulo 29 squared, and the corresponding 29-primary cokernel torsion is Z/29.

The characteristic-zero source is Ilin–Ono, *Ternary polynomial relations and periods of cusp forms*, https://arxiv.org/abs/2610.19309. That paper already proves positive-characteristic reconstruction failure; this packet gives the separately stated classifications, lower-degree constructions, arithmetic obstruction, and finite exceptional-prime computation. It makes no external-priority claim.

## Reading order

1. `TERNARY_FROBENIUS_RECONSTRUCTION.md`: definitions, supported classification, universal families, exact census, and limits.
2. `CYCLIC_NORM_LIFT_OBSTRUCTIONS.md`: all-precision supported theorem, Smith corollary, finite full-cyclic lift tests, and finite-logarithm bridge.
3. `EXCEPTIONAL_PRIME_29.md`: explicit 23-alternant exception and its interior lift obstruction.
4. `review/PERIOD_PROOF_REVIEW.md`: independent proof review, exact target hashes, and reproduction instructions.

## Verification

From this directory, the standard-library checks are:

    python3 -B verify_periods_independent.py
    python3 -B verify_lift_independent.py
    python3 -B verify_exception_29.py

The independent review directory contains six additional verifiers using SymPy and NumPy. Normal and Python -O controls are retained. Production code is kept separate from the review implementations. `SHA256SUMS` records the published payload's file hashes, excluding itself.

## Open boundaries

The two construction formulas are upper bounds, not a general minimum formula. The all-prime classification of earlier exceptional kernels remains open. Full cyclic lifting outside the supported lattice is established only for the explicitly tested cases. No cusp-form representation for the new rational quadratic family, characteristic-zero counterexample, Jacobian-conjecture result, or cellular-automaton theorem is claimed.
