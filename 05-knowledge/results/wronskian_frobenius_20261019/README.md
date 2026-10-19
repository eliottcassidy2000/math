# Weighted Wronskians and scalar Frobenius lifts

2026-10-19. Scoped proofs, exact controls, and an independent review.

## Main completed result

For every prime p and precision m>=1, a polynomial or Laurent polynomial f over F_p has an integer Laurent-polynomial lift A satisfying

    A(X)^p = A(X^p) modulo p^m

if and only if f=0 or all support exponents of f lie in one residue class modulo p^(m-1). A degree-preserving lift is explicit. For nonmonomial f, the p-adic content of its support differences gives its exact last compatible precision.

Read [SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md](SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md), especially Theorem3.1. It also proves the exact late-time p-grid of every scalar additive cellular-automaton kernel. Normalized nonmonomial kernels always have the full precision depth m-1.

## Supporting mechanisms

- [WRONSKIAN_FROBENIUS_TRANSFER.md](WRONSKIAN_FROBENIUS_TRANSFER.md): weighted/nonuniform Hermite rank, Frobenius graph jet losses, sorted determinant valuation, sharp scalar templates, divided-power/Hasse duality, and a precise prime-race projection boundary.
- [FINITE_LOGARITHM_BRIDGE.md](FINITE_LOGARITHM_BRIDGE.md): the same carry specializes to a classical finite logarithm and a Fermat quotient; nontrivial F_p roots are exactly double, with explicit order-three symmetry. These elementary identities retain classical attribution and do not imply a Wieferich-prime finiteness theorem.
- [ROOT_PROOF_REVIEW.md](ROOT_PROOF_REVIEW.md): independent proof review of the Wronskian and scalar-classification statements, plus an import-independent exhaustive lift check.

## Reproduce

All scripts use only Python's standard library. Run from this directory:

    python3 verify_wronskian_frobenius.py > verification.out
    python3 -O verify_wronskian_frobenius.py > verification.optimized.out
    cmp verification.out verification.optimized.out
    python3 verify_template_depth.py > template_depth_verification.out
    python3 -O verify_template_depth.py > template_depth_verification.optimized.out
    cmp template_depth_verification.out template_depth_verification.optimized.out
    python3 independent_lift_check.py > independent_lift_check.json
    python3 -O independent_lift_check.py > independent_lift_check_optimized.json
    cmp independent_lift_check.json independent_lift_check_optimized.json

The first producer replay passes 118,920 exact gates, including 24 nonzero same-minor checks and positive valuation controls. The second passes 118,286 exact gates, including every one of 83,085 coefficient lifts in seven declared finite boxes. The independent path checks all 120,490 coefficient lifts in four declared boxes. These universes overlap; they are verification workloads, not a claim of that many distinct mathematical examples. Normal and optimized outputs agree byte-for-byte on every path.

The finite tests are controls for the proofs; they are not substitutes for the universal arguments. Scalar commutativity, coefficient precision, original rule support, seed content, and quotient aliasing remain explicit. The Rule30 and noncommuting-matrix hostiles block unlicensed extensions. No major-conjecture closure is claimed.

`SHA256SUMS` pins the frozen packet. The producer worked outside the checkout; repository publication is handled separately.
