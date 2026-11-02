# Session139: incidence, arithmetic scope, and sharp stability

**Checkpoint 1 — independently accepted, 2026-11-02.**

## Current result

For a weighted tournament A, let `V=sum_{i<j} Aij Aji` and
`Delta=H(A)-det(I+A)-abs(det A)`. The optimal coefficient in
`Delta>=K_n V` is exactly, writing `R=floor(n/2)`,

```text
kappa_(2R)   = (4/R)[1-(3/4)^R-(1/4)^R],
kappa_(2R+1) = (4/R)[1-(3/4)^R].
```

This packet proves `K_n=kappa_n` for **every n>=748**. The proof reuses one
wedge budget and retains the binary baseline in a square completion; it does
not scan the missing orders. Positive-variance equality is exactly maximal
ordinal fair pairs, with one singleton in odd order. Inherited results cover
n=3..9; at n=3 every weighted triangle has `Delta=V`.
**Orders 10..747 remain OPEN.** The sufficient threshold 748 is not claimed
optimal, and failure of this certificate at 747 is not a counterexample.

A separate exact theorem identifies **weighted order five** as the first
spectral reconstruction failure for H and Delta. An orthogonally conjugate
rational pair has identical full zeta and all-one rooted resolvent but
`H(B)-H(A)=Delta(B)-Delta(A)=1/3456`. The missing scalar is disjoint
edge/complement-triangle incidence `M23`; adding it restores H at order five.
The binary first-blindness boundary remains six. No RH/GRH claim follows.

## Routes

- [Synthesis and honest frontier](SYNTHESIS.md)
- [Concept board and connection contracts](BOARD.md)
- [Verification and reproducible controls](VERIFICATION.md)
- [Frozen checkpoint scope](CHECKPOINT1.md)
- [Shared-budget proof](stability/PROOF_N748.md) and [independent audit](independent_threshold/AUDIT.md)
- [Weighted zeta boundary](riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md)
- [Historical recovery and corrections](riemann_five_eleven/HISTORICAL_RECOVERY.md)
- [Independent bridge audit](independent/AUDIT.md)
- [Inherited Session138 packet](../session138_stability_and_actual_contexts_20261102/README.md)

The all-moduli/root-cocycle continuation is not promoted by this checkpoint.
The historical recovery's residue-ring statements are algebraic, not ordered
weighted-tournament inequalities. Spectral blindness refutes exact recovery,
not every possible spectral sufficient bound. Unrestricted tournament spectrum,
general strong-quotient scalar zero-rigidity, and LRC(14) remain OPEN.
