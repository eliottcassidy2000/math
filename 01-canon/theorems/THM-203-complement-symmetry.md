# THM-203: Complement Symmetry — H(T̄) = H(T)

**Status:** VERIFIED (exhaustive n ≤ 7)
**Source:** opus-2026-03-14-S89

## Statement

For every tournament T on n vertices, H(T̄) = H(T), where T̄ is the complement tournament (all arcs reversed).

## Proof

If v₁ → v₂ → ⋯ → vₙ is a Hamiltonian path in T, then vₙ → vₙ₋₁ → ⋯ → v₁ is a Hamiltonian path in T̄. This reversal bijection preserves the count.

## Consequence: Odd Walsh Levels Vanish

In the Walsh-Hadamard expansion H(T) = Σ_S ĉ_S χ_S(T):

- ĉ_S for T̄ picks up factor (-1)^|S|
- Since H(T̄) = H(T): ĉ_S((-1)^|S| - 1) = 0
- For |S| odd: ĉ_S = 0

**Therefore E_{2k+1} = 0 for all k ≥ 0.**

Only even Fourier levels contribute to Var(H).
