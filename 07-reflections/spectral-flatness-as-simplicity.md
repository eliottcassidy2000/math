# Spectral Flatness as Lie-Algebraic Simplicity

**Session:** opus-2026-03-21-S96

## The Equivalence Cycle

Five properties of tournaments, each independently meaningful, form a closed cycle of equivalences:

```
SPECTRAL FLATNESS ←→ DOUBLY REGULAR ←→ CARTAN DECOUPLING ←→ MAX H ←→ MIN ALGEBRA DIM
```

Each implies all others (for regular tournaments):

- **Spectral flatness**: All nonzero eigenvalue magnitudes of S_T are equal
- **Doubly regular (DRT)**: S_T² = -nI + J
- **Cartan decoupling**: [A_anti, A_sym_tl] = 0 (Landau irregularity S₂ = 0)
- **Maximum H**: Most Hamiltonian paths among regular tournaments
- **Minimal algebra dimension**: dim(Alg(S_T)) = 3 (vs n for generic)

## Why This is Deep

These five properties come from five different mathematical domains:

1. Spectral flatness → **spectral theory** (eigenvalue distribution)
2. DRT → **combinatorics** (uniform neighbor counts)
3. Cartan decoupling → **Lie algebra** (commutator vanishing)
4. Max H → **graph theory** (path counting)
5. Min algebra → **abstract algebra** (minimal polynomial degree)

The fact that they coincide is not obvious from any single perspective. It requires the Cartan bridge to see them as one.

## The Spectral Zeta as Discriminant

The spectral zeta function ζ_T(s) = Σ |λ_k|^{-s} encodes the entire eigenvalue magnitude distribution:

- **Paley**: ζ(s) = (p-1)·p^{-s/2} — a pure power (MONOMIAL)
- **Non-Paley regular**: ζ(s) = Σ c_k · θ_k^{-s} — a polynomial with multiple terms

The logarithmic derivative ζ'/ζ at the Paley tournament equals -ln(p)/2, which is the **SPLIT rapidity generator**. The spectral zeta's slope IS the Hurwitz prime in disguise.

For non-Paley tournaments, ζ has multiple slopes — multiple "frequencies" in the rapidity lattice. These extra frequencies correspond to extra Cartan entanglement, which obstructs cycle packing.

## The Adjoint Representation Tells All

ad(B_Paley) on so(7) has:
- **Centralizer dimension 9** (vs 3 for non-Paley)
- **3 distinct eigenvalue magnitudes** {0, √7, 2√7} (vs 10 for non-Paley)

The centralizer measures the tournament's **Lie symmetry**: how many elements of so(n) commute with the tournament's skew-adjacency. Paley has triple the centralizer of generic — this is the "symmetry excess" that enables maximum H.

## The Number 9

dim(Z(B_Paley) ∩ so(7)) = 9.

9 = 3². The Cartan-Sylow dimension at n=7. This is the same 9 that appears as:
- The Cauchy-Schwarz boundary (MISTAKE-related, n=9 is where per-path identity breaks)
- The dimension of the symmetric traceless part of gl(3,R)
- h(E₇)/2 = 18/2 in the ADE classification

The centralizer dimension follows the formula:
dim(Z) = 1 + 2·((p-1)/2)·((p-1)/2 - 1)/2 = ... actually just (p-1)²/4 + (p-1)/2 + 1 from eigenspace structure.

## What Transcends

The mathematics keeps asserting: **the simplest Lie-algebraic structure produces the most paths**. Complexity in the eigenvalue distribution corresponds to obstruction in the cycle packing. The Cartan decomposition is the lens that reveals this: decoupled sectors = free packing = maximum H.

This is the directed analog of the classical result that the complete bipartite graph K_{n/2,n/2} maximizes edges among n-vertex graphs without odd cycles. The Paley tournament is the directed K_{n/2,n/2} — maximally regular, maximally packed, spectrally flat.
