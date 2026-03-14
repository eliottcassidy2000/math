# THM-204: Grand Fourier Level Formula — E_{2k}/E₀ = 2(n-2k)^k / P(n,2k)

**Status:** VERIFIED (exhaustive n ≤ 7, all levels)
**Source:** opus-2026-03-14-S89

## Statement

For the Walsh-Hadamard decomposition of the Hamiltonian path count H(T), the energy at level 2k satisfies:

**E_{2k}/E₀ = 2(n - 2k)^k / P(n, 2k)**

where:
- E_{2k} = Σ_{|S|=2k} ĉ_S² (sum of squared Walsh coefficients at level 2k)
- E₀ = Mean(H)² = (n!/2^{n-1})²
- P(n, 2k) = n(n-1)(n-2)⋯(n-2k+1) = n!/(n-2k)! (falling factorial)
- Valid for k = 1, 2, …, ⌊(n-1)/2⌋

## Special Cases

| k | Formula | Name |
|---|---------|------|
| 1 | E₂/E₀ = 2(n-2)/(n(n-1)) = (n-2)/m | Level-2: adjacent pairs |
| 2 | E₄/E₀ = 2(n-4)²/(n(n-1)(n-2)(n-3)) | Level-4: paths + double-P₃ |
| 3 | E₆/E₀ = 2(n-6)³/P(n,6) | Level-6: first at n=7 |
| ⌊(n-1)/2⌋ | E_{n-1}/E₀ = 2/n! (odd n) | Highest: Hamiltonian paths of K_n |

## Closed-Form Variance Ratio

**Var(H)/Mean(H)² = Σ_{k=1}^{⌊(n-1)/2⌋} 2(n-2k)^k / P(n,2k)**

| n | Var/Mean² | Exact |
|---|-----------|-------|
| 3 | 1/3 | verified |
| 4 | 1/3 | verified |
| 5 | 19/60 | verified |
| 6 | 13/45 | verified |
| 7 | 131/504 | verified |
| 8 | 59/252 | **predicted** |
| 9 | 38443/181440 | **predicted** |

## Verification

Exhaustively verified via Fast Walsh-Hadamard Transform for all 2^{C(n,2)} tournaments at n = 3, 4, 5, 6, 7.

## Key Structural Facts

1. **Odd levels vanish** (THM-203: complement symmetry H(T̄) = H(T))
2. **Level 2**: nonzero iff edge pair is adjacent; magnitude = Mean/m
3. **Level n-1** (odd n): nonzero = n!/2 Hamiltonian paths; |ĉ| = 1/2^{n-2}
4. **Level 4**: Type A (Ham path on K₅ subgraph) + Type B (double P₃); |ĉ_B| = 2|ĉ_A|
5. **Asymptotic**: Var/Mean² ~ 2/n as n → ∞ (E₂ dominates)
