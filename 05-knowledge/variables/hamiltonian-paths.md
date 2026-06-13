# Variable: H(T) — Hamiltonian Path Count

**Symbol:** `H(T)` or just `H`
**Type:** positive odd integer (Redei's theorem)
**Defined in:** `01-canon/definitions.md`; THM-001

## Definition
H(T) = number of Hamiltonian paths in tournament T on n vertices.

## Values at small n

| n | Min H | Max H | # tournaments | Source |
|---|-------|-------|---------------|--------|
| 1 | 1 | 1 | 1 | trivial |
| 2 | 1 | 1 | 1 | trivial |
| 3 | 1 | 3 | 2 | `04-computation/tournament_lib.py` |
| 4 | 1 | 5 | 4 | exhaustive |
| 5 | 1 | 15 | 12 | exhaustive; OEIS A038375 |
| 6 | 1 | 45 | 56 | exhaustive; OEIS A038375 |
| 7 | 1 | 189 | 456 | exhaustive; max at Paley T_7; OEIS A038375 |
| 8 | 1 | 661 | 6880 | exhaustive; SC maximizer; OEIS A038375 |
| 9 | 1 | 3357 | - | sampled; OEIS A038375 |

## Key equations

- **OCF (THM-002):** H(T) = I(Omega(T), 2)
- **Redei (THM-001):** H(T) is odd for all tournaments T
- **Trace formula:** H(T) = tr(M) for odd n; H = sum of all M[a,b] for even n
- **W-polynomial:** H(T) = W(T, 1/2) * 2^{n-1} [evaluation at r=1/2]
- **Forward polynomial:** H(T) = F(T, 1) = sum_k F_k
- **Complement invariance:** H(T) = H(T^op) [THM-030, from even-r structure of W]
- **Deletion-contraction:** H(T) = H(T\e) + H(T/e) for any edge e [verified, 0 violations]
- **Alpha decomposition:** H(T) = 1 + 2*(t3 + t5 + t7 + ...) + 4*(bc33 + ...) + ... = 1 + 2*sum(alpha_k)

## Computation methods
- **Brute force:** O(n!) via permutation enumeration
- **DP (Held-Karp):** O(n^2 * 2^n) via `ham_path_count_dp()`
- **Transfer matrix:** O(n * 2^n) via diagonal computation

## Relationships
- H(T) = tr(M) links to [transfer-matrix.md](transfer-matrix.md)
- H(T) = I(Omega(T), 2) links to [omega-graph.md](omega-graph.md) and [independence-poly.md](independence-poly.md)
- H decomposes into [fourier-coefficients.md](fourier-coefficients.md) via W-polynomial
- Max H achieved by Paley tournaments (see [cycle-counts.md](cycle-counts.md))

## Tags
#hamiltonian-paths #core #OCF #Redei #trace #W-polynomial #complement-invariance
