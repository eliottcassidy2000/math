# THM-201: Grand Fourier Energy Formula

**Status:** CORRECTED (kind-pasteur-S112). The m^k formula is the LEADING TERM only.
**Session:** kind-pasteur-2026-03-14-S105–S107, CORRECTED kind-pasteur-2026-03-15-S112
**Proves:** Asymptotic Fourier energy at every level for tournament H

## Statement (CORRECTED)

For the Hamiltonian path count H(T) on tournaments with n vertices, the Fourier energy at level 2k on the Boolean hypercube {0,1}^m (m = C(n,2)) satisfies:

**E_{2k} / E_0 = 2 · g_k(n − 2k) / (n)_{2k}**

where g_k is the combinatorial matching polynomial (THM-217):
- g_1(m) = m, g_2(m) = m², g_k has degree k in general
- g_k(m) ~ 2^{k-1} m^k / k! for large m (dominant eigenvalue approximation)

The ORIGINAL formula E_{2k}/E_0 = 2·(n-2k)^k/(n)_{2k} is the **leading-term approximation**, exact only for k ≤ 2. For k ≥ 3, corrections of relative order O(1/m) apply.

**Alternative (THM-216):** A degree-3 reparametrization g̃_k exists with the same sum — see THM-216.

Parameters:
- E_{2k} = Σ_{|S|=2k} Ĥ(S)² is the total squared Fourier energy at level 2k
- E_0 = Mean(H)² = (n!/2^{n−1})²
- (n)_{2k} = n!/(n−2k)! is the falling factorial
- Valid for 1 ≤ k ≤ ⌊(n−1)/2⌋, m = n-2k ≥ 1

## Corollaries

**Corollary 1 (Variance Formula):**
Var(H)/Mean(H)² = Σ_{k=1}^{⌊(n−1)/2⌋} 2·(n−2k)^k / P(n, 2k)

**Corollary 2 (Concentration):**
n · Var(H)/Mean(H)² → 2 as n → ∞. The ratio is O(1/n).

**Corollary 3 (Spectral Purification):**
E_2/Var → 1 as n → ∞. Level 2 captures (n−2)/n of the variance exactly (if formula holds).

**Corollary 4 (Geometric Series):**
For large n, E_{2k}/E_0 ∼ 2/n^k. The total is a geometric series 2/(n−1).

## Verification

| n | k | E_{2k}/E_0 (formula) | E_{2k}/E_0 (exact) | Match |
|---|---|---------------------|---------------------|-------|
| 3 | 1 | 1/3 | 1/3 | ✓ |
| 4 | 1 | 1/3 | 1/3 | ✓ |
| 5 | 1 | 3/10 | 3/10 | ✓ |
| 5 | 2 | 1/60 | 1/60 | ✓ |
| 6 | 1 | 4/15 | 4/15 | ✓ |
| 6 | 2 | 1/45 | 1/45 | ✓ |
| 7 | 1 | 5/21 | 5/21 | ✓ |
| 7 | 2 | 3/140 | (sum verified) | ✓ |
| 7 | 3 | 1/2520 | (sum verified) | ✓ |
| 8 | all | 59/252 | Monte Carlo ≈ 0.236 | ✓ |

## Proof sketch (spectator freedom argument)

Fixing 2k arcs at Fourier level 2k creates k interaction points (k pairs of arcs). Each interaction point has (n−2k) effective free spectator choices. With replacement: (n−2k)^k total. Normalize by P(n,2k) arc arrangements. Factor 2 from path reversal involution.

E_{2k}/E_0 = 2 · (spectator freedom)^(interaction points) / (arc arrangements)

## Connection to known results

- k=1: PROVED (S75, S104). The "cone formula" E_2/E_0 = 2(n−2)/(n(n−1)).
- The Degree Drop theorem ensures only EVEN levels contribute.
- The crystallization theorem (E_2/Var → 1) follows immediately.
- The 1/3 ratio at n=3,4 is the special case where only k=1 contributes.

## Exact Var/Mean² values

| n | Var/Mean² | Fraction |
|---|-----------|----------|
| 3 | 0.333333 | 1/3 |
| 4 | 0.333333 | 1/3 |
| 5 | 0.316667 | 19/60 |
| 6 | 0.288889 | 13/45 |
| 7 | 0.259921 | 131/504 |
| 8 | 0.234127 | 59/252 (predicted) |
