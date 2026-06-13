# THM-292: Exact Formula for Σ H — Total Hamiltonian Path Count Over All Tilings

**Status:** PROVED (algebraic, verified n=3,...,7 exhaustively, n=8,9,10 via formula)
**Filed by:** opus-2026-04-04-S5
**Depends on:** THM-284, OCF, tiling model

## Statement

For a tournament on n vertices with fixed base path, the total Hamiltonian path count
summed over all 2^m tilings (m = C(n-1,2)) is:

  **Σ_t H(T(t)) = Σ_{σ succession-free} 2^{m - (n-1-bp(σ))}**

where:
- The sum is over all succession-free permutations σ of {0,...,n-1}
  (OEIS A000255: permutations with no i such that σ(i+1) = σ(i)+1)
- bp(σ) = number of base-path arcs used (descents by exactly 1)
- m = C(n-1, 2)

Equivalently: Σ H = Σ_{bp=0}^{n-1} N(n, bp) · 2^{C(n-1,2)-n+1+bp}

where N(n, bp) = number of succession-free permutations of [n] with exactly bp
descents-by-1.

## Proof

Each labeled permutation σ = (v₁,...,v_n) is a valid Hamiltonian path in EXACTLY those
tournaments where every arc v_i → v_{i+1} is present. Three cases:
1. **Base-path arc** (v_i = v_{i+1}+1): always present. No tile constraint.
2. **Tile arc** (|v_i - v_{i+1}| ≥ 2): constrains one tile to a specific direction.
3. **Reverse base-path arc** (v_{i+1} = v_i+1): NEVER present in any tiling.

If σ has any reverse base-path arc, it's valid in ZERO tilings.
Otherwise, σ has bp base-path arcs and (n-1-bp) tile arcs, each constraining a
DISTINCT tile. The remaining m - (n-1-bp) tiles are free (2 choices each).

So: |{tilings where σ is HP}| = 2^{m-(n-1-bp)} if σ is succession-free, 0 otherwise.

Summing: Σ_t H(T(t)) = Σ_{succession-free σ} 2^{m-n+1+bp(σ)}. ∎

## Computed Values

| n | m | Σ H | E[H] | A000255(n) |
|---|---|-----|------|------------|
| 3 | 1 | 4 | 2 | 3 |
| 4 | 3 | 32 | 4 | 11 |
| 5 | 6 | 632 | 9.875 | 53 |
| 6 | 10 | 29,696 | 29 | 309 |
| 7 | 15 | 3,251,200 | 99.21875 | 2,119 |
| 8 | 21 | 815,136,768 | 388.6875 | 16,687 |
| 9 | 28 | 461,027,409,920 | 1717.46 | 148,329 |
| 10 | 36 | 580,881,441,882,112 | 8452.9375 | 1,468,457 |

## New OEIS Sequences

**Sequence 1: Σ H(n)** = Total Hamiltonian path count over all tilings at n.
4, 32, 632, 29696, 3251200, 815136768, 461027409920, 580881441882112
(NOT found in OEIS as of 2026-04-04)

**Sequence 2: 2^m · E[H]** = Same as Σ H (integer by construction).

**Known sequence:** The succession-free permutation count IS OEIS A000255.

## The BP Distribution N(n, bp)

| n\bp | 0 | 1 | 2 | 3 | 4 | 5 | 6 | Total |
|------|---|---|---|---|---|---|---|-------|
| 3 | 0 | 2 | 1 | | | | | 3 |
| 4 | 2 | 5 | 3 | 1 | | | | 11 |
| 5 | 14 | 20 | 14 | 4 | 1 | | | 53 |
| 6 | 90 | 115 | 72 | 26 | 5 | 1 | | 309 |
| 7 | 646 | 790 | 467 | 168 | 41 | 6 | 1 | 2119 |

This is a JOINT distribution of successions-free permutations by their descent-by-1 count.
It may correspond to a triangle in OEIS (e.g., A046739 or related).

## See Also
- THM-284 (linear coefficient = 2^(skip-1))
- THM-291 (Mode B recursion)
- OEIS A000255 (succession-free permutations)
- Scripts: sum_H_formula.py, sequences_exact.py
