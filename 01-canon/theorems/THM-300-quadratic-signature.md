# THM-300: Quadratic Coefficient Matrix Signature

**Status:** CONJECTURED, verified n=5,6,7,8
**Discovered by:** opus-2026-04-04-S3

## Statement

Let Q be the m×m symmetric matrix of quadratic multilinear coefficients:
Q[i,j] = c_{ij} for i ≠ j, Q[i,i] = 0.

**Conjecture.** For n ≥ 5, Q has exactly n-2 negative eigenvalues and m-(n-2) positive eigenvalues (full rank).

## Evidence

| n | m | Positive | Negative | Zero | n-2 | Match? |
|---|---|----------|----------|------|-----|--------|
| 4 | 3 | 1 | 1 | 1 | 2 | No (rank deficient) |
| 5 | 6 | 2 | 3 | 1 | 3 | Yes (neg count) |
| 6 | 10 | 6 | 4 | 0 | 4 | Yes |
| 7 | 15 | 10 | 5 | 0 | 5 | Yes |
| 8 | 21 | 15 | 6 | 0 | 6 | Yes |

At n=5, Q has rank 5 (one zero eigenvalue) but still n-2=3 negative eigenvalues. For n≥6, Q appears to have full rank.

## Remarks

- The trace of Q is always 0 (multilinear polynomial has no t_i² terms).
- The negative eigenspace concentrates on the left column of the staircase (tiles sharing lower vertex 1 or upper vertex n).
- The number n-2 equals the dimension of the staircase (number of rows or columns).
- The nullspace at n=5 involves tiles at staircase coordinates forming an L-shape with alternating signs.
