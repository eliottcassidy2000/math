# THM-217: Binomial Truncation of g_k

**Status:** VERIFIED (computationally, k=1..10)
**Session:** opus-2026-03-15-S89c
**Depends on:** THM-216

## Statement

For all k ≥ 1, the polynomial g_k(m) has support on at most {C(m,0), C(m,1), C(m,2), C(m,3)} in the binomial basis. Specifically:

$$g_k(m) = g_k(0) \cdot \binom{m-1}{2} + m + 2(k-1) \cdot \binom{m}{2} + 2a_k \cdot \binom{m}{3}$$

where:
- g_k(0) is the polynomial extrapolation to m=0 (nonzero for k ≥ 4)
- a_k is the leading coefficient of 3·g_k(m) = a_k·m³ + b_k·m² + c_k·m + d_k
- The identity 1 - m + C(m,2) = C(m-1,2) is used

## Consequences

1. **Degree bound**: g_k(m) has degree ≤ 3 for all k (degree exactly 3 for k ≥ 3)
2. **Two free parameters**: The entire polynomial is determined by (g_k(0), a_k)
3. **Universal boundary**: g_k(1) = 1 and g_k(2) = 2k hold automatically
4. **Binomial cancellation**: In the naive formula g_k^naive(m) = Σ C(k-1,r-1)·C(m,r)·2^{r-1}, the C(m,r) coefficient for r ≥ 4 is C(k-1,r-1)·2^{r-1}. The true coefficient is EXACTLY ZERO. The correction kills these terms completely.

## Binomial basis decomposition

| k | Δ⁰ = g_k(0) | Δ¹ | Δ² | Δ³ = 2a_k |
|---|---|---|---|---|
| 1 | 0 | 1 | 0 | 0 |
| 2 | 0 | 1 | 2 | 0 |
| 3 | 0 | 1 | 4 | 4 |
| 4 | -8 | 9 | -2 | 20 |
| 5 | -592 | 593 | -584 | 776 |
| 6 | -114320 | 114321 | -114310 | 139320 |
| 7 | -33338240 | 33338241 | -33338228 | 39652540 |
| 8 | -12475185560 | 12475185561 | -12475185546 | 14619453484 |
| 9 | -5629549881808 | 5629549881809 | -5629549881792 | 6525375440480 |
| 10 | -2973062116837472 | 2973062116837473 | -2973062116837454 | 3415543797850416 |

Pattern: Δ¹ = 1 - Δ⁰, Δ² = Δ⁰ + 2(k-1)

## Relationship to cluster correlations

The naive formula assumes independent cluster weights: each tiling with r clusters of sizes s_1,...,s_r gets weight ∏ h_{s_i}(n) where h_s(n) = 2/(n)_{2s}. This gives C(m,r) coefficient = C(k-1,r-1)·2^{r-1}.

The true weight includes correlation corrections between clusters sharing a permutation. These corrections EXACTLY cancel all r ≥ 4 contributions. The mechanism is NOT that individual 4-cluster expectations are zero (they're nonzero), but that the sum over all placements cancels.

Key correlation facts (verified for small n):
- Two separated clusters (gap ≥ 1): ratio = n(n-1)/((n-2)(n-3)), independent of gap
- Adjacent clusters (gap = 0): different ratio n(n-1)/(2(n-2)(n-3))
- Gap independence holds for all cluster size combinations tested

## Open questions

1. Prove the truncation analytically (why do r ≥ 4 terms cancel?)
2. Find recurrence or GF for a_k: 2, 10, 388, 69660, 19826270, 7309726742, 3262687720240, 1707771898925208
3. Find recurrence or GF for g_k(0): 0, -8, -592, -114320, -33338240, -12475185560, -5629549881808, -2973062116837472
4. Both sequences are new to OEIS
5. Ratio a_k/a_{k-1} grows roughly linearly: 5, 38.8, 179.5, 284.6, 368.7, 446.3, 523.4
6. Ratio g_k(0)/a_k slowly approaches ~-2

## Evidence

- Exact polynomial match for n=3..21 (gk_verify_89c.py)
- Binomial basis computed via forward differences (gk_corrections_89c.py)
- 4-cluster cancellation verified via brute-force expectation (gk_cumulant_89c.py)
