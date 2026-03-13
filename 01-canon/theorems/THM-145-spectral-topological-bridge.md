# THM-145: Spectral-Topological Bridge (Conjecture)

**Status:** CONJECTURE (verified p=7, p=11, p=13 exhaustive)
**Author:** opus-2026-03-12-S69
**Dependencies:** THM-125 (constant symbol matrix), HYP-562 (integer ESF)

## Statement

For a circulant tournament C_p^S on prime p vertices, define:
- Q_k = |S_hat(k)|² (Fourier magnitude squared, k=1,...,m where m=(p-1)/2)
- e_j = e_j(Q_1,...,Q_m) = j-th elementary symmetric function of Q_k
- Ω_m = per-eigenspace dimension of the m-th omega space (GLMY path complex)

**Claim:** The tuple (e_2, e_3, ..., e_m) UNIQUELY DETERMINES the entire Ω profile
(Ω_0, Ω_1, ..., Ω_{p-1}), and hence all topological invariants (Betti numbers,
Euler characteristic).

Note: e_0 = 1, e_1 = m(p-m)/2 are universal (same for all S).

## Verified Data

### p = 7 (exhaustive, 8 orientations, 2 orbits):
| (e_2, e_3) | Ω profile | chi_per |
|-------------|-----------|---------|
| (5, 1) | [1,3,6,11,14,9,2] | 0 |
| (12, 8) | [1,3,6,9,9,6,3] | 1 |

6 orientations map to (5,1) and 2 to (12,8). Each gives unique Ω.

### p = 11 (exhaustive, 32 orientations, 4 orbits):
| (e_2, e_3, e_4, e_5) | Ω profile | chi_per |
|-----------------------|-----------|---------|
| (35,28,9,1) | [1,5,20,74,224,522,926,1222,1115,611,148] | 0 |
| (68,127,97,23) | [1,5,20,70,200,439,711,827,648,301,64] | 2 |
| (79,171,130,23) | [1,5,20,70,201,430,620,596,384,151,26] | 0 |
| (90,270,405,243) | [1,5,20,70,205,460,700,690,450,180,30] | 1 |

Each e_j pattern uniquely determines Ω.

### p = 13 (exhaustive, 64 orientations, 6 orbits):
| (e_2, e_3, e_4, e_5, e_6) | Ω profile | chi_per |
|---------------------------|-----------|---------|
| (70,84,45,11,1) | [1,6,30,140,556,1823,4983,11241,19730,24780,20344,9642,1989] | 1 |
| (122,305,357,180,27) | [1,6,30,135,517,1622,4129,8436,13292,15078,11316,4988,977] | -3 |
| (148,448,604,362,79) | [1,6,30,135,520,1627,3915,6965,9012,8325,5157,1906,321] | -8 |
| (161,552,812,375,53) | [1,6,30,135,524,1675,4176,7707,10295,9681,6081,2298,399] | 4 |
| (161,552,851,570,131) | [1,6,30,135,524,1663,4001,6876,8409,7266,4237,1493,238] | 1 |
| (174,721,1566,1701,729) | [1,6,30,135,528,1707,4245,7503,9177,7770,4374,1488,237] | -17 |

6 distinct e_j patterns → 6 distinct Ω profiles. Note e_2=161 appears twice but
e_4 differs (812 vs 851), confirming higher ESFs are needed.

## Consequences

1. **Q-spectrum encodes topology**: The Fourier magnitudes of S completely
   determine the GLMY path homology of C_p^S.

2. **Spectral-topological duality**: Flat Q (small Σ Q_k²) → large |chi_per|
   (more homology). Peaked Q (large Σ Q_k²) → small |chi_per|.

3. **Algebraic bridge**: Since e_j are integers (HYP-562), the map
   Z^{m-1} → Z^p : (e_2,...,e_m) ↦ (Ω_0,...,Ω_{p-1}) is well-defined.
   Finding this map explicitly would give a FORMULA for Betti numbers
   in terms of Fourier data.

4. **Potential proof strategy**: If we can express Ω_m as a polynomial
   (or rational function) of e_j, combined with the Morgan-Voyce
   structure for Interval (e_j = C(m+j,2j)), we could prove
   chi_per(Interval) = 0 at p ≡ 3 mod 4 algebraically.

## Open Questions

1. What is the explicit map (e_2,...,e_m) → (Ω_3,...,Ω_{p-1})?
   Note: NOT a linear function of e_j (checked at p=7).

2. Does the map extend to ALL primes? (Verified at p=7,11,13; need p=17.)

3. Is the map a polynomial? A rational function? A piecewise formula?

## Related
- THM-125: Constant symbol matrix (enables per-eigenspace decomposition)
- THM-144: Eigenspace boundary rank anomaly
- HYP-562: Integer ESF
- HYP-579: Spectral-topological duality
