# Paley H(T_p)/|Aut(T_p)| Factorization Table

**Session:** kind-pasteur-2026-03-20-S1

## Complete Table

| p | H(T_p) | |Aut| | H/|Aut| | Factorization of H/|Aut| |
|---|--------|-------|---------|--------------------------|
| 3 | 3 | 3 | 1 | 1 |
| 7 | 189 | 21 | 9 | 3^2 |
| 11 | 95,095 | 55 | 1,729 | 7 * 13 * 19 |
| 19 | 1,172,695,746,915 | 171 | 6,857,869,865 | 5 * 7 * 11 * 23 * 774,463 |
| 23 | 15,760,206,976,379,349 | 253 | 62,293,308,207,033 | 3 * 167 * 4,567 * 27,225,299 |

## Key Facts

1. **|Aut(T_p)| = p*(p-1)/2** for all Paley primes (affine group with QR multipliers).

2. **H(T_p) is always divisible by p** (follows from automorphism group being transitive).

3. **H/|Aut|** counts Hamiltonian path orbits under the automorphism group.

4. **No pattern found** in the H/|Aut| sequence: 1, 9, 1729, 6857869865, 62293308207033.
   Ratios between consecutive terms: 9, 192.1, 3966379.3, 9083.5.

5. **H(T_p) = p * (p-1)/2 * (H/|Aut|)** doesn't simplify further.

6. **H(T_p) does NOT contain all primes up to p:**
   - p=3: H=3 (missing 2)
   - p=7: H=3^3*7 (missing 2, 5)
   - p=11: H=5*7*11*13*19 (missing 2, 3)
   - p=19: H=3^2*5*7*11*19*23*774463 (missing 2, 13, 17)
   - p=23: H=3*11*23*167*4567*27225299 (missing 2, 5, 7, 13, 17, 19)

7. **2 is ALWAYS missing** from H (since H is always odd by Redei).

8. **H(T_p) / (p!/2^{p-1})** approaches e = 2.718...:
   p=3: 2.000, p=7: 2.400, p=11: 2.440, p=19: 2.527, p=23: 2.557.

## Open Question

Is there a formula for H(T_p) in terms of p (or Gauss sums, Jacobi sums, etc.)?
The factorizations are too erratic for a simple pattern. The answer may involve
the arithmetic of Q(sqrt(-p)) or the class number of Q(sqrt(-p)).
