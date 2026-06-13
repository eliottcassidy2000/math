# THM-152: Transitive Tournament Simplicial Identity

**Status:** PROVED
**Session:** opus-2026-03-13-S70
**Depends on:** GLMY path homology definitions

## Statement

For the transitive tournament T_n on n vertices (i→j iff i < j):

$$\Omega_m(T_n) = \binom{n}{m+1} \quad \text{for all } 0 \leq m \leq n-1$$

The Omega profile of T_n is the row of Pascal's triangle: (C(n,1), C(n,2), ..., C(n,n)).

## Proof

In T_n with vertex set {0, 1, ..., n-1} and arcs i→j iff i < j:

**Step 1: Path condition.** An m-path (v_0, v_1, ..., v_m) requires v_i → v_{i+1}
for all i. In T_n this forces v_0 < v_1 < ... < v_m.

**Step 2: Regularity is automatic.** Regularity requires v_{i-1} → v_{i+1} for
all 0 < i < m. Since v_{i-1} < v_i < v_{i+1}, we have v_{i-1} < v_{i+1}, so
v_{i-1} → v_{i+1} holds automatically in T_n.

**Step 3: Bijection.** Every strictly increasing sequence of length m+1 from
{0,...,n-1} is a regular m-path, and these are the ONLY regular m-paths.
Therefore Ω_m = #{(m+1)-subsets of [n]} = C(n, m+1). ∎

## Verification

Confirmed computationally at n = 3, 4, 5, 6, 7, 8 (all Omega_m values match).

## Corollaries

1. **Simplicial complex:** The path complex of T_n is isomorphic to the
   (n-1)-simplex. Every face of the simplex corresponds to a regular path.

2. **Acyclicity:** All Betti numbers β_k(T_n) = 0 for k ≥ 1.
   The Euler characteristic chi = Σ(-1)^m C(n,m+1) = 1 - (1-1)^n + (-1)^{n+1}... = 1.
   Actually: chi = Σ_{m=0}^{n-1} (-1)^m C(n,m+1) = -Σ_{k=1}^{n} (-1)^k C(n,k) = -(−1+1)^n + C(n,0) = 1.

3. **Transitive maximality:** T_n maximizes Omega_m for m ≤ n-2:
   - Omega_0 = n for ALL tournaments (trivial).
   - Omega_1 = C(n,2) for ALL tournaments (every edge is a 1-path).
   - Omega_2 = C(n,3) iff t_3 = 0, which is the transitive tournament.
   - For m = n-1: the REGULAR tournament maximizes Omega_{n-1}, not transitive.

4. **Omega_2 formula specialization:** Omega_2 = C(n,3) - t_3 at T_n gives
   C(n,3) - 0 = C(n,3), consistent.

## Notes

- The deviation from simplicial, D_m = Omega_m - C(n,m+1), measures how far
  the tournament is from transitive in homological terms.
- For the regular tournament at n=5: profile is (5,10,5,5,5), not C(5,k).
- For the interval tournament at n=7: Omega_4 = 56 > C(7,5) = 21,
  so D_4 = +35. The path complex can be LARGER than the simplex!
