# Variable: F(T,x) — Forward-edge (ascent) polynomial

**Symbol:** `F(T,x)` or `F_k(T)`
**Type:** polynomial in x of degree n-1 with nonneg integer coefficients
**Defined in:** opus-S37 (Worpitzky deep dive)

## Definition
F(T,x) = sum_{k=0}^{n-1} F_k x^k

where F_k = #{Hamiltonian paths P with exactly k ascents}.
An "ascent" at position i means P_i < P_{i+1} in the natural vertex labeling.

## Key properties

| Property | Status | Source |
|----------|--------|--------|
| F(T,1) = H(T) | trivial | definition |
| F_k(T) = F_{n-1-k}(T^op) | PROVED | path reversal, worpitzky_F_at_2.py |
| Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)} | PROVED | each perm is HP of 2^extra tournaments |
| F is always unimodal | CONJECTURED | 100% through n=8, HYP-204 |
| F is almost always log-concave | CONJECTURED | 1020/1024 at n=5, 100% at n=6-8, HYP-205 |
| Roots ~89% real | OBSERVED | stable rate n=5-8, HYP-206 |
| F is NOT palindromic in general | PROVED | 0/64 at n=4 |
| F(-1) is always odd | TRIVIAL | F(-1) ≡ H (mod 2) |

## Worpitzky expansion
F(T,x) = sum_k w_k C(x+k, n-1)

where w_k are always integers.
- w_{n-1} = F_0 = #{fully decreasing HPs}
- w_{n-2} = H - n*F_0
- w_0 = (-1)^{n-1} F(-1)

## Generating function
sum_{m>=0} F(T,m) x^m = W^*(T,x) / (1-x)^n

where W^*(x) = sum_k w_{n-1-k} x^k is the reversed Worpitzky polynomial.

## Relationship to W-polynomial (THM-K)
W(T,r) = (r - 1/2)^{n-1} * F(T, (2r+1)/(2r-1))

## Relationship to Eulerian polynomial
Average F(T,x) = A_n(x) / 2^{n-1}.
Decomposition: F(T,x) = (H/n!) * A_n(x) + R(T,x) where R(T,1) = 0.
R is NOT palindromic in general.

## Special cases
- Transitive tournament: F = [0,...,0,1], H = 1
- Max-H at n=5: F = [1,8,6,0,0], H = 15
- Regular tournament (n=5): F = [1,10,4,0,0], H = 15

## Scripts
- worpitzky_connection.py, worpitzky_deep.py (initial exploration)
- worpitzky_restricted_eulerian.py (Sum_T result, residual analysis)
- worpitzky_ehrhart.py (Ehrhart connection, Worpitzky structure)
- worpitzky_F_at_2.py (complement duality, GF)
- worpitzky_roots.py (real-rootedness, log-concavity, unimodality)
- worpitzky_exceptions.py (exceptions analysis, parity)
- worpitzky_n6_test.py (sampling at n=6,7,8)

## Tags
#F-polynomial #ascent #Worpitzky #unimodality #log-concavity
