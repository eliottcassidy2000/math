# THM-219: NUD-Poisson-Euler Connection

**Status:** VERIFIED (computationally, n=3..27) + HEURISTIC PROOF
**Session:** opus-2026-03-15-S89c
**Depends on:** THM-216, THM-217

## Statement

Let NUD(n) = number of permutations of [n] with no unit descent (σ(i) - σ(i+1) = 1 never occurs), and let W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)} where adj1(σ) counts unit ascents (positions where σ(i+1) = σ(i) + 1). Then:

### Part A: NUD identification
NUD(n) = A000255(n-1) with recurrence NUD(n) = (n-1)·NUD(n-1) + (n-2)·NUD(n-2), NUD(1)=NUD(2)=1.
EGF: Σ NUD(n) x^n/n! = exp(-x)/(1-x)².

### Part B: Poisson limit
In a uniformly random NUD permutation of [n], the number of unit ascents adj1(σ) converges in distribution to Poisson(1) as n → ∞.

### Part C: W/NUD ratio
W(n)/NUD(n) → e as n → ∞, where e = Euler's number.

### Part D: CV² asymptotics
CV²(n) = W(n)/n! - 1 ~ 2/n as n → ∞. More precisely, CV²(n) · n → 2.

## Proof sketch

**Part A:** Well-known. NUD perms avoid consecutive pattern 21. EGF is standard (Goulden-Jackson/Elizalde-Noy framework for consecutive pattern avoidance).

**Part B (heuristic):** The indicators I_i = 1_{σ(i+1)=σ(i)+1} for i=1,...,n-1 each have marginal probability ≈ 1/n in a random NUD perm. There are n-1 indicators, so E[adj1] ≈ 1 + 1/n → 1. The indicators are weakly dependent (each depends only on local values), so by the Chen-Stein Poisson approximation theorem, adj1 → Poisson(1).

Computational verification: the distribution of adj1 converges to Poisson(1) for n=8,10,12 with ratios P(adj1=k)/P_Poisson(k) → 1.

**Part C:** Follows from Part B: E[2^X] where X ~ Poisson(1) = Σ 2^k e^{-1}/k! = e^{-1} · e^2 = e.

**Part D:** CV²(n) = W(n)/n! - 1. Since NUD(n)/n! → 1/e (from EGF: exp(-x)/(1-x)² at x=1 doesn't converge, but NUD(n)/n! ~ 1/e by standard asymptotics of A000255), and W(n)/NUD(n) → e, we get W(n)/n! → 1. The correction: NUD(n)/n! ≈ 1/e + c₁/n + c₂/n², and the 2^adj1 weighting adds a correction. The leading behavior CV² ~ 2/n follows from the contribution of g_1(n-2) = n-2 term: 2(n-2)/n(n-1) = 2/n - 6/n² + ...

## Data

| n | W(n)/NUD(n) | CV²(n) | CV²·n |
|---|---|---|---|
| 5 | 2.981 | 0.317 | 1.583 |
| 10 | 2.947 | 0.193 | 1.927 |
| 15 | 2.883 | 0.131 | 1.972 |
| 20 | 2.846 | 0.099 | 1.985 |
| 25 | 2.822 | 0.080 | 1.991 |

### Part E: Hertzsprung connection
Define N(n,j) = #{σ ∈ NUD(n) : adj1(σ) = j}. Then:
- Row sums: Σ_j N(n,j) = NUD(n) = A000255(n-1)
- Column j=0: N(n,0) = A002464(n) (Hertzsprung numbers — perms where |σ(i+1)-σ(i)| ≠ 1 for all i)
- W(n) = Σ_j 2^j N(n,j)

The Hertzsprung column satisfies:
N(n,0) = (n+1)·N(n-1,0) - (n-2)·N(n-2,0) - (n-5)·N(n-3,0) + (n-3)·N(n-4,0)

This recurrence does NOT extend to j > 0; the bivariate recurrence for N(n,j) is unknown.

## Related sequences

- NUD(n): A000255 (offset by 1)
- N(n,0): A002464 (Hertzsprung numbers)
- W(n): NEW — not in OEIS as of 2026-03-15
- N(n,j) triangle: NEW — not in OEIS as of 2026-03-15
- W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904, 46963358, ... for n=1,2,3,...

## Open questions

1. Prove Part B rigorously (Chen-Stein bounds on the indicator coupling)
2. Find the EGF of W(n) explicitly (bivariate EGF marking adj1 in NUD perms)
3. Determine the next correction: CV²(n) = 2/n + c/n^α + ... (what is c and α?)
4. Submit W(n) to OEIS
5. Submit N(n,j) triangle to OEIS
6. Find bivariate recurrence for N(n,j)
