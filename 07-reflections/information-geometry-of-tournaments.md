# Information Geometry of Tournament Space

**Session**: kind-pasteur-2026-03-22-S20y

## The Three Frameworks

Tournament space {0,1}^m (where m = C(n,2)) is simultaneously:
1. A **hypercube** in combinatorics (vertices = tournaments, edges = arc flips)
2. A **statistical manifold** in information geometry (Boltzmann family on H)
3. A **Morse landscape** in topology (sublevel/superlevel set filtration)

These three perspectives are unified by the Walsh-Fourier transform.

## Walsh-Fourier: Only Even Orders

The most surprising finding: H has **exactly zero odd-order Walsh coefficients** at both n=5 and n=6. The energy distribution:

| Order | n=5 (% of Var) | n=6 (% of Var) | Interpretation |
|-------|----------------|-----------------|----------------|
| 0     | -- (mean)      | -- (mean)       | E[H] |
| 2     | 94.7%          | 92.3%           | Pairwise arc interactions |
| 4     | 5.3%           | 7.7%            | Quadruple arc interactions |
| 6     | 0%             | 0%              | Absent |
| odd   | 0%             | 0%              | **Exactly zero** |

**Why odd = 0**: H(T) = H(T^complement) by path reversal. The complement map on {0,1}^m sends x_j -> 1-x_j, which acts as (-1)^{|S|} on Walsh coefficient f_hat(S). For H to be complement-invariant, all odd-order coefficients must vanish.

**Why order 2 dominates**: The OCR (97% of H variation explained by scores) is the information-geometric reflection of order-2 dominance. Scores are the marginal arc biases -- they capture pairwise interactions. The 3% residual (the OCR gap) comes from order-4 interactions.

**Connection to Stadler-Reidys**: In the combinatorial landscapes framework, H is a "quasi-elementary" landscape: the elementary landscape test gives R^2 = 0.977 (n=5) and R^2 = 0.990 (n=6), improving with n. The average neighbor function is NEARLY linear in H. This means gradient ascent is NEARLY optimal -- the landscape has low "ruggedness."

## Fisher Information: Maximum at Uniformity

The Boltzmann family P_beta(T) = exp(beta*H(T))/Z(beta) is a 1-parameter exponential family. Its Fisher information I(beta) = Var_beta[H].

Key finding: **max Fisher is at beta=0** (the uniform distribution). This means:
- The H-landscape has NO sharp phase transition in temperature
- The distribution over H values changes MOST rapidly near uniformity
- This is the opposite of what happens in mean-field models (where max Fisher is at the critical temperature)

The reason: the H-landscape is "globally convex" in the information-geometric sense. There are no deep valleys creating barriers between states. The single-basin property at n=5 (and near-single-basin at n=6) confirms this.

## The Hessian and Morse Theory

The discrete Hessian Hess_H[j,k] = H(x^j^k) - H(x^j) - H(x^k) + H(x) captures the curvature of H at each tournament.

**At n=5 local maxima (H=15)**:
- Signature: (7-, 0, 3+) -- 7 negative eigenvalues
- Trace = -42 (strong negative curvature overall)
- The negative eigenvalues are LARGE (up to -17.7)

**At n=6 global max (H=45)**:
- Signature: (8-, 0, 7+) -- Morse index 8
- Trace = -108
- Both large negative AND large positive eigenvalues

**At n=6 secondary max (H=37)**:
- Signature: (7-, 0, 8+) -- Morse index 7
- **Different parity of Morse index!**
- det < 0 (odd index) vs det > 0 (even index) for H=45

The Morse index parity distinguishes the global from the secondary peak. In Morse-theoretic Euler characteristic:
- chi = sum_p (-1)^{index(p)} over all critical points
- H=45 contributes +1 (even index 8)
- H=37 contributes -1 (odd index 7)

This is consistent with the hypercube having Euler characteristic 0 (2^m vertices, m*2^{m-1} edges, etc. alternate to give 0 for all m >= 1).

## Sublevel Persistence: Trivially Simply Connected

The sublevel filtration {T : H(T) <= h} has exactly ONE connected component from the first tournament added (at H=1) all the way to the last. No persistent homology features at all -- Betti_0 = 1 at every threshold.

This means: **the H-landscape has no topological barriers at any threshold.** You can always walk from any tournament to any other while staying below any common upper bound on H (or equivalently, above any common lower bound).

Combined with the Walsh-Fourier near-elementarity, this gives:
- The landscape is **smooth** (quasi-elementary)
- The landscape is **simply connected** (trivial persistence)
- Gradient ascent works because there are no topological obstructions

## Connection to the Zonotopal Framework

Kolesnik-Sanchez (DCG 2024) showed that the score map pi: {0,1}^m -> Pi_n (the permutohedron) is a zonotopal projection. The fiber pi^{-1}(s) is the set of tournaments with score s.

In our Walsh-Fourier language:
- The **order-0** component is the global mean E[H]
- The **order-2** component is captured by scores (the projection to the permutohedron)
- The **order-4** component is the **within-fiber** variation (tournaments with the same score but different H)

This gives a fiber bundle structure:
```
{0,1}^m  --pi-->  Pi_n
   |                |
   H                E[H|score]  (order-0 + order-2)
   |                |
   L               Var[H|score] (order-4 = OCR gap)
```

The OCR = Var(E[H|score])/Var(H) = (order-2 energy)/(order-2 + order-4 energy) = 92-95%.

## Ollivier-Ricci Curvature

For the bare hypercube, Ollivier-Ricci curvature is kappa = 2/m for every edge (Lin-Lu-Yau). This is POSITIVE, meaning the hypercube is "positively curved" (distances contract under random walk).

The H-weighted version (where we bias the random walk toward high-H neighbors) would give HIGHER curvature near local maxima (where neighbors converge) and LOWER curvature near saddle points (where neighbors diverge). Computing this requires the full Wasserstein distance, which we can implement for small n.

The Eidi-Jost (2019) extension handles directed graphs, so we can also compute Ricci curvature of tournaments THEMSELVES (not just of the flip graph).

## Applications

1. **Tournament optimization as natural gradient**: The IGO framework (Ollivier et al., JMLR 2020) gives a principled algorithm for finding high-H tournaments. Put a Boltzmann distribution over {0,1}^m, update beta via natural gradient. The quasi-elementary landscape means this converges FAST.

2. **Compression of tournament sequences**: The Walsh-Fourier spectrum tells us that order-2 captures 92% of H-variation. Any compression scheme that preserves scores preserves most of the information. This validates the arc-flip compressor from S163.

3. **Ranking stability**: The Hessian eigenvalues at each tournament measure the sensitivity of H to arc flips. Large negative eigenvalues = vulnerable directions. The Morse index = number of "danger" directions. Tournaments with low Morse index are more stable under perturbation.

4. **Network curvature analysis**: The Ollivier-Ricci approach extends to any tournament (not just the flip graph). Positive curvature = community structure. This connects to the arborescence analysis from S20w (trees = positive curvature, cycles = negative curvature).

## Open Questions

1. **Does the even-only Walsh pattern persist for all n?** (Conjectured yes, by complement invariance.)
2. **Is the order-4 fraction monotone in n?** (5.3% at n=5, 7.7% at n=6 -- increasing?)
3. **What is the exact elementary landscape eigenvalue?** The effective eigenvalue lambda = m*(1-b) where b is the regression slope. At n=5: lambda = 10*(1-0.579) = 4.21. At n=6: lambda = 15*(1-0.713) = 4.31. Are these approaching a limit?
4. **Can we prove the landscape is "smooth" for all n?** (R^2 -> 1 as n -> infinity?)
5. **The Hessian parity phenomenon**: H=45 has even Morse index, H=37 has odd. Is there a pattern relating Morse index parity to the SC/score structure?

## Key References

- Stadler & Reidys (2002): Combinatorial Landscapes
- Kolesnik & Sanchez (DCG 2024): Geometry of Random Tournaments
- Ollivier (JMLR 2020): Information-Geometric Optimization
- Beerenwinkel et al. (2007): Geometric epistasis on Boolean hypercube
- Lin, Lu & Yau (2011): Ricci curvature of graphs
- Eidi & Jost (2019): Ollivier-Ricci curvature for directed hypergraphs
- Engstrom (2009): Morse-Fourier on simplicial complexes
- De Visser, Park & Krug (2022): Discrete Morse on fitness landscapes
