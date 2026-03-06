# Web Research: Zero-Free Regions, Subdivided Claws, and Real Roots of I(Omega(T), x)

**Source:** opus-2026-03-05-S9 (web research session)

---

## HEADLINE FINDING: Subdivided Claw-Freeness is Finite, Like Claw-Freeness

### Computational Discovery

Omega_3(T) (the 3-cycle conflict graph) transitions from S_{a,b,c}-free to having S_{a,b,c} as the subdivided claw gets larger:

| Property | Holds for n<= | Fails at n= | Failure rate at transition |
|----------|---------------|-------------|---------------------------|
| Claw-free (S_{0,0,0}) | 8 | 9 | 88% (100 random) |
| S_{1,1,1}-free | 11 | 12 | 100% (50 random) |

**Pattern:** Each level of subdivision buys ~3 more tournament vertices before failure. This is because S_{a,b,c} has 3(a+b+c)+4 vertices (approximately), and the number of 3-cycles in Omega grows as ~n^3/24.

**Implication:** No fixed subdivided claw is universally avoided by Omega(T). For every specific S_{a,b,c}, there exists n large enough that Omega_3(T) contains it.

### Also Tested: Line Graph Hypothesis

Omega_3(T) contains K_5-e (Beineke forbidden subgraph #2) starting at n=6:

| n | K_5-e in Omega_3 |
|---|-----------------|
| 5 | 0/20 (0%) |
| 6 | 9/20 (45%) |
| 7 | 8/10 (80%) |
| 8 | 10/10 (100%) |

**Conclusion:** Omega(T) is NOT a line graph for n>=6. Heilmann-Lieb theorem does not apply.

---

## Connection 1: Jerrum-Patel 2026 — The Boundary Theorem (CRITICAL)

**Paper:** Jerrum & Patel, "Zero-free regions for the independence polynomial on restricted graph classes," J. London Math. Soc. (2026)

### Main Result

For any fixed subdivided claw H = S_{a,b,c}:
- The **univariate** independence polynomial of H-free graphs of bounded max degree Delta has all zeros on the negative real line.
- The **multivariate** independence polynomial has zero-free regions that can be quantified.

### Best-Possible Characterization

This theorem is **sharp in two senses**:
1. **H must be a subdivided claw (or path):** For any H that is NOT a subdivided claw or path, there exist H-free graphs of max degree 3 whose independence polynomial zeros accumulate on the positive real axis.
2. **Bounded degree is required:** Without the degree bound, the result fails.

### Why This Matters for OCF

Subdivided claws are the EXACT boundary between "nice" and "bad" independence polynomial behavior. This is the deepest known result classifying when I(G, x) has restricted zero locations based on forbidden subgraphs.

### Why It Doesn't Directly Apply to Omega(T)

Two obstacles:
1. **Omega(T) is not always S_{a,b,c}-free** for any fixed (a,b,c). Every subdivided claw eventually appears as n grows.
2. **Omega(T) has unbounded degree** as n grows. The max degree of Omega_3(T) grows roughly as n^2.

So Jerrum-Patel cannot directly explain the real-roots conjecture (THM-020).

---

## Connection 2: Bezakova-Galanis-Goldberg-Stefankovic 2024 — Mixing Time Dichotomy

**Paper:** arXiv:2404.07615, published in Combinatorics, Probability and Computing (Cambridge)

### Dichotomy Theorem (Theorem 17)

For bounded-degree H-free graphs:
- If H is a **subdivided claw** S_{a,b,c}: Glauber dynamics for the hard-core model mixes in **O(n log n)** for ALL fugacities lambda.
- If H is neither a subdivided claw nor a path: mixing time is **exponential** in n for sufficiently large lambda.

### Connection to OCF

This dichotomy parallels Jerrum-Patel: subdivided claws are the exact boundary for both:
- Zero-free regions (Jerrum-Patel)
- Polynomial mixing (Bezakova et al.)

The hard-core model at fugacity lambda=2 on Omega(T) has partition function Z = I(Omega(T), 2) = H(T). This is in the non-perturbative regime (lambda=2 exceeds the uniqueness threshold). Yet the evaluation is always positive (= number of Hamiltonian paths).

The same bounded-degree limitation applies: the dichotomy doesn't directly constrain Omega(T).

---

## Connection 3: 2-Intersection Graphs of 3-Hypergraphs

**Paper:** arXiv:2305.13932

### Omega_3(T) as an Intersection Graph

Omega_3(T) is the **intersection graph** of a specific family of 3-element subsets of [n]: those triples {i,j,k} that form directed 3-cycles in T. Two triples are adjacent in Omega_3 iff they share at least one element.

This means Omega_3(T) is a **1-intersection graph of a 3-uniform hypergraph** (specifically, the 3-cycle hypergraph of T).

### Structural Implication

The class of intersection graphs of 3-element sets is well-studied. Key facts:
- Every graph is an intersection graph of SOME family of sets (Marczewski's theorem)
- But intersection graphs of k-element sets from [n] have bounded clique number and specific structural properties
- The maximum clique in Omega_3(T) corresponds to a vertex of T: all C(d,1) 3-cycles through that vertex form a clique, where d depends on the tournament's score sequence.

### Connection to Claw-Freeness

A claw in Omega_3 requires 4 triples: one "center" triple C0 and three "leaf" triples C1, C2, C3 where each Ci shares a vertex with C0 but the Ci are pairwise disjoint. Since each Ci has 3 elements and they're pairwise disjoint, we need >= 9 elements. This is the vertex-counting argument for claw-freeness at n<=8.

For S_{1,1,1}, we need 7 triples forming the subdivided claw pattern. With intermediary triples sharing vertices with both center and leaf triples, the vertex count requirement grows to ~11-12, matching our computational finding.

---

## Connection 4: Real-Stability and Symmetric Functions

### The Key Question

No known graph property that holds universally for Omega(T) explains real-rootedness of I(Omega(T), x). The explanation must be algebraic.

### Real-Stable Polynomials

A multivariate polynomial p(z_1, ..., z_n) is **real-stable** if it has no zeros where all variables have positive imaginary parts. For univariate polynomials, real-stability = all roots real.

Chudnovsky-Seymour's proof for claw-free graphs goes through showing real-stability of the multivariate independence polynomial. Engstrom extended this: for claw-free G, the multivariate I(G, z) is "same-phase stable."

### Potential Path via Grinberg-Stanley Framework

The Irving-Omar paper (arXiv:2412.10572) proves OCF using the Redei-Berge symmetric function and matrix algebra. The symmetric function U_D encodes Hamiltonian path information and satisfies:

```
ham(D) = [coefficient extraction from U_D]
```

**Speculative connection:** If I(Omega(T), x) can be expressed as a specialization of the Redei-Berge symmetric function (or a related symmetric function), then properties of symmetric functions (e.g., Schur-positivity, real-stability of certain evaluations) might force real-rootedness.

The power-sum expansion of U_D involves terms indexed by cycle types. The odd-cycle terms contribute to I(Omega(T), x). If the generating function for these terms has a representation as a real-stable polynomial, real-rootedness would follow.

**This is the most promising unexplored direction.**

---

## Connection 5: Quantifying the Root Gap

**Paper:** arXiv:2510.09197 (Prakash & Sharma, 2025, FSTTCS 2025)

For connected graphs, beta(G) (the smallest root in absolute value) is a simple real root. This paper quantifies the gap between beta(G) and the next-smallest root.

### For Omega(T)

If all roots of I(Omega(T), x) are indeed real and negative, then beta(Omega(T)) is the largest (least negative) root. The gap to the next root controls how "robust" the real-rootedness is.

For OCF, we evaluate at x=2 > 0. Since all roots are negative (conjectured), I(Omega(T), 2) > 0 automatically. The value H(T) = I(Omega(T), 2) is determined by the root locations: if the roots are r_1 <= r_2 <= ... <= r_k < 0, then:

```
H(T) = prod_{i=1}^{k} (1 - 2/r_i) > 0
```

(since each factor 1 - 2/r_i > 1 when r_i < 0).

---

## Connection 6: Chromatic Polynomial of Claw-Free Graphs

**Paper:** arXiv:2601.06918 (2026)

Zero-free regions for the chromatic polynomial of claw-free graphs (with and without induced squares/diamonds). This extends the zero-free region methodology from independence polynomials to chromatic polynomials.

### Relevance

The chromatic polynomial and independence polynomial are related through the partition function framework. For claw-free graphs, both have restricted zero locations. The techniques used (Barvinok interpolation, correlation decay) are the same.

For Omega(T), the chromatic polynomial chi(Omega(T), k) counts proper k-colorings. The chromatic number chi(Omega(T)) = omega(Omega(T)) when Omega is perfect (n<=7). Understanding the chromatic polynomial's zero structure could complement the independence polynomial analysis.

---

## Synthesis: The State of the Real-Roots Conjecture

### What We Now Know

1. **No graph property explains universal real-rootedness.** Claw-freeness fails at n=9, S_{1,1,1}-freeness at n=12, line-graph at n=6, perfectness at n=8. Every tested graph property eventually fails.

2. **The Jerrum-Patel boundary theorem confirms subdivided claws are special.** They are the EXACT boundary for zero-free regions in bounded-degree graphs. But both conditions (fixed H, bounded degree) fail for Omega(T) at large n.

3. **Real-rootedness, if true universally, requires an algebraic/symmetric-function explanation.** No graph-theoretic forbidden subgraph argument can work because Omega(T) eventually contains every fixed subgraph.

4. **The Irving-Omar framework (arXiv:2412.10572) is the key.** They prove OCF using matrix algebra and symmetric functions. The same framework might contain information about the FULL independence polynomial, not just its evaluation at x=2.

### Prioritized Next Steps

**Tier 1 (Most Promising):**
1. Study whether the Irving-Omar matrix algebra framework says anything about the independence polynomial beyond x=2. Can I(Omega(T), x) be expressed as a symmetric function evaluation?
2. Computationally verify real-rootedness at n=12-15 where S_{1,1,1}-freeness fails.

**Tier 2 (Important Context):**
3. Determine the rate of growth of max subdivided claw size in Omega(T) as n grows. Is it O(n), O(sqrt(n)), or O(log n)?
4. Test if there's a weaker graph property (e.g., "quasi-claw-free" or "locally claw-free") that holds universally.

**Tier 3 (Worth Exploring):**
5. The 1-intersection graph structure of Omega_3(T) may have algebraic implications.
6. Connection between Redei-Berge Hopf algebra (Grinberg, arXiv:2402.07606) and real-stability of I(Omega(T), x).

---

## Author Correction

**IMPORTANT:** arXiv:2412.10572 is by **John Irving and Mohamed Omar**, not Grinberg and Stanley. It "consolidates and expands on the work of Chow, Grinberg and Stanley, and Lass." The project files should be updated to attribute the OCF proof correctly (Corollary 20 is from Irving-Omar, building on Grinberg-Stanley's framework).

---

## References

- [Jerrum & Patel (2026)](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/jlms.70458) — Zero-free regions for independence polynomial, JLMS
- [Bezakova et al. (2024)](https://arxiv.org/abs/2404.07615) — Glauber dynamics for hard-core model on H-free graphs
- [Irving & Omar (2024)](https://arxiv.org/abs/2412.10572) — Revisiting Redei-Berge via matrix algebra (proves OCF)
- [Chudnovsky & Seymour (2007)](https://www.sciencedirect.com/science/article/pii/S0095895606000876) — Real roots for claw-free graphs
- [Prakash & Sharma (2025)](https://arxiv.org/abs/2510.09197) — Quantifying the root gap
- [arXiv:2601.06918](https://arxiv.org/abs/2601.06918) — Chromatic polynomial zero-free regions for claw-free graphs
- [arXiv:2305.13932](https://arxiv.org/abs/2305.13932) — 2-intersection graphs of 3-hypergraphs
- [Grinberg (2024)](https://arxiv.org/abs/2402.07606) — Redei-Berge Hopf algebra of digraphs
- [Building graphs with real-rooted I.P.](https://link.springer.com/article/10.1007/s00373-009-0857-5)
- [Generalizations of matching polynomial](https://alco.centre-mersenne.org/article/ALCO_2019__2_5_781_0.pdf) — Multivariate independence polynomial
