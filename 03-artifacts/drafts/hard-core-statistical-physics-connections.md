# Hard-Core Model and Statistical Physics Connections to H(T) = I(Omega(T), 2)

**Author:** kind-pasteur-2026-03-05-S13
**Date:** 2026-03-05
**Purpose:** Web research summary on connections between the OCF identity H(T) = I(Omega(T), 2) and statistical physics / combinatorial theory.

---

## 1. The Hard-Core Lattice Gas at Fugacity lambda=2

### Background

The **hard-core lattice gas model** (or hard-core model) on a graph G = (V, E) with fugacity lambda > 0 assigns to each independent set I of G the weight lambda^|I|. The **partition function** is:

    Z_G(lambda) = sum_{I independent} lambda^|I| = I(G, lambda)

This is precisely the independence polynomial I(G, x) evaluated at x = lambda.

**Our identity H(T) = I(Omega(T), 2) therefore says: the number of Hamiltonian paths in tournament T equals the partition function of the hard-core model on the odd-cycle conflict graph Omega(T) at fugacity lambda = 2.**

### Is lambda = 2 special?

**Short answer: lambda = 2 is ABOVE the uniqueness threshold for essentially all nontrivial graphs, placing us firmly in the non-perturbative / phase-transition regime.**

The **uniqueness threshold** (tree threshold) for the hard-core model on graphs of maximum degree Delta is:

    lambda_c(Delta) = (Delta - 1)^(Delta-1) / (Delta - 2)^Delta

Specific values:
- Delta = 3: lambda_c = 4/1 = 4.0
- Delta = 4: lambda_c = 27/16 = 1.6875
- Delta = 5: lambda_c = 256/243 ~ 1.053
- Delta = 6: lambda_c = 3125/4096 ~ 0.763
- Delta >= 7: lambda_c < 1

**For Delta >= 5, lambda = 2 exceeds lambda_c.** Since Omega(T) typically has maximum degree growing with n, our fugacity lambda = 2 is deep in the non-perturbative regime for all but the smallest tournaments.

This means:
1. Standard perturbative tools (cluster expansion, polymer models) do NOT converge at lambda = 2 for Omega(T).
2. Approximate counting of I(Omega(T), 2) via MCMC is #P-hard in general (Sly, 2010; Sly-Sun, 2012).
3. The fact that H(T) gives the EXACT value of this partition function is remarkable -- it provides an exact combinatorial evaluation of a partition function in a regime where no general efficient algorithm exists.

### Implications

The non-perturbative nature means that the OCF identity is not a consequence of any low-fugacity approximation or cluster expansion. It is an exact algebraic identity that holds at a "difficult" value of the fugacity. This is consistent with the project's earlier finding (kind-pasteur/opus-S5) that lambda = 2 is above all convergence thresholds.

**Reference:** Sly, "Computational Transition at the Uniqueness Threshold," arXiv:1005.5584 (2010).

---

## 2. I(G, 2): Combinatorial Interpretation

The evaluation I(G, 2) = sum_{k >= 0} alpha_k * 2^k has a clean combinatorial interpretation:

**I(G, 2) counts all "2-labeled independent sets" of G.** That is, for each independent set S in G, we assign one of 2 labels (say, "red" or "blue") to each vertex in S. Then I(G, 2) is the total number of such labeled independent sets.

Equivalently: I(G, 2) = |{(S, f) : S independent in G, f: S -> {1,2}}|.

For our setting: **H(T) counts the number of pairs (S, f) where S is a set of vertex-disjoint odd directed cycles in T and f assigns one of two colors to each cycle in S.** The empty collection contributes 1 (the "transitive" base case).

This "2-coloring of independent cycle collections" interpretation was noted in T046 and may be the key to a bijective proof. Each Hamiltonian path of T should correspond to a unique such pair (S, f), but the bijection (if it exists) appears non-trivial and may require a global construction (e.g., Lindstrom-Gessel-Viennot style).

---

## 3. The Lovasz Local Lemma (LLL) Connection

### Scott-Sokal Theorem (2005)

The foundational paper connecting statistical mechanics to the LLL is:

**Scott & Sokal, "The Repulsive Lattice Gas, the Independent-Set Polynomial, and the Lovasz Local Lemma," J. Stat. Phys. 118, 1151-1261 (2005). arXiv:cond-mat/0309352.**

**Main theorem:** The conclusion of the LLL holds for dependency graph G and probabilities {p_x}_{x in V} if and only if the multivariate independence polynomial Z_G(p_1, ..., p_n) is nonvanishing in the polydisc |z_x| <= p_x.

This means: the LLL is equivalent to the nonvanishing of the hard-core partition function in a certain region of the complex plane.

### Implications for our setting

The LLL framework is most useful for establishing that the independence polynomial is NONZERO (positive) at certain values. Since I(Omega(T), 2) = H(T) >= 1 (by Redei's theorem, every tournament has at least one Hamiltonian path), we already know I(Omega(T), 2) > 0. The LLL connection doesn't directly help us compute the value, but it does confirm that the partition function is well-defined and positive.

The Shearer criterion (the optimal version of LLL) gives the exact boundary of the zero-free region: the independence polynomial of a graph of max degree Delta is nonzero for |lambda| < lambda_S(Delta) = Delta^Delta / (Delta+1)^(Delta+1). For lambda = 2, this is violated whenever Delta >= 4 (lambda_S(4) = 256/3125 ~ 0.082... wait, that's for the Shearer zero-free DISC, not the real-positive threshold). On the positive real axis, the situation is better: the Shearer bound gives a zero-free interval [0, lambda_S(Delta)) but lambda = 2 exceeds this for moderate Delta.

**Key insight:** The fact that I(Omega(T), 2) > 0 always is NOT guaranteed by LLL/Shearer bounds for large n. It is guaranteed by the combinatorial identity H(T) >= 1. This means the OCF identity provides information that goes beyond what LLL-type bounds can deduce.

---

## 4. Shearer's Bound

Shearer's bound gives the optimal zero-free disc for the independence polynomial of bounded-degree graphs:

    |lambda| < lambda_S(Delta) = (Delta-1)^(Delta-1) / Delta^Delta ~ 1/(e*Delta)

For lambda = 2 on the positive real axis:
- Delta = 1: lambda_S = 1/1 = 1.0 (but max deg 1 means matching, I(G,2) trivially positive)
- Delta = 2: lambda_S = 1/4 = 0.25
- General: lambda_S(Delta) ~ 1/(e*Delta)

**lambda = 2 is far beyond the Shearer zero-free disc for any nontrivial maximum degree.** This means Shearer's bound cannot be used to establish positivity of I(Omega(T), 2), let alone compute it. The positivity comes from the combinatorial identity, not from any analytic bound.

**Reference:** Shearer, "On a problem of Spencer," Combinatorica 5 (1985), 241-245. See also Scott-Sokal (2005) for the optimal version.

---

## 5. Scott-Sokal Theory: Zeros and the Hard-Core Model

### Zero-free regions

Scott and Sokal (2005) provide the definitive treatment of zero-free regions for the independence polynomial. Their key results:

1. **Dobrushin-Shearer criterion:** I(G, z) is nonzero whenever |z_v| <= (Delta-1)^(Delta-1)/Delta^Delta for all v (where Delta = max degree + 1). This is the optimal bound for the class of all graphs of bounded degree.

2. **Lee-Yang interpretation:** Zero-free activity domains correspond to parameter regimes WITHOUT phase transitions. Phase transitions occur precisely at the boundary of the zero-free region.

3. **Connection to convergence:** The cluster expansion (Mayer expansion) converges in exactly the zero-free region. Outside this region, the expansion diverges.

### What this says about lambda = 2

Since lambda = 2 is outside the zero-free region for moderate-degree graphs, the independence polynomial CAN have zeros near lambda = 2 for general graphs. The fact that I(Omega(T), 2) > 0 always is a SPECIAL PROPERTY of the conflict graph Omega(T), not a general property of independence polynomials at lambda = 2.

This makes the structure of Omega(T) crucial. What is special about conflict graphs of tournaments?

### Claw-free graphs and real roots

**Chudnovsky-Seymour Theorem (2007):** If G is claw-free (no induced K_{1,3}), then I(G, x) has ALL REAL ROOTS (all roots are negative real numbers).

**Reference:** Chudnovsky & Seymour, "The roots of the independence polynomial of a clawfree graph," J. Combin. Theory Ser. B 97 (2007), 350-357.

This extends the Heilmann-Lieb theorem (1972) which proved the same for line graphs (matching polynomials).

**For our setting:** Omega(T) is trivially claw-free for n <= 8 (vertex counting argument: 3 pairwise vertex-disjoint odd cycles + 1 touching all three requires >= 9 vertices). It FAILS at n = 9 (90% of random tournaments have claws in Omega(T); see opus-S7). So the Chudnovsky-Seymour real-roots result applies only for small n.

**When Omega(T) IS claw-free:** All roots of I(Omega(T), x) are real and negative, which implies:
- I(Omega(T), x) > 0 for all x >= 0
- The coefficients alpha_k are log-concave and unimodal
- I(Omega(T), 2) > 0 is guaranteed (but this is already known from H(T) >= 1)

**Jerrum-Patel (2026):** Recent paper extending zero-free regions to H-free graph classes beyond claw-free, including subdivided claws. Published in J. London Math. Soc., 2026.

**Reference:** Jerrum & Patel, "Zero-free regions for the independence polynomial on restricted graph classes," J. London Math. Soc. (2026). DOI: 10.1112/jlms.70458.

---

## 6. Conflict Graphs from Combinatorial Structures

### What is known about conflict graphs of directed cycles?

The **induced odd cycle packing number** iocp(G) is the maximum number of vertex-disjoint odd cycles in G. This is precisely the independence number alpha(Omega(T)) when we restrict to odd cycles.

The conflict graph Omega(T) has special structure:
- **Vertices** are directed odd cycles of T (of lengths 3, 5, 7, ..., up to n if n is odd)
- **Edges** connect cycles sharing at least one vertex
- For small n (n <= 7), Omega(T) is always **perfect** (no odd holes or antiholes of length >= 5)
- At n = 8, Omega(T) can be **imperfect** (C5 subgraphs appear in ~54% of tournaments)
- For n <= 8, Omega(T) is always **claw-free** (trivially, by vertex counting)
- At n = 9, claw-freeness fails (90% of random tournaments)

### Interval graph structure

For the specific tournament T_full_n (full tiling tournament), Omega(T_full_n) is an **interval graph** -- all odd cycles are consecutive intervals, and two cycles conflict iff their intervals overlap. This is the strongest possible structure and explains why I(Omega(T_full_n), 2) satisfies the Tribonacci recurrence.

For general tournaments, Omega(T) is NOT an interval graph. But it may belong to other restricted classes (e.g., chordal in some cases).

### Cycle packing in digraphs

The study of vertex-disjoint cycle packings in digraphs is a well-developed area. Noga Alon proved that every digraph with minimum outdegree >= r contains vertex-disjoint cycles, with explicit bounds. The Erdos-Posa property for directed cycles (and its failure for odd directed cycles) is relevant but has not been directly connected to the independence polynomial framework.

---

## 7. The Grinberg-Stanley Proof and Its Implications

### The proof

The OCF identity H(T) = I(Omega(T), 2) is proved by combining:

1. **Grinberg & Stanley, "The Redei-Berge symmetric function of a directed graph," arXiv:2307.05569 (2023):** For a tournament D on [n], U_D can be written as a polynomial in p_1, 2p_3, 2p_5, ... with nonneg integer coefficients. The power-sum expansion encodes odd cycle structure.

2. **Irving & Omar, "Revisiting The Redei-Berge Symmetric Functions via Matrix Algebra," arXiv:2412.10572 (2024):** Corollary 20 states:

       ham(D-bar) = sum_{sigma in S(D), all cycles odd} 2^{psi(sigma)}

   where psi(sigma) = number of nontrivial cycles of sigma, S(D) = permutations whose nontrivial cycles are all D-cycles.

For tournaments: D-bar = D^op (the complement of a tournament is its converse/opposite), and ham(D^op) = ham(D) by path reversal. The RHS counts pairs (collection of vertex-disjoint odd D-cycles, 2-coloring of cycles) = I(Omega(D), 2).

### Proof technique

The proof uses matrix algebra: adjacency matrices, generating functions for walks, permanent-determinant formulas, character theory, and Jacobi/Sylvester determinant identities. It is entirely different from our combinatorial approach (vertex deletion, arc flips, swap involution).

### What the proof does NOT tell us

- No bijective or combinatorial interpretation of WHY H(T) = I(Omega(T), 2)
- No connection to statistical mechanics or the hard-core model
- No explanation of the special role of lambda = 2
- No structural insight into Omega(T) that would explain positivity

---

## 8. Summary of Key Findings

| Topic | Finding | Relevance to OCF |
|-------|---------|-------------------|
| Hard-core at lambda=2 | Non-perturbative regime; above uniqueness threshold for Delta >= 5 | Cluster expansion fails; OCF is an exact non-perturbative identity |
| I(G,2) interpretation | Counts 2-labeled independent sets (S, f: S -> {1,2}) | Each Ham path corresponds to a colored independent cycle collection |
| LLL/Shearer | Zero-free disc too small to reach lambda=2 for Delta >= 2 | LLL cannot prove I(Omega(T), 2) > 0; positivity comes from H(T) >= 1 |
| Scott-Sokal | Comprehensive theory of zero-free regions = no phase transition | lambda=2 is in the phase transition regime for general graphs |
| Chudnovsky-Seymour | Claw-free => all roots real and negative | Applies to Omega(T) only for n <= 8; fails at n >= 9 |
| Jerrum-Patel 2026 | Extended zero-free regions for H-free classes | May apply to restricted tournament families |
| Conflict graph structure | Omega(T) is perfect for n <= 7, claw-free for n <= 8, interval for T_full | Special structure explains positivity and tractability for small n |
| Grinberg-Stanley proof | Matrix algebra + symmetric functions | Completely different technique; no stat mech connection |

### Most important insight

**The OCF identity H(T) = I(Omega(T), 2) is a rare example of an EXACT evaluation of a hard-core partition function in the non-perturbative regime.** Most results in statistical mechanics focus on approximate counting or zero-free regions. The tournament structure provides enough algebraic structure (via the Redei-Berge symmetric function) to give an exact formula.

This suggests that tournament-theoretic methods could potentially be used to study hard-core models on other specially structured graphs, or conversely, that statistical mechanics techniques beyond the perturbative regime could yield new tournament results.

---

## References

1. **Scott & Sokal (2005).** "The Repulsive Lattice Gas, the Independent-Set Polynomial, and the Lovasz Local Lemma." J. Stat. Phys. 118, 1151-1261. arXiv:cond-mat/0309352.

2. **Grinberg & Stanley (2023).** "The Redei-Berge symmetric function of a directed graph." arXiv:2307.05569.

3. **Irving & Omar (2024).** "Revisiting The Redei-Berge Symmetric Functions via Matrix Algebra." arXiv:2412.10572.

4. **Chudnovsky & Seymour (2007).** "The roots of the independence polynomial of a clawfree graph." J. Combin. Theory Ser. B 97, 350-357.

5. **Sly (2010).** "Computational Transition at the Uniqueness Threshold." arXiv:1005.5584.

6. **Shearer (1985).** "On a problem of Spencer." Combinatorica 5, 241-245.

7. **Heilmann & Lieb (1972).** "Theory of monomer-dimer systems." Comm. Math. Phys. 25, 190-232.

8. **Sokal (2001).** "A Personal List of Unsolved Problems Concerning Lattice Gases and Antiferromagnetic Potts Models." arXiv:cond-mat/0004231.

9. **Jerrum & Patel (2026).** "Zero-free regions for the independence polynomial on restricted graph classes." J. London Math. Soc. DOI: 10.1112/jlms.70458.

10. **Bencs & Buys (2025).** "Optimal zero-free regions for the independence polynomial of bounded degree hypergraphs." Random Structures & Algorithms. DOI: 10.1002/rsa.70018.

11. **Harvey, Srivastava & Vondrak (2018).** "Computing the Independence Polynomial: from the Tree Threshold down to the Roots." arXiv:1608.02282.
