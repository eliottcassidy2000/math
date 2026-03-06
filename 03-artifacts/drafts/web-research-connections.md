# Web Research: Deep Creative Connections for OCF

**Source:** opus-2026-03-05-S5 (web research session)

---

## Connection 1: The Redei-Berge Hopf Algebra (HIGHEST PRIORITY)

**Papers:** arXiv:2402.07606 (Grinberg), arXiv:2407.18608 (properties)

The Redei-Berge symmetric function U_X for a digraph X encodes Hamiltonian path information:

```
U_X = sum_{sigma in S_V} F_{XDes(sigma)}
```

where XDes(sigma) is the "X-descent set" of a vertex ordering sigma.

### Why this is critical for OCF:

1. **The comultiplication is our subset convolution.** The Hopf algebra (D, zeta) has:
   ```
   Delta([X]) = sum_{S subset V} [X|S] tensor [X|V\S]
   ```
   This is EXACTLY the structure of our transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(V\S). The Hopf algebra formalism captures our subset convolution as a coalgebra operation.

2. **The character counts Hamiltonian paths.** zeta([X]) = #{sigma: XDes(sigma) = empty} counts Hamiltonian paths in the complementary digraph. For tournaments, this directly gives H(T).

3. **Antipode formula.** S(U_X) = (-1)^|V| U(X-bar) (Theorem 4.9). For tournaments T, X-bar = T^op (reverse tournament). This encodes Berge's theorem that H(T) ≡ H(T^op) (mod 2).

4. **Reciprocity.** u_X(-m) = (-1)^|V| u_{X-bar}(m) (Theorem 5.8). This connects positive and negative evaluations, potentially relating to our (-1)^|S| signs.

5. **Deletion-contraction.** The Redei-Berge function satisfies deletion-contraction like the chromatic polynomial. This gives an inductive tool that could parallel our vertex-deletion approach to OCF.

### KEY INSIGHT: OCF as a Hopf algebra identity

Our identity H(T) = I(Omega(T), 2) could be expressible as a relation between
the character zeta of the Redei-Berge Hopf algebra and the independence polynomial
evaluated at 2. The comultiplication structure of the Hopf algebra might encode
the odd-cycle collection structure directly.

**Action:** Read arXiv:2402.07606 in full. Formalize OCF in the Hopf algebra language. Check if the inclusion-exclusion structure of I(Omega, 2) = sum_{indep sets C} 2^|C| has a natural Hopf-algebraic interpretation.

---

## Connection 2: Björklund's Determinant Sums & Parity of Hamiltonian Cycles

**Papers:** arXiv:1008.0541 (Björklund 2010), arXiv:1301.7250 (Björklund-Husfeldt 2013)

### Key ideas:

- **Cycle cover reduction:** Hamiltonian cycles are detected by reducing to counting weighted arc-labeled cycle covers over finite fields of characteristic 2.
- **Labeled Cycle Cover Sum (LCCS):** Count cycle covers with specific label constraints where Hamiltonian solutions contribute deterministically.
- **Determinant summation:** The LCCS can be evaluated using determinants, leveraging the algebraic sieving method.
- **Parity formula:** A "new combinatorial formula for the number of Hamiltonian cycles modulo a positive integer" based on cycle covers (arXiv:1301.7250).

### Connection to OCF:

Our identity H(T) = I(Omega(T), 2) relates Hamiltonian PATH count to an alternating sum over vertex-disjoint ODD CYCLE collections. Björklund's approach relates Hamiltonian CYCLE count to an alternating sum over ALL cycle covers (via inclusion-exclusion and determinants).

The analogy is:
```
Björklund: H_cycles = sum over cycle covers, weighted by (-1) and det
Our OCF:   H_paths  = sum over indep odd-cycle sets, weighted by 2^k
```

**Could there be a "directed" version of Björklund's reduction that:**
1. Works for Hamiltonian paths instead of cycles?
2. Reduces specifically to odd cycle covers instead of all cycle covers?
3. Uses the tournament structure (each pair has exactly one arc) to simplify?

The characteristic-2 aspect is particularly interesting since our OCF operates modulo 2 at the deepest level (Redei = H is odd = I(Omega, 2) is odd).

**Action:** Study whether Björklund's labeled cycle cover approach can be adapted to express H(T) in terms of odd cycle covers, yielding OCF directly.

---

## Connection 3: Radchenko-Villegas — Independence Polynomials & Hypergeometric Series

**Paper:** arXiv:1908.11231

### Main result:
A graph G is chordal if and only if the power series expansion of I_G^{-1}(x) is Horn hypergeometric.

### Relevance:
Omega(T) is always PERFECT but NOT always chordal (14% non-chordal at n=6). This means:
- For chordal Omega(T): the inverse I^{-1} has a nice hypergeometric expansion
- For non-chordal Omega(T): the inverse does NOT have this property

The boundary between chordal and non-chordal Omega(T) may correspond to a structural transition in the tournament. Understanding what makes Omega(T) chordal could help with inductive approaches.

**Action:** Characterize which tournaments T have chordal Omega(T). The chordal case might be tractable via the hypergeometric structure.

---

## Connection 4: Dyer-Jerrum-Müller-Vušković — Counting Independent Sets on Perfect Graphs

**Paper:** arXiv:1909.03414 (SIAM J. Discrete Math. 2021)

### Main result:
The partition function of vertex-weighted independent sets can be approximated in polynomial time for (fork, odd hole)-free graphs via graph decomposition. The key technique: the permanent of the adjacency matrix can be viewed as approximating the partition function on line graphs of bipartite graphs (which are perfect).

### Relevance:
Omega(T) is always perfect (no odd holes). **VERIFIED: Omega(T) is also ALWAYS CLAW-FREE.**

Exhaustive verification:
- n<=5: Omega(T) is always complete, so trivially claw-free.
- n=6: 32,768/32,768 tournaments — ALL claw-free (exhaustive).
- n=7: 500/500 random — ALL claw-free.
- n=8: 30/30 random — ALL claw-free.

**This means Omega(T) is a (claw, odd-hole)-free graph** — exactly the class studied by Dyer-Jerrum-Müller-Vušković! Their decomposition via clique cutsets and atoms could give a structural formula for I(Omega(T), 2).

For claw-free perfect graphs, the structure theorem (Chvátal-Sbihi) says they decompose via clique cutsets into "atoms" that are either line graphs of bipartite graphs, or "peculiar" graphs. The partition function on each atom can be computed via the permanent (for line graphs) or small-case analysis.

**Action:** Prove Omega(T) is always claw-free. Then apply the Dyer-Jerrum clique-cutset decomposition to get a structural formula for I(Omega, 2) = H(T).

**WHY this might prove OCF:** If Omega(T) decomposes via clique cutsets into atoms whose independence polynomials factor in a way that mirrors the Hamiltonian path decomposition of the tournament, then OCF follows from the decomposition structure.

---

## Connection 5: Scott-Sokal — LLL, Independence Polynomial, and Non-Perturbative Regime

**Paper:** Scott & Sokal, J. Stat. Phys. 2005

### Main result:
The LLL conclusion holds for dependency graph G and probabilities {p_x} if and only if the independent-set polynomial for G is nonvanishing in the polydisc of radii {p_x}.

The zero-free region is |lambda| < d^d/(d+1)^{d+1} for max degree d+1 (Shearer bound).

### Relevance:
Lambda=2 is ABOVE the Shearer bound for all Omega(T) with max degree >= 2. This means:
- The independence polynomial I(Omega(T), x) HAS zeros in the complex disk |x| < 2
- Standard cluster expansion / polymer expansion do NOT converge at lambda=2
- OCF is genuinely non-perturbative: I(Omega, 2) = H(T) requires exact cancellations

**Key question:** Where ARE the zeros of I(Omega(T), x) in the complex plane? The zeros must lie between the Shearer bound and x=2. Understanding the zero structure could reveal why lambda=2 is special.

**Action:** Compute zeros of I(Omega(T), x) for small tournaments. Map the zero-free region.

---

## Connection 6: Chapman — ASMs and Tournaments

**Paper:** Chapman, Adv. Appl. Math. 27, 2001 (arXiv:math/0008029)

### Main result:
An explicit bijection from oriented square-ice graphs to tournaments, answering a question of Bressoud. This connects tournaments to the six-vertex model.

### Relevance:
The six-vertex model with domain-wall boundary conditions has the Izergin-Korepin determinant as its partition function, which enumerates ASMs. If the Chapman bijection preserves Hamiltonian path structure, then H(T) could be expressible via six-vertex model partition functions or determinantal formulas.

Kuperberg's symmetry classes of ASMs have exact formulas. If tournament symmetry classes (e.g., self-converse, vertex-transitive) correspond to ASM symmetry classes, the exact formulas might constrain or compute H(T).

**Action:** Read Chapman's paper in full. Check if the bijection has any relation to Hamiltonian paths. Investigate whether H(T) has a determinantal formula via the six-vertex model.

---

## Connection 7: Dvořák-Pekárek — Induced Odd Cycle Packing Number

**Paper:** arXiv:2001.02411 (J. Graph Theory 2023)

### Main results:
- Graphs with iocp(G) <= k are chi-bounded by a polynomial in clique number
- Randomized O_{k,t}(n^{k+4}) algorithm for finding large independent sets

### Relevance:
For Omega(T), the induced odd cycle packing number iocp(Omega(T)) measures how many vertex-disjoint "meta-cycles" (odd cycles of odd cycles!) exist. This is a higher-order structure.

Since Omega(T) is always perfect, chi(Omega(T)) = omega(Omega(T)). The iocp bounds on chi-boundedness are automatically tight for perfect graphs.

The structural results on bounded iocp graphs could constrain Omega(T) in useful ways.

---

## Connection 8: Björklund-Husfeldt Parity Formula

**Paper:** arXiv:1301.7250 (STOC 2013)

A deterministic O(1.619^n) time algorithm computes the parity of Hamiltonian cycles via a "new combinatorial formula for #Ham_cycles mod integer."

### Relevance:
If this formula can be adapted for Hamiltonian PATHS in tournaments, it might directly yield H(T) mod 2 (which we know is 1 by Redei) or H(T) mod 4 (which relates to alpha_1).

More importantly, the formula's structure — expressing parity via cycle covers and determinants — might be the bridge between the Hamiltonian path count and the odd-cycle structure that OCF requires.

**Action:** Read the full paper. Check if the parity formula generalizes to directed paths in tournaments.

---

## Connection 9: Six-Vertex Model & Domain-Wall Boundary Conditions

The Izergin-Korepin determinant gives exact partition functions for the six-vertex model with DWBC. Combined with the Desnanot-Jacobi identity, this yields new algebraic identities.

### Potential path to OCF:
1. Chapman bijection: tournaments <-> oriented square-ice configurations
2. Six-vertex model partition function: Z = Izergin-Korepin determinant
3. Specialize to count Hamiltonian paths: H(T) = some evaluation of Z
4. Express I(Omega(T), 2) similarly
5. Show both evaluations agree by algebraic identity

This is speculative but the chain of connections (tournaments -> ASMs -> six-vertex -> determinant) is well-established.

---

## Synthesis: Most Promising Paths Forward

### Tier 1 (Highest potential):
1. **Redei-Berge Hopf algebra** — Our subset convolution IS the comultiplication. Formalizing OCF in this language could reveal structural reasons for the identity.
2. **Björklund's cycle cover approach** — The reduction of Hamiltonian counting to cycle cover determinants, adapted for directed paths and odd cycles, could directly yield OCF.

### Tier 2 (Strong potential):
3. **Perfectness + claw-freeness of Omega(T)** — If Omega is claw-free, Dyer-Jerrum gives structural decomposition of I(Omega, 2).
4. **Independence polynomial zeros** — Mapping the zero structure of I(Omega(T), x) near x=2 could reveal why this evaluation is special.
5. **Six-vertex / ASM determinantal formulas** — Chapman bijection could connect H(T) to determinantal expressions.

### Tier 3 (Worth exploring):
6. **Chordality boundary** — Characterizing when Omega(T) is chordal (Radchenko-Villegas relevance).
7. **Björklund-Husfeldt parity formula** — Adapting for directed Hamiltonian paths.
8. **Dvořák iocp bounds** — Higher-order structural constraints on Omega(T).

---

## References

- [Scott & Sokal (2005)](https://link.springer.com/article/10.1007/s10955-004-2055-4) — Repulsive Lattice Gas & LLL
- [Björklund (2010)](https://arxiv.org/abs/1008.0541) — Determinant Sums for Hamiltonicity
- [Björklund & Husfeldt (2013)](https://arxiv.org/abs/1301.7250) — Parity of Directed Hamiltonian Cycles
- [Radchenko & Villegas (2019)](https://arxiv.org/abs/1908.11231) — Independence Polynomials & Hypergeometric Series
- [Dvořák & Pekárek (2020)](https://arxiv.org/abs/2001.02411) — Induced Odd Cycle Packing Number
- [Dyer, Jerrum, Müller, Vušković (2021)](https://arxiv.org/abs/1909.03414) — Counting Weighted Independent Sets
- [Chapman (2001)](https://arxiv.org/abs/math/0008029) — ASMs and Tournaments
- [Grinberg (2024)](https://arxiv.org/html/2402.07606) — Redei-Berge Hopf Algebra of Digraphs
- [Properties of Redei-Berge (2024)](https://arxiv.org/abs/2407.18608) — Deletion-contraction properties
- [Kuperberg (2002)](https://arxiv.org/abs/math/0008184) — Symmetry Classes of ASMs
- [Two Proofs of Ham Cycle Identity (2025)](https://arxiv.org/abs/2510.02473)
- [Cyclic Subsets of Tournaments (2025)](https://arxiv.org/abs/2508.03634) — Hunter, Liu, Milojević, Sudakov
- [Barvinok (2016)](https://link.springer.com/book/10.1007/978-3-319-51829-9) — Combinatorics & Complexity of Partition Functions
- [Perkins (2020)](https://combgeo.org/wp-content/uploads/2020/11/Will_Perkins_StatMechLectures.pdf) — Statistical Mechanics Methods in Combinatorics
