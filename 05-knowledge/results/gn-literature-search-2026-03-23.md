# Literature Search: Is G_n (Tournament Isomorphism Class Graph) Known?

**Date:** 2026-03-23
**Searcher:** opus-2026-03-23
**Object:** G_n = graph on tournament isomorphism classes, edges = single arc flip crossing orbit boundaries

---

## Executive Summary

**G_n itself — the arc-flip graph on tournament ISOMORPHISM CLASSES — appears to be a new object.** No paper was found that studies this exact quotient. However, several closely related structures exist in the literature, and understanding how G_n differs from each is illuminating.

---

## Search 1: "Cayley graph of symmetric group on tournaments"

**Verdict: G_n is NOT a Cayley graph.**

A Cayley graph Cay(G, S) has vertices = group elements and edges from a generating set. G_n has vertices = S_n-orbits of tournaments (not group elements). The connection is:

- The **full tournament flip graph** (on labeled tournaments) IS a quotient of the hypercube Q_{C(n,2)} — the binary cube on C(n,2) coordinates, where each coordinate is an arc. G_n is the further quotient Q_{C(n,2)} / S_n, quotienting by the vertex-relabeling action.

- The paper [Isomorphism Classes of Vertex-Transitive Tournaments](https://arxiv.org/abs/2301.09993) studies iso classes of vertex-transitive tournaments specifically but does NOT construct a flip graph on them.

**Key insight:** G_n = (Hamming graph on binary staircase fillings) / S_n. This is a "quotient of a hypercube by a group action on coordinates" — a known construction framework but not specifically studied for tournaments.

---

## Search 2: "Flip graph tournaments / arc reversal reconfiguration"

**Verdict: The interchange graph is the CLOSEST known object, but it operates WITHIN a fixed score sequence, not across all tournaments.**

### The Interchange Graph (Brualdi-Li, 1980s; Ryser)

- **Definition:** Fix a score sequence s. The interchange graph G(s) has:
  - Vertices = all tournaments with score sequence s
  - Edges = pairs differing by a cyclic triangle reversal
- **Ryser's theorem:** G(s) is always connected. Diameter ≤ (n-1)(n-2)/2.
- **Key difference from G_n:** The interchange graph operates on LABELED tournaments within a fixed score sequence. G_n operates on ISOMORPHISM CLASSES across ALL score sequences, with arbitrary single-arc flips (not just triangle reversals).

### Coxeter Interchange Graphs (2023)

- Paper: [Coxeter interchange graphs](https://arxiv.org/abs/2312.04532)
- Extends interchange graphs to "Coxeter tournaments" on signed graphs (including collaborative/solitaire games)
- These are regular graphs; degree determined by distances in Coxeter permutahedra
- Paper: [Random walks on Coxeter interchange graphs](https://arxiv.org/abs/2401.17210) proves rapid mixing for sampling

**Summary:** The interchange graph is the labeled, fixed-score-sequence cousin of G_n. Nobody has taken the quotient by isomorphism AND allowed cross-score-sequence transitions.

---

## Search 3: "Moduli space tournaments"

**Verdict: No moduli space interpretation found for tournaments specifically.**

The search returned general material on GIT quotients and moduli spaces but nothing tournament-specific. However, G_n could naturally be viewed as the 0-skeleton of a moduli space:

- Tournaments = binary fillings of staircase δ_{n-1}
- S_n acts by row/column permutation
- The orbifold [Q_{C(n,2)} / S_n] would be the moduli space
- G_n is its 1-skeleton (adjacency of orbits under single-bit flips)

This framing exists implicitly in the literature on group actions on binary strings but has not been specialized to tournaments.

---

## Search 4: "Johnson graph quotient / orbital graph on pairs"

**Verdict: G_n is related to but distinct from Johnson and orbital graphs.**

- The Johnson graph J(n,k) has vertices = k-subsets of [n], edges = pairs sharing k-1 elements.
- G_n's vertices are orbits of S_n on {0,1}^{C(n,2)}, not k-subsets.
- The **double-coset graph** construction Cay(G, KzK)/K (where (G,K) is DCIS) is structurally analogous: take a Cayley graph and quotient by a subgroup. For G_n, one would need G = Z_2^{C(n,2)} (the binary cube group), K = S_n embedded via its action on pairs, and the generating set = single-coordinate flips.

**Key reference:** [Symmetric Graphs and their Quotients](https://arxiv.org/pdf/1306.4798) discusses quotients of symmetric graphs by normal subgroups and block systems. G_n fits this framework but isn't an example from the paper.

---

## Search 5: "Schur graph / orbital graph on pairs"

**Verdict: The Schur-Weyl graph is a different but spiritually related object.**

- The [Schur-Weyl graph](https://link.springer.com/article/10.1134/S0016266321030023) is a graded graph arising from the RSK correspondence and Schur-Weyl duality. It connects representations, not tournaments.
- **Suborbital graphs** of S_n acting on r-element subsets are studied (e.g., [SUBORBITAL GRAPHS OF THE SYMMETRIC GROUP](https://www.ajol.info/index.php/jagst/article/view/112798/102539)) but these act on subsets of [n], not on binary labelings of pairs.
- For tournaments: S_n acts on pairs {i,j}, inducing an action on {0,1}^{C(n,2)}. The orbitals of this action give a coherent configuration. G_n is the quotient of the Hamming graph by this action — this is an "orbital quotient" but doesn't seem to have been named or studied.

---

## Search 6: "Expander graph with bounded degree and exponential vertices"

**Verdict: G_n belongs to a known regime but may have unusual properties.**

G_n has:
- |V| ~ 2^{C(n,2)} / n! vertices (exponential in n^2)
- Average degree ~ C(n,2) × avg_class_size / |V| = O(n^2) (polynomial in n)
- This gives a VERY sparse graph: degree/vertices → 0 exponentially

Known expander families with similar profile:
- **Cayley graphs of S_n** with transposition generators: n! vertices, degree O(n)
- **Cayley graphs of SL_2(F_p)**: p^3 vertices, bounded degree — these are Ramanujan
- **Quotients of hypercubes** by group actions generally

The question of whether G_n is an expander is open and interesting. The very high connectivity correlation (corr(class_size, degree) ≈ +0.96) suggests large classes are well-connected, but tiny classes (high symmetry) may be bottlenecks.

---

## Search 7: "Association scheme on tournaments"

**Verdict: Coherent configurations on tournaments exist, but not on tournament ISO CLASSES.**

- [Schurian tournaments](https://arxiv.org/abs/1108.5645): A schurian tournament is one whose coherent configuration is schurian (= arises from a permutation group). Polynomial-time isomorphism testing exists for these.
- **Doubly regular tournaments** give non-schurian coherent configurations (smallest: 15 vertices).
- The coherent configuration of S_n acting on tournaments gives an association scheme on the SET of labeled tournaments. Quotienting by S_n collapses this to G_n, but whether G_n itself carries an association scheme structure is unstudied.

**Potential lead:** If G_n has an association scheme structure, its eigenvalues would decompose into idempotent matrices, connecting to the adjacency algebra. The observed eigenvalue ±1 (from complement symmetry) would be one relation class.

---

## Search 8: "Bruhat order on tournaments"

**Verdict: YES — there is a Bruhat order on tournaments, and it is closely related to G_n's H-gradient DAG.**

### Key paper: [Bruhat order of tournaments](https://www.sciencedirect.com/science/article/pii/S0024379514003681) (Brualdi-Fritscher, 2014)

- **Definition:** Two Bruhat orders on tournaments are defined. In the first, T ≤ T' if T can be obtained from T' by a sequence of "source-to-sink path reversals" (reversing arcs along a directed path from source to sink of a transitive subtournament).
- **Key property:** Each such reversal DECREASES the score of the source by 1 and INCREASES the score of the sink by 1 (interior scores unchanged).
- **Cover relation:** T covers T' iff they differ by a SINGLE arc reversal that changes the score sequence in a specific way.
- **For minimum score vectors:** The Bruhat order is isomorphic to the Boolean lattice 2^[k] for some k.

**Connection to G_n:** The H-gradient DAG on G_n (where edges point from lower H to higher H) is an ISOMORPHISM-CLASS-LEVEL QUOTIENT of the Bruhat order. The Bruhat order operates on labeled tournaments; G_n's H-gradient DAG operates on iso classes. The observation that "29/30 edges at n=5 increase H" means the Bruhat structure is almost perfectly preserved at the quotient level.

---

## Search 9 (bonus): "Seidel switching classes"

**Verdict: Seidel switching on GRAPHS is the undirected analogue of what G_n does for tournaments.**

- **Seidel switching:** For a graph G and subset S of vertices, switch all edges/non-edges between S and V\S. Two graphs are switching-equivalent if related by a switching.
- The **switching class graph** (graph on switching equivalence classes, with adjacency from single-vertex switches) is the undirected analogue of the tournament arc-flip graph on iso classes.
- Recent paper: [On identity Seidel switches](https://arxiv.org/abs/2601.04530) studies switches that preserve a graph up to isomorphism — these form a 2-group.
- **Key analogy:** Tournament arc-flip on iso classes : Seidel vertex-switching on iso classes :: tournaments : graphs

---

## Search 10 (bonus): "Quotientopes"

**Verdict: G_n may be the skeleton of a quotientope.**

- **Quotientopes** (Pilaud-Santos) are polytopes whose skeleta are cover graphs of lattice quotients of the weak Bruhat order on S_n.
- The family includes: permutahedra, associahedra, hypercubes, and ~2^{2^n} others.
- Every quotientope skeleton has a Hamilton path (proved by greedy algorithm).
- G_n could potentially be realized as a quotientope if the H-gradient poset is a lattice quotient of the weak order. This would be a remarkable structural result.

**Key reference:** [Combinatorial generation via permutation languages II: Lattice congruences](https://arxiv.org/pdf/1911.12078v1) and the thesis [Shard polytopes and quotientopes](https://hal.science/tel-03499594/).

---

## Conclusion: What IS New About G_n

G_n combines three features that have been studied SEPARATELY but never TOGETHER:

| Feature | Known analogue | How G_n differs |
|---------|---------------|-----------------|
| Arc-flip adjacency | Interchange graph (Brualdi-Li) | Interchange graph fixes score seq; G_n allows all transitions |
| Quotient by S_n | Burnside counting of iso classes | Burnside counts vertices; G_n also tracks edges |
| H-gradient DAG | Bruhat order on tournaments | Bruhat order on labeled T's; G_n quotients by isomorphism |
| Binary cube quotient | Necklaces, switching classes | Those quotient by cyclic/dihedral group; G_n quotients by S_n on pairs |
| Edge coloring (blue/black) | SC-preserving structure | Novel: no analogue in switching or interchange literature |

**G_n is genuinely new as a named, studied object.** The closest existing constructions are:
1. **Interchange graph** (within fixed score sequence, labeled tournaments)
2. **Bruhat order** (labeled tournaments, partial order not graph)
3. **Seidel switching class graph** (undirected graphs, not tournaments)
4. **Quotientope skeletons** (quotients of weak order, not binary cube)

A paper defining and studying G_n systematically would be novel and could connect to all four of these existing threads.

---

## Recommended Next Steps

1. **Check if G_n is a lattice quotient of the weak Bruhat order** — if yes, it's a quotientope
2. **Compute whether G_n carries an association scheme** — check if the adjacency matrices for different "types" of transitions commute
3. **Determine if G_n is an expander** — compute spectral gap at n=5,6
4. **Write up G_n as a formal definition** with the interchange graph and Bruhat order connections as context
5. **Investigate the "modular flip graph" framing** — G_n = flip graph of staircase fillings modulo the mapping class group (= S_n for staircases)
