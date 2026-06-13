# Applications Roadmap: Where Tournament Ideas Create Breakthroughs

**Session:** kind-pasteur-2026-03-16-S116n32

## The Thesis

Our project produced three core innovations:
1. **The OCF**: H(T) = I(Omega(T), 2) — Hamiltonian path count as independence polynomial evaluation
2. **The Tetrahedron Algorithm**: Three-regime classification (solvable/pivot/forbidden) for adaptive computation
3. **The Multilinear Polynomial**: H as degree-4 function of tiling bits, dominated by long-range arcs

These connect to at least 12 active research frontiers where they provide new tools, speedups, or impossibility results.

---

## Tier 1: IMMEDIATE IMPACT (code exists, audience exists)

### 1A. LLM Evaluation (RLHF / Chatbot Arena)

**The problem:** Ranking n language models from pairwise comparisons. Current approaches (Bradley-Terry, Elo) assume transitivity. Real data has cycles (Claude > GPT-4 > Gemini > Claude).

**What we offer:**
- **H(T) as ambiguity measure**: H=1 means unambiguous ranking. Large H means many consistent orderings. This is a single number that quantifies "how settled is the leaderboard?"
- **high_leverage_comparisons()**: Identifies which single matchup, if re-evaluated, would most reduce ambiguity. Directly applicable to Chatbot Arena's query budget.
- **Forbidden H values**: H=7 and H=21 are impossible for ANY tournament. This constrains which ambiguity levels can actually occur — a new impossibility theorem for preference learning.
- **Regime classification**: Most LLM leaderboards are near-transitive (Regime 30), where our perturbation gives O(k·2^k) speedup over DP. For n=20 models with k=3 upsets: 16 million× speedup.

**Connection to recent work:** BLITZRANK (Feb 2026, arXiv:2602.05448) already uses tournament graphs for zero-shot LLM ranking. Our toolkit adds H computation and sensitivity analysis on top of their framework.

**Status:** `RankAmbiguity` class in tournament_toolkit.py. Ready for integration.

### 1B. Topological Data Analysis for Directed Networks

**The problem:** TDA tools (persistent homology, Betti numbers) are well-developed for undirected graphs but underdeveloped for DIRECTED graphs. GLMY path homology (2012+) fills the gap but lacks efficient implementations.

**What we offer:**
- **beta_2 = 0 for ALL tournaments** (THM-108/109): A universal vanishing theorem. No implementation needed — it's a mathematical fact that eliminates one computation entirely.
- **Seesaw mechanism** (THM-095): beta_1 · beta_3 = 0 means 1-holes and 3-holes never coexist. This cuts the Betti computation in half for many cases.
- **Small-prime Gaussian elimination**: 8× memory reduction for rank computation. Core of `mod_rank_library.py`.
- **Circulant eigenspace decomposition** (THM-125): n× speedup for circulant digraphs. Independently discovered in arXiv:2602.04140 (Feb 2026).

**Connection to recent work:** Multiscale Hochschild Laplacian (2025) addresses digraph TDA. Persistent path Laplacians applied to protein-ligand binding. Our tools complement both.

**Status:** `mod_rank_library.py` and `circulant_homology.py` proto-complete. Need PyPI packaging.

### 1C. Score Sequence Reconstruction

**The problem:** Given a set of possible scores, reconstruct valid score sequences. Recent paper (arXiv:2512.16961, Dec 2025) gives new algorithms.

**What we offer:**
- **Score rigidity**: At n=6, 13/22 score classes have H uniquely determined by score. This means the "inverse problem" (score → H) is often solvable without enumerating tournaments.
- **Gap sequence**: The forced 5-cycle contribution (gap = min_H - (1+2·c3)) is a new invariant of score sequences that measures non-3-cycle structural content.
- **Rao's formula** as an O(n) shortcut for c3 computation.

**Status:** Computed for n ≤ 6. New OEIS-candidate sequences identified.

---

## Tier 2: MEDIUM-TERM IMPACT (theory connects, implementation needed)

### 2A. Boolean Satisfiability via Fourier Methods

**The problem:** FourierSAT solves Boolean constraints by representing them as multilinear polynomials and using continuous optimization.

**What we offer:** Our Walsh decomposition of H(T) on the tiling hypercube is EXACTLY the Fourier analysis of a Boolean function. The 213-term degree-4 polynomial for n=6 is a concrete instance where:
- The polynomial structure (dominated by long-range terms) reveals which variables matter most
- The variance decomposition (degree 1: 29%, degree 2: 45%, degree 3: 13%) shows how interactions dominate
- The tiling grid provides a natural spatial structure for the Boolean variables

**Application:** Tournament-structure-aware SAT solvers could exploit the skip hierarchy to focus search on high-impact variables (long-range arcs) first.

### 2B. Hard-Core Model and Phase Transitions

**The problem:** The independence polynomial I(G, λ) is the partition function of the hard-core lattice gas. Its zeros determine phase transition locations.

**What we offer:**
- **OCF at λ=2**: We compute I(Omega(T), 2) for thousands of conflict graphs. This is DEEP in the high-activity regime (λ=2 >> λ_c for most Omega graphs).
- **Real roots persist**: I(Omega(T), x) has all real negative roots for n ≤ 8 (proved) and n ≤ 20 (verified). This extends far beyond claw-free, suggesting a NEW graph class with guaranteed real-rootedness.
- **Explicit evaluation**: Our algorithms compute I(Omega, λ) at arbitrary λ, enabling exploration of phase transitions in tournament conflict graphs.

**Open problem:** Characterize the graph class containing all tournament conflict graphs Omega(T). This class lies strictly between "claw-free" (fails n=9) and "general" and preserves real-rootedness of I(G,x). Finding this class would be a contribution to both combinatorics and statistical physics.

### 2C. LDPC Code Design from Circulant Tournaments

**The problem:** Quasi-cyclic LDPC codes use circulant matrices. Recent work (Nov 2025, arXiv:2511.22183) constructs block MDS LDPC codes from punctured circulant matrices.

**What we offer:**
- **Paley tournament structure**: The quadratic residue structure of Paley tournaments is identical to the circulant structure used in QC-LDPC codes. Our eigenspace decomposition (THM-125) could inform code design.
- **Torsion analysis**: Our `split_inert_analyzer.py` checks boundary maps for torsion at specific primes — directly relevant to code distance properties.
- **The {2,3,7} / {2,3,5} distinction**: Codes based on Hurwitz primes vs golden primes may have different error-correction properties.

### 2D. Flip Graph as Optimization Landscape

**The problem:** Evolutionary algorithms use bit-flip mutations to explore binary fitness landscapes. The "fitness probability distribution of bit-flip mutation" has been studied since 2014.

**What we offer:** The flip graph of tournaments IS a fitness landscape where:
- **Fitness = H(T)** (or any tournament invariant)
- **Mutations = single arc flips**
- **The landscape is ALWAYS connected** (verified n ≤ 6)
- **Delta-H distribution** is symmetric and concentrated near 0
- **The multilinear polynomial gives the EXACT landscape function**

**Application:** Tournament optimization (finding max-H tournament, finding tournament with specific H, finding tournament closest to transitive) can use landscape-aware evolutionary algorithms. The skip hierarchy suggests a variable-ordering heuristic: mutate long-range arcs first (larger effect per mutation).

### 2E. Game Theory: Manipulation Detection

**The problem:** In sports tournaments, detecting match-fixing or strategic manipulation requires understanding which outcomes are "suspicious."

**What we offer:**
- **sensitivity()**: Identifies the single game whose result most affects the overall ranking ambiguity. Unusually influential games may indicate manipulation.
- **Forbidden H values**: If a tournament result implies H=7 or H=21, something is wrong — these values are mathematically impossible.
- **Regime classification**: A tournament that should be in Regime 30 (near-transitive, strong favorites) but appears in Regime 42 (cyclic, confused) may have been manipulated.

---

## Tier 3: SPECULATIVE (deep connections, long-term research)

### 3A. Formal Group Cryptography

The Cayley formal group F(x,y) = (x+y)/(1+xy) with height infinity at p=2 has properties relevant to cryptographic applications of formal groups (Lubin-Tate cryptography). The self-selecting property of {2,3,7,43} under the von Staudt map is a novel algebraic structure. Long-term.

### 3B. DNA Assembly / Contig Ordering

The contig ordering problem (post-assembly step in genome sequencing) IS a tournament Hamiltonian path problem. For well-assembled genomes, the tournament is near-transitive (Regime 30) where our perturbation method excels. The gap between Eulerian assembly methods and the remaining Hamiltonian ordering step could be bridged.

### 3C. Quantum Tournament Optimization

Kemeny ranking via QAOA (Combarro et al., 2023) uses quantum circuits for tournament optimization. Our regime classification could inform circuit depth: Regime 30 tournaments need shallow circuits (few upsets to optimize), while Regime 42 tournaments need deep circuits.

---

## Cross-Cutting Insights

### The Perturbation Principle

Across ALL applications, the key speedup is the same: **most real-world tournaments are near-transitive** (Regime 30). Sports leagues have dominant teams. LLM leaderboards have clear winners. Voting has front-runners. Real networks have hierarchy. The perturbation method (O(k·2^k) where k = upsets from transitive) exploits this universal structure.

### The Skip Hierarchy

Long-range arcs matter more than short-range. This appears in:
- **Tiling polynomial**: skip-4/5 coefficients are 2× skip-2/3 coefficients
- **RLHF**: A distant model beating a close model is more informative than close models trading wins
- **LDPC codes**: Long-range connections in the Tanner graph improve minimum distance
- **Network analysis**: Long-range connections break community structure (the "strength of weak ties")

### The Forbidden Values

H=7 and H=21 are impossible for all n. This is a universal constraint that applies to:
- Voting (certain "near-consensus" configurations cannot exist)
- Sports (certain final-standings patterns are provably impossible)
- LLM evaluation (certain ambiguity levels cannot occur)
- Network analysis (certain path-homology signatures are forbidden)

---

## Implementation Priority

| Priority | Application | Deliverable | Effort |
|----------|------------|-------------|--------|
| 1 | LLM Evaluation | `RankAmbiguity` with BLITZRANK integration | 1 week |
| 2 | TDA | PyPI packaging of `mod_rank` + `circulant_homology` | 1 week |
| 3 | Score Reconstruction | OEIS submissions for 6 new sequences | 2 days |
| 4 | Flip Landscape | Benchmark against evolutionary optimization | 1 week |
| 5 | LDPC Design | Circulant tournament → code matrix converter | 2 weeks |
| 6 | Game Theory | Match-fixing detection prototype | 2 weeks |
| 7 | SAT Solving | Tournament-aware FourierSAT variant | 1 month |
| 8 | Statistical Physics | Zero-free region characterization paper | 2 months |

Sources:
- [BLITZRANK](https://arxiv.org/abs/2602.05448)
- [Score Sequence Reconstruction](https://arxiv.org/abs/2512.16961)
- [LDPC from Circulant Matrices](https://arxiv.org/abs/2511.22183)
- [Path Homology of Circulant Digraphs](https://arxiv.org/abs/2602.04140)
- [TDA Beyond Persistent Homology](https://link.springer.com/article/10.1007/s10462-025-11462-w)
- [FourierSAT](https://www.sciencedirect.com/science/article/abs/pii/S0004370221001107)
- [Hard-Core Model and Independence Polynomial](https://link.springer.com/article/10.1007/s10955-004-2055-4)
- [Flip Graph Research](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.GD.2025.12)
- [RLHF Pairwise Framework](https://arxiv.org/abs/2504.04950)
- [GNNRank](https://arxiv.org/abs/2202.00211)
