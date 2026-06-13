# The Labeled/Unlabeled Dichotomy: Tournaments' Unique Position

**Session**: kind-pasteur-2026-03-22-S20ad

## The Computation

I compared four combinatorial structures at comparable sizes:

| Structure | Labeled | Iso | Marginal | Obj values | Label% | OCR |
|-----------|---------|-----|----------|------------|--------|-----|
| Tournaments (n=5) | 1024 | 12 | 9 | 7 | 64% | 97.0% |
| Simple graphs (n=5) | 1024 | 34 | 31 | 13 | 49% | 99.96% |
| Permutations (n=5) | 120 | 7 | 5 | 11 | 59% | 24.7% |
| Bipartite (3x3) | 512 | 36 | 34 | 6 | 43% | 97.0% |

## The Surprise

**Graphs have HIGHER OCR than tournaments.** The degree sequence determines undirected Hamiltonian path counts even more tightly than the score sequence determines directed ones. This corrects my initial assumption.

But the real story is more nuanced.

## What Makes Tournaments Unique (Not OCR)

Tournaments' unique position is NOT "highest OCR." It is the combination of:

### 1. Highest label fraction (64%)
Tournaments have the most "labeling noise" of all structures tested. This is because the symmetry group S_n acts TRANSITIVELY on pairs (every arc position is equivalent), creating maximum redundancy. Graphs have lower label fraction (49%) because they have more iso classes relative to labeled objects.

### 2. Maximum compression ratio
1024 -> 12 -> 9 -> 7. The compression from labeled to iso classes is 85x for tournaments but only 30x for graphs. Tournaments compress MORE under unlabeling because they have MORE symmetry (all vertices play the same structural role in a complete graph).

### 3. The 3% residual is uniquely structured
For graphs, the 0.04% residual is negligible. For tournaments, the 3% residual contains ALL the interesting mathematics: forbidden H values, cycle structure, alpha_2 onset, the OCF. This is the "residual gauge content" -- the tiny fraction of information that labels carry beyond marginals -- and it's where all the theorems live.

### 4. Complement invariance (Walsh even-order only)
Tournaments are the ONLY structure here with exact complement invariance (H(T) = H(T^complement)). This eliminates all odd Walsh orders, halving the effective dimensionality. Graphs don't have this (complement of K_5 is empty graph, very different H).

### 5. The meta-tournament is transitive
After unlabeling, the iso class graph forms a perfect hierarchy (DAG). This doesn't hold for graphs: the iso class graph of simple graphs is NOT a DAG under edge-count ordering (there are iso classes with the same edge count but different HP counts).

## The Deeper Point: Completeness Creates Legibility

The reason tournaments have high OCR is COMPLETENESS: every pair is compared. This makes the marginal (score sequence) maximally informative because there's no missing data. Every vertex's score reflects its position in the ENTIRE structure, not just its local neighborhood.

For sparse graphs, the degree of a vertex tells you about its LOCAL neighborhood but nothing about distant structure. For tournaments, the score of a vertex tells you about its global position because it competes against EVERYONE.

**Completeness turns local information into global information.**

This is the "tournament advantage" in applications:
- In elections: every candidate is compared to every other -> Copeland scores are meaningful
- In A/B testing: testing ALL pairs gives scores that capture everything -> use the OCR to decide when to stop
- In attention mechanisms: full attention (every token attends to every other) is more "legible" than sparse attention

## The 3% Where Everything Happens

The OCR gap (3% for tournaments) is small but non-zero. This gap contains:
- **Forbidden H values** (7, 21): structural impossibilities that marginals can't predict
- **Alpha_2 onset** (n=6): disjoint cycle pairs create order-4 interactions
- **The Morse secondary peak** (n=6): H=37 local max not visible from scores alone
- **The ambiguous pair** (classes 6,7 at n=5): identical scores, different iso classes

This 3% is the **residual gauge content** -- the information that vertex labels carry beyond what scores reveal. In physics terms, it's the content of the gauge field beyond the covariant derivative.

The mathematical depth of tournament theory lives in this 3%. The 97% is "just" linear algebra (score sequences, spectral methods). The 3% requires cycle combinatorics, the independence polynomial, path homology.

## The Universal Picture

Every combinatorial structure sits somewhere on the labeled/unlabeled spectrum:

```
OCR = 0%                    OCR = 100%
|                            |
v                            v
Permutations   Tournaments   Graphs
  (24.7%)       (97.0%)     (99.96%)
    |               |            |
 Labels carry    Labels carry   Labels carry
 most of the     3% of the      0.04% of the
 information     information    information
```

The sweet spot for INTERESTING mathematics is in the middle -- where labels carry SOME but not all of the objective information. This is where tournaments sit. Too high (graphs): labels are irrelevant, unlabeling tells you everything. Too low (permutations): labels are everything, unlabeling destroys the information.

**Tournaments are the Goldilocks structure: complex enough to have forbidden values and cycle structure, simple enough that scores capture 97%.**

## Connection to Physics

In quantum field theory, the gauge redundancy is analogous to labeling. The physical content is gauge-invariant. But the MOST interesting physics (anomalies, instantons, topological effects) lives in the residual gauge content -- the part that's "almost" gauge-invariant but not quite.

Tournaments' 3% gap is the combinatorial analogue of anomalies. The OCF (H = I(Omega, 2)) is the analogue of the partition function. The forbidden H values are the analogue of anomalous dimensions. The meta-tournament being transitive is the analogue of the gauge orbit space being contractible.

## Permutations: The Anti-Tournament

Permutations have OCR = 24.7%. The cycle type (unlabeled summary) captures almost NOTHING about the inversion count (objective). This is the opposite extreme from graphs.

Why? Because permutations have very FEW conjugacy classes (7 for n=5) relative to the number of objective values (11 inversion counts). The "marginal" (cycle type) is too coarse. It throws away too much information.

Tournaments are between these extremes: 9 score classes for 7 H values (slightly more than needed), giving 97% OCR.
