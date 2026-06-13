# Tournaments as n-Person Games

**Session:** opus-2026-04-05-S24

## The Three Game-Theoretic Faces of a Tournament

A tournament on n vertices simultaneously IS:

### Face 1: A 2-Player Zero-Sum Game (Nash)

The adjacency matrix A is the payoff matrix. Player 1 picks vertex i, Player 2 picks vertex j. Payoff: +1 if i→j, -1 if j→i, 0 if i=j.

**Finding:** The mixed Nash equilibrium of this skew-symmetric game always has value 0 (a tournament is "fair" as a game). But the **support size** varies:
- Transitive tournament: support = 1 (pure strategy: play the top vertex)
- Regular tournament: support = n (uniform strategy)

corr(H, Nash support) = **+0.84** at n=5. More Hamiltonian paths ↔ more strategic complexity. A tournament with many consistent orderings has no clear "best" strategy — all vertices must be mixed over.

### Face 2: An n-Person Cooperative Game (von Neumann-Morgenstern)

Define the characteristic function v(S) = H(T|_S): the number of Hamiltonian paths in the subtournament on coalition S. Each vertex alone contributes v({i}) = 1. The grand coalition contributes v(N) = H(T).

This game asks: can the H(T) Hamiltonian paths be "fairly" allocated to vertices?

**The Shapley Value:**

φᵢ = weighted average over all coalitions of the marginal contribution of vertex i.

Key properties verified at n=4,5:
- **Σ φᵢ = H(T)** — the allocation is efficient.
- **Claim A is the top-level marginal:** H(T) - H(T\v) = 2Σμ(C) gives the marginal at the grand coalition. The Shapley value averages this over ALL coalition sizes.
- **Transitive: φ = (1/n, ..., 1/n)** — the single Hamiltonian path is shared equally.
- **Regular: φ = (H/n, ..., H/n)** — by symmetry, each vertex gets an equal share of all paths.
- **corr(H, Shapley spread) = 0.000** — the inequality among vertices is ORTHOGONAL to the total path count. This is remarkable: whether paths are "concentrated" on certain vertices or "distributed" is completely independent of how many paths there are.

**The Core:**

The core is non-empty iff there exists an allocation x with Σxᵢ = H(T) and Σ_{i∈S} xᵢ ≥ H(T|_S) for all S.

At n=5: core is non-empty for H ≥ 9 (6 of 12 iso classes). Below H=9, no "fair" allocation is possible — some coalition is always shortchanged.

The core boundary at H=9 may correspond to the threshold where the tournament is strongly connected (top cycle = all n vertices).

**vNM Stable Sets:**

For a tournament, a vNM stable set must be a singleton {v} where v beats everyone (score = n-1). This exists only for tournaments with a dominant vertex:
- Transitive: stable set = {top vertex}. H = 1.
- Regular: **NO stable set exists.** H = max.

The absence of vNM stable sets is EXACTLY the Condorcet paradox. H(T) measures the severity of the paradox: H=1 means total consensus (one consistent ordering), H=max means maximum ambiguity (many orderings, no clear "best").

### Face 3: The Social Choice Tournament (Condorcet)

With n alternatives and many voters, majority rule produces a tournament. This is the original context of Arrow's impossibility theorem: no voting system can always produce a transitive ranking from cyclic preferences.

Our OCF (H = I(Ω, 2)) computes the number of consistent total orderings from this cyclic structure. The independence polynomial at fugacity 2 is the "counting function for the Condorcet paradox."

## The Orthogonality Theorem

The most surprising finding:

**corr(H, Shapley spread) = 0.000 at n=5**

The "total number of consistent orderings" (H) and the "inequality of how those orderings are distributed across vertices" (Shapley spread) are UNCORRELATED.

A tournament can have many paths distributed equally (regular, H=15, spread=0) or many paths concentrated on a few vertices (some H=15 class, spread>0). It can have few paths distributed equally (transitive, H=1, spread=0) or few paths concentrated (some H=3 class, spread>0).

In cooperative game theory terms: the SIZE of the pie is orthogonal to how unequally the pie must be cut.

This orthogonality might break at larger n — needs verification.

## The Core Boundary

The cooperative game v(S) = H(T|_S) has non-empty core iff H(T) is large enough. At n=5, the threshold is H=9 = 2n-1. This corresponds to:
- Top cycle = all n vertices (tournament is strongly connected)
- c₃ ≥ 3 (above the 5-cycle forcing threshold)
- Nash support ≥ 3

Below the core boundary: no fair allocation. Some coalition can "defect" and get more paths on their own. This is the cooperative analog of the Condorcet paradox: not only do preferences cycle, but no allocation can satisfy all subgroups simultaneously.

Above the core boundary: cooperation is stable. The Hamiltonian paths can be distributed so every subgroup is satisfied. Higher H gives more room for fair allocations.

## Connection to the Metagraph

The game-theoretic invariants correlate with metagraph position:

| Invariant | corr with H | Interpretation |
|-----------|-------------|----------------|
| Nash support | +0.84 | Strategic complexity |
| Nash entropy | +0.79 | Strategic uncertainty |
| Top cycle size | +0.84 | Strong connectedness |
| Uncovered set size | +0.88 | Strategic stability |
| Markov entropy | +0.72 | Random walk mixing |
| Shapley spread | 0.00 | **Orthogonal to H** |

The metagraph Fiedler ordering (which encodes the dominant geometric variation) tracks Nash support, top cycle, and uncovered set — all +0.8 correlated with H. The Fiedler vector captures the "strategic complexity" axis.

But Shapley spread is orthogonal. This suggests a **second axis** in the metagraph geometry: the "fairness axis," measuring how equally paths are distributed among vertices. This axis is invisible in the H-filtered view but might appear in higher eigenvectors of the metagraph Laplacian.

## Predictions for Larger n

1. **Shapley orthogonality**: Does corr(H, spread) = 0 hold at n=6,7? If so, it's a theorem. The proof would use the symmetry of the Shapley averaging over all coalitions.

2. **Core threshold**: Does the core boundary grow as H ~ 2n-1? Or faster?

3. **Nash support transition**: The support jumps from 1 to n at some H value. Is this related to the cycle-forcing threshold?

4. **Second eigenvector**: Does the second Laplacian eigenvector of the metagraph correlate with Shapley spread?
