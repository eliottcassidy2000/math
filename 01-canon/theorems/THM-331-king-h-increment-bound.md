---
theorem: THM-331
name: King H-Increment Lower Bound
status: PROVED
session: opus-2026-05-27-S1
verified: exhaustively n=3..6 (0 violations)
depends_on: THM-070 (Claim A), THM-019 (King theorem / max-degree reachability)
---

## Statement

Let T be a tournament on n ≥ 3 vertices. Let Q be any vertex with **maximum outdegree** d⁺(Q). Define:
- **Court**: N⁺(Q) = {a : Q → a}, |Court| = d⁺(Q)
- **Rivals**: N⁻(Q) = {b : b → Q}, |Rivals| = n − 1 − d⁺(Q)

Then:

**H(T) − H(T−Q) ≥ 2 · |N⁻(Q)|**

with H(T−Q) = number of Hamiltonian paths in the tournament T with vertex Q removed.

## Proof

By Claim A (proved, THM-070; Grinberg-Stanley):
$$H(T) - H(T-Q) = 2 \sum_{C \ni Q,\, C \text{ odd}} \mu(C) \geq 2 \cdot \#\{\text{directed odd cycles through } Q\}$$

since each μ(C) ≥ 1. It therefore suffices to show that the number of directed odd cycles through Q is at least |N⁻(Q)|.

**The King theorem step.** Q has maximum outdegree, so Q is a king: every vertex is reachable from Q in ≤ 2 steps (classical result, see e.g. Moon 1966). In particular, for every rival b ∈ N⁻(Q) (b beats Q), since b is not in N⁺(Q), it must be reachable via some court member a: Q → a → b. This gives a directed **3-cycle Q → a → b → Q** (since b → Q already). So every rival b is in at least one 3-cycle through Q.

Therefore: #{directed odd cycles through Q} ≥ #{directed 3-cycles through Q} ≥ |N⁻(Q)|.

Combining: H(T) − H(T−Q) ≥ 2|N⁻(Q)|. □

## Tight Cases

The bound is tight (H(T)−H(T−Q) = 2|N⁻(Q)|) in the following cases (exhaustively classified for n ≤ 6):

| n | Tight+SC | Tight+non-SC |
|---|----------|--------------|
| 3 | 6 | 6 |
| 4 | 24 | 56 |
| 5 | 0 | 560 |
| 6 | 0 | 7584 |

**Key finding (see THM-334):** For n ≥ 5, **tight ↔ non-strongly-connected.** No tight SC case exists at n=5 or n=6.

**Tightness at n=3:** All (T,Q) pairs are tight. At n=3 every 3-cycle gives delta=2=2×1, and every transitive tournament gives delta=0=2×0.

**Tightness at n=4:** Tight SC cases have score (1,1,2,2), d⁺(Q)=2, |rivals|=1, delta=2. Here T−Q is the 3-cycle (H=3), H(T)=5=3+2.

## Corollaries

**Corollary 1 (Dominant vertex).** If d⁺(Q) = n−1 (Q beats all others), then |N⁻(Q)| = 0 and Q appears only at position 0 in every HP (see THM-332). Hence H(T) = H(T−Q) and the bound gives H(T)−H(T−Q) ≥ 0 (trivially tight).

**Corollary 2 (Regular lower bound).** For a regular tournament (d⁺(Q) = (n−1)/2, odd n):
$$H(T) - H(T-Q) \geq 2 \cdot \frac{n-1}{2} = n-1$$

This gives the recursive inequality H(T) ≥ H(T−Q) + (n−1). Unrolling: H(regular n-tournament) ≥ 1 + (n−1) + (n−2) + ... = C(n,2). (Very weak; actual growth is much faster.)

**Corollary 3 (Q-P symmetry).** The dual statement: let P have minimum outdegree. For all v in T: v can reach P in ≤ 2 steps. So any vertex v (including Q) is in at least one directed 3-cycle through P, giving H(T) − H(T−P) ≥ 2|N⁺(P)|.

## Upper Bound on 3-Cycles Through Q

From computation: c₃(Q) ≤ |Court| × |Rivals|. Verified n=3..6, 0 failures.

Proof sketch: each 3-cycle Q→a→b→Q uses one court member a and one rival b. The total is bounded by the number of pairs (a,b) ∈ Court × Rivals, i.e., d⁺(Q) × |N⁻(Q)|. □

Full range: c₃(Q) / |rivals| ∈ [1, n−2] at n=5,6.

## Excess Analysis for SC Tournaments

For strongly connected T:

| n | Min excess (delta − 2|rivals|) for SC |
|---|--------------------------------------|
| 3 | 0 |
| 4 | 0 |
| 5 | **2** (score (1,1,2,3,3), deg_Q=3, 1 rival, H=9) |
| 6 | **4** (score (1,1,2,3,4,4), deg_Q=4, 1 rival, H=15) |

The minimum excess for SC grows with n. Conjecture: min SC excess = 2⌊n/2⌋ − 2 for n ≥ 5.
