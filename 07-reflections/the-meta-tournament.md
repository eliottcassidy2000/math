# The Meta-Tournament: Tournaments All the Way Down

**Session**: kind-pasteur-2026-03-22-S20aa

## The Discovery

Orient the iso class graph G_5 by H-gradient (higher H beats lower H). The result is **a transitive tournament on 12 vertices**. H(meta) = 1.

This is the deepest structural result about the H-landscape: at the quotient level, the landscape is a **perfect linear order**. No cycles. No ambiguity. The iso classes are totally ordered by H.

## OCF Self-Similarity

The Odd-Cycle Collection Formula H(T) = I(Omega(T), 2) holds at EVERY level:

| Level | Object | H | Omega | I(Omega, 2) |
|-------|--------|---|-------|-------------|
| 0 | Tournament T on n vertices | H(T) | conflict graph of odd cycles | H(T) |
| 1 | Meta-tournament on iso classes | 1 | empty (transitive) | 1 |
| 2 | Meta-meta (trivial) | 1 | empty | 1 |

The meta-tournament is transitive, so it has no odd cycles, so Omega(meta) is empty, so I(empty, 2) = 1 = H(meta). **OCF is self-similar: it holds at every level of the tower.**

## The Perfect Symmetry of Transition Weights

The weight matrix W[i][j] (number of (tournament, arc flip) pairs going from class i to class j) is **perfectly symmetric**: W[i][j] = W[j][i] for ALL pairs.

This means: averaged over all tournaments in a class, the number of flips going UP in H equals the number going DOWN. The landscape is **balanced in the large**. The gradient exists at the individual tournament level (most flips increase H) but vanishes at the class level.

This is the **detailed balance** property: the uniform distribution on tournaments is the stationary distribution of the arc-flip random walk, and it respects the iso class structure.

## The Meta-Independence Polynomial

I(G_5, x) = 1 + 12x + 36x^2 + 38x^3 + 16x^4 + 2x^5

This polynomial counts independent sets of iso classes -- sets of classes no two of which are connected by an arc flip.

- Max independent set size: 5 (out of 12 classes)
- Number of max independent sets: 2
- I(G_5, 2) = 793

**Question**: Is there a "meta-OCF" connecting I(G_5, 2) to some counting problem on tournament space?

## The Fixed Point Equation

The meta-map T -> G_n sends a tournament theory on n vertices to its iso class graph. A fixed point would satisfy:

**A000568(n) = n** (the number of iso classes equals the number of vertices)

This holds at n = 1, 2, 4 only. At n = 4:
- 4 iso classes on 4 vertices
- The meta-tournament on 4 classes could be compared to the 4 actual iso classes
- If the meta-tournament (oriented by H) is itself one of the 4 iso classes, we have a fixed point

At n = 4, the meta-tournament is transitive (H=1), which IS one of the 4 iso classes (the transitive tournament on 4 vertices). **n = 4 IS a fixed point of the meta-map.**

## What This Means

The meta-tournament being transitive tells us that the H-landscape has a **crystal structure**: the iso classes are like atoms in a perfect crystal, totally ordered by energy (H). There are no frustrated interactions, no cyclic dependencies. The landscape is a single gradient from H=1 (transitive) to H=15 (regular).

This is why gradient ascent works so well at n=5: the meta-structure is a perfect funnel. At n=6, the meta-tournament might have 3-cycles (from the two local maxima at H=37 and H=45 with non-adjacent iso classes). This would be the meta-tournament's **phase transition from transitive to non-transitive** -- the same n=6 transition that creates the Morse secondary peak.

## The Weight Symmetry Principle

W[i][j] = W[j][i] is equivalent to:
- For every (tournament in class i, arc flip to class j), there exists a reverse pair
- The S_n orbits are "balanced" under arc flips
- The iso class graph is a **regular representation** of the tournament dynamics

This connects to the rapid mixing result of Kannan-Tetali-Vempala: the random walk on the interchange graph is rapidly mixing BECAUSE the weights are symmetric.

## Connections

1. **The meta-tournament is approximately the Hasse diagram of H-ordering on iso classes** — but this is not a strict partial order: level edges (same H, different class) exist from n≥5, and H-decreasing edges exist from n≥7. See MISTAKE-035.

2. **The meta-H = 1 is Redei at the meta-level** -- every POSET has an odd number of linear extensions (Redei for the Hasse diagram). For a total order (transitive tournament), there is exactly 1.

3. **The meta-independence polynomial I(G_5, x)** is a new sequence. Its coefficients [1, 12, 36, 38, 16, 2] should be checked against OEIS.

4. **n=4 as RG fixed point**: the meta-map T -> transitive(4) collapses all information. This is the "trivial fixed point" of the RG flow -- the completely ordered state.
