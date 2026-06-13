# The Dual Systems: Interaction vs Self-Loops

*opus-2026-03-23-S261*

## Two Systems, Same Counts

The Royle et al. identity creates two parallel systems with identical cardinalities:

**System A (Interacting):** Even Graphs + Odd Graphs
- |Even| = A000568(n), |Odd| = A000088(n) - A000568(n)
- They INTERACT: edge toggles can cross the even/odd boundary
- The EO boundary has 36 edges at n=5 (49% of all graph meta-graph edges)

**System B (Self-Contained):** Tournaments (no shadow)
- |Tournaments| = A000568(n) — same count as even graphs!
- They DON'T interact with anything: arc reversal stays in tournament-land
- Instead, they have SELF-LOOPS: 7 of 12 classes at n=5

## The Complementarity Principle

**Interaction and self-loops are complementary:**

| Property | System A (Even+Odd) | System B (Tournaments) |
|----------|--------------------|-----------------------|
| Can leave your type? | YES (even→odd via toggle) | NO (tournament stays tournament) |
| Self-loops? | NO (toggle always changes class) | YES (reversal can preserve class) |
| Boundary edges | 36 (EO interaction) | 0 (no boundary) |
| Internal structure | EE=9, OO=29 | E=30, self-loops=7 |

The energy that goes into BOUNDARY INTERACTION in System A goes into SELF-INTERACTION in System B.

## Why This Happens

Both systems are quotients of the same m-cube Q_m by S_n, but with different **semantics**:
- System A: bit k = 1 means edge {i,j} is **present** (undirected)
- System B: bit k = 1 means arc goes **from i to j** (directed)

The Burnside sum is IDENTICAL: (1/n!) Σ_g 2^{c(g_E)}. Each 2^{c(g_E)} term contributes to BOTH counts. But the flip operation differs:
- In System A: toggle bit k = add/remove edge (can create new automorphisms)
- In System B: flip bit k = reverse arc (always gives valid tournament)

This semantic difference creates:
- **Self-loops in B**: reversing an arc can give an isomorphic tournament
- **No self-loops in A**: toggling an edge ALWAYS changes the degree sequence

## The Shadow

The user identifies the "tournament shadow" as the anti-tournament: configurations forbidden by the antisymmetry constraint A[i][j] + A[j][i] = 1. These correspond to the even-order Burnside terms — the same terms that count odd graphs.

|Shadow| = |Odd graphs| = A000088 - A000568.

The shadow doesn't exist as actual objects. It's the COUNT of what's forbidden. This is the particle-antiparticle analogy:
- Even graphs = particles
- Odd graphs = antiparticles
- Tournaments = stable particles (self-contained)
- Anti-tournaments = virtual antiparticles (counted by |Odd|)

## The Gauge Analogy

The two systems represent different **gauges** of the same underlying structure:
- **Graph gauge**: interaction energy lives on the EO boundary
- **Tournament gauge**: interaction energy lives in self-loops
- **Gauge-invariant quantity**: the total Burnside sum Z = A000088(n)

Changing gauge (reinterpreting bits as undirected vs directed) moves energy between boundary interaction and self-interaction, leaving the total count invariant.

## Open Questions

1. Is there a quantitative conservation law relating EO edges and self-loop count?
2. Can we define a "partition function" that unifies both systems?
3. Does the Burnside ring multiplication (kind-pasteur S20dy) reveal the gauge structure?
4. At large n, both systems become "free": EO → 0 (almost all graphs even) and self-loops → 0 (almost all arcs non-neutral). What governs the rate of approach?
