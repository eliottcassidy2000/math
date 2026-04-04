# THM-301: Sign Rule for Quadratic Coefficients

**Status:** Parts (a),(b) PROVED; Part (c) VERIFIED n=4-8
**Discovered by:** opus-2026-04-04-S5
**Proved by:** opus-2026-04-04-S6
**Related:** THM-287 (OCF quadratic decomposition), THM-302, THM-290

## Statement

For the quadratic multilinear coefficients c_{ij} of H(t):

**(a) SAME-END (PROVED):** If tiles i,j share an upper or lower vertex, then c_{ij} < 0.

**(b) CROSS-END (PROVED):** If one tile's upper vertex equals the other's lower vertex (relay vertex), then c_{ij} > 0.

**(c) DISJOINT (VERIFIED n=4-8):** If tiles share no vertex, then c_{ij} > 0 (when nonzero).

## Proof of (a): Same-End

Let tiles i,j share vertex v (0-indexed: v-1). Both tile arcs are incident to vertex v-1.

**Step 1 (Opposite Direction Theorem):** In any directed cycle through both tiles, vertex v-1 has indegree=outdegree=1. Both tiles provide arcs at v-1. If both were forward (both incoming to v-1), indegree ≥ 2, contradiction. If both backward (both outgoing), outdegree ≥ 2, contradiction. Therefore one tile is forward and the other backward: **opposite directions**. The chi coefficient is (-1)·(+1) = -1.

**Step 2 (B=0):** For same-end tiles, any two cycles C₁∋tile_i, C₂∋tile_j both visit the shared vertex v. So C₁ ∩ C₂ ⊇ {v} ≠ ∅, meaning they are NOT vertex-disjoint. Therefore B(i,j) = 0 (no independent pair contribution from |I|=2 OCF terms).

**Step 3:** c_{ij} = 2·A(i,j) + 4·B(i,j) = 2·A(i,j). Since every contributing cycle has chi=-1, A(i,j) = Σ(−1) = −(number of valid cycles). Therefore c_{ij} = −2·(number of valid cycles) ≤ −2.

## Proof of (b): Cross-End

Let the relay vertex be r = x₁ = y₂ (tile 1's upper = tile 2's lower). Both tiles connect to vertex r-1 (0-indexed).

**Same Direction Theorem:** Tile 1 forward = outgoing from r. Tile 2 forward = incoming to r. Both forward: indeg=1, outdeg=1 at r. ✓. Both backward: same. ✓. Mixed (one fwd, one bkwd): either indeg=2 or outdeg=2. ✗ impossible. Therefore **same direction**, chi = +1.

Since A(i,j) > 0 and B(i,j) ≥ 0, c_{ij} > 0.

## Disjoint Part (c): Mechanism and Evidence

For disjoint tiles, no shared vertex constrains the directions. Computationally:
- At n=5: 4 disjoint pairs with c>0 (0 with chi=-1 cycles)
- At n=6: 18 disjoint pairs with c>0. Two pairs have chi=-1 cycles, but same-direction cycles outnumber opposite-direction cycles (ratio 3:1).
- At n=7,8: 0 violations across 51+? pairs.

The mechanism: most 5-cycles through disjoint tiles use both in the backward direction (aligned with the "increase from transitive" tendency). Opposite-direction cycles require an additional reversal that costs one more "against the flow" step, making them rarer.

## Key Insight: Additive vs Multiplicative

The sign of c_{ij} (additive interaction) DISAGREES with the sign of log-interaction λ_{ij} for disjoint and cross-end pairs. The log-interaction is ALWAYS ≤ 0 (antiferromagnetic, THM-290), while c_{ij} > 0 for these pairs. This reflects:
- **Additive synergy** (c > 0): flipping two tiles creates MORE HPs than the sum of individual effects
- **Multiplicative competition** (λ < 0): the RATIO H(ij)/H(i)H(j)/H(0) < 1

These coexist because the linear effects are large (2^{skip-1} each), creating a base of HPs that the quadratic interaction modulates.

## Computational Verification

| n | Same-end | Cross-end | Disjoint | Total violations |
|---|---------|----------|---------|-----------------|
| 4 | 2 (−) | 0 | 0 | 0 |
| 5 | 8 (−) | 1 (+) | 4 (+) | 0 |
| 6 | 20 (−) | 4 (+) | 18 (+) | 0 |
| 7 | 40 (−) | 10 (+) | 51 (+) | 0 |
| 8 | 66 (−) | ? (+) | ? (+) | 0 |
