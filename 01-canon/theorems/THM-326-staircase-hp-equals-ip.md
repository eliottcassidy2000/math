---
id: THM-326
title: Universal Identity H(T) = I(Omega(T), 2)
status: VERIFIED (n=2..6, 36,866 tournaments; staircase k=2..6)
session: opus-2026-05-27-S6
tags: [independence-polynomial, hamiltonian-paths, OCF, universal]
related: [THM-318, THM-319, THM-323, THM-324, THM-325, CONJ-001, MISTAKE-004]
---

## Statement

For any tournament $T$ on $n$ vertices:

$$H(T) = I(\Omega(T), 2)$$

where:
- $H(T)$ = number of Hamiltonian paths in $T$
- $\Omega(T)$ = conflict graph on **all distinct directed odd cycles** of $T$ (two cycles conflict iff they share a vertex)
- $I(G, x) = \sum_k \alpha_k x^k$ = independence polynomial (plain, no $\mu$-weights)

## Verification

- **Exhaustive n=2..6**: 36,866 labeled tournaments (2+8+64+1024+32768), **0 mismatches**
- **Staircase T_k, k=2..6**: verified individually with explicit computation

| k | n | H(T_k) | I(Ω(T_k), 2) |
|---|---|---------|---------------|
| 2 | 4  | 5      | 5             |
| 3 | 6  | 29     | 29            |
| 4 | 8  | 233    | 233           |
| 5 | 10 | 2489   | 2489          |
| 6 | 12 | 33773  | 33773         |

## Critical Implementation Note

**Ω(T) must contain each distinct DIRECTED odd cycle as a separate vertex**, even if two directed cycles have the same vertex set. Deduplicating by vertex set (using frozensets) undercounts Ω vertices and breaks the identity.

Example (n=5 tournament): a 5-cycle vertex set {0,1,2,3,4} can have 3 distinct directed Hamiltonian cycles. These are 3 separate Ω vertices, all pairwise conflicting. This gives $\alpha_1 = 7$ (4 triangles + 3 five-cycles), all pairwise conflicting, so $I(\Omega, 2) = 1 + 14 = 15 = H(T)$.

If vertex-set-deduplicated: $\alpha_1 = 5$ (4 triangles + 1 five-cycle set), $I = 11 \neq H = 15$. **Wrong.**

## Consequence for A038375

$$a(n) = \max_T H(T) = \max_T I(\Omega(T), 2)$$

The maximum Hamiltonian path count equals the maximum of the independence polynomial evaluated at 2. This gives a representation-theoretic perspective on A038375: the optimal tournament maximizes the total weighted count of vertex-disjoint odd cycle collections.

## Connection to OCF

This identity is a consequence of the Odd Cycle Collection Formula (OCF): $H(T)$ = sum over all weighted cycle-path covers of $T$. With weights $p_1=1$, $p_{2j+1}=2$, $p_{2j}=0$ (single vertices contribute 1, odd cycles contribute 2, paths of length ≥2 and even cycles contribute 0), the formula reduces to summing $2^{|S|}$ over all collections $S$ of vertex-disjoint directed odd cycles — i.e., $I(\Omega(T), 2)$.

## Prior history

MISTAKE-004 (retracted) had an incorrect "refutation" of this identity based on mu-weighted computation. The identity was re-confirmed by opus-2026-03-05-S2 with exhaustive n≤6 verification. The current session (opus-2026-05-27-S6) independently re-verified and found a related bug (vertex-set deduplication) that could cause false refutations.
