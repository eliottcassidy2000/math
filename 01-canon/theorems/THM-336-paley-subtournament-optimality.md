---
theorem: THM-336
name: Paley Sub-Tournament Optimality
status: CONJECTURED (verified p=7, p=11; strong lower bounds p=19, p=23)
session: opus-2026-05-27-S7
verified: computationally p=7,11 (exact via exhaustive search); p=19,23 (lower bounds only)
depends_on: THM-332, OCF
---

## Statement

For prime p ≡ 3 (mod 4), let QR_p denote the Paley tournament on p vertices
(i → j iff (j−i) mod p is a quadratic residue mod p).

**Conjecture:** For all such p,
- H(QR_p) = a(p) — the Paley tournament maximizes HP count over all p-vertex tournaments
- H(QR_p − v) = a(p−1) — the vertex-deleted Paley tournament maximizes HP count over all (p−1)-vertex tournaments

where a(n) = A038375(n) = max Hamiltonian paths in any n-vertex tournament.

## Verified Values

| p  | a(p−1) = H(QR_p − v)     | a(p) = H(QR_p)            | c(p)                   | Status          |
|----|--------------------------|---------------------------|------------------------|-----------------|
| 7  | 45                       | 189                       | 72                     | CONFIRMED exact |
| 11 | 15,745                   | 95,095                    | 39,675                 | CONFIRMED exact |
| 19 | 117,266,659,317          | 1,172,695,746,915         | 527,714,543,799        | lower bound     |
| 23 | 1,313,333,107,451,805    | 15,760,206,976,379,349    | 7,223,436,934,463,772  | lower bound     |
| 31 | (requires ~8 GB RAM)     | (requires ~8 GB RAM)      |                        | open            |

"CONFIRMED exact" means the value equals a(n) as verified by exhaustive search (a038375 solver).
"Lower bound" means H(QR_p) ≥ a(p) and H(QR_p − v) ≥ a(p−1).

## Immediate Corollaries

1. **Vertex-transitivity leverage**: QR_p is vertex-transitive (Z_p acts), so all H(QR_p − v) are equal.
2. **Dominant-vertex version**: By THM-332, if QR_p has a dominant vertex, then H(QR_p) = H(QR_p − dom). But QR_p is regular (no dominant vertex), so this doesn't directly apply.
3. **H-increment**: H(QR_p) = H(QR_p − v) + 2·c(p), where c(p) = #{directed odd simple cycles through any fixed vertex of QR_p}.

## Asymptotic Formula (new, session S7)

**Empirical law:** c(p)/a(p−1) → (p−1)/4 as p → ∞.

| p  | c(p)/a(p−1)  | (p−1)/4 |
|----|-------------|---------|
| 7  | 1.6000      | 1.5000  |
| 11 | 2.5198      | 2.5000  |
| 19 | 4.5001      | 4.5000  |
| 23 | 5.5001      | 5.5000  |

Error decays as O(1/p²) empirically. This gives:
- **a(p)/a(p−1) → (p+1)/2** as p → ∞ (Paley growth rate)

Exact formula: a(p) = a(p−1) · (p+1)/2 + ε(p) where ε(p) = 2·a(p−1)·δ(p)
and δ(p) = c(p)/a(p−1) − (p−1)/4 → 0.

## Why QR_p is Optimal (heuristic)

QR_p has three properties that likely make it extremal:
1. **Regularity**: all vertices have equal outdegree (p−1)/2. By THM-335, the Q-P gap = 0 → near the top of the H-hierarchy for its score sequence.
2. **Strong connectivity**: QR_p is strongly connected (well-known; follows from THM-333 applied to any spanning Hamiltonian cycle).
3. **Cycle richness**: c(p) = (p−1)/4 · a(p−1) + O(1/p) — as many odd cycles per vertex as possible given regularity.

The combination of regularity + SC + maximum cycle density seems to force optimality. No proof.

## Computation Method

H(QR_p) is computed via the **circulant-reduced bitmask DP** (see `04-computation/paley_hp_counter.c`):

Standard DP: O(2^p · p) memory. For p=23: 1.5 GB.
Circulant DP: O(2^{p−1}) memory (factor 2p savings). For p=23: 33 MB.

Key insight: QR_p is Z_p-invariant. Normalize current vertex = 0. State = bitmask of visited offsets. All p starting vertices are equivalent → initialize weight = p.

**CRITICAL**: must process states in order of path length (popcount), not mask value, since new_mask_hi can be < old mask_hi.

## Connection to Sequence Extensions

With p=31 (next Paley prime ≡ 3 mod 4): a(30) = H(QR_31 − v), a(31) = H(QR_31).
Required: 2^30 × 8 = 8 GB RAM and ~150s compute.

Subsequent primes: p = 31, 43, 47, 59, 67, 71, 79, 83, ...
These give new terms a(30), a(31), a(42), a(43), ... in A038375.
