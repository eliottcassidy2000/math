# THM-098: Pfaffian-Betti Separation at n=6

**Status:** PROVED (exhaustive computation, 32768 tournaments)
**Proved by:** kind-pasteur-2026-03-08-S40 (exhaustive verification), opus-2026-03-07-S46e (initial discovery)

## Statement

For any tournament T on 6 vertices with skew-adjacency matrix S = A - A^T:

1. If β₁(T) > 0 (C-phase), then |Pf(S)| ∈ {1, 3}
2. If β₃(T) > 0 (S-phase), then |Pf(S)| ∈ {7, 9}
3. The sets {1,3} and {7,9} are DISJOINT, so |Pf(S)| perfectly separates C from S

Equivalently, using skew eigenvalues ±ia, ±ib, ±ic (a ≤ b ≤ c, a²+b²+c²=15):

4. If β₁(T) > 0, then spectral gap c - a ≥ 2.922
5. If β₃(T) > 0, then spectral gap c - a ≤ 1.646
6. The gap ranges [2.922, 3.464] and [1.268, 1.646] are disjoint

## Key Data

There are exactly 5 skew eigenvalue triples for n=6 tournaments (on the sphere a²+b²+c²=15):

| (a, b, c) | |Pf|=abc | gap | Achievable phases |
|---|---|---|---|
| (0.146, 2.103, 3.249) | 1 | 3.10 | P, C |
| (0.268, 1.000, 3.732) | 1 | 3.46 | P, C |
| (0.505, 1.732, 3.427) | 3 | 2.92 | P, C |
| (0.727, 2.236, 3.078) | 5 | 2.35 | P only |
| (1.000, 2.646, 2.646) | 7 | 1.65 | P, S |
| (1.732, 1.732, 3.000) | 9 | 1.27 | P, S |

|Pf| = 5 is a **phase boundary** — only P-phase (contractible) appears there.

## Phase Distribution (n=6, exhaustive)

- P (contractible): 27648 (84.4%)
- C (β₁>0, circle-like): 4800 (14.6%)
- S (β₃>0, sphere-like): 320 (0.98%)

## Mechanism

The skew-adjacency matrix S of any tournament on n vertices satisfies:
- tr(S²) = -n(n-1) (constant)
- This constrains eigenvalue sum of squares: Σ aᵢ² = n(n-1)/2

At n=6: a²+b²+c² = 15 is a fixed sphere.

- **C-phase** (β₁>0) requires ANISOTROPIC eigenvalues: one dominant, others small.
  This makes abc = |Pf| small and gap c-a large.
- **S-phase** (β₃>0) requires NEAR-DEGENERATE eigenvalues: all similar.
  This makes abc = |Pf| large and gap c-a small.

The Pfaffian abc measures the "uniformity" of the eigenvalue distribution on the constraint sphere.
By AM-GM, abc is maximized when a=b=c=√5, giving |Pf| ≈ 11.18 (unreachable since |Pf| must be odd integer).

## Extension to n=7, n=8

At n=7 (odd): det(S)=0, so |Pf| is not available. But the spectral gap of S still separates phases:
- S-phase: mean gap 1.82, max 2.59
- C-phase: mean gap 3.47, min 2.40
- P-phase: intermediate (mean 2.62)

At n=8: |Pf| separation is NOT perfect but strongly correlated.
- C-phase has |Pf| ∈ {1,3,5} (small Pfaffian)
- S-phase spans full |Pf| range but enriched at large values
- |Pf| is monotonically negatively correlated with spectral gap

## Connection to OCF

The S-phase tournaments have the highest H(T) on average (mean 28 vs 18 for P at n=6),
consistent with near-degenerate eigenvalues corresponding to regular/near-regular tournaments.
The Paley tournament (all eigenvalues equal, conference matrix) is the extreme S-phase endpoint.

## Scripts

- `04-computation/pfaffian_betti_n7.py`
- `04-computation/pfaffian_topology_deep.py`
- `04-computation/pfaffian_betti_mechanism.py`
- `04-computation/spectral_betti_gap.py`
- `04-computation/spectral_topology_n8.py`

## Open Questions

1. Prove the separation algebraically (WHY does anisotropy force β₁>0?)
2. At what n does the Pfaffian separation break down completely?
3. Is the spectral gap a universal separator at all n? (Appears to break at n≥8)
4. What is the precise c₃ constraint for S-phase? (c₃ ∈ {2, 8} at n=6)
