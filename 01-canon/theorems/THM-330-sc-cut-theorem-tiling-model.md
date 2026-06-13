---
theorem: THM-330
name: Strong Connectivity Cut Theorem (Tiling Model)
status: PROVED
session: opus-2026-05-27-S1
verified: exhaustively n=3..7 (0 errors)
---

## Statement

In the tiling model on n vertices (base path n-1 → n-2 → ... → 0, with tiles at non-consecutive pairs), a tiling t is **strongly connected** if and only if for every cut index k ∈ {1, ..., n-1}, at least one tile (x, y) with x ≥ k, y < k, x ≥ y+2 has bit value 1 (upward orientation: y beats x).

Equivalently: the tournament T_t is NOT strongly connected if and only if there exists k such that **all** tiles crossing cut k are in the default downward orientation (x beats y).

## Proof

**Setting up the cut structure.** In any directed graph, a tournament T is not strongly connected iff there exists a non-empty proper subset S ⊂ V such that all arcs between S and V\S go from S to V\S (S is a dominating set).

**Base path forces upward-closed dominating sets.** The base path contributes arcs k → k-1 for k = 1,...,n-1. Suppose S dominates V\S. If k ∈ V\S and k-1 ∈ S, the arc k → k-1 goes from V\S to S, contradicting domination. Hence: if any vertex k-1 ∈ S, then k ∈ S (propagating upward). By induction: S must be an **upward-closed interval** S = {k, k+1, ..., n-1} for some k ∈ {1,...,n-1}.

**Cut k is the bottleneck.** With S = {k,...,n-1}, the base path arc k → k-1 goes from S to V\S (consistent). All remaining arcs between S and V\S are the **tiles crossing cut k**: tiles (x,y) with x ∈ S (x ≥ k), y ∈ V\S (y < k), x ≥ y+2 (non-consecutive). S dominates iff ALL these tiles are in the downward orientation.

**Both directions.** T is not SC ↔ ∃ such S ↔ ∃ k with all cut-k tiles downward. Contrapositive: T is SC ↔ every cut k has at least one upward tile. □

## Cut Sizes

The number of tiles crossing cut k is:

f(n, k) = (k−1) + k(n−1−k)

**Minimum:** f(n,1) = f(n, n-1) = n-2. These are the **minimum cuts**, achieved at the endpoints. Cut k=1 consists of tiles {(x,0): x ≥ 2} (the horizontal leg of the staircase). Cut k=n-1 consists of tiles {(n-1, y): y ≤ n-3} (the vertical leg).

**Maximum:** f(n, ⌊n/2⌋) ≈ n²/4.

## Corollaries

**Corollary 1 (Apex SC Theorem — see THM-333).** If the apex tile (n-1, 0) has bit=1, all n-1 cuts are satisfied, so T is SC.

**Corollary 2 (Geometric interpretation).** Strong connectivity requires that neither the source (vertex n-1) nor the sink (vertex 0) is "absolute" in the tiling sense. Specifically:
- Cut k=1 satisfied ↔ vertex 0 is NOT an absolute sink (has at least one out-tile).
- Cut k=n-1 satisfied ↔ vertex n-1 is NOT an absolute source (at least one vertex beats it).

**Corollary 3 (SC fraction).** The number of strongly connected tilings grows rapidly toward 2^m. Computed values:

| n | m=C(n-1,2) | Total tilings | SC tilings | Non-SC |
|---|-----------|---------------|------------|--------|
| 3 | 1         | 2             | 1          | 1      |
| 4 | 3         | 8             | 5          | 3      |
| 5 | 6         | 64            | 50         | 14     |
| 6 | 10        | 1024          | 903        | 121    |
| 7 | 15        | 32768         | 30773      | 1995   |
| 8 | 21        | 2097152       | 2032504    | 64648  |

Non-SC tiling counts: 1, 3, 14, 121, 1995, 64648 (see HYP-1739 for pattern).

## Relation to Staircase Triangle

The three sides of the staircase triangle (see `07-reflections/everything-is-the-triangle.md`) correspond to the cut structure:
- **Vertical leg** (source column): tiles at cut k=n-1 — the minimum cut involving vertex n-1.
- **Horizontal leg** (sink row): tiles at cut k=1 — the minimum cut involving vertex 0.
- **Hypotenuse** (apex region): the apex tile (n-1,0) crosses ALL cuts simultaneously.

Strong connectivity requires the triangle to be "permeable" in both the horizontal and vertical directions, at every intermediate level.
