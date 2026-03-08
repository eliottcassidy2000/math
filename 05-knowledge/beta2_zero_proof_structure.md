# β₂ = 0 for All Tournaments: Proof Structure

## Statement
**Theorem**: For every tournament T on n ≥ 3 vertices, β₂(T) = 0 in GLMY path homology.

## Proof Strategy (Induction via LES)

### Setup
For a vertex v in tournament T, the long exact sequence (LES) of the pair (T, T\v) gives:
```
... → H₂(T\v) → H₂(T) →^{j*} H₂(T,T\v) →^δ H₁(T\v) → ...
```

### Base Cases
- n = 3: dim(Ω₂) ≤ 1, dim(Ω₃) = 0 for all 3-tournaments. Verified β₂ = 0.
- n = 4: Verified exhaustively (64 tournaments, all have β₂ = 0).

### Inductive Step
**Given**: β₂(S) = 0 for all tournaments S on < n vertices.
**Want**: β₂(T) = 0 for any n-tournament T.

Choose an interior vertex v (1 ≤ d⁺(v) ≤ n-2). Such v exists for n ≥ 3.

By induction: β₂(T\v) = 0, so the LES gives:
```
0 = H₂(T\v) → H₂(T) →^{j*} H₂(T,T\v) →^δ H₁(T\v) → ...
```
Since H₂(T\v) = 0, j* is injective: H₂(T) ↪ H₂(T,T\v).

**KEY CLAIM (HYP-258)**: There exists an interior v such that H₂(T,T\v) = 0.

If this claim holds, then H₂(T) ↪ 0, so H₂(T) = 0. QED.

## Computational Evidence for HYP-258

| n | Method | Result |
|---|--------|--------|
| 5 | Exhaustive (1024) | ✓ ALL have interior v with h₂_rel = 0 |
| 6 | Exhaustive (32768) | ✓ ALL have interior v with h₂_rel = 0 |
| 7 | Sampled (500) | ✓ ALL have interior v with h₂_rel = 0 |

## Universal Bound (HYP-262)

**Σ_v h₂(T,T\v) ≤ 3** for all tournaments T (verified n = 5,6,7).

Since every tournament on n ≥ 7 vertices has at least n-2 ≥ 5 > 3 interior vertices,
pigeonhole guarantees some interior v with h₂_rel = 0.

For n = 4,5,6: verified exhaustively.

## Structural Findings

### Dimension universality (HYP-257)
When nonzero: h₂(T,T\v) = 1 and β₁(T\v) = 1. The connecting map δ is a 1D → 1D linear map.

### Necessary condition
h₂(T,T\v) > 0 ⟹ β₁(T\v) > 0. Equivalently, β₁(T\v) = 0 ⟹ h₂(T,T\v) = 0.

### β₁(T) > 0 case (HYP-263)
If β₁(T) > 0, then Σ_v h₂(T,T\v) = 0 (all vertices give trivial relative H₂).
So β₂(T) = 0 is immediate when T itself has nontrivial H₁.

### Source/sink failure (HYP-260)
δ-injectivity fails ONLY at source (d⁺ = n-1) or sink (d⁺ = 0) vertices.
Interior vertices always give injective δ (when δ is defined).
Mechanism: source v appears only at position 0 in paths, limiting boundary diversity.

## What Remains for Full Proof

1. **Prove HYP-262** (Σ h₂_rel ≤ 3) algebraically, OR
2. **Prove HYP-258** directly (∃ interior v with h₂_rel = 0), OR
3. **Prove HYP-259** (δ injective for all interior v)

Any ONE of these, combined with induction, completes the proof.

## Files
- `04-computation/beta2_h2rel_zero_vertex.py` — Tests HYP-258
- `04-computation/beta2_h2rel_sum.py` — Tests HYP-262
- `04-computation/beta2_delta_dimensions.py` — Tests HYP-257
- `04-computation/beta2_delta_sourcesink.py` — Tests HYP-259/260
- `04-computation/beta2_h2rel_which_vertex.py` — Correlation analysis

Author: opus-2026-03-08-S49
