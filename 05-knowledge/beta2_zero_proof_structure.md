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
| 8 | Sampled (200) | ✓ ALL have interior v with h₂_rel = 0 |

## Universal Bound (HYP-262)

**Σ_v h₂(T,T\v) ≤ 3** for all tournaments T (verified n = 5,6,7,8).

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

### H₁-killing reformulation (NEW)
**Σ h₂_rel = Σ_v dim(ker(i*_v))** where i*_v: H₁(T\v) → H₁(T) is the inclusion map.

h₂_rel(T,T\v) = 1 iff vertex v "fills" a cycle in T\v — the H₁ class of T\v becomes
a boundary when v-paths are available. The sum counts how many vertices are "critical fillers."

### Key identity
For any tournament T (assuming β₂(T\v) = 0 by induction):
```
n·β₂(T) = Σ_v h₂_rel - Σ_v β₁(T\v) + Σ_v rk(i*_v)
```

### Cone contraction (partial)
The cone operator c_v (prepending v) gives a chain contraction ∂₃∘c_v + c_v∘∂₂ = id
for tournaments where v is a source. For non-source vertices, the cone fails.
At n=5: single-vertex cone works for 600/1024 tournaments; 424 need other mechanisms.

### Worst case analysis
At n=5 with Σ=3: unique isomorphism class (120 tournaments, all relabelings).
- Score sequence [1,2,2,2,3], t₃=4, β₁=0
- 3 critical vertices (cyc=2, scores {1,2,3}); 2 non-critical (cyc=3, scores {2,2})
At n=6 with Σ=3: t₃ ∈ {7,8}, critical vertices have lower cycle counts.

### Cycle count correlation
At n=5: h₂_rel=1 requires exactly t₃(T\v) = 2 (since β₁ = 1 iff t₃ = 2 at n=4).
At n=6: h₂_rel=1 requires β₁(T\v) > 0, threshold varies.
Σ h₂_rel correlates with t₃ but is NOT determined by it.

### β₁(T) vs t₃ relationship
β₁(T) > 0 does NOT correspond to t₃ > 0. A single 3-cycle CAN be a boundary.
- n=4: β₁ = 1 ⟺ t₃ = 2. (t₃=0,1 give β₁=0)
- n=5: β₁ depends on more than t₃. (t₃=3,4 split between β₁=0,1)

### Relative Ω₂ structure
- dim(Ω₂^rel) = d₂(T) - d₂(T\v), varies by vertex
- d₁^rel = n-1 always (universal), d₀^rel = 1 always
- dim(Ω₂) ≠ #TT in general (NT combinations can contribute to Ω₂)
- BUT: relative H₂ generators consist only of TT paths
- ALL relative H₂ generators have ∂₂ mapping to non-v 1-paths

## What Remains for Full Proof

1. **Prove HYP-262** (Σ h₂_rel ≤ 3) algebraically, OR
2. **Prove HYP-258** directly (∃ interior v with h₂_rel = 0), OR
3. **Prove HYP-259** (δ injective for all interior v)

Any ONE of these, combined with induction, completes the proof.

### Most promising approaches
- **H₁-killing bound**: Show that at most 3 vertices can simultaneously be "critical fillers"
  of independent cycle classes. This would prove HYP-262.
- **Cycle counting**: The constraint Σ cyc(v) = 3t₃ limits how many vertices can have
  the right cycle count for h₂_rel = 1, but needs tournament-specific refinement.
- **Direct existence for large n**: For n ≥ 7, the average h₂_rel ≤ 3/n → 0, so a
  probabilistic or pigeonhole argument might work if combined with HYP-262.

## Files
- `04-computation/beta2_h2rel_zero_vertex.py` — Tests HYP-258
- `04-computation/beta2_h2rel_sum.py` — Tests HYP-262
- `04-computation/beta2_h2rel_sum_n8.py` — Tests HYP-262 at n=8
- `04-computation/beta2_delta_dimensions.py` — Tests HYP-257
- `04-computation/beta2_delta_sourcesink.py` — Tests HYP-259/260
- `04-computation/beta2_h2rel_which_vertex.py` — Correlation analysis
- `04-computation/beta2_sum_structure.py` — Algebraic structure of Σ h₂_rel
- `04-computation/beta2_sum_vs_t3.py` — t₃ correlation
- `04-computation/beta2_euler_char.py` — Euler characteristic analysis
- `04-computation/beta2_vpaths_analysis.py` — v-path structure
- `04-computation/beta2_h1_killing.py` — H₁-killing mechanism
- `04-computation/beta2_cone_contraction.py` — Cone operator test
- `04-computation/beta2_cone_allvertex.py` — Cone on Z₂ for all (T,v)
- `04-computation/beta2_h2rel_combinatorial.py` — Combinatorial characterization
- `04-computation/beta2_cycle_count_condition.py` — Cycle count conditions
- `04-computation/beta2_beta1_deletion.py` — β₁ deletion existence
- `04-computation/beta2_h1_survival.py` — H₁ survival under inclusion

Author: opus-2026-03-08-S49
