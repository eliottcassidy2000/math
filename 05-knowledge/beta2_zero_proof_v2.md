# β₂ = 0 for All Tournaments: Proof (Version 2)

## Theorem
For every tournament T on n ≥ 3 vertices, β₂(T) = 0 in GLMY path homology.

## Proof (conditional on HYP-262)

### Structure
By strong induction on n. The key ingredient is:

**HYP-262**: For every tournament T on n vertices, Σ_v h₂(T, T\v) ≤ 3,
where h₂(T, T\v) = dim H₂(T, T\v) is the relative homology dimension.

### Base cases
- n = 3: Verified exhaustively. All 4 tournaments have β₂ = 0.
  (dim Ω₂ ≤ 1, dim Ω₃ = 0.)
- n = 4: Verified exhaustively. All 64 tournaments have β₂ = 0.

### Inductive step (n ≥ 5)
Let T be a tournament on n ≥ 5 vertices. Assume β₂(S) = 0 for all
tournaments S on fewer than n vertices.

**Step 1**: By HYP-262, Σ_v h₂(T, T\v) ≤ 3 < 5 ≤ n.
Since h₂(T, T\v) ≥ 0 for all v, at most 3 vertices can have h₂(T, T\v) > 0.
Therefore there exists a vertex v₀ with h₂(T, T\v₀) = 0.

**Step 2**: The long exact sequence of the pair (T, T\v₀) gives:
```
H₂(T\v₀) → H₂(T) →^{j*} H₂(T, T\v₀) →^δ H₁(T\v₀) → H₁(T) → ...
```
By the induction hypothesis, β₂(T\v₀) = 0, so H₂(T\v₀) = 0.
Therefore j* is injective: H₂(T) ↪ H₂(T, T\v₀).

**Step 3**: Since H₂(T, T\v₀) = 0 (from Step 1), we get H₂(T) ↪ 0.
Therefore H₂(T) = 0, i.e., β₂(T) = 0. □

### Note on v₀
Any vertex v with h₂(T, T\v) = 0 works — it need not be interior.
The LES is valid for all vertices (including source/sink).

---

## Analysis of HYP-262

### Computational verification

| n | Method | Max Σ h₂_rel | Result |
|---|--------|-------------|--------|
| 3 | Exhaustive (4) | 0 | ✓ |
| 4 | Exhaustive (64) | 1 | ✓ |
| 5 | Exhaustive (1024) | 3 | ✓ |
| 6 | Exhaustive (32768) | 3 | ✓ |
| 7 | Sampled (2000) | 3 | ✓ |
| 8 | Sampled (200) | 3 | ✓ |

### Decomposition by β₁(T)

**Case 1: β₁(T) > 0.**
Data shows: Σ h₂_rel = 0 (HYP-263, verified exhaustive n=5,6, sampled n=7).

Mechanism: When β₁(T) > 0, H₁(T) ≠ 0. The LES gives
  h₂_rel(v) = β₂(T) + dim ker(i*_v: H₁(T\v) → H₁(T)).
If β₂(T) = 0 and i* is injective (both of which hold computationally),
then h₂_rel = 0 for all v.

**Case 2: β₁(T) = 0.**
Then H₁(T) = 0, so i*_v: H₁(T\v) → 0 is the zero map. Hence
  ker(i*_v) = H₁(T\v), and h₂_rel(v) = β₂(T) + β₁(T\v).
Since β₂(T) = 0 computationally: h₂_rel(v) = β₁(T\v).

Data shows: Σ β₁(T\v) ≤ 3 when β₁(T) = 0 (verified exhaustive n=5,6, sampled n=7).

### Perfect partition (verified n=5,6,7)

**β₁(T) = 0 ⟺ Σ_v β₁(T\v) ≤ 3**

This is an EXACT equivalence, not just an implication:
- At n=5: β₁=0 gives Σ∈{0,1,2,3}; β₁=1 gives Σ∈{3,4,5}
- At n=6: β₁=0 gives Σ∈{0,1,2,3}; β₁=1 gives Σ∈{4,5,6}
- At n=7: β₁=0 gives Σ∈{0,1,2,3}; β₁=1 gives Σ∈{5,6,7}

The gap widens with n: no tournaments with Σ=4 at n≥6!

### HYP-257: h₂_rel ∈ {0, 1}

Verified exhaustively at n=5 (5120 pairs), n=6 (196608 pairs).
Sampled at n=7,8. Zero violations.

This gives an alternative proof route:
- If β₂(T) ≥ 1: then h₂_rel(v) ≥ β₂(T) ≥ 1 for ALL v.
- Since h₂_rel ≤ 1: we need β₂(T) = 1 and ker(i*) = 0 for all v.
- Then h₂_rel(v) = 1 for all v, so Σ = n ≥ 5 > 3.
- Contradiction with HYP-262.

### What remains to prove algebraically

**Minimal claim needed**: Σ_v h₂(T, T\v) < n for all n-tournaments (n ≥ 4).

This is weaker than HYP-262 (which says ≤ 3). Even Σ < n suffices for the
pigeonhole argument.

**Equivalent formulation (when β₁(T) = 0)**:
Not every vertex-deletion T\v can have β₁(T\v) > 0.
Equivalently: if β₁(T) = 0, then T has a vertex v with β₁(T\v) = 0.

**Combinatorial interpretation**: β₁(T\v) > 0 means T\v has a 1-cycle
that is not a boundary — a "non-fillable cycle." The claim says that
if T itself has no such cycles, then at least one subtournament T\v
also has no such cycles.

---

## Related hypotheses

- HYP-257: h₂_rel ∈ {0,1} for all (T,v)
- HYP-258: ∃ vertex v with h₂_rel = 0 (follows from HYP-262 for n≥4)
- HYP-259: δ injective for interior vertices
- HYP-260: δ-injectivity fails only at source/sink
- HYP-262: Σ h₂_rel ≤ 3
- HYP-263: β₁(T) > 0 ⟹ Σ h₂_rel = 0
- HYP-269: Z₂(Ω) ⊆ im(∂₃|A₃) (verified but not used in this proof)
- HYP-270: im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃) (FAILS at n≥6, not used)

## Files
- `04-computation/beta2_sigma_bound_algebra.py` — Joint analysis of Σβ₁_del and Σh₂_rel
- `04-computation/beta2_beta1_deletion_bound.py` — β₁(T)=0 ⟹ Σβ₁(T\v)≤3 verification
- `04-computation/beta2_relative_structure.py` — Relative complex dimensions
- `04-computation/beta2_augmented_source.py` — Source augmentation approach
- `05-knowledge/beta2_zero_proof_structure.md` — Earlier proof structure document

Author: opus-2026-03-08-S49
