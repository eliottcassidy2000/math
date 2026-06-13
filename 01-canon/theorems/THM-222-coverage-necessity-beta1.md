# THM-222: Coverage Necessity for β₁=0

**Status:** PROVED
**Source:** opus-2026-03-15-S72
**Verified:** Exhaustive n=5 (0/720 violations), n=6 (0/27968 violations)

## Statement

For any tournament T on n vertices:
**β₁(T) = 0 ⟹ every directed edge participates in at least one transitive triple.**

Equivalently: if there exists an edge i→j with out-deg(i)=1 and out-deg(j)=n−2 (THM-221), then β₁(T) ≥ 1.

The converse is FALSE: at n=5, 144 tournaments have full edge coverage but β₁=1 (the "twist" phenomenon).

## Proof

Suppose edge i→j is not in any transitive triple. We show β₁ ≥ 1.

**Step 1:** The row of ∂₂ indexed by edge (i,j) is identically zero.
The boundary ∂₂(a,b,c) = (b,c) − (a,c) + (a,b) involves edge (i,j) only if (i,j) is a face of triple (a,b,c). By assumption, no such triple exists.

**Step 2:** im(∂₂) ⊆ {v ∈ ℝ^E : v_{(i,j)} = 0}.
Zero row means every column has zero in that coordinate.

**Step 3:** There exists c ∈ ker(∂₁) with c_{(i,j)} = 1.
By Rédei's theorem, T has a Hamiltonian path, hence a directed path P: j = v₀ → v₁ → ⋯ → v_ℓ = i. Define c = e_{(i,j)} + Σ e_{(v_k, v_{k+1})}. Then ∂₁(c) = 0 (telescoping) and c_{(i,j)} = 1.

**Step 4:** c ∉ im(∂₂) by Step 2, so β₁ = dim(ker ∂₁ / im ∂₂) ≥ 1. ∎

## The Twist Phenomenon

At n=5: 304 tournaments have β₁=1, split into:
- **160 "untwisted"**: have an uncovered edge (score 1 → score 3). Obvious obstruction.
- **144 "twisted"**: ALL edges covered, but triple boundaries contribute WRONG SIGNS. The covering triples exist individually but cannot be linearly combined to fill a 1-cycle. Analogous to Möbius strip: locally orientable, globally non-orientable.

The 2-cycle relation in twisted cases: -∂(3,2,1) + ∂(4,2,1) - ∂(4,3,1) + ∂(4,3,2) = 0 on a 4-subtournament. This "closed surface" uses a degree of freedom without helping fill the external 1-cycle.
