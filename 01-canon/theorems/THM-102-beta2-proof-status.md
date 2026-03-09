# THM-102: β₂ = 0 for Tournaments — Proof Status

**Status:** COMPUTATIONALLY VERIFIED, ALGEBRAIC PROOF INCOMPLETE
**Filed by:** opus-2026-03-08-S43

## Statement

**Conjecture (THM-100):** For every tournament T on n vertices, β₂(T) = 0 in GLMY path homology.

## Computational Evidence

| n | Method | Total | β₂=0 | Failures |
|---|--------|-------|------|----------|
| 3 | exhaustive | 8 | 8 | 0 |
| 4 | exhaustive | 64 | 64 | 0 |
| 5 | exhaustive | 1024 | 1024 | 0 |
| 6 | exhaustive | 32768 | 32768 | 0 |
| 7 | sampled | ~5000 | ~5000 | 0 |
| 8 | sampled | ~5000 | ~5000 | 0 |
| 9 | sampled | ~1000 | ~1000 | 0 |
| 10 | sampled | 50 | 50 | 0 |

**Additionally**: Paley T_7 and T_11 have β₂=0. rank(d₃) = ker(d₂) EXACTLY for all tested.

## Structural Discoveries

### 1. DT+Cancellation Filling (THM-101)
The space Z₂ = ker(∂₂|Ω₂) is spanned by boundaries from:
- **DT 4-paths**: (a,b,c,d) with a→c AND b→d (doubly-transitive)
- **Cancellation pairs**: differences of 3-paths sharing a bad face

Verified exhaustive at n≤6 with 0 failures.

### 2. DT Face Structure
A DT path (a,b,c,d) has face pattern:
- (b,c,d): always TT (b→d is DT condition)
- (a,b,c): always TT (a→c is DT condition)
- (a,c,d): TT iff a→d, NT iff d→a
- (a,b,d): TT iff a→d, NT iff d→a

**Key insight**: The edge a↔d is FREE (not constrained by DT). When d→a,
DT boundaries include NT 2-path components, allowing DT boundaries to
reach the NT part of Z₂.

### 3. DT-Only Deficit
- n=5: DT alone fills Z₂ for ALL 1024 tournaments (100%)
- n=6: DT alone fails for 960/32768 (2.9%), deficit always exactly 1
  - All failures: score (1,2,2,3,3,4) or (2,2,2,3,3,3)
  - |DT| = 9, dim(Z₂) = 10, rk(∂₃|DT → Z₂) = 9

### 4. Ω₂ Structure
- dim(Ω₂) = |TT triples| + dim(NT cancellation space)
- Z₂ has NT components (784/1024 at n=5, 95% at n=6)
- DT boundaries reach NT components through the free a↔d edge

### 5. Edge Removal
- Removing one edge from a tournament can create β₂=1 (4/512 at n=5)
- The unfillable 2-cycle always involves BOTH endpoints of the removed edge
- Completing any oriented graph to a tournament kills β₂ (n=4 exhaustive)
- No bidirected edge addition creates β₂>0 from a tournament (n=4)

### 6. Relative Homology
- H₂(T, T\v) = 0 for all T, v: verified exhaustive n≤5, sampled n=6
- β₁ can increase under vertex deletion (840/5120 at n=5)
- But ∃ vertex v with β₁(T\v) ≤ β₁(T) always (n=5 exhaustive)

## Proof Approaches

### A. Direct DT+Cancellation (Most Developed)
**Idea**: Show im(∂₃|Ω₃) ⊇ Z₂ by decomposing any 2-cycle into DT boundaries plus cancellation elements.

**Status**: Computationally verified n≤6. Algebraic proof needs:
- Why DT boundaries span a large enough subspace of Z₂
- How the free a↔d edge provides enough NT coverage
- Why the deficit is always fillable by cancellation elements

### B. Vertex Deletion Induction
**Idea**: Induction on n using the long exact sequence:
  0 = H₂(T\v) → H₂(T) → H₂(T,T\v) → H₁(T\v) → H₁(T)

**Status**: Reduces to proving H₂(T,T\v) = 0, which is verified computationally.
The relative complex R_*(v) = Ω_*(T)/Ω_*(T\v) needs:
- dim(R₂) = 5-7, dim(R₃) = 2-9 at n=5
- ker(∂₂^rel) small (0-2), rk(∂₃^rel) ≥ ker(∂₂^rel)
- Need to show this rank inequality algebraically

### C. Rank Formula
**Idea**: Prove rk(∂₂) + rk(∂₃) = dim(Ω₂) directly.

**Status**: Verified exhaustive n=5. Equivalent to β₂=0 but gives no independent insight.

### D. Discrete Morse Theory
**Idea**: Construct acyclic matching on Ω₃ → Ω₂ that eliminates all critical 2-cells.

**Status**: Not yet attempted. Kozlov's theory (arXiv:cs/0504090) applies to free chain complexes. Tournaments are NOT transitive, so existing matchings don't work. Need a tournament-specific construction.

### E. Spectral Sequence
**Idea**: Use filtration of Ω₂ by "bad face count" or vertex support.

**Status**: Caputi-Menara (arXiv:2503.06722) proves results for TRANSITIVE tournaments only. General case open.

### F. Arc-Flip Invariance (NEW — Most Promising, S41)
**Idea**: Prove β₂ = 0 by showing β₂ is invariant under arc flips, then using β₂(T_trans) = 0.

**Discovery** (kind-pasteur-S41): Under ANY arc flip u→v to v→u in a tournament:
  δ(dim Z₂) = δ(rk ∂₃) EXACTLY.

This means β₂ = dim(Z₂) - rk(∂₃) is an **arc-flip invariant**.

**Proof strategy**:
1. β₂(T_trans) = 0 (transitive tournament = simplex, contractible) ✓
2. Any tournament is reachable from T_trans by arc flips ✓
3. β₂ is invariant under arc flips (HYP-233, to prove)
4. Therefore β₂ = 0 for all tournaments.

**Equivalent reformulation**: δ(rk ∂₃) = δ(dim Ω₂) + δ(β₁).
Key structural facts:
- dim(Z₁) = (n-1)(n-2)/2 is CONSTANT for all tournaments
- δ(rk ∂₂) = -δ(β₁) (follows from constant dim(Z₁))
- δ(dim Ω₂) = p_{vu} - p_{uv} (analytic formula, PROVED)
- δ(β₁) is NON-LOCAL (depends on global tournament structure)

**Status**: Verified exhaustive n=5 (10240 flips), sampled n=6 (15000), n=7 (2500), n=8 (500). 0 mismatches. N=6 exhaustive in progress.

**Structural insight**: dim(Ω₂) = |TT triples| + |NT cancellation dimensions|. Under arc flip, TT triples and NT cancellation pairs rebalance to maintain β₂ = 0.

### G. Cone-from-T' Construction (NEW — S42)
**Idea**: For each vertex v, swap cycles z at v can be filled by
  w = Σ α_{abc} [(v,a,b,c) + (a,b,c,v)]
where (a,b,c) ranges over allowed 2-paths in T' = T\{v}.

**Key properties**:
- T'-internal faces cancel: d₃(v,a,b,c) + d₃(a,b,c,v) has zero (a,b,c) component
- The B matrix B·α = z is always solvable (unfiltered: ALL T' paths)
- Filling is automatically in Ω₃ at n=5,6 (100%). Breaks slightly at n≥7 (~93-98%)
- Rank surplus grows with n: min(rank(B) - swap_dim) = 2(n=5), 6(n=7), 11(n=8), 15(n=9)

**Status**:
- Unfiltered single-vertex: 500/500 at n=7,8; 200/200 at n=9. Zero failures.
- Filtered (v→a, c→v only): fails 1/1000 at n=8 (too restrictive)
- Multi-vertex (all v): ALWAYS works, including n=8 filtered failure case

**Proof needs**: Show rank(B_unfiltered) ≥ swap_dim algebraically.

**Scripts**: `beta2_filtered_cone.py`, `beta2_cone_failure.py`, `beta2_unfiltered_large.py`, `beta2_omega_membership.py`

## Dead Ends for Algebraic Proof
- TT subcomplex: not exact (S42)
- A_* projected complex: ∂∂ ≠ 0 (S42)
- Flag/simplicial complex: path H₂=0 even when simplicial H₂≠0 (S41)
- Extension lemma: 50% of TT triples at n=4 have NO DT extension (S43)
- Z₂ ⊆ span(TT): FALSE, Z₂ has NT components (S43)
- Euler characteristic: χ ≠ 1-β₁ when β₃>0 (S43)
- Naive front/back cone: sign reversal prevents direct formula (S42)
- Single-vertex filtered cone: fails at n=8 (too few valid T' paths) (S42)
- LES via i_* injectivity: i_*: H_1(T\v)→H_1(T) rarely injective (S42)

## Key Open Question

**Why does tournament completeness force β₂ = 0?**

The essential mechanism is: for every pair (u,v), exactly one of u→v or v→u holds. This means:
1. Every 4-vertex subtournament has β₂=0 (finitely many cases)
2. DT paths exist abundantly (tournament completeness creates transitive structure)
3. The free edge a↔d in DT paths provides NT boundary components
4. Edge removal breaks this by creating a "gap" that 2-cycles can hide in

The proof should formalize why there are "enough" DT paths and cancellation elements to fill every 2-cycle, using tournament completeness as the key ingredient.

## Scripts
- `04-computation/beta2_dt_cancel_filling.py` — main DT+cancel verification
- `04-computation/beta2_section_construction.py` — section s: Z₂ → Ω₃
- `04-computation/beta2_dt_deficit_analysis.py` — deficit structure at n=6
- `04-computation/beta2_twin_obstruction.py` — completing to tournament kills β₂
- `04-computation/beta2_edge_removal_anatomy.py` — unfillable cycle anatomy
- `04-computation/beta2_boundary_structure.py` — Z₂ TT/NT decomposition
- `04-computation/beta2_omega2_relations.py` — Ω₂ relations, DT face NT-ness
- `04-computation/beta2_extension_lemma.py` — TT extension to DT
- `04-computation/beta2_local_to_global.py` — 4-vertex subset filling
- `04-computation/beta2_ses_proof.py` — relative homology computation
- `04-computation/beta2_relative_correct.py` — β₁ monotonicity
- `04-computation/beta2_arcflip_exactness.py` — arc-flip invariance discovery (S41)
- `04-computation/beta2_arcflip_invariance.py` — verification at n=7,8 (S41)
- `04-computation/beta2_arcflip_mechanism.py` — algebraic mechanism analysis (S41)
- `04-computation/beta2_arcflip_example.py` — concrete examples (S41)
- `04-computation/beta2_arcflip_exhaustive_n6.py` — exhaustive n=6 verification (S41)

## See Also
- THM-100 (original conjecture)
- THM-101 (DT+cancel filling theorem)
- HYP-207 through HYP-232 (related hypotheses)
- HYP-233 (arc-flip invariance)
- HYP-234 (proof via arc-flip invariance)
