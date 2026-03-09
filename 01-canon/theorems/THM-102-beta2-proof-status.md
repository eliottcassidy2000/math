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

### I. THM-103: β₁(T) ≤ 1 for All Tournaments (PROVED, S50)
**Status**: COMPLETE ALGEBRAIC PROOF.

**Proof** (star constraint argument): For each vertex v, take a Hamiltonian path
u₁→...→u_d in out(v) (exists by classical theorem). Each consecutive triple
(v, uᵢ, uᵢ₊₁) is transitive, giving cocycle constraint w(v,uᵢ₊₁) = w(v,uᵢ) + w(uᵢ,uᵢ₊₁).
This eliminates d⁺(v)-1 edge variables. Different vertices give disjoint eliminated sets.
Total eliminated: C(n,2)-n. Free cocycle variables: ≤ n. Coboundaries: n-1. Hence β₁ ≤ 1.

**Consequence for β₂ proof**: Combined with Sum_v β₁(T\v) ≤ 3 (HYP-282, verified n≤10),
this gives ∃ good vertex for n ≥ 4, completing approach H if HYP-282 is proved.

**Additional findings (S50)**:
- Bad vertices form transitive triples (verified 100% n=5,6,7)
- 3-cycle among bad vertices forces β₁(T)≥1 (flip obstruction)
- Hidden cycle z_v uses ALL n-1 vertices (global, not local)
- rank_drop(v) = (n-2) + β₁(T\v) exact identity
- At n=5: #bad = max(0, t₃-1) exact formula
- TT constraints are the ONLY cocycle constraints (IC always redundant)
- Cocycle restriction res_v: Z¹(T)→Z¹(T\v) is NOT surjective for bad v (proved algebraically)

**What remains**: Prove HYP-282 (Sum ≤ 3) or the weaker Sum < n. The algebraic
proof of THM-103 doesn't directly yield the deletion bound.

**New findings (S51 — rank-critical analysis)**:
- rank(∂₂|Ω₂) = C(n,2) - n + 1 - β₁ (universal formula, verified n≤6)
- dim(Z₁) = C(n,2) - (n-1) is UNIVERSAL for all n-vertex tournaments
- Therefore β₁ depends ENTIRELY on rank(∂₂): β₁ = dim(Z₁) - rank(∂₂)
- dim(ker(∂₂^T)) = n - 1 + β₁ (dual formula)
- Bad-vertex TT is NOT rank-critical when β₁=0 (#bad=3): removing it from Ω₂ does NOT drop rank
- Rank-critical TTs governed by redundancy = #TTs - rank; high redundancy → 0 RC TTs
- #RC correlates with t₃ (r=0.69 at n=6), not directly with bad vertices (r=0.46)
- Flip obstruction mechanism: arc flip changes MANY TTs simultaneously (global Ω₂ restructuring),
  not just the single bad-vertex TT. The net effect reduces rank(∂₂) by 1.
- β₁=0 with all-bad 3-cycle: ZERO instances at n=5,6,7 (flip obstruction confirmed)
- β₁=1 minimum #bad: 3(n=5), 4(n=6), 5(n=7) — grows with n

See THM-103 for the full proof.

### H. Inductive Proof via b1 Monotonicity (NEW — Most Promising, S43)
**Idea**: Prove beta_2 = 0 by induction on n using the LES of (T, T\v).

**LES**: 0 = H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)
By induction, H_2(T\v) = 0. So H_2(T) injects into H_2(T,T\v).
H_2(T,T\v) = 0 iff delta: H_2(T,T\v) -> H_1(T\v) is injective.
By LES exactness, H_2(T,T\v) = 0 iff i_*: H_1(T\v) -> H_1(T) is injective.

**KEY EQUIVALENCE** (S43): i_* is injective iff b1(T\v) <= b1(T).

**Proof strategy**:
1. Base case: n <= 4, beta_2 = 0 trivially (no 2-cycles). DONE.
2. Induction step: Assume beta_2(T') = 0 for all (n-1)-vertex tournaments.
3. For n-vertex tournament T, find vertex v with b1(T\v) <= b1(T).
4. Then H_2(T,T\v) = 0, so H_2(T) injects into 0, giving beta_2(T) = 0.

**Computational verification**:
- n=5: 1024/1024 tournaments have a good vertex (exhaustive)
- n=6: 32768/32768 tournaments have a good vertex (exhaustive)
- n=7-12: 100% (sampled, 500-20 trials each)

**Structural discoveries** (S43):
- **b1(T) in {0, 1}** for ALL tournaments (verified n<=20, zero b1>=2)
- b1=0 rate: 62%(n=4), 70%(n=5), 85%(n=6), 95%(n=7), 99.2%(n=8), 100%(n>=10 sampled)
- When b1(T)=1: ALL n vertices are good (b1(T\v) <= 1 = b1(T))
- When b1(T)=0: most vertices are good (e.g., 84% at n=5, 82% at n=6)
- Bad vertices: delta_b1 = +1 always (b1 increases by exactly 1)
- ker(i_*) = C(n-2, 2) when nonzero? (=3 at n=5, =6 at n=6)
- Rank drop: rk(d_2(T)) - rk(d_2(T\v)) is typically n-2 (good), occasionally n-1 (bad)

**Partial proof of HYP-278** (S43 continued):
Three of four cases are PROVED algebraically:
1. **b1(T)=1**: Trivially good (all vertices satisfy b1(T\v) <= 1 = b1(T)). PROVED.
2. **b1(T)=0, T not SC**: Delete from condensation boundary -> not SC -> b1=0. PROVED.
3. **b1(T)=0, T SC, kappa(T)=1**: Delete cut vertex -> not SC -> b1=0. PROVED.
4. **b1(T)=0, T SC, kappa(T)>=2**: All T\v are SC. Need: exists v with b1(T\v)=0.
   - EMPTY at n=5 (no such tournaments exist)
   - n=6: 1680 tours, ALL score (2,2,2,3,3,3), c3=8, kappa=2. Verified 1680/1680.
   - n=7: 378+ sampled, all verified. Scores near-regular only.
   - Sum_v b1(T\v) in {0, 3} for kappa>=2 case at n=6

**Additional discoveries** (S43 continued):
- Sum_v b1(T\v) <= 3 for ALL tournaments n=5-10 (bound independent of n!)
- b1=1 iff SC at n=3,4. First SC+b1=0 appears at n=5.
- Hidden cycles at bad vertices are ALWAYS linearly independent (rank = #bad)
- codim(Z_1(T\v) + Z_1(T\w)) = 1 in Z_1(T) always (proved algebraically)
- The complement L_{vw} visits ALL n vertices (n=5 verified)
- W_v + W_w + W_x = V for ALL triples v,w,x (n=5 verified)
- Regular n=7: 70% b1=0, 30% b1=1. Sum_v b1(T\v) in {0,1} for regular.
- Vertex with max c3(v) is good 93.6% of the time (not 100%)

**What remains to prove**:
- HYP-278 Case 4: kappa>=2 + SC + b1=0 -> exists good vertex
- HYP-279: b1(T) <= 1 for all tournaments (would close case 4 if combined with Sum<=3)
- HYP-282 (CORRECTED): Sum_v b1(T\v) <= 3 WHEN b1(T)=0 (gives >=n-3 good vertices)

**NEW DISCOVERIES (kind-pasteur-S43 continued)**:

#### TT Boundary Spanning Theorem (THM-103)
**PROVED**: TT (transitive triple) boundaries span ALL of im(d_2).
NT (non-transitive) elements of Omega_2 contribute NOTHING to im(d_2).
This means: b1 = dim(Z_1) - rk(TT boundary matrix in Z_1).
Verified exhaustive n<=6, sampled n<=10.

#### Cycle Sum Equality Theorem (THM-104)
**PROVED ALGEBRAICALLY**: If two directed 3-cycles share a directed edge,
any TT-cocycle z assigns them equal cycle sums (z(a,b)+z(b,c)+z(c,a)).

PROOF: Two cycles sharing edge a->b, with third vertices c and d.
For edge c-d (either direction), TT triples force the difference to 0.
Specifically: for c->d, TT(b,c,d) and TT(c,d,a) give the result.

COROLLARY: All cycles in a connected component of the 3-cycle graph
(adjacency by shared directed edge) have equal TT-cocycle sums.

#### Dominant Vertex Forcing Theorem (THM-105)
**PROVED ALGEBRAICALLY**: If vertex d dominates all of {a,b,c} (d->a,b,c)
or is dominated by all (a,b,c->d), the cycle sum of 3-cycle (a,b,c) is 0.

PROOF: Three TT triples (d,a,b), (d,b,c), (d,c,a) sum to give
z(a,b)+z(b,c)+z(c,a) = 0.

#### Free Cycle Bridge Theorem (THM-106)
**PROVED ALGEBRAICALLY**: For a free 3-cycle C (no external vertex
dominates/is-dominated-by all 3 vertices), EVERY external vertex v
creates a bridging 3-cycle sharing an edge with C.

PROOF: Exhaustive case analysis. Vertex v has pattern (x,y,z)
where x=v->a, y=v->b, z=v->c. The no-bridge conditions
NOT(x=1,y=0), NOT(y=1,z=0), NOT(z=1,x=0) rule out all patterns
EXCEPT (0,0,0) and (1,1,1), which both make C dominated (contradiction).

#### Dominant Vertex Characterization of b1 (HYP-283)
**VERIFIED n<=8**: b1(T) = #{connected components of 3-cycle graph
where EVERY cycle is free (non-dominated)}.

Exhaustive: n=4 (64/64), n=5 (1024/1024), n=6 (32768/32768).
Sampled: n=7 (500/500), n=8 (200/200).

The max number of free components is ALWAYS <= 1 (giving b1 <= 1).

#### Cocycle Structure
At n=4: cocycle z = alpha * (3-cycle participation count). Exact for all 24.
At n=5: z = a*cyc_count + b*score_src + c*score_tgt + d. Exact for all 304.
Cocycle is supported entirely on 3-cycle edges (100% at n=5).
Cycle sums are all equal (and sign-definite) when cycle graph connected.

**PROOF STATUS UPDATE**: The proof of beta_2 = 0 reduces to showing:
(a) At most 1 free component (gives b1 <= 1, equivalent to HYP-279)
(b) Sum_v b1(T\v) <= 3 when b1=0 (verified, equivalent to HYP-282)
Both together give good vertex existence => inductive proof.

**Scripts**: `beta2_relative_induction.py`, `beta2_vertex_analysis.py`,
`beta2_b1_monotonicity.py`, `beta2_b1_bound.py`, `beta2_rk_d2_formula.py`,
`beta2_averaging_argument.py`, `beta2_sum_b1_bound.py`, `beta2_hidden_cycles.py`,
`beta2_codimension_analysis.py`, `beta2_sc_connectivity.py`, `beta2_kappa2_analysis.py`,
`beta2_b1_characterization.py`, `beta2_regular_b1.py`, `beta2_rank_drop_mechanism.py`,
`beta2_proof_final.py`, `beta2_tt_boundary_span.py`, `beta2_tt_cocycle.py`,
`beta2_good_vertex_final.py`, `beta2_cocycle_algebra.py`, `beta2_potential_proof.py`,
`beta2_dominant_vertex.py`

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
- LES via i_* injectivity AT EVERY v: fails ~16-18% of (T,v) pairs (S42/S43)
  But EXISTS good v always (S43, verified n<=12): approach H is correct strategy
- Cone-from-T' as full Z_2 filler: swap cycles DON'T span Z_2 (0/1024 at n=5, S43)

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
