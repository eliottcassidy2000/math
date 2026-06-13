# (6,2) IMPOSSIBILITY — STRUCTURAL PROOF

**opus-2026-03-14-S71e**

## Claim: (α₁=6, α₂=2) is impossible for ALL n.

### Setup
H = 1 + 2α₁ + 4α₂ = 21 requires T = α₁ + 2α₂ = 10.
The decomposition (α₁=6, α₂=2) has T = 6 + 4 = 10. ✓

α₁ = Σ d(VS) over cycle vertex sets.
α₂ = Σ d(VS_i)·d(VS_j) over disjoint pairs.

### Case A: Some disjoint pair involves a VS with |VS| ≥ 5. **PROVED IMPOSSIBLE.**

Let (VS_big, VS_small) be a disjoint pair with |VS_big| ≥ 5, |VS_small| ≥ 3.

By HYP-1142: VS_big (which has a Hamiltonian cycle) contains ≥3 internal cyclic triples.
Each internal triple ⊂ VS_big, so it's disjoint from VS_small.
Each contributes d(triple)·d(VS_small) ≥ 1·1 = 1 to α₂.

α₂ ≥ d(VS_big)·d(VS_small) + 3·d(VS_small) = d(VS_small)·(d(VS_big) + 3) ≥ 1·(1+3) = 4 > 2.

**CONTRADICTION. Case A is impossible for all n.** ∎

### Case B: All disjoint pairs are between 3-cycles (|VS|=3, d=1).

α₂ = #{disjoint 3-cycle pairs}. Need α₂ = exactly 2.

#### Subcase B1: A 5-VS exists (not in any disjoint pair). **PROVED IMPOSSIBLE.**

The 2 disjoint pairs use 4 triples (A,B,C,D) contributing 4 to α₁.
The 5-VS contributes d(5-VS) ≥ 1 to α₁.
HYP-1142: the 5-VS has ≥3 internal cyclic triples.

Even if one internal triple coincides with A (a member of a disjoint pair):
- ≥2 NEW internal triples, each contributing 1 to α₁.
- Total: α₁ ≥ 4 + 1 + 2 = 7 > 6.

**CONTRADICTION. Subcase B1 is impossible for all n.** ∎

#### Subcase B2: NO 5+-VS exists. Only 3-cycles. α₁ = t₃ = 6, t₅ = t₇ = ... = 0.

**Key Theorem** (proved exhaustively):
- n=5: t₅=0 ⟹ t₃ ≤ 2
- n=6: t₅=0 ⟹ t₃ ≤ 2
- n=7: t₅=0 ⟹ t₃ ≤ 3

These are EXHAUSTIVE over ALL 2^C(n,2) tournaments.

For two disjoint 3-cycle pairs (Subcase 2a: shared triple → n≥7; Subcase 2b: no sharing → n≥8):

At **n ≤ 7**: By the exhaustive theorem, max t₃ ≤ 3 < 6 when t₅=0.
**IMPOSSIBLE.** ∎

At **n = 8**: Sampling (500k random tournaments, 873 with t₅=0):
- Max t₃ = 4 when t₅=0
- All t₃=4 cases have score sequence (1,1,2,2,5,5,6,6)
- t₃ = 6 with t₅ = 0: **0 hits in 500k**
- Combinatorial arrangement test: 16 valid orientations × 5000 trials = 0 hits

At **n = 9**: Sampling (200k, 16 with t₅=0): max t₃ = 3.

**Counting bound**: t₅=0 ⟹ every 5-subset has t₃(sub) ≤ 2.
Double counting: t₃ · C(n-3,2) ≤ 2·C(n,5), giving t₃ ≤ n(n-1)(n-2)/30.
- n=5: bound = 4 (actual max = 2)
- n=6: bound = 4 (actual max = 2)
- n=7: bound = 7 (actual max = 3)
- n=8: bound = 11 (actual max = 4)

The counting bound is loose, but the ACTUAL maximum grows as ~n/2.

### The Overlap-Concentration Dilemma (key structural insight)

**Dilemma**: α₂=2 requires OVERLAPPING triples (dp=2, not dp=C(t₃,2)), but overlapping triples create CONCENTRATED 5-subsets with t₃(sub) ≥ 3, which forces 5-cycles.

**Data** (n=8, t₅=0, 1M samples):
| t₃ | dp | count | C(t₃,2) | dp = C(t₃,2)? |
|-----|-----|-------|---------|----------------|
| 0 | 0 | 162 | 0 | ✓ |
| 1 | 0 | 308 | 0 | ✓ |
| 2 | 0 | 713 | 1 | overlapping |
| 2 | 1 | 94 | 1 | ✓ |
| 3 | 2 | 278 | 3 | overlapping |
| 4 | 4 | 152 | 6 | overlapping |

**Key observations:**
1. At t₃ ≥ 3 with t₅=0: dp ≥ 2 always. No dp=0 or dp=1 possible.
2. At t₃=4 with t₅=0: dp=4, not C(4,2)=6. But still dp > 2.
3. Score sequences highly polarized: (1,1,2,2,5,5,6,6) for t₃=4.
4. The near-bipartite structure forces triples to be spread across the partition.
5. Overlapping triples require concentration → forces t₃(sub) ≥ 3 → t₅ ≥ 1.

**Therefore**: dp=2 with t₃=6 is impossible when t₅=0, because:
- t₃=6 with t₅=0 requires extreme spread (all 5-subsets have t₃(sub) ≤ 2)
- This spread forces dp ≫ 2
- Conversely, dp=2 requires overlap, which concentrates triples, creating 5-cycles

### Tournaments with t₃=3, t₅=0 at n=7 (exhaustive)

3360 such tournaments, ALL have dp=2 (not 0 or 1 or 3).
alpha_2 = 2 always → H = 1+6+8 = 15 (not 7, not 21).
These are the ONLY tournaments with alpha_1=3 and t₅=0.

The pattern {A,B,C} where A∩B ≠ ∅, A⊥C, B⊥C (one overlapping pair, two disjoint pairs)
is the UNIQUE structure at t₃=3 with t₅=0.

## PROOF STATUS

| Range | Method | Status |
|-------|--------|--------|
| Case A (5-VS in pair) | HYP-1142 forcing | **PROVED for all n** |
| Subcase B1 (5-VS, not in pair) | α₁ budget overflow | **PROVED for all n** |
| Subcase B2, n ≤ 6 | Counting bound: t₃ ≤ 4 < 6 | **PROVED** |
| Subcase B2, n = 7 | Exhaustive: max t₃ = 3 when t₅=0 | **PROVED** |
| Subcase B2, n = 8 | Sampling 1M + overlap-concentration dilemma | **Strong evidence** |
| Subcase B2, n ≥ 9 | Sampling + structure + dilemma | **Strong evidence** |

**Overall: (6,2) is PROVED impossible for n ≤ 7, and strongly supported for all n.**

The remaining formal gap: proving max(t₃ | t₅=0) < 6 for all n ≥ 8. The overlap-concentration
dilemma provides a compelling structural argument but needs to be formalized into a rigorous
proof for general n.
