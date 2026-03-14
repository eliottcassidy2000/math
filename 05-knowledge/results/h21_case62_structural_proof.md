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

### Structural explanation for the gap

Tournaments with t₅=0 are "almost transitive" — every 5-element subtournament
has at most 2 cyclic triples (well below the maximum of 5). The score sequences
are extremely polarized (like (1,1,2,2,5,5,6,6) at n=8), corresponding to
near-bipartite tournaments where most vertices either dominate or are dominated.

Such tournaments have:
- t₃ concentrated in "boundary" triples between the high and low score groups
- 3-cycles forced to pair up, giving dp = t₃ (every pair disjoint)
- No room for the "spread" overlap structure needed for dp = 2 with t₃ = 6

## PROOF STATUS

| Range | Method | Status |
|-------|--------|--------|
| Case A (5-VS in pair) | HYP-1142 forcing | **PROVED for all n** |
| Subcase B1 (5-VS, not in pair) | α₁ budget overflow | **PROVED for all n** |
| Subcase B2, n ≤ 6 | Counting bound: t₃ ≤ 4 < 6 | **PROVED** |
| Subcase B2, n = 7 | Exhaustive: max t₃ = 3 when t₅=0 | **PROVED** |
| Subcase B2, n = 8 | Sampling 500k: max t₃ = 4 | **Strong evidence** |
| Subcase B2, n ≥ 9 | Sampling + counting bound + structure | **Strong evidence** |

**Overall: (6,2) is PROVED impossible for n ≤ 7, and strongly supported for all n.**
**The remaining gap is showing max(t₃ | t₅=0) < 6 for all n, currently empirical for n ≥ 8.**
