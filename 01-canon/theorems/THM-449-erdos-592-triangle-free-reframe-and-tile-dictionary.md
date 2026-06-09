# THM-449: Erdős 592 — the triangle-free reframe and the tile/shuffle dictionary

**Status:** STUB — reserved by mac-mini-2026-06-09-S1 (content this session)
**Source:** mac-mini-2026-06-09-S1

## Scope (what this file will contain when filled)

Erdős Problem 592 ($1000, OPEN): characterise countable ordinals α such that α→(α,3)².

Planned content, each part marked PROVED only when a proof is written here:
1. (folklore, proof to be included) α→(α,3)² ⟺ every triangle-free graph G on vertex set α
   has an independent set of order type α.
2. (dictionary, to prove) Tournament form: identifying a 2-coloring of [α]² with a tournament
   on α (red = arrow agrees with ordinal order), α→(α,3)² says: every tournament on α contains
   an order-faithful transitive subtournament on a set of type α, or an order-reversed
   transitive triple.
3. (Galvin's observation, proof to be included) If α is additively decomposable then the
   bipartite split coloring witnesses α↛(α,3)². Z₂-grading view.
4. (identification, to verify computationally) Larson's interaction-scheme forms on
   U={(a,b):a<b<ω} ≅ ω² = the staircase tile-pair geometry classes of this repo.

## Known results recorded (literature, not ours)
Specker: ω²→(ω²,m)², ω^n↛(ω^n,3)² for 3≤n<ω. Chang/Milner: ω^ω→(ω^ω,m)².
Galvin–Larson: β≥3 with property ⟹ β=ω^γ. Schipperus: ω^(ω^γ)→(ω^(ω^γ),3)² for γ a sum of
≤2 indecomposables; ↛ for ≥4. OPEN: exactly 3 summands (smallest ω^(ω³)).
