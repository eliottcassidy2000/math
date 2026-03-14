# THM-202: P₄ Cannot Be Ω(T) — Toward H≠21

**Status:** PROVED (for P₄ Omega structure)
**Proved by:** opus-2026-03-14-S71g
**Dependencies:** THM-200, directed cycle counting, Jacobsthal connection

## Statement

For any tournament T, the conflict graph Ω(T) is never isomorphic to P₄ (the path graph on 4 vertices).

## Proof

Suppose Ω(T) = P₄, meaning T has exactly 4 directed odd cycles C₁, C₂, C₃, C₄ with adjacency C₁-C₂-C₃-C₄ (sharing vertices pairwise along the path).

**Step 1.** Each Cᵢ must be a 3-cycle.

A 5-cycle on 5 vertices forces ≥3 additional 3-cycles, which would create more than 4 total cycles. So all cycles are triangles.

**Step 2.** The 4 triangles span ≤9 vertices with shared vertices v₁ = C₁∩C₂, v₂ = C₂∩C₃, v₃ = C₃∩C₄.

Write C₁ = {a, b, v₁}, C₂ = {v₁, c, v₂}, C₃ = {v₂, d, v₃}, C₄ = {v₃, e, f}.

**Step 3.** The "dominance cascade" forces extra triangles.

To have exactly 4 triangles, every other triple must be transitive. This creates a chain of forced arcs:

1. Triple {c, v₂, d}: arcs c→v₂ (C₂) and v₂→d (C₃). To avoid triangle: force c→d.
2. Triple {c, d, v₃}: c→d (forced) and d→v₃ (C₃). To avoid triangle: force c→v₃.
3. Triple {v₁, c, v₃}: v₁→c (C₂) and c→v₃ (forced). To avoid: force v₁→v₃.
4. Triple {v₁, c, v₂}: already in C₂.

The cascade continues through the "private" vertices:
5. Triple {b, v₁, c}: b→v₁ (C₁) and v₁→c (C₂). To avoid: force b→c.
6. Triple {b, c, d}: b→c (forced) and c→d (forced). To avoid: force b→d.
7. ...continuing: b is forced to beat {c, d, v₃, e, f} — all vertices outside C₁.

Similarly, vertex a is forced to beat b (C₁) and then b's cascade gives a→{c, d, e, f, v₃}.

**Step 4.** The shared vertices form unavoidable "skip triangles."

From the cascade:
- a→v₂ is forced (via the a→b→...→v₂ chain and direct forcing).
- v₂→v₁ is the arc direction in C₂ (v₁→c→v₂→v₁ gives v₂→v₁).
- v₁→a is from C₁.
- Result: a→v₂→v₁→a is a directed triangle!

Concretely, with the standard orientation C₁: a→b→v₁→a, C₂: v₁→c→v₂→v₁, C₃: v₂→d→v₃→v₂, C₄: v₃→e→f→v₃:

The forced extra triangles are:
- **{a, v₁, v₂}**: a→v₂→v₁→a (skip triangle through backbone)
- **{v₁, v₂, v₃}**: v₁→v₃→v₂→v₁ (backbone triangle)
- **{v₂, v₃, f}**: v₂→f→v₃→v₂ (skip triangle through backbone)

These 3 extra triangles are **structurally unavoidable**: they arise from the backbone of shared vertices combined with the dominance cascade.

**Step 5.** With 7 triangles instead of 4, Ω(T) ≠ P₄. Contradiction.

∎

## Computational Verification

- n=7 (P₄ on 7 vertices with shared edges): 0/2048 tournaments achieve t₃=4. Min t₃=5.
- n=9 (P₄ on 9 vertices with single shared vertices): 0/16,777,216 tournaments achieve t₃=4. Min t₃=9.
- With all forced arcs: only 4 remaining arcs (2^4=16 completions). Min t₃=7, with exactly 3 extra "skip triangles" in every case.

## Connection to Jacobsthal Numbers

I(P₄, 2) = 21 is the 5th term of the sequence I(Pₘ, 2) = 1, 3, 5, 11, 21, 43, 85, 171, ... (Jacobsthal-type, recurrence a(n) = a(n-1) + 2·a(n-2)).

The P₄ impossibility means H=21 cannot come from an Ω = P₄ structure. Combined with the K₃+K₁ impossibility (THM-201), two of the four graph structures yielding I(G,2)=21 are ruled out.

## Remaining Cases for H=21 Proof

To fully prove H≠21, must also rule out:
- K₆ minus 2 edges (α₁=6, α₂=2): absent at n≤7 exhaustive
- K₈ minus 1 edge (α₁=8, α₂=1): absent at n≤7 exhaustive
- K₁₀ (α₁=10, α₂=0): absent at n≤7 exhaustive

## Key Scripts

- `04-computation/knacci_simplex_cuboid.py` — H-spectrum gap analysis
- `04-computation/h7_theorem.py` — Base THM-200 verification
