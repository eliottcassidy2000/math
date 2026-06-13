# THM-151: Explicit Formula for Omega_3

**Status:** VERIFIED (n=5 exhaustive)
**Session:** opus-2026-03-13-S70
**Depends on:** GLMY path homology definitions

## Statement

For a tournament T on n vertices:

$$\Omega_3(T) = \binom{n}{4} - \sum_{v} \left[ c_3(T|_{N^{-}(v)}) + c_3(T|_{N^{+}(v)}) \right]$$

where:
- N^-(v) = in-neighborhood of v, N^+(v) = out-neighborhood of v
- c_3(G) = number of 3-cycles in the subtournament induced on vertex set G
- The sum is over all vertices v ∈ V(T)

## Proof sketch

**Step 1:** Omega_3 counts regular 3-paths (v₀,v₁,v₂,v₃) with:
- Path: v₀→v₁→v₂→v₃
- Regularity: v₀→v₂ and v₁→v₃

**Step 2:** For each 4-vertex subset {a,b,c,d}, omega3_local = #regular 3-paths using
exactly these 4 vertices. By exhaustive computation over all 4 tournament types:

| Score sequence | c₃ (internal) | omega3_local | Type |
|---|---|---|---|
| (0,1,2,3) | 0 | 1 | Transitive |
| (1,1,2,2) | 2 | 1 | 4-cycle |
| (0,2,2,2) | 1 | 0 | Sink + 3-cycle |
| (1,1,1,3) | 1 | 0 | Source + 3-cycle |

**Observation:** omega3_local = 1 iff c₃ is even; = 0 iff c₃ is odd.

At n=4, c₃ ∈ {0, 1, 2} only, so "odd" means c₃ = 1.

**Step 3:** Omega_3 = Σ_{4-subsets S} omega3_local(S) = C(n,4) - #{S : c₃(T|_S) = 1}.

**Step 4:** A 4-subset S with c₃(T|_S) = 1 is either type (0,2,2,2) or (1,1,1,3):
- Type (0,2,2,2): vertex v dominated by all 3 others, which form a 3-cycle.
  The 3-cycle is a subset of N^-(v). Count: c₃(T|_{N^-(v)}) for each v.
- Type (1,1,1,3): vertex v dominating all 3, which form a 3-cycle.
  Count: c₃(T|_{N^+(v)}) for each v.

Each such 4-subset is counted exactly once (unique dominated/dominating vertex).

## Verification

- n=5, exhaustive (1024 tournaments): all match
- Omega_3 ∈ {3, 5} at n=5 → B ∈ {0, 2}
- At n=6: Omega_3 ∈ {9, 10, 11, 12, 15}, C(6,4)=15, B ∈ {0, 3, 4, 5, 6}
- Transitive tournament: all neighborhoods transitive → B=0 → Omega_3 = C(n,4)

## Corollaries

1. **Walsh degree 4:** The formula involves c₃ of neighborhoods, which depends on
   4-vertex substructure. In Walsh-Fourier terms, Omega_3 is a degree-4 polynomial
   in the edge variables (verified: only nonzero at Walsh degrees {0, 4}).

2. **Orthogonality to t₃:** corr(Omega_3, t₃) = 0 because Omega_3 depends on the
   DISTRIBUTION of 3-cycles among neighborhoods, not the global count t₃.
   Two tournaments with the same t₃ can have different neighborhood 3-cycle distributions.

3. **Comparison with Omega_2:** Omega_2 = C(n,3) - t₃ depends on GLOBAL 3-cycle count
   (Walsh degree 2). Omega_3 depends on LOCAL 3-cycle distribution (Walsh degree 4).
   They are Walsh-orthogonal.

## Notes

- This formula generalizes: Omega_m likely has a similar structure involving
  local counts of (m-1)-vertex cycle patterns.
- The even/odd parity of c₃ determining omega3_local is reminiscent of
  the sign in simplicial homology.
