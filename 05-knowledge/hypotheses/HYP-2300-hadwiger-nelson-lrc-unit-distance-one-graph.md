# HYP-2300 — Hadwiger–Nelson = Lonely Runner = unit distance: one graph, three questions, bound by χ·α ≥ n

**Session:** claudebox-2026-06-04-S634. **Prompt (user):** pursue the Hadwiger–Nelson / LRC / unit-distance
unification; see how they are the same underlying thing; how insights of one are the keys to the other. **Threads:**
HYP-2299 (colorings), HYP-2235 (unit/CM), HYP-2215 (Delsarte LP), HYP-2195 (lonely measure).

## One object: the unit-distance graph G of a homogeneous space
All three problems are questions about the **unit-distance graph** `G` (`x∼y ⟺ d(x,y)=1`) of a homogeneous space
(the plane for unit distance / HN; the circle/torus for LRC):
| question about `G` | the problem |
|---|---|
| max **edges** on `n` points | **unit distance** |
| **chromatic number** `χ(G)` of the whole space | **Hadwiger–Nelson** |
| max **independent set** `α(G)` (the 1-avoiding / lonely set) | **Lonely Runner** (the lonely set, `p₀`) |

## The unifying inequality χ·α ≥ n (formalized + verified)
A proper `k`-coloring partitions the `n` vertices into `k` independent sets, so the largest is an independent set `S`
with `k·|S| ≥ n`. Hence **`χ(G)·α(G) ≥ n`** — Hadwiger–Nelson's `χ`, the Lonely Runner's `α`, and the unit-distance
graph in ONE inequality. **Formalized** `Math/Combinatorics/HadwigerNelsonBridge.lean` `exists_indepSet_card_mul_ge`
(a `k`-coloring yields an independent set with `k·|S| ≥ n`). **Verified:** Moser spindle `χ=4, α=2` (`8≥7`); triangular
patches. So a **dense** unit-distance graph (small `α` — loneliness *hard*) **forces a large chromatic number** (HN),
and a **large** lonely set (LRC easy) **caps** the chromatic number below. The three problems are three corners of one
graph, all driven by its edge-density.

## The measurable bridge: χ_m = 1/(1-avoiding density), and p₀ is the same quantity
HN's measurable chromatic number `χ_m(ℝ²) ≥ 1/m₁`, where `m₁` = the max density of a **1-avoiding set** (no two points
at distance 1; `m₁ ∈ [0.229, 0.259]`, ⟹ `χ_m ≥ 5`, Falconer). The LRC **lonely set IS a 1-avoiding (forbidden-avoiding)
set** on the circle, and `p₀` (the lonely measure, HYP-2195) is **the m₁-analogue**. So *the 1-avoiding density of the
plane and the lonely measure of the circle are the same quantity* on two homogeneous spaces — and both are bounded by
the **same Delsarte/Fourier LP** (HYP-2215; Witsenhausen–Oler on the HN side, my covering-depth LP on the LRC side).

## The keys (insight transfer)
- **HN → LRC:** the Witsenhausen–Oler / Delsarte Fourier LP that bounds the 1-avoiding density is the same LP that
  bounds `p₀` — the autocorrelation of the forbidden region. Import the HN spectral bounds to LRC.
- **LRC → HN:** the resonance / additive-chain collapse (where `p₀=0`, HYP-2195) = where the 1-avoiding set is forced
  small = the unit-distance clusters (Moser-spindle rigidity) that force `χ` up. The LRC's "why loneliness fails" is
  the HN's "why χ is high."
- **unit distance → both:** the extremal configs. The **triangular lattice** (densest unit distances, 3-colorable =
  hexagonal, `δ=1/6` chord) is the **tame corner** of all three; the **CM-tower / non-lattice** are the **records**
  (the UD disproof, de Grey `χ≥5`) — both win by *leaving the lattice* (S631).
- **the π/3 thread:** the triangular UDG has `χ=3` (hexagonal 3-coloring), the densest unit distances, the `δ=1/6`
  chord; the cube-root / `Φ₃` is the shared `3` — the lattice's tame value on all three problems.
- **tie/coloring monotonicity (S633):** `G ≤ H ⟹ χ(G) ≤ χ(H)`; adding unit ties simultaneously raises `χ` (HN),
  shrinks `α` (LRC), and counts edges (UD) — one monotone family.

## Synthesis
Unit distance, Hadwiger–Nelson, and the Lonely Runner are three questions about one unit-distance graph on a
homogeneous space — max edges, chromatic number, max independent set — locked together by `χ·α ≥ n` and by the shared
Delsarte LP. The lattice (hexagonal, `χ=3`, `δ=1/6`) is their common tame extremum; the records on all three come from
non-lattice algebra (CM tower / spindle rigidity). The Lonely Runner's lonely set is the circle's 1-avoiding set; its
`p₀` is the plane's avoiding density; the perspective key (free vs fixed) is the lattice-vs-algebra axis on all three.

## Formalized (math-lean, sorry-free) — `Math/Combinatorics/HadwigerNelsonBridge.lean`
`exists_indepSet_card_mul_ge` (the `χ·α ≥ n` bridge).

## Open
- The Delsarte LP transfer made explicit: derive an LRC `p₀` bound from a Witsenhausen-style HN/1-avoiding LP (or vice
  versa) on matched homogeneous spaces.
- The resonance ↔ chromatic-cluster correspondence (additive chains ↔ the spindle-forcing subgraphs).
- Formalize `α(G) ≥ n/χ(G)` as the named corollary (and the measurable/fractional version).
