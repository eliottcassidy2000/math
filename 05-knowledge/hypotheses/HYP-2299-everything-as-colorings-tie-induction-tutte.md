# HYP-2299 — Everything as graph colorings: tie induction, the Tutte unification, and Hadwiger–Nelson

**Session:** claudebox-2026-06-04-S633. **Prompt (user):** pursue tie induction; see everything as graph colorings
(some nodes, some edges, some both); creatively reframe as many problems as possible. **Threads:** HYP-2290 (tie-graphs),
HYP-2265 (3-cycle/π-3), HYP-2250 (deletion-contraction), HYP-2235 (unit distance).

## The unifying object: the Tutte polynomial / deletion–contraction
**Tie induction = deletion–contraction**, the Tutte recursion `T(G) = T(G−e) + T(G/e)`. Every invariant of the arc is a
specialization: the **chromatic polynomial** (proper colorings), the **independence polynomial** (`H`, my partition
function), flow/reliability. "Some nodes, some edges, some both" = the **0/1/2-cells**: nodes (points/runners/numbers),
edges (distances/arcs/relations), both = **2-cells = 3-cycles / unit-triangles = the resonance** (the filled triangle).
The resonance lives at the 2-cells — and the 2-cell triangle is the minimal **odd cycle**, the obstruction to
2-coloring, the cube-root `π/3` (S265).

## Tie induction for colorings (formalized)
A **tie is a removed edge**; removing edges (`G ≤ H`) only makes coloring easier — a proper coloring of `H` restricts
to `G`. So **`χ(G) ≤ χ(H)`: adding ties only lowers the chromatic number** — the deletion side of deletion–contraction.
Formalized `Math/Combinatorics/ChromaticTieInduction.lean`: `Coloring.ofLE` (coloring restricts to a subgraph),
`Colorable.ofLE`, `chromaticNumber_mono`.

## The coloring reframe of every problem
| problem | the coloring |
|---|---|
| tournament `H` | color = vertex *order* (Ham path); 3-cycles = odd cycles = chromatic obstruction (`χ≥3`) |
| **unit distance** | color the PLANE so unit-apart points differ = **Hadwiger–Nelson**, `χ(ℝ²)∈{5,6,7}`; the UDG is the tie-graph — maximizing unit distances forces `χ` up |
| LRC | color the clock by nearest runner / by depth; lonely = an uncovered color; covering = a coloring |
| Collatz | color by **parity** (the 2-coloring) = the parity vector; even/odd = bipartite |
| Goldbach/Lemoine | even = unordered node-pair, odd = ordered arc = edge 2-coloring (S630) |
| forbidden `H={7,21}` | 3-cycle = `K₃` = `χ=3` (odd cycle); resonance = non-bipartite = chromatic obstruction |

## The 3-cycle = odd cycle = resonance = chromatic obstruction (verified)
The triangle `K₃` (the resonance atom, the 2-cell) has `χ=3`; even cycles are 2-colorable (bipartite), odd cycles
`χ=3`. So **resonance = odd cycle = non-bipartite = the chromatic obstruction**, the same `π/3` cube-root (3 colors).
The tight/collapse configs are where the odd cycles (resonances) concentrate.

## Hadwiger–Nelson: the unit-distance chromatic number (verified)
- **Moser spindle (n=7): `χ = 4`** (forces `χ(ℝ²) ≥ 4`; de Grey 2018: `≥5`) — the non-lattice rigid 4-chromatic UDG.
- **Triangular UDG: `χ = 3`** (verified n=7,12,19) — the **hexagonal 3-coloring**, the `π/3` cube-root structure (3
  colors = the 3 of the hexagonal lattice). So the lattice's tie-graph is 3-colorable; forcing `χ` higher needs the
  non-lattice rigidity (Moser/de Grey), the same lattice-vs-algebra gap as the disproof (S631). Maximizing the
  tie-graph (unit distances) and maximizing its chromatic number are the two extremal faces.

## Synthesis
"Everything as a coloring" = everything is a Tutte/deletion–contraction invariant of a graph whose 2-cells (triangles)
are the resonance. Tie induction (chromatic monotonicity, deletion) is the inductive engine; the chromatic obstruction
is the odd-cycle/3-cycle/π-3 resonance; the unit-distance chromatic number is Hadwiger–Nelson, where the lattice's
3-coloring (hexagonal) is beaten by non-lattice rigidity — the coloring shadow of the grid disproof.

## Formalized (math-lean, sorry-free) — `Math/Combinatorics/ChromaticTieInduction.lean`
`Coloring.ofLE`, `Colorable.ofLE`, `chromaticNumber_mono`.

## Open
- The contraction side: formalize the full deletion–contraction / a chromatic-polynomial recursion (the Tutte engine).
- The unit-distance chromatic number ≥ 4 (Moser spindle) in Lean (a specific UDG, `χ ≥ 4`).
- LRC depth-coloring as a proper coloring / the covering as a chromatic statement.
