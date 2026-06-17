---
id: THM-519
title: The Alcuin number of the OCF conflict graph — the boat-tax +1 is exactly the ideal-gas regime (Alcuin(Ω)=τ(Ω)+1 ⟺ Ω edgeless ⟺ H=3^{α₁}); the Kuratowski threshold K₅ = five mutually-overlapping odd cycles; and Ω Hamiltonian-path ⟺ Ω connected ⟺ H does not factor over odd-cycle-overlap components
status: PROVED in the α(Ω)=1 regime for ALL n (Ω=K_m; large-boat ⟺ m=1 by Lemma 4.3; Kuratowski K₅ ⟺ m≥5); VERIFIED exhaustively over all tournament iso classes n=3,4,5,6 (0 mismatches: large-boat ⟺ Ω edgeless ⟺ H=3^{α₁}; α=1⟹Ω=K_m & planar⟺m≤4) + Paley T₇/regular T₅. The α(Ω)≥2 small-boat direction is verified n≤6, conjectured all n (HYP-2550).
source: kind-pasteur-2026-06-16-S1
depends_on:
  - THM-002   # OCF: H = I(Ω,2); the conflict graph Ω of directed odd cycles
  - THM-517   # H ≤ 3^{α₁} with equality iff Ω edgeless (the ideal-gas bound) — the centerpiece pivot
  - THM-029   # Ω=K₃ non-realizable (the realizable-Ω cone; bears on G↦T_G)
related:
  - HYP-2550  # general-n: large-boat ⟺ ideal gas (the α≥2 small-boat direction)
  - HYP-2551  # "Ω(T) planar" is a hereditary tournament property; Robertson-Seymour obstruction
  - HYP-2552  # the G↦T_G tiling map and conflict-graph realizability
  - OPEN-Q-106
  - T828
  - reflection: the-alcuin-boat-tax-is-the-ideal-gas-and-kuratowski-counts-overlapping-odd-cycles-kps
external: Csorba, Hurkens, Woeginger, "The Alcuin number of a graph and its connections to the vertex cover number", SIAM J. Discrete Math. 24(3) (2010) 757–769 (jstor 41642576).
---

# THM-519 — the Alcuin number of the conflict graph

**Background (Csorba–Hurkens–Woeginger 2010).** For any graph `G`, the **Alcuin number**
`Alcuin(G)` (min boat capacity for a feasible river-crossing schedule where each unattended
bank must induce a *stable/independent* set) satisfies `τ(G) ≤ Alcuin(G) ≤ τ(G)+1` (their
Observation 2.1; `τ` = vertex cover number). So every graph is **small-boat** (`Alcuin=τ`)
or **large-boat** (`Alcuin=τ+1`). Their Lemma 4.3: **two distinct maximum independent sets ⟹
small-boat.** Their structure theorem 3.1 is the NP-certificate. For the conflict graph
`Ω(T)` (vertices = directed odd cycles of `T`, edge = sharing a `T`-vertex): `H(T)=I(Ω,2)`,
`α(Ω)=ν_odd` = max vertex-disjoint odd-cycle packing `≤⌊n/3⌋`, `τ(Ω)=|V(Ω)|−α(Ω)` (Gallai)
`= (#odd cycles) − ν_odd`.

## A. The boat-tax +1 is exactly the ideal gas (centerpiece)

> **`Alcuin(Ω(T)) = τ(Ω)+1` (large-boat) ⟺ `Ω(T)` is edgeless with `≥1` vertex ⟺ all odd
> cycles of `T` are pairwise vertex-disjoint ⟺ `H(T)=3^{α₁}` (the ideal-gas / free-board
> maximum, THM-517, with `α₁≥1`). Otherwise `Alcuin(Ω)=τ(Ω)`.**

VERIFIED 0 mismatches over all iso classes `n=3..6` (11 large-boat classes, all edgeless `Ω`,
all `H=3^{α₁}`). **Conceptual reason:** an edgeless `Ω` (`m≥1` conflict-free items) has `τ=0`
but still needs **one boat seat** to ferry a nonempty independent set across — that lone seat
*is* the `+1`. When `Ω` has an edge, the vertex cover the ferryman already carries supplies
that seat for free (small-boat). So **the Alcuin `+1` is the "boat tax" for transporting a
conflict-free set, and it is levied exactly in the OCF ideal-gas regime** `H=3^{α₁}` — where
the lattice gas has no interactions (no excluded volume `|E(Ω)|`, THM-517). The `α(Ω)=1`
sub-case is fully general: `α=1 ⟺ Ω=K_m`; `m=1 ⟹` large-boat, `m≥2 ⟹` small-boat (the `m`
singleton max-independent-sets trigger Lemma 4.3).

## B. Kuratowski: K₅ = five mutually-overlapping odd cycles

When `α(Ω)=1` (all odd cycles pairwise overlap), `Ω = K_m` (`m=#odd cycles`). `K_m` is planar
iff `m≤4`. So:

> **The first non-planar conflict graph is `K₅` — five pairwise-overlapping odd cycles —
> first realized at `n=5`, `H=11` (`I(K₅,2)=1+2·5=11`).** `H=9` (`K₄`) is the last planar one
> in this regime.

The Kuratowski obstruction `K₅` is literally "5 odd cycles that pairwise share a vertex."
VERIFIED: 0 anomalies over 45 `α=1` classes (`Ω=K_m`, planar ⟺ `m≤4`); non-planar onset at
`n=5,H=11`; 43/56 classes have non-planar `Ω` at `n=6`; Paley `T₇` `Ω` (80 vtx) non-planar.

## C. Hamiltonicity of Ω ⟺ overlap-connectivity ⟺ non-factoring H

> **`Ω(T)` has a Hamiltonian path ⟺ `Ω` is connected** (verified `n≤6` + Paley `T₇`; forced
> when `α=1` since `Ω=K_m`). And **`Ω` disconnected ⟺ `H` factors:** `I(Ω₁⊔Ω₂,x)=I(Ω₁,x)I(Ω₂,x)`,
> so `H=∏ H(component)`. Hence a **Hamiltonian path in `Ω` certifies that `H` does NOT factor
> over odd-cycle-overlap components** (the cyclic structure is overlap-irreducible). A
> **Hamiltonian cycle** is the stronger "no overlap-bridge" condition (2-connected tour of all
> odd cycles). The maximally non-Hamiltonian `Ω` (edgeless) is exactly the large-boat /
> ideal-gas extreme of (A), where `H=3^{α₁}=∏_{α₁} 3` is fully factored into single-cycle blocks.

## Scope / honesty

PROVED for all `n` in the `α(Ω)=1` regime (the dominant small-`n` case): `Ω=K_m`, large-boat
⟺ `m=1`, Kuratowski `K₅` ⟺ `m≥5`. VERIFIED exhaustively `n=3..6` (all iso classes) + Paley
`T₇`/regular `T₅` for the full statements including `α≥2`. The `α(Ω)≥2` small-boat direction
(Ω-with-an-edge ⟹ small-boat, even with a unique maximum packing) is verified `n≤6` and
conjectured for all `n` (HYP-2550). The Alcuin/structure-theorem machinery is
Csorba–Hurkens–Woeginger's; the contributions here are the **conflict-graph specialization**,
the **`+1` = ideal-gas identification** (via THM-517), the **Kuratowski-K₅ = 5-overlapping-odd-cycles**
threshold, and the **Ham-path ⟺ connected ⟺ non-factoring-H** correspondence. The `G↦T_G`
tiling map (edge⟹transitive arc, non-edge⟹flip) does NOT satisfy `Ω(T_G)≅G` (e.g. `C₅↦` a
tournament with `Ω=K₅`, `P₅↦K₄`); most graphs are not conflict graphs (THM-029), see HYP-2552.
