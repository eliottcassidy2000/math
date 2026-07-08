---
source: kind-pasteur-2026-07-07-S76
status: THEOREM (THM-658) + a proved sandwich + a well-verified conjecture. A cross-domain
  convergence: the circular-coloring "linearization defect" coincides exactly with the
  Haralambis mu>M separation locus.
tags:
  - lonely-runner
  - LRC14
  - circular-chromatic-number
  - distance-graphs
  - homomorphism-ladder
  - Haralambis
  - cross-domain
---

# The linearization defect is the μ>M locus

**kind-pasteur-2026-07-07-S76 (THM-658).** Three research strands in this project turned out
to be the same fact seen from three sides. Writing it down because the coincidence is the
content.

## The three strands

1. **opus-S141's homomorphism ladder.** LRC reformulated as graph theory:
   `LRC(14) ⟹ GRAPH-14 (χ_c ≤ 14 for 13-generator distance graphs) ⟹ MOTZKIN-14 (μ ≥ 1/14)`.
   The named open question: is `χ_c(G_S) = 1/M(S)` *identically* (linearization)? If yes, LRC
   *is* graph theory; if no, the "linearization defect `1/M − χ_c`" is a graph-theoretic
   location for the moat.

2. **opus-S145's THM-652.** At the Goddyn–Wong tight instance `GW = {1..11,13,24}`:
   `χ_f = 13`, `χ = 14 = 1/M`, and `χ_c ∈ (13,14]` left open — the "decisive rung question."

3. **Haralambis 1977's `μ > M` phenomenon.** Some distance sets have Motzkin density `μ`
   strictly above the loneliness `M` (a "separation"); `{1,3,4,5}` is the smallest, GW is a
   `|S|=13` example, Lucas `{1,3,4,7}` another.

## The single fact

For any finite `S`, two universal bounds sandwich the circular chromatic number:

> `1/μ(S) = χ_f(G_S) ≤ χ_c(G_S) ≤ 1/M(S)`.

The left is vertex-transitivity (`χ_f = 1/`independence-ratio). The right is the **linear
coloring** `c(x) = a·x mod N` where `a/N` is the loneliness witness — the rotation coloring,
the direct image of an LRC witness in the graph world. Since `μ ≥ M` always, the sandwich is
nonempty, and it collapses precisely on the `μ = M` locus:

> **`χ_c(G_S) = 1/M(S)  ⟺  μ(S) = M(S)`.**

The `μ = M ⟹ χ_c = 1/M` half is **proved** (squeeze). The `μ > M ⟹ χ_c < 1/M` half is a
**conjecture**, verified on 11 instances with zero counterexamples — and at GW I built the
explicit witness (THM-658): a quasi-periodic `(27,2)`-coloring with two alternating
color-increments `{9,16}`, holding an effective gap `2/27 > 1/14 = M`, something **no single
rotation can do**. So `χ_c(G_GW) ≤ 13.5 < 14`.

## Why the three strands are one

- **Strand 1's linearization defect = strand 3's separation.** `χ_c < 1/M` (defect) happens
  exactly when `μ > M` (separation). The graph-theoretic "moat location" opus asked for *is*
  the Haralambis locus. The defect is not mysterious — it is the density slack `μ − M` given
  a coloring interpretation.
- **Strand 2 (GW) is the flagship instance.** GW is a `μ > M` set (`μ = 1/13 > M = 1/14`), so
  by the law it must have `χ_c < 1/M = 14`. The rung question was never open in spirit — it
  was forced by the separation, once the sandwich is seen.
- **The mechanism is linearization vs. non-linearization.** `1/M` is the best *single-frequency*
  (rotation) gap. `χ_c` allows *variable-speed* colorings. When `μ = M` there is no slack to
  exploit and rotation is optimal (`χ_c = 1/M`); when `μ > M` the extra independent-set
  density is exactly the room a variable-speed winding needs to beat every rotation.

## What it means for LRC(14)

A quiet but useful negative: **the graph-coloring reformulation cannot prove LRC(14) via
linearization.** `GRAPH-14` (`χ_c ≤ 14`) is a valid *consequence* of LRC but is strictly
weaker at every `μ > M` tight instance (GW witnesses a `0.5` gap), so `χ_c = 1/M` — the
identity that would make "LRC = graph theory" — is false off the `μ = M` locus. The density
floor (`μ_{1/7}`) work is on the right object; the circular-chromatic surrogate throws away
the tight instances. This redirects the graph route from "prove LRC" to "understand the
`μ > M` locus," which is now the *same* problem as the linearization defect.

## Ledger

- Proved: the sandwich `1/μ ≤ χ_c ≤ 1/M`; the easy half `μ = M ⟹ χ_c = 1/M`; the GW witness
  `χ_c(G_GW) ≤ 27/2 < 14` (THM-658, verified certificate).
- Conjectured (11/11, 0 counterexamples): `χ_c = 1/M ⟺ μ = M`; equivalently the defect half
  `μ > M ⟹ χ_c < 1/M`.
- Cross-domain: the circular-coloring linearization defect (opus-S141) = the Haralambis
  `μ > M` separation locus = the room a variable-speed coloring needs. One fact, three names.
- Open: the defect half in general (a construction from a `μ > M` witness to a sub-`1/M`
  coloring — GW's two-speed winding is the template); the exact `χ_c(G_GW) ∈ (13, 13.5]`.
- Files: `lrc_chic_gw_quasiperiodic`, `lrc_chic_gw_sat`, `lrc_chic_linearization_locus`
  `_kps_S76`; THM-658.
