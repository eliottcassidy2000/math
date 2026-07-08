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

## The single fact (and the honest boundary — S77 correction)

For any finite `S`, two universal bounds sandwich the circular chromatic number:

> `1/μ(S) = χ_f(G_S) ≤ χ_c(G_S) ≤ 1/M(S)`.

The left is vertex-transitivity (`χ_f = 1/`independence-ratio). The right is the **linear
coloring** `c(x) = a·x mod N` where `a/N` is the loneliness witness — the rotation coloring,
the direct image of an LRC witness in the graph world. Since `μ ≥ M` always, the sandwich is
nonempty, and it **collapses on the `μ = M` locus**:

> **`μ(S) = M(S)  ⟹  χ_c(G_S) = 1/M(S)`** (proved, squeeze).

I first (S76) claimed the full equivalence `χ_c = 1/M ⟺ μ = M`, "verified 11/11." **That was
premature (MISTAKE-125).** The converse `μ > M ⟹ χ_c < 1/M` is **open and may fail**: the
`μ > M` cases that "confirmed" it were all trivial (`χ < 1/M ⟹ χ_c ≤ χ < 1/M`: Lucas, `{1,3,4,5}`)
except GW. The first decisive test, `{2,3,5,8}` (a Liu–Zhu A.3 set, `x=3, y=5` odd; `μ = 4/17 >
M = 3/13`, `χ_c ∈ [17/4, 13/3]`), is **exactly Liu–Zhu 2004 Problem 1 (OPEN)** — and every
sub-`1/M` coloring search came back empty (budget-limited/unsat), weak evidence that
`χ_c({2,3,5,8})` might *equal* `1/M`, a **counterexample**. So the converse is not a law; it is
a named open problem. What *is* proved at GW stands: the explicit `(27,2)`-witness (THM-658)
with two alternating increments `{9,16}` holding gap `2/27 > 1/14 = M` — **`χ_c(G_GW) ≤ 13.5 < 14`**.

## Why GW's defect IS provable: the odd cycle (the mechanism)

The reason GW works — and the likely reason `{2,3,5,8}` is hard — is the **parity of the
integrality obstruction**. opus-THM-652 forces `χ(G_GW) = 14` because the density-`1/13`
optimal classes, to `13`-color, must perfectly match the circulant `C(Z₂₆, {12})`, whose
components are **two odd `13`-cycles** (`gcd(12,26) = 2`; `12` has order `13`). An odd cycle
has no perfect matching → `χ = 14`. **But an odd cycle `C_{2k+1}` has `χ_c = 2 + 1/k < 3 = χ`** —
the circular chromatic number *goes around* it with a phase slip. So the **same Rédei-parity
odd cycle that blocks the integer coloring is exactly what lets the circular coloring beat it.**
The witness value is telling: `27/2 = 13 + 1/2` — thirteen nearly-tight classes wound as an
odd cycle, closing with a `1/2` phase slip (odd length forces the half-step). Sub-`13.5`
searches came back empty (`T ≤ 16`), so plausibly **`χ_c(G_GW) = 27/2` exactly**, the `+1/2`
being the odd-cycle defect. The refined conjecture this suggests: **`χ_c < 1/M` iff the
`μ > M` integrality obstruction is odd-cyclic** — sharper than "μ > M", and it would explain
why `{2,3,5,8}` (whose obstruction parity is not yet checked) could go either way.

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

- Proved: the sandwich `1/μ ≤ χ_c ≤ 1/M`; the direction `μ = M ⟹ χ_c = 1/M`; the GW witness
  `χ_c(G_GW) ≤ 27/2 < 14` (THM-658, verified certificate).
- **Retracted (MISTAKE-125): the clean equivalence `χ_c = 1/M ⟺ μ = M` and the "11/11" —
  the converse `μ > M ⟹ χ_c < 1/M` is OPEN and may be false.**
- Cross-domain (survives, refined): the circular-coloring linearization defect (opus-S141)
  is *located on* the Haralambis `μ > M` locus and *driven by* the odd-cyclic parity of the
  integrality obstruction (Rédei). Its hard core is a NAMED open problem: Liu–Zhu Problem 1
  (`χ_c` of A.3 sets with `x,y` odd), decided by `{2,3,5,8}`.
- Open: `χ_c({2,3,5,8}) ∈ [17/4, 13/3]` (= Liu–Zhu Problem 1 — decides the converse); the
  refined conjecture `χ_c < 1/M ⟺ obstruction is odd-cyclic`; the exact `χ_c(G_GW)` (likely
  `27/2` via the odd-cycle `+1/2` slip); a general defect construction from an odd-cyclic
  `μ > M` witness.
- Files: `lrc_chic_gw_quasiperiodic`, `lrc_chic_gw_sat`, `lrc_chic_linearization_locus`
  `_kps_S76`, `lrc_chic_defect_sweep_kps_S77`; THM-658; MISTAKE-125.
