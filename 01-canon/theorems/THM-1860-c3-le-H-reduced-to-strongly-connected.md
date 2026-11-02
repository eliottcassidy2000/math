---
id: THM-1860
title: "Triangle count bound as an inherited OCF corollary; historical SCC reduction"
status: >
  PROVED as an immediate corollary of the inherited odd-cycle collection
  formula: H(T)>=1+2*c3(T), hence c3(T)<=H(T), for every tournament.
  The Session137 status repair corrects an obsolete open-base description;
  it is not a new theorem or a new independent proof of OCF.
  The historical SCC laws and formalized arithmetic kernel remain valid.
source: kind-pasteur-2026-07-21-S128c135 (owner: work the full WOWII loop for a long session)
depends_on:
  - THM-1845    # the tournament WOWII conjecture generator that proposed c₃ ≤ H
related: [THM-1830, THM-1580]
lean: 04-computation/lean/TournamentH7/TournamentH7/SumLeProd.lean
script: 04-computation/wowii_candidates_largen_kps_S128c135.py (+ .out), /tmp scbase check
---

# THM-1860 — current status correction

**Inherited OCF corollary, not open.** In `H(T)=sum_C 2^|C|` over vertex-disjoint odd-cycle collections, the empty collection contributes1 and each single directed triangle contributes2. All other terms are nonnegative. Therefore

    H(T)>=1+2*c3(T)>=c3(T).

This is already supplied by [the canonical OCF](THM-002-ocf.md); the 2026-11-02 status repair does not claim a new result. The earlier conditional reduction below is preserved as history, with its obsolete open-base language superseded by this paragraph. Its SCC product and sum laws and arithmetic kernel remain valid.

## Historical reduction and provenance


The WOWII generator (THM-1845) proposed `c₃ ≤ H` (3-cycles ≤ Hamiltonian paths). It survived
every test; here it is reduced to a clean base, with the arithmetic kernel formalized.

## The reduction (proved)

Decompose `T` into strongly-connected components `C₁ ⋗ C₂ ⋗ ⋯ ⋗ C_r` in the transitive
condensation order (every vertex of `Cᵢ` beats every vertex of `Cⱼ` for `i < j`).

- **`H(T) = ∏ᵢ H(Cᵢ)`.** A Hamiltonian path, once it leaves `Cᵢ`, cannot return (no back-arcs),
  so it visits all of `C₁`, then all of `C₂`, …; it is exactly a concatenation of a Hamiltonian
  path of each `Cᵢ`.
- **`c₃(T) = Σᵢ c₃(Cᵢ)`.** A directed 3-cycle needs a directed cycle, so all three vertices lie
  in one component.

Singleton components have `(c₃, H) = (0, 1)`; non-trivial components are strongly connected with
`H(Cᵢ) ≥ 2` (in fact the minimum over strongly-connected tournaments is `3, 5, 9, 15, 25` for
`n = 3..7`, by Moon's vertex-pancyclicity). So, given the per-component base `c₃(Cᵢ) ≤ H(Cᵢ)`,

> `c₃(T) = Σᵢ c₃(Cᵢ) ≤ Σᵢ H(Cᵢ) ≤ ∏ᵢ H(Cᵢ) = H(T)`,

the middle inequality being **sum ≤ product for a list of naturals each ≥ 2** (singletons drop
from both the sum and the product). ∎ (modulo the base)

## The residual base

All that remains is **`c₃ ≤ H` for strongly-connected tournaments** — verified exhaustively
`n ≤ 5`, on tens of thousands of strongly-connected samples at `n = 6, 7`, and to `n = 23` on
the high-`c₃` families (Paley, rotational, near-regular), with zero violations. It is where the
inequality has the most room: a strongly-connected tournament is vertex-pancyclic (Moon), so it
has many Hamiltonian paths while `c₃ = O(n³)`.

## The formalized kernel

`SumLeProd.sum_le_prod` (Lean 4 / Mathlib, sorry-free): a list of naturals each `≥ 2` has
`sum ≤ prod`. This is the exact arithmetic step of the reduction, and it is the WOWII "formalize"
discipline applied to the tournament graffiti — the step the fleet's directed-WOWII work
(klein-S397) and prototypes had not yet done.

## The loop, closed

The generator (THM-1845) proposed the pair `c₃ ≤ H` and `H ≤ 2^{n−2}c₃+1`. Pushing both past the
`n ≤ 7` filter:

- **`H ≤ 2^{n−2}c₃+1` is refuted at `n = 10`** — a WOWII *off-by-scale* bound, tight at the
  transitive tournament (`H = 1 = 2^{n−2}·0+1`) but beaten once `H ∼ n!/2ⁿ` outruns `2ⁿ·c₃`
  (near-regular-10: `8767 > 8193`; Paley-11: `95095 > 28161`).
- **`c₃ ≤ H` survives** and reduces to the strongly-connected base above.

That is the whole WOWII loop on a self-proposed pair: **generate → filter → push past the filter
→ one refuted with an explicit witness, one proved-modulo-a-clean-base and formalized.**

## Named next

- **The strongly-connected base `c₃ ≤ H`.** Likely routes: an injection from 3-cycles into
  Hamiltonian paths through a fixed median order, or a Moon-style pancyclicity count giving
  `H ≥ ` (something `≥ c₃`).
- Feed the refuted `H ≤ 2^{n−2}c₃+1` back to the generator as a *lower* bound search: what is the
  correct `H`-vs-`c₃` envelope?
