---
id: THM-1860
title: "c₃ ≤ H (3-CYCLES ≤ HAMILTONIAN PATHS) IS PROVED MODULO THE STRONGLY-CONNECTED BASE, VIA THE SCC DECOMPOSITION — the strongest survivor of the WOWII tournament conjecture loop (THM-1845), and its arithmetic kernel is formalized in Lean. For every tournament T, the number of directed 3-cycles c₃ is at most the number of directed Hamiltonian paths H. (0) STATUS: verified exhaustively for all tournaments n ≤ 6, on 95661+ strongly-connected samples at n = 6,7, and to n = 23 on high-c₃ structured families (Paley, rotational, near-regular) — 0 violations. (1) THE SCC REDUCTION, PROVED: order the strongly-connected components C₁ ⋗ ⋯ ⋗ C_r transitively; a Hamiltonian path must traverse each Cᵢ completely before the next (no back-arcs), so H(T) = ∏ᵢ H(Cᵢ), while every 3-cycle lies inside one component, so c₃(T) = Σᵢ c₃(Cᵢ). Singletons contribute (c₃,H) = (0,1); non-trivial components are strongly connected with H(Cᵢ) ≥ 3 (verified min H over strongly-connected tournaments = 3,5,9,15,25 for n = 3..7). Hence, given the per-component base c₃(Cᵢ) ≤ H(Cᵢ), Σᵢ c₃(Cᵢ) ≤ Σᵢ H(Cᵢ) ≤ ∏ᵢ H(Cᵢ) — the last step being SUM ≤ PRODUCT for naturals each ≥ 2 (the trivial components drop out of both). (2) SO c₃ ≤ H FOLLOWS FROM THE STRONGLY-CONNECTED CASE c₃ ≤ H, which is the only remaining piece — verified n ≤ 7 and to n = 23, and where H is genuinely large (Moon vertex-pancyclicity gives many Hamiltonian paths). (3) FORMALIZED: the arithmetic kernel 'a list of naturals each ≥ 2 has sum ≤ product' is proved sorry-free in Lean 4 / Mathlib (SumLeProd.lean). (4) THIS IS THE WOWII LOOP CLOSING: the generator (THM-1845) proposed c₃ ≤ H; pushing past the n ≤ 7 filter refuted the sibling H ≤ 2^{n−2}c₃+1 (breaks at n = 10: near-regular-10 has H = 8767 > 8193, Paley-11 has 95095 > 28161 — a WOWII off-by-scale bound) while c₃ ≤ H survived to n = 23 and then reduced to a clean strongly-connected base"
status: >
  (0) VERIFIED as stated (exhaustive n ≤ 6; strongly-connected samples n = 6,7; structured
  families to n = 23; H ≤ 2^{n−2}c₃+1 refuted at n = 10 exactly).
  (1) PROVED: H(T) = ∏ H(Cᵢ) and c₃(T) = Σ c₃(Cᵢ) are the standard SCC-order facts (a
  Hamiltonian path is a concatenation of per-component Hamiltonian paths; 3-cycles are
  intra-component); the sum ≤ product step is proved and formalized.  H(Cᵢ) ≥ 2 for
  non-trivial components is Moon (strongly connected ⟹ Hamiltonian ⟹ ≥ n paths); verified
  min = 3,5,9,15,25.
  (2) The residual STRONGLY-CONNECTED base c₃ ≤ H is NOT proved — it is the whole remaining
  content, strongly evidenced (n ≤ 7 exhaustive/sample, n ≤ 23 structured, 0 violations) but
  open.  So this is 'c₃ ≤ H reduced to strongly-connected', not 'c₃ ≤ H proved'.
  (3) Lean: SumLeProd.sum_le_prod, sorry-free, builds under Mathlib v4.30.0.
  A candidate theorem with a proved reduction + formalized kernel; honest about the open base.
source: kind-pasteur-2026-07-21-S128c135 (owner: work the full WOWII loop for a long session)
depends_on:
  - THM-1845    # the tournament WOWII conjecture generator that proposed c₃ ≤ H
related: [THM-1830, THM-1580]
lean: 04-computation/lean/TournamentH7/TournamentH7/SumLeProd.lean
script: 04-computation/wowii_candidates_largen_kps_S128c135.py (+ .out), /tmp scbase check
---

# THM-1860 — c₃ ≤ H, reduced to the strongly-connected case

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
