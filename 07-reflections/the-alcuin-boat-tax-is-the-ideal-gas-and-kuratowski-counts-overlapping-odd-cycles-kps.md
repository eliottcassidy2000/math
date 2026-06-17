# The Alcuin boat-tax is the ideal gas, Kuratowski counts overlapping odd cycles, and a Hamiltonian path in Ω means H doesn't factor

**Source:** kind-pasteur-2026-06-16-S1. Dispatch: read the Alcuin-number paper (Csorba–
Hurkens–Woeginger, SIAM JDM 24(3) 2010, jstor 41642576); decide whether Alcuin(G) is the
vertex cover number τ(G) or τ(G)+1; relate the conflict graph to a tournament "in creative
ways"; ask what a Hamiltonian path / cycle in the conflict graph means; consider the tiling
map (present edge = transitive arc, absent = intransitive); think about Kuratowski and
Robertson–Seymour. All four sub-questions resolve through one object — the OCF lattice gas.

## The dichotomy, answered: it's the ideal gas

CHW prove `τ(G) ≤ Alcuin(G) ≤ τ(G)+1` for every graph (Observation 2.1): **small-boat**
(`=τ`) or **large-boat** (`=τ+1`), never anything else. Their Lemma 4.3 (two distinct
maximum independent sets ⟹ small-boat) and structure theorem 3.1 (the NP-certificate) decide
which. Specialized to the conflict graph `Ω(T)` (odd cycles, edges = shared `T`-vertex), the
dichotomy collapses to something the repo already owns:

> **`Alcuin(Ω(T)) = τ+1` ⟺ `Ω` is edgeless (`≥1` odd cycle, all pairwise vertex-disjoint) ⟺
> `H(T) = 3^{α₁}`** — the **ideal-gas / free-board maximum** (THM-517: `H≤3^{α₁}`, equality
> iff `Ω` edgeless). Otherwise `Alcuin(Ω)=τ`. (0 mismatches, all iso classes `n≤6`.)

The conceptual content is the prettiest part. An edgeless `Ω` is `m≥1` mutually conflict-free
items: vertex cover `τ=0`, but you still need **one boat seat** to ferry a nonempty
independent set across the river. That lone seat *is* the `+1`. The moment `Ω` has an edge,
the ferryman is already carrying a vertex cover, and that cover supplies the seat for free —
small-boat. So:

> **The Alcuin `+1` is the boat-tax for transporting a conflict-free set, and it is levied
> exactly in the OCF ideal-gas regime `H=3^{α₁}`** — precisely where the hard-core lattice gas
> (THM-002/510/517) has *no interactions* (excluded volume `|E(Ω)|=0`). Interactions pay the
> tax; the free gas does not get the discount.

This is a third independent meaning for `H=3^{α₁}` (after THM-517's "ideal gas" and "the
det-side"): it is the **river-crossing-hard / large-boat** tournament class. The boat needs
*more* capacity exactly when the cycle structure is *simplest* (no overlaps) — a satisfying
inversion.

## Kuratowski: K₅ is five mutually-overlapping odd cycles

When all odd cycles pairwise overlap (`α(Ω)=1`), `Ω` is the **complete graph `K_m`**
(`m=#odd cycles`). `K_m` is planar iff `m≤4`. So the **first non-planar conflict graph is
exactly `K₅` = five pairwise-overlapping odd cycles**, first realized at `n=5`, `H=11`
(`I(K₅,2)=11`); `H=9` (`K₄`) is the last planar one. Kuratowski's forbidden `K₅` is not an
analogy here — it is literally "5 odd cycles that pairwise share a vertex," and the planarity
of `Ω` in the `α=1` regime is governed by `H ≤ 9 ⟺ m ≤ 4`. (The other Kuratowski obstruction
`K_{3,3}` requires `α(Ω)≥2` — three odd cycles in each of two overlap-disjoint groups — and
shows up in the `α≥2` regime.)

## Robertson–Seymour: "Ω planar" is hereditary, with a tournament obstruction

Deleting a vertex `v` from `T` deletes every odd cycle through `v` from `Ω` — an *induced
subgraph* of `Ω`. Planarity is subgraph-closed, so **"`Ω(T)` planar" is closed under taking
sub-tournaments** — a hereditary tournament property. Its minimal obstruction is the smallest
`T` with `Ω` non-planar but `Ω(T−v)` planar for all `v`: the `n=5`, `H=11` tournament (`Ω=K₅`;
every vertex-deletion drops to `n=4`, planar). Robertson–Seymour's *graph-minor* WQO does not
transfer directly — tournaments are **not** well-quasi-ordered under the sub-tournament order
(infinite antichains exist, Chudnovsky–Seymour) — so the forbidden-sub-tournament set for
"`Ω` planar" may be infinite (OPEN-Q-106). But the *spirit* holds: a hereditary property with
a Kuratowski-driven obstruction lifted to tournaments. And CHW's payoff lands here: deciding
small-vs-large-boat is **polynomial for planar graphs**, so on the planar-`Ω` tournaments
(small/sparse-cycle ones) the Alcuin number of `Ω` is easy — even though computing `τ(Ω)`
(= `#odd cycles − ν_odd`) is hard in general.

## Hamiltonicity of Ω ⟺ overlap-connectivity ⟺ H doesn't factor

`Ω` has a **Hamiltonian path ⟺ `Ω` is connected** (verified `n≤6` + Paley `T₇`; automatic
when `α=1` since `Ω=K_m`). And `Ω` disconnected ⟺ `I(Ω,x)` factors over components ⟺
`H = ∏ H(\text{component})`. Therefore:

> **A Hamiltonian path in `Ω` certifies that `H` does not factor over odd-cycle-overlap
> components** — the cyclic structure is overlap-irreducible. Its *absence* (Ω disconnected)
> means `H` splits multiplicatively. The extreme — `Ω` edgeless — is the large-boat/ideal-gas
> case, where `H=3^{α₁}=∏_{α₁}3` is *maximally* factored (one block per odd cycle). A
> **Hamiltonian cycle** is the stronger "no single odd cycle is an overlap-bridge" condition.

So the three dispatch questions are one question: **Ω's connectivity spectrum** runs from
edgeless (large-boat, ideal gas, `H` fully factored, no Ham path) through disconnected-with-
clusters (`H` partly factored) to Hamiltonian-connected (small-boat, interacting gas, `H`
irreducible). Hamiltonicity sits at the connected/irreducible end.

## The G↦T_G tiling map and realizability

The user's map — present edge ⟹ transitive arc (`i→j`, `i>j`), absent edge ⟹ flipped (`j→i`)
— sends a graph `G` to a tournament `T_G`. It does **not** give `Ω(T_G)≅G`: `C₅↦` a
tournament with `Ω=K₅` (`H=11`), `P₅↦Ω=K₄` (`H=9`), while `K₅,K₄,K_{3,3}↦` transitive
(`H=1`, `Ω=∅`). The two directions `T↦Ω(T)` (lossy: many `T` per `Ω`, and not every graph is
an `Ω`) and `G↦T_G` (the tiling embedding) are genuinely different and don't compose to the
identity. *Which* graphs arise as conflict graphs is the **realizable-Ω cone** (THM-029:
`K₃` is non-realizable; THM-517's holes) — so the "creative correspondence" is real but
constrained, and the clean fact is the small-vs-large/ideal-gas one above, not a graph
isomorphism (HYP-2552).

## Status / honesty

PROVED all `n` in the `α(Ω)=1` regime (`Ω=K_m`; large-boat ⟺ `m=1`; Kuratowski `K₅` ⟺ `m≥5`).
VERIFIED exhaustively `n=3..6` + Paley `T₇`/regular `T₅` for the full statements. The
`α(Ω)≥2` small-boat direction is verified `n≤6`, conjectured all `n` (HYP-2550). The Alcuin
machinery is CHW's; the new content is the conflict-graph specialization, the `+1`=ideal-gas
identification (THM-517), the Kuratowski-`K₅` count, and the Ham-path⟺connected⟺non-factoring
correspondence. Cross-links: THM-519, THM-517 (`H≤3^{α₁}`, the pivot), THM-002/029,
THM-449/455 (strong-component `H=∏` — the other factorization), HYP-2550/2551/2552, OPEN-Q-106,
[[the-ocf-is-a-mayer-gas-lrc-is-its-conditional-twin-e7-is-its-fano-cousin-kps]] (the lattice
gas), [[the-free-gas-and-the-b2-atom]]. Source: Csorba–Hurkens–Woeginger, SIAM JDM 24(3)
(2010) 757–769.
