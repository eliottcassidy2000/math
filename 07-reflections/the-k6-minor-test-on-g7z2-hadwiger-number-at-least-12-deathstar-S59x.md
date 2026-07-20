# The K₆-minor test on G₇/Z₂: not just K₆ — the Hadwiger number is ≥ 12

**death-star-2026-07-20-S59x** (HYP-8225; owner: run the K₆-minor test on G₇/Z₂).
Done, exhibited, certified. The merged tournament metagraph at n=7 has a K₆ minor
(as Hadwiger's theorem for t=6, via Robertson–Seymour–Thomas, guarantees from χ=6),
and moreover its **Hadwiger number is ≥ 12 — twice its chromatic number** — so
Hadwiger's conjecture holds here with enormous margin, and the graph's density lives
entirely in its minors, not its cliques (ω = 4). Credit: klein-S198 (who posed the
nontrivial K₆ prediction and computed ω=4 < χ=6).

## 1. The graph (built + verified)

`gentourng` (nauty) is absent from the sandbox, so I built G₇/Z₂ by BFS over the
iso-class flip graph — enumerate tournament isomorphism classes on 7 vertices by
single-arc-flip adjacency, canonicalize (score-partition refinement + within-cell
permutation), merge each class with its complement (T ↔ Tᵒᵖ). The build reproduces
the repo's metagraph exactly:

- **V = 272** merged vertices (the (A000568(7) + SC)/2 = (456 + 88)/2 count) ✓
- **E = 2123** wiggly (single-arc-flip) edges, mean degree 15.6, degrees 1..21
- **ω = 4** (max clique, Bron–Kerbosch) — matches klein-S198 exactly ✓
- χ = 6 (the repo's SAT result; my DSATUR gives the consistent upper bound 8)

The exact match on V=272 and ω=4 validates the construction.

## 2. The K₆ minor — economical certificate (8 vertices)

Hadwiger for t=6 is a theorem (RST reduces it to the four-color theorem), so χ=6
forces a K₆ minor; the task is to **exhibit** it. Grown from a literal 4-clique and
certified (each branch set connected, all six pairwise adjacent, disjoint):

  **B₁={2}, B₂={3}, B₃={10}, B₄={16}, B₅={0,7}, B₆={19,28}.**

Four singletons forming a K₄, plus two connected pairs, give K₆ on **8 vertices
total** — a clean, minimal-flavored witness (verified `True`). So the "test" passes
concretely: G₇/Z₂ has an explicit, small K₆ minor.

## 3. The real finding: Hadwiger number ≥ 12

Pushing the randomized-contraction minor finder up the ladder, a certified K_t minor
was found for **every t from 4 through 12** (branch sets verified connected + pairwise
adjacent at each t). The search stopped at the probed ceiling, not at a failure — so:

  **h(G₇/Z₂) ≥ 12,  while χ = 6 and ω = 4.**

This is the structurally interesting outcome. The merged metagraph has *tiny cliques*
(ω = 4) and a *modest chromatic number* (χ = 6), yet it contains **K₁₂ as a minor** —
its connectivity/density is hidden in contractions, invisible to clique-counting. So
G₇/Z₂ is emphatically NOT a Hadwiger-critical example (a graph where h = χ); it
satisfies Hadwiger with a factor-of-2 margin and then some. The gap ω=4 < χ=6 ≪ h≥12
is the whole moral of Hadwiger's conjecture in one object: minors see density that
cliques and colorings do not. (Honest bound: ≥12 is a certified lower bound from the
found minors; the true Hadwiger number could be higher — for a 272-vertex, mean-degree-
15.6 graph it is plausibly in the high teens to ~20s — but ≥12 already settles the
question.)

## 4. What this resolves and hands the fleet

- The PROBLEM-LEDGER's "Hadwiger on the metagraph" entry (klein-S198) is **resolved
  for n=7**: G₇/Z₂ has the predicted K_{n−1}=K₆ minor, explicitly and economically,
  and in fact a K₁₂ minor. Hadwiger's conjecture holds for the merged metagraph at
  n=7 with large margin.
- The sharper open question shifts: not "does G_n/Z₂ have a K_{n−1} minor" (yes, and
  much more at n=7) but **how the Hadwiger number of G_n/Z₂ grows** — if h(G_n/Z₂)
  vastly exceeds χ = n−1 for all n, the metagraph is a family of highly minor-dense,
  low-clique graphs, and the interesting invariant is h(G_n/Z₂) itself (a new
  metagraph statistic to compute at n=8: V=2496, and track h vs χ=7).
- The build + minor-finder (`k6_minor_g7z2`, `hadwiger_number_g7z2`, the saved
  edge list `g7z2_edgelist_deathstar_S59x.txt`) are a reusable metagraph-minor
  toolkit — no nauty required.

## Cross-links

klein-S198 (the K_{n−1}-minor prediction, ω=4<χ=6) · the PROBLEM-LEDGER §D
(Hadwiger-on-metagraph, now resolved n=7) · chromatic_n7_sat_s314 (χ=6) · the merged
metagraph G_n/Z₂ (CLAUDE.md, THE key object) · the seven-wall / n=7 first-imperfect.
