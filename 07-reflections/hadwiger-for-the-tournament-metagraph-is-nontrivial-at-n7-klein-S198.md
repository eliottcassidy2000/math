# Hadwiger's conjecture is nontrivial for the tournament metagraph at n ≥ 7 — a K_{n−1}-minor prediction

**klein-2026-07-09-S198.** A genuine "consideration" of the Hadwiger conjecture on the project's
central object, the merged tournament metagraph `G_n/Z_2`.

## Why it is nontrivial here (and where)

Hadwiger's conjecture: **every graph `G` has a `K_{χ(G)}` minor** (equivalently, `K_t`-minor-free `⟹`
`(t−1)`-colorable). For a **perfect** graph it is trivial — `χ = ω`, so a `K_χ` *clique* (subgraph, hence
minor) exists. So Hadwiger only bites where `ω < χ`.

The merged metagraph `G_n/Z_2` (arc-reversal iso-class graph, complement-factored) has (opus-S314,
`chromatic-number-synthesis.md`):

| n | V | ω (clique) | χ | perfect? |
|---|---|-----------|---|----------|
| 4 | 3 | 3 | 3 | yes (K_3) |
| 5 | 10 | 4 | 4 | yes |
| 6 | 34 | 5 | 5 | yes (chordal) |
| **7** | **272** | **4** | **6** | **NO — ω < χ** |

At `n = 7` the graph **ceases to be perfect** (`ω = 4 < χ = 6`; it has odd holes). This is exactly where
Hadwiger becomes a real question — and it is a question about **the KEY object of the project** (CLAUDE.md).
Since the canonical conjecture `χ(G_n/Z_2) = n−1` has `χ` growing linearly while `ω` stays bounded (`≈ 4–5`),
the gap `χ − ω` grows, and Hadwiger predicts, for every `n ≥ 7`:

> **`G_n/Z_2` contains a `K_{n−1}` minor** — `n−1` disjoint connected branch-sets, pairwise adjacent —
> even though its largest *clique* is only `≈ 4–5`.

At `n = 7` the concrete, checkable claim is: **`G_7/Z_2` (272 vertices, `ω = 4`) has a `K_6` minor.**

## The candidate minor structure — the score ladder / SC spine

`G_n/Z_2` is a flip graph: an edge reverses ONE arc, which changes the **score sequence** by L¹-distance
`2` (two scores shift `±1`). The natural `(n−1)`-fold structure that a `K_{n−1}` minor would ride on is the
**score/`H`-gradient spine** (CLAUDE.md's *principal line*): the chain from the **transitive** class
(`H = 1`, score sequence `0,1,…,n−1`) through the self-complementary backbone to the **regular** class
(all scores `(n−1)/2`). The `n−1` distinct score-levels are the natural branch-set candidates; the SC-SC
"spine" edges + SC-NS "ribs" supply the cross-adjacencies a `K_{n−1}` minor needs. So the prediction is not
just numerical: **the minor, if it exists as Hadwiger says, should be the score-ladder contracted along
the SC spine.** Verifying this would identify the chromatic obstruction `χ = n−1` with a concrete
topological (minor) witness — turning the coloring bound into a *structural* theorem.

## The two readings of "Hadwiger" — and both are live in this project

1. **Hadwiger (minors), above** — nontrivial on `G_n/Z_2` for `n ≥ 7`; a clean dual-mandate lead.
2. **Hadwiger–Nelson (unit-distance / distance-graph chromatic number)** — already a live thread here
   (`alternating-group-graph-…`, `chromatic-number-of-the-plane-…`), and the natural bridge to the **LRC**
   side: the Lonely Runner is a distance/circulant-graph covering statement, and its "forbidden-distance"
   colorings are Hadwiger–Nelson-flavored. The two Hadwigers meet at the same place the project already
   lives — *the chromatic number of a Cayley/relation graph*.

## Concrete backlog lead

- **Build `G_7/Z_2`** (2²¹ tournaments → 456 iso-classes → 272 merged; adjacency = single-arc-reversal
  between classes) and **run a `K_6`-minor test** (contraction heuristic; `K_6` is fixed-size). Expected:
  a `K_6` minor exists (Hadwiger holds); the interest is its structure vs. the score-ladder/SC-spine
  prediction. If, against expectation, none is found, that is a Hadwiger counterexample — extraordinary,
  so verify triply.
- **Pattern check** `n = 8`: `χ = 7`, `ω = ?` (compute), Hadwiger `⟹ K_7` minor; does `χ − ω` keep growing,
  and does the score-ladder minor scale?

The methodological pairing with the same session's Mertens note is deliberate: **Mertens** (a *disproved*
belief about cancellation) steers the LRC near-resonance sum toward the cancellation-free route; **Hadwiger**
(a *believed, open* structural conjecture) is the lens under which the metagraph's chromatic number
`χ = n−1` should acquire a minor-witness. One warns against an assumption; the other proposes a structure.
