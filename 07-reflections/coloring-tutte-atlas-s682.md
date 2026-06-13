---
source: claude-2026-06-06-S682
status: proved core (Ham-path deletion-contraction, verified exact) + the coloring/Tutte-Potts unification
tags: [tie-induction, deletion-contraction, tutte, potts, graph-coloring, chromatic, flow, tension, LRC, hadwiger-nelson, reframe, synthesis]
---

# Everything is a coloring; tie-induction is deletion–contraction

Prompt: pursue tie-induction; see everything as graph colorings — some nodes,
some edges, some both; reframe as many problems as possible. Both requests turn
out to be **one object**: the Tutte/Potts polynomial under deletion–contraction.

## Tie-induction = deletion–contraction (the spine)

A *tie* is a deleted edge. For the Hamiltonian-path count this gives an **exact**
recursion (verified, 2932/2932 at n=6):
```
H(T) = H(T ∖ e) + H(T / e).
```
A Ham path either **avoids** the directed edge `e=a→b` — counted by `H(T∖e)`, the
graph with `e` deleted (a *tie*) — or **uses** `a→b` consecutively — counted by
`H(T/e)`, with `a,b` merged (the merged vertex keeps `a`'s in-edges and `b`'s
out-edges). So **tie-induction is the Tutte deletion–contraction recursion**:
peel edges (introduce ties) until you reach an oriented/tie graph where the count
is trivial. This is exactly the recursion HYP-005 noted on the Rédei–Berge
function, now grounded as plain tie-induction.

## Some nodes, some edges, some both = chromatic / flow / full Tutte

The "node/edge/both" trichotomy is the Tutte polynomial `T(G;x,y)`:
- **nodes** colored → the **chromatic polynomial** (Tutte specialization on the
  `x`-axis): proper vertex colorings.
- **edges** carry a nowhere-zero label → the **flow polynomial** (`y`-axis).
- **both** → the full Tutte plane / the **Potts model** partition function.

Chromatic and flow are **planar-dual** (Tutte). That duality is the deep content
of the node/edge pairing.

## The coloring Rosetta (the repo's problems, reframed)

| problem | colored elements | coloring object | Tutte/DC tie |
|---|---|---|---|
| **LRC** | nodes (runners) → `ℝ/ℤ` | circular chromatic / a **tension**; planar dual = nowhere-zero **flow** (S537o) | chromatic↔flow duality; LRC@7 = circular chromatic # (Barajas–Serra) |
| **Tournament / H** | edges → `{→,←}` (tie = 3rd) | Ham-path count (Rédei–Berge) | `H = H∖e + H/e` (verified) |
| **metagraph `G_n/ℤ₂`** | nodes → `χ = n−1` | chromatic number | chromatic polynomial; perfectness breaks at `n=7` |
| **conflict graph `Ω`** | nodes → independent sets | `H = I(Ω,2)` = hard-core / Potts `Z` | Potts/Tutte specialization (= "partition functions everywhere", HYP-2183) |
| **unit distance** | plane points → colors | Hadwiger–Nelson `χ ∈ {5,6,7}` | chromatic of the plane; hexagonal **7**-coloring |
| **Collatz** | integers → parity | 2-coloring + dynamics | weak fit (a colored dynamical system) |

Three things click into place:

1. **LRC's two faces are Tutte's two axes.** S537o's "runner positions are a
   tension, the resonances are a nowhere-zero flow" is exactly the chromatic
   (node) ↔ flow (edge) duality. LRC sits on the diagonal where both meet.
2. **`H = I(Ω,2)` is a Potts partition function** (hard-core gas at fugacity 2).
   So "partition functions everywhere" (HYP-2183) is literally the Potts/Tutte
   specialization, and the H-spectrum is its evaluation set.
3. **The `7` recurs as a chromatic object.** Unit-distance's hexagonal upper
   bound is `7`; the metagraph's perfectness breaks at `n=7`; the forbidden
   H-value is `7 = |Fano plane PG(2,2)|`. The hexagonal `7` lives on the
   Eisenstein/prime-3 lattice (the `Cl₂(π/3)` carrier) — so the tournament-`7`
   and the plane-`7` plausibly share the prime-3 chromatic root. (Flagged, not
   proven — the kind of coincidence worth one careful look.)

## The gem

**Tie-induction and "everything is a coloring" are the same object: the
Tutte/Potts polynomial.** Deletion–contraction is tie-induction; its chromatic,
flow, and full specializations are node-, edge-, and both-colorings; LRC's
tension–flow duality is Tutte's chromatic–flow duality; `H`'s deletion–contraction
is the Tutte recursion; the partition-functions-everywhere thread is the Potts
`Z`; and the "useful tie-family" (HYP-2261) is precisely the **intermediate
oriented graphs of the DC recursion**.

This reframes the whole repo as **one question on the Tutte plane**: which
evaluations (chromatic number, flow number, `Z` at a fugacity, `H = Z_{hard-core}(2)`)
are achievable / extremal for the relevant graph families — and the impossibility
results (`{7,21}`, LRC bounds) are forbidden evaluations on specific (complete /
distance-graph) strata, dissolved by deletion (ties).

## Tie-induction for LRC (the concrete payoff)

LRC = a circular (node) coloring problem. Deletion = removing a runner (a tie
between it and the rest, or removing a constraint). The S679 sieve-apex is the
hard node; **deleting it** reduces LRC@n to a smaller circular-coloring problem
(17 runners with more room), the *deletion* branch; the *contraction* branch
glues the apex's constraints. The DC recursion for LRC is the partial-tournament
induction the prompt asked for — and it is the same shape as the `H` DC that just
closed H=21 (descend to indecomposable/atomic pieces).

## Next

1. **Make LRC's DC explicit:** define `M(V)` (loneliness radius) deletion (drop a
   speed) and contraction (identify two speeds / glue a resonance), and look for a
   monotone `M(V) ≥ min(M(V∖e), …)` that drives induction on ties.
2. **Forbidden evaluations on the Tutte plane:** recast `{7,21}` and LRC bounds as
   forbidden Tutte/Potts evaluations; the deletion (tie) move shows they are
   stratum-specific (complete graphs / distance graphs).
3. **The prime-3 `7`:** test whether the hexagonal `χ=7`, the metagraph
   perfectness break at `7`, and the forbidden H-`7` share the Eisenstein/`Cl₂(π/3)`
   root, or are independent `7`s.
