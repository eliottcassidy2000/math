# HYP-2250 — Partition functions everywhere: the deletion-contraction recursion = the transfinite game value (Hamkins)

**Session:** claudebox-2026-06-03-S625. **Prompt (user):** Hamkins on transfinite game values / infinite Go; the
"n+2 recursion" connections; make math improvements; come to see partition functions everywhere. **Threads:** HYP-2245
(H=I(Ω,2) partition function), HYP-2200 (depth-GF), HYP-2180 (Collatz altitude = iterated-log), HYP-2215 (Delsarte LP).

## The one recursion, three semirings
The **deletion-contraction recursion** of the independence polynomial — the hard-core partition function —
`I(Ω,x) = I(Ω∖v,x) + x·I(Ω∖N[v],x)` is the **universal partition-function transfer**. For a path/chain conflict
graph it is the **n−2 Fibonacci recursion** `I_n = I_{n−1} + x·I_{n−2}` → **Jacobsthal** `1,3,5,11,21,43,85,…` at
`x=2` (verified; `(2^{n+1}−(−1)^{n+1})/3`). Read over different semirings it is every "value-by-recursion" object:

| semiring | recursion | object |
|---|---|---|
| `(+, ×)` ordinary | `Z = Σ x^{size}` | independence polynomial / H / depth-GF / unit-distance count |
| `(min, +1)` / `(sup,+1)` tropical | `val = min/sup over moves (+1)` | **Hamkins's transfinite game value** (infinite Go/chess) — the well-founded rank; shortest path |
| ordinal | well-founded recursion | the open-game ordinal value; the Collatz altitude |

**So "partition functions everywhere" is literal:** the count, the game value, and the rank are the SAME recursion
over `(×,+)`, `(min,+1)`, and the ordinals.

## Hamkins → our problems
- **Collatz is an open game.** Position `n`, move = shortcut step, won = reach 1. Its **game value = the well-founded
  ordinal rank = the altitude tower** (HYP-2180): floor-1 (`log n` steps) and floor-2 (`loglog n` epochs) say the
  ordinal is `ω`-ish, not finite. Collatz conjecture = every position has a countable ordinal value = the game is
  open (no infinite play) = the recursion is **well-founded**. (Game-value non-termination = a cycle/divergence.)
- **LRC / unit distance / H** are the `(×,+)` face of the same recursion: the covering-depth GF, the unit-distance
  count, the tournament `H = I(Ω,2)` (HYP-2245) all obey deletion-contraction; their obstructions (collapse `p₀=0`,
  forbidden `H = 7·3^k`) are spectrum constraints on the partition function.
- **The Jacobsthal 21 vs the forbidden H 21.** `21 = I(P_6, 2)` (a path conflict graph) is achievable as an
  independence-polynomial value, but `H = 21` is FORBIDDEN for tournaments (THM-079): not every conflict graph `Ω`
  arises from a tournament's 3-cycle structure. The tournament constraint carves the gaps out of the full
  partition-function spectrum.

## The improvement
The deletion-contraction recursion is the engine to formalize/exploit: it gives (a) the free-vertex `(1+x)` transfer
⟹ the `3^k` baseline (last session, `indepPoly_edgeless`), (b) the component multiplicativity ⟹ the `7·3^k`
propagation (THM-079), (c) the Fibonacci/Jacobsthal chains, and (d) — over the ordinal semiring — the Collatz/open-game
rank. **One recursion to formalize, four problems served.** (The free-vertex case is formalized; the general
deletion-contraction `I(Ω) = I(Ω∖v) + x·I(Ω∖N[v])` is the next Lean target.)

## Formalized (math-lean, sorry-free, prior — `Math/Tournaments/IndependencePolynomial.lean`)
`indepPoly` (= the partition function), `indepPoly_edgeless`/`indepPoly_two_edgeless` (the free-vertex/`3^k`
baseline = the recursion's base case), `indepPoly_le_edgeless`/`indepPoly_two_le` (`H ≤ 3^k`).

## Open
- Formalize the general deletion-contraction recursion `I(Ω,x) = I(Ω∖v,x) + x·I(Ω∖N[v],x)` (attempted S625; the
  free-vertex `(1+x)` case is clean, the closed-neighborhood bijection is the work).
- The Collatz ordinal game value as a rigorous object: altitude (HYP-2180) = rank; well-foundedness = the conjecture.
- The tropical/ordinal ↔ ordinary semiring transfer as a single theorem (the "partition functions everywhere" law).
