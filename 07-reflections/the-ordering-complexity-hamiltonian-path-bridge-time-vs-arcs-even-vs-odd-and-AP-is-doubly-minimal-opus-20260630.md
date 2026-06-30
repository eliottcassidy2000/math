# The ordering-complexity ↔ Hamiltonian-path bridge (pushed, with a self-correction): both O(S) and H(T) count realized linear orders of n elements with a marked observer — the LRC realizes them by TIME (the runners' snapshot family, a closed walk on the permutohedron, EVEN because time-reversal t↔1−t closes it, O=Φ(n−1)=Farey for the AP) and the tournament by ARCS (the Ham paths, consistency with the orientation, ODD by Rédei) — so the duality is TIME↔ARCS, EVEN↔ODD, Farey↔OCF, family↔consistency; there is NO clean O=H (the parity blocks it), but the EXTREMAL MINIMA match: the AP is doubly MINIMAL (min M = covering-min, min O = fewest orders — correcting my earlier "max O / doubly extremal" error) ↔ the TRANSITIVE tournament (min H=1)

*opus-2026-06-30. Owner: push the ordering-complexity ↔ Hamiltonian-path bridge. Pushed — and it forced a
self-correction (the AP MINIMIZES O, not maximizes). The bridge is structural (time vs arcs), not a value
equality; the parity blocks O=H, but the minima match (AP ↔ transitive).*

## Self-correction first: the AP is doubly MINIMAL (not extremal)
Earlier I claimed the AP *maximizes* `O` (doubly extremal: min M, max O). **Wrong.** Recomputed `O(S)` over
`(n−1)`-sets of distinct positive integers (n=5):
| set | `O` |
|---|---|
| **AP `{1,2,3,4}`** | **6 = Φ(4) (the MIN)** |
| `{1,2,3,5}` | 10 |
| `{1,2,4,5}` | 10 |
| `{1,3,5,7}` | 14 |
> **The AP MINIMIZES `O`** — it has the smallest max speed `n−1`, hence the smallest Farey sequence `F_{n−1}`,
> hence `O = Φ(n−1)`, the minimum (verified: min `O` over `4`-subsets of `[1..8]` is `6` at the AP). So **the
> AP is doubly MINIMAL: min `M` (covering-min, hardest to escape) and min `O` (fewest realized orders) — the
> tightest AND simplest LRC config.** (My "doubly extremal" reflection is corrected to "doubly minimal.")

## The bridge: time vs arcs (the honest structure)
Both invariants count **realized linear orders of `n` elements with a marked observer**, by different
mechanisms:
| | `O(S)` (LRC) | `H(T)` (tournament) |
|---|---|---|
| realizes orders by | **TIME** — the runners move, snapshots as `t` sweeps | **ARCS** — orders consistent with the orientation |
| structure | a **closed walk on the permutohedron** (the snapshot cycle) | a **subset** of orders (the Ham paths) |
| parity | **EVEN** — time-reversal `t↔1−t` reverses each order, so the set is reversal-closed | **ODD** — Rédei (`H=I(Ω,2)=1+2·…`); reverse of a Ham path is `R(T)`'s, not `T`'s |
| value (AP / transitive) | `Φ(n−1)` (Farey length) | `1` (transitive) |
| number-theory | **Farey / totient summatory** `Φ` | **OCF** `I(Ω,2)` |
> **The duality is TIME ↔ ARCS:** the LRC unfolds orders dynamically (the runners, the snapshot walk); the
> tournament fixes them statically (the arcs, the consistency). The **observer** (the marked origin / vertex)
> is the fixed point of both. This is the deep bridge — the LRC is the *dynamic* order-realizer, the
> tournament the *static* one.

## No clean O=H — and why (the parity is the obstruction)
> `O(S)` is **EVEN** (the snapshot family is closed under time-reversal `t↔1−t`, which reverses each order);
> `H(T)` is **ODD** (Rédei). So **no tournament has `H = O`** — the parities never match. The reversal
> symmetry that the LRC HAS (time-reversal = the complement, `c→−c` from the inhomogeneous reframe) and the
> tournament LACKS (the reverse of a Ham path lives in `R(T)`) is exactly the even/odd obstruction. The
> tournament's irreducible `+1` (Rédei) is the parity defect; the LRC's reversal-pairing is the evenness.

## What DOES match: the minima (AP ↔ transitive)
- **AP** = min `O` = min `M` = the simplest LRC config (fewest orders, hardest to escape).
- **Transitive** = min `H` = `1` = the simplest tournament (one order, the linear order).
> So the bridge's clean correspondence is at the BOTTOM: **AP ↔ transitive**, both the minimal/simplest in
> their order-count. (And both are the "ordered" configs — the AP equally-spaced, the transitive linearly
> ordered; both time-reversal/complement self-symmetric = SC.) The matched minima, not the values, are the
> bridge's anchor.

## Status
- **Corrected (opus, honest):** the AP MINIMIZES `O` (not maximizes); doubly MINIMAL (min `M`, min `O`); my
  earlier `new-invariants` reflection's "doubly extremal / max O" is wrong and is superseded here.
- **The bridge (pushed):** `O` (TIME, even, Farey, permutohedron walk, family) vs `H` (ARCS, odd, OCF,
  consistency); both count observer-marked linear orders; no `O=H` (parity); minima match (AP ↔ transitive).
- **The duality:** TIME ↔ ARCS, EVEN ↔ ODD, Farey ↔ OCF — the LRC is the dynamic order-realizer, the
  tournament the static one, sharing the observer.
- **Open thread:** the LRC snapshot walk's structure (the allowable sequence for generic speeds; coincident
  crossings for the AP) and whether the *directed*/reversal-quotient `O/2` matches any tournament invariant.

Related: new-invariants-…ordering-complexity (CORRECTED here: min not max O), the-observer-on-the-tournament-
side (`H`=Ham paths), the-inhomogeneous-lrc-complement-reframe (time-reversal = complement = the parity),
reconciliation-…IS-the-cusp (AP ↔ transitive, SC); A002088 (Φ), Rédei; OPEN-Q-039/108.
