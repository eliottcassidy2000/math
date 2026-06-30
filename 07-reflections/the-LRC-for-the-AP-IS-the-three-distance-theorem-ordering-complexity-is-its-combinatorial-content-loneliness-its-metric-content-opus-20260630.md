# Chasing the open thread: O/2 is NOT a tournament invariant (it is #Farey fractions in (0,½], a Farey count — honest negative, the parity obstruction is removed but the count stays Farey not OCF) — but the chase opened a real unification: the LRC for the AP {1,…,n−1} IS the THREE-DISTANCE (Sós–Steinhaus) theorem on {0,t,2t,…,(n−1)t}, with the ORDERING COMPLEXITY O=Φ(n−1) = #three-distance REGIMES (its COMBINATORIAL content) and the LONELINESS M=1/n at t=1/n = the optimal min-gap (its METRIC content); three-distance is the SHARED ENGINE already living in the project as the construction gaps {1,n,2n}, and the loneliness integral L(S) converges to ≈0.27

*opus-2026-06-30. Owner: chase the open thread (O/2 ↔ tournament invariant?) and any others. The asked thread
is a clean negative; chasing it surfaced that the LRC's snapshots are Sós permutations, so the whole LRC(AP)
is the three-distance theorem — combinatorics (O) and metric (M) both fall out of it.*

## The asked thread (O/2): honest NEGATIVE
The reversal-quotient `O/2` (undirected snapshot orders) is **not** a tournament invariant:
| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| `O/2` | 2 | 3 | 5 | 6 | 9 | 11 | 14 |
| score-seqs | 2 | 4 | 9 | 22 | 59 | 167 | … |
| tournaments A000568 | 4 | 12 | 56 | 456 | … | | |
> `O/2 = Φ(n−1)/2 = #Farey fractions in (0, ½]` — a **Farey count**. The `n=4` coincidence (`O/2 = 2 =`
> score-seqs) breaks at `n=5` (`3 ≠ 4`). So even with the parity obstruction removed (`O/2` is mixed-parity),
> the LRC's order-count is genuinely **Farey, not OCF**. The bridge to the tournament stays *structural* (time
> vs arcs), not a count-identity. **Negative recorded** — the LRC counts orders the Farey way; the tournament
> counts them the OCF way; they are different functions.

## The unification the chase opened: LRC(AP) = the three-distance theorem
The LRC for the AP is exactly Sós–Steinhaus on the orbit `{0, t, 2t, …, (n−1)t} mod 1`:
- **Three-distance verified:** every snapshot has **≤ 3 distinct gaps**, and (Steinhaus) the largest gap = the
  sum of the other two (checked n=6,8,10 at generic `t`).
- **`O = Φ(n−1)` = the COMBINATORIAL content** — the number of **three-distance regimes**: the Farey intervals
  of `t` on which the gap-pattern (hence the ordering) is constant. The orderings change exactly at the Farey
  fractions where the three-gap structure rearranges. So the ordering complexity *is* the count of
  Sós–Steinhaus regimes.
- **`M = 1/n` at `t = 1/n` = the METRIC content** — the loneliness `max_t min_{k≥1} ‖kt‖` is the optimal
  min-gap, attained at the **equally-spaced** `t = 1/n` (all gaps `1/n`); verified n=6,8,10,12 (optimal `t=1/n`
  exactly). The covering-min is the three-distance theorem's metric extremum.
> **So the LRC(AP) is the three-distance theorem wearing two hats:** its *combinatorics* is the ordering
> complexity `O = Φ(n−1)` (the regimes), its *metric* is the loneliness `M = 1/n` (the optimal gap). The Farey
> sequence indexes the regimes; the equally-spaced point realizes the extremum.

## Three-distance is the SHARED ENGINE (it was already here)
Three-distance already lives in the project as the **construction gaps `{1, n, 2n}`** (mac-mini's
three-distance synthesis, HYP-3702 taxonomy) — the gap-set of the non-extremal covering construction. The
chase shows three-distance *also* generates:
- the **ordering complexity** (the regimes, `O = Φ(n−1)`),
- the **covering-min** (the metric extremum, `M = 1/n`).
> So the same Sós–Steinhaus engine produces the construction's static gap-set `{1,n,2n}`, the LRC's dynamic
> ordering count, AND the covmin. Three-distance is not a side-fact — it is the **engine of the AP side** of
> the LRC, the time-domain counterpart of the tournament's OCF engine. (The bridge restated: **time-side =
> three-distance/Sós; arc-side = OCF/Rédei**.)

## Other thread (loneliness integral): converges to ≈0.27
`L(AP_n) = ∫ M_c dc`: `0.313, 0.288, 0.277, 0.271, 0.267` for `n = 6,10,14,18,22` — slowly decreasing toward
a constant `≈ 0.26–0.27` (limit not yet pinned; candidate `1/4`? still `> 0.25` at n=22). The **mean-observer
loneliness** is `Θ(1)` while the worst-observer `M_0 = 1/n → 0` — the LRC difficulty stays localized at `c=0`.
(Open: the exact limit of `L`.)

## Status
- **Negative (opus, honest):** `O/2` is a Farey count (`#Farey in (0,½]`), NOT a tournament invariant; the
  LRC↔tournament bridge is structural (time vs arcs), not a count-identity.
- **Unification (opus, verified):** LRC(AP) = the three-distance theorem; `O = Φ(n−1)` = #three-distance
  regimes (combinatorial); `M = 1/n` at `t=1/n` = optimal min-gap (metric); Steinhaus largest-gap law holds.
- **Connection:** three-distance is the shared engine — construction gaps `{1,n,2n}` + ordering complexity +
  covmin all from Sós–Steinhaus; the time-side engine dual to the tournament's OCF.
- **Open:** the exact limit of the loneliness integral `L` (≈0.27); whether the three-distance ↔ continued-
  fraction structure of `t` refines the regime count.

Related: the-ordering-complexity-hamiltonian-path-bridge (time vs arcs; O/2 thread resolved here),
new-invariants-…ordering-complexity (`O=Φ(n−1)`), SECOND-CORRECTION-…AP-scaled (`M=1/n`); mac-mini HYP-3702
(three-distance taxonomy), the construction `{1,n,2n}`; A002088 (Φ); OPEN-Q-108.
