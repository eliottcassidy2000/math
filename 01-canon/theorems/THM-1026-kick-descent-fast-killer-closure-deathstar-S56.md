# THM-1026 — The kick descent (renumbered from 1015; klein first-pushed THM-1015 "interval survival"): the two-scale fast-killer stratum, and the reversal of the near-equal obstruction (death-star-2026-07-18-S56)

**Status:**
- **Lemma K (PROVED, trivially):** a single kick that lands the whole killer block in the safe band
  certifies `M(S) ≥ 1/14` (explicit rational witness).
- **Theorem (PROVED, `j ≤ 6` sufficiently-fast killers):** the kick always exists, so `M(S) ≥ 1/14`.
  Union bound: each fast killer blocks density `1/7` of the kick range; `j ≤ 6 ⟹` union `< 1`.
- **The near-equal reversal (VERIFIED):** near-equal killers' bad-kick-sets **overlap**, shrinking the
  union — so the *clustered-killer* case that defeats THM-1011's (BG-K) gluing is **closed** by the
  kick. The full THM-1011 obstruction battery (`[257,258,263]`, `[300,301,302]`, `[500,502,505,509]`,
  …) is closed by explicit witnesses.
- **Apex-7:** `j = 7 ⟹` union `→ 1` — the wall appears exactly here.

Advances the named-next open step of THM-1011 (the near-equal route). Source HYP-7305; the kick base is
THM-1010/THM-1002 (`M ≥ ρM'/(ρ+1)`), here iterated over a whole block instead of one element. Scripts:
`04-computation/lrc_cluster_lift_deathstar_S56.py`, `lrc_collapse_measure_deathstar_S56.py` (+`.out`).
WLOG positive speeds.

> **Convergence note.** klein's THM-1015 (interval survival, first-pushed 09:56) proves the same clustered closure `Σ1/k < L_P(7−r)`, `r<7`, WITHOUT a spacing hypothesis — subsuming the union-bound half here. This file's distinct content is the **kick-space view**, which yields the `q(K) = meas(F)` merge (THM-1027) and the near-equal reversal; kept as an independent derivation.

---

## Setting and Lemma K

`S = P ⊔ K`: the **core** `P` (`c = |P|` slow speeds) and the **killer block** `K` (`j = |K|` fast
speeds, `min K > 13·max P`). Let `μ = M(P) ≥ 1/(c+1)` (LRC(c+1), settled), attained at `t₀`.

**Lemma K.** If there is a kick `s` with `|s| ≤ (μ − 1/14)/max(P)` such that `‖k(t₀+s)‖ ≥ 1/14` for
every `k ∈ K`, then `M(S) ≥ 1/14`.
*Proof.* For `p ∈ P`: `‖p(t₀+s)‖ ≥ ‖p t₀‖ − |p·s| ≥ μ − max(P)·|s| ≥ μ − (μ−1/14) = 1/14`. So
`t = t₀+s` is `1/14`-lonely for all of `S`. ∎ (Explicit rational witness — a per-family certificate.)

## The existence theorem (`j ≤ 6`, fast killers)

**Theorem.** Let `s_max = (μ − 1/14)/max(P)`. The **bad-kick-set** of killer `k`,
`B_k = {s ∈ [−s_max, s_max] : ‖k(t₀+s)‖ < 1/14}`, is a union of intervals of width `1/(7k)` spaced
`1/k`, so `meas(B_k) ≤ 2s_max/7 + 2/(7k)`. Hence
```
meas(⋃_{k∈K} B_k) ≤ j·(2s_max/7) + (2/7)·Σ_{k∈K} 1/k.
```
If `j ≤ 6` and the killers are fast enough that `(2/7)Σ1/k < 2s_max(1 − j/7)` — i.e.
`Σ_{k∈K} 1/k < (μ − 1/14)(7 − j)/max(P)` — then `⋃B_k` does not cover `[−s_max, s_max]`, so a good
kick exists and `M(S) ≥ 1/14`.

The threshold `j = 7` is exact: at `j = 7` the density bound is `2s_max·(7/7) = ` the whole range —
seven fast killers can, in principle, tile the kick interval. **This is the apex-7 wall, seen in
kick-space.**

## The near-equal reversal (the point)

THM-1011 found that its two-block gluing (BG-K) **fails on near-equal killer blocks** — near-equal
killers have nearly coincident *bad sets*, producing long runs and large H-oscillation that the gluing
cannot absorb. The kick descent turns this on its head: near-equal killers have **overlapping
bad-kick-sets**, so `⋃B_k` is *smaller* than the union bound, not larger. Measured (core `{1..10}`,
kick range):

| `K` | union `meas(⋃B_k)` | `Σ meas(B_k)` | good kick? |
|---|---|---|---|
| `[257,258,263]` (2.3% spread) | 0.376 | 0.429 | ✓ |
| `[300,301,302]` (tight) | 0.372 | 0.481 (collapsed) | ✓ |
| `[200,300,450]` (spread) | 0.287 | 0.409 | ✓ |

So **near-equality is an asset for the kick, not an obstruction** — the exact reversal of THM-1011's
finding for the gluing route. Every family in the THM-1011 clustered-killer battery (and 2–4 killer
clusters, near-equal and spread, over cores `{1..8}…{1..11}`) is closed by an explicit witness
(`cluster_lift` .out), including families where `Σ1/k` exceeds the clean threshold above — there the
union bound is loose (short kick range / few periods) but the overlap makes the actual union `< 1`.

## Honest scope

- **Uniformly PROVED:** `j ≤ 6` sufficiently-fast killers (`Σ1/k < (μ−1/14)(7−j)/max(P)`) — closes the
  spread and the well-separated two-scale strata.
- **Closed per-family (rigorous certificates):** the entire THM-1011 near-equal obstruction battery,
  and `j = 7,8` near-equal/spread clusters over small cores — via the kick witness search, a decision
  procedure that outputs an explicit `t` with `min_v‖vt‖ ≥ 1/14`.
- **Residual:** a *uniform* theorem for moderately-fast near-equal blocks needs the moiré/overlap bound
  (the bad-kick-sets' correlation), not the union bound; and `j = 7` is the apex-7 wall. The genuine
  hard case remains the **single-scale comparable** family (no fast/slow split, so no block to kick) —
  where the kick has no room and `M = 1/14` is possible (the AP).

**Net.** The two-scale renormalization is now *exact on the fast-killer stratum*: a whole near-equal
killer block is lifted into the band by one kick, closing the clustered case THM-1011 left open. What
remains of "covering ⟹ M > 1/14" is the single-scale comparable core — the apex-7 wall, `j = 7` in
kick-space — which no peeling/kicking reaches, consistent with the whole session's boundary.

→ THM-1010/1011/1014, THM-1002/1000, THM-726/735, THM-995 (IX/X), HYP-7305, HYP-3901.
