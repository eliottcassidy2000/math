---
id: THM-678
title: The multi-detuned counting dispatch — for v = g·H ∪ D (|D| = d detuned), a good branch (u+c)/g exists whenever Σ_i N_i/q_i < 1 (N_i = orbit points of coordinate i in the open 1/7 danger arc, q_i = g/gcd(δ_i, g)); corollaries d = 2 dispatched UNLESS q_1 = q_2 = 2, d = 3 dispatched when all q_i ≥ 8; the d = 2 residual is EXACTLY the congruent half-harmonic pair (δ_1 ≡ δ_2 ≡ g/2 mod g), whose mod-2g lift is the 2-adic descent — the branch-side face of the descent-burden/half-sum-moduli structure
status: PROVED (counting argument, below; elementary). Machine-verified: the dispatch fires on adversarial d = 2/3 batteries exactly on the predicted region; the q=q=2 residual instances verified to need the lift; the blocking/supply duality measured (companion .out).
source: monad-explorer-2026-07-09-S8 (HYP-5787) — the S7 handoff (d ≥ 2 generalization).
depends_on: []   # LRC(≤13) citation for the harmonic time u, as in THM-668
related:
  - THM-668   # the d = 1 case (branch pigeonhole; formalized in LRCDetunedDispatch.lean)
  - THM-676 (mac-mini, descent-burden) / THM-671-B5 + THM-673 (klein, quintic Bonferroni + aggregated supply)  # the duality this connects (below)
  - kps-S18 LRCClusterGcd  # the mechanism's contrapositive ancestor
---

# THM-678 — the multi-detuned counting dispatch

## Statement

Let `v = g·H ∪ D` be 13 nonzero speeds: `g ≥ 2`, `H` the harmonic multipliers
(`|H| = 13 − d`), `D = {δ_1, …, δ_d}` with `g ∤ δ_i`. Let `u` be an LRC(14 − d) time for
`H` (clearance `1/(14−d) ≥ 1/13`, citation). For the branch times `τ_c = (u + c)/g`,
`c = 0, …, g−1`:

- every harmonic clearance equals `‖m·u‖ ≥ 1/(14−d)` at EVERY branch (integer shift);
- detuned coordinate `i` runs over an equally-spaced orbit of `q_i = g/gcd(δ_i, g) ≥ 2`
  values (each hit `g/q_i` times); let `N_i` = the number of its orbit values inside the
  open danger arc `(−1/14, 1/14)`.

> **If `Σ_{i=1}^{d} N_i / q_i < 1`, some branch `c*` has ALL detuned clearances ≥ 1/14,
> hence `M(v) ≥ 1/(14−d) > 1/14` and `v` is lonely.**

**Orbit-count bounds.** `N_i ≤ 1` when `q_i ≤ 7` (spacing `1/q_i ≥ 1/7` — an open arc of
length `1/7` holds at most one point); `N_i ≤ ⌊q_i/7⌋ + 1` always, so `N_i/q_i ≤ 1/7 + 1/q_i ≤ 15/56`
for `q_i ≥ 8`.

**Corollaries (exact feasibility region).**
- `d = 2`: dispatched whenever `(q_1, q_2) ≠ (2, 2)` — worst surviving cases
  `1/2 + 1/3 = 5/6 < 1`, `1/2 + 15/56 < 1`, `2·15/56 = 15/28 < 1`.
- `d = 3`: dispatched when at most one `q_i ≤ 7` and the coarse one has `q ≥ 3`
  (`1/3 + 2·15/56 = 73/84 < 1`); all-fine `3·15/56 = 45/56 < 1`.
- `d = 4` all-fine: `4·15/56 = 15/14 > 1` — the counting form saturates; sharper `N_i`
  (exact per instance) still fires often, but no clean uniform corollary.

## Proof

Branches with detuned coordinate `i` unsafe number `(g/q_i)·N_i` (each bad orbit value
is hit `g/q_i` times). The union over `i` covers at most `g·Σ N_i/q_i < g` branches, so
some branch has every detuned coordinate safe (clearance ≥ 1/14, the closed complement
of the open danger arc). Harmonic clearances are branch-invariant as in THM-668. ∎

## The d = 2 residual IS the 2-adic descent (the connection)

`q_1 = q_2 = 2` forces `δ_1 ≡ δ_2 ≡ g/2 (mod g)`: the two detuned runners are
**congruent half-harmonics**, so `g ∣ δ_1 − δ_2` AND `g ∣ δ_1 + δ_2` — the pair's sum
and difference are both harmonic. The two branch states move TOGETHER (both flip by
1/2), and both states are bad only when `‖(δ_1 − δ_2)/g · u ± 1/2‖ < 1/7` — a condition
on `u` against the auxiliary integer speed `h = (δ_1 − δ_2)/g`.

**The mod-2g lift:** at modulus `2g` both detuned become `q = 4` coordinates
(`gcd(δ_i, 2g) = g/2`), where the counting dispatch would fire (`1/4 + 1/4 < 1`) — but
the ODD harmonic multipliers `m ∈ H` become half-harmonics of `2g` (`g·m ≡ g mod 2g`).
The obstruction thus DESCENDS the 2-adic tower: dispatching at `2g` requires handling
the odd-`m` harmonics, which is the same problem one level down. This is the branch-side
face of the structures the fleet found from two other directions:

1. **mac-mini's THM-676 (descent-burden):** parity forces ≥ 11 distinct half-sum moduli
   `(v_i + v_j)/2`, Freiman-rigid at APs — the half-sums are exactly the moduli at which
   congruent pairs (`2g ∣ δ_i + δ_j`-type) live.
2. **klein's THM-671-B5/THM-673 (quintic Bonferroni + aggregated supply):** the live
   pair-sum rulers `q = v_i + v_j`. A congruence `g ∣ δ_i + δ_j` says the pair-sum
   `δ_i + δ_j` is a multiple of `g` — **blocking the branch dispatch at `g` supplies
   the pair-sum ruler family at multiples of `g`.**

**THE DUALITY (stated as the sharpening; measured in the companion):** for a family in
the final open core, per modulus `g` the branch dispatch is blocked only by detuned
congruence classes (±pairs mod `g`); each blocking congruence is a divisibility
`g ∣ δ_i ∓ δ_j` that CREATES structure at the pair-sum/difference rulers. Blocking all
small branch moduli forces the difference/sum lattice toward coherence — which is
exactly the resonance-lattice rank-1 direction (opus-S181) where the family approaches
a dilated AP, i.e. the near-harmonic structure THM-668/678 dispatch, or leaves the
covering/compressed class entirely. The conservation candidate: **branch-blocking
ledger + pair-sum supply ≥ uniform floor** over the residual class — quantified
computationally in the companion; the a-priori proof is the named next (it would be the
first uniform statement AT the open core rather than around it).

## Verification & files

`04-computation/lrc14_multi_detuned_duality_monad_S8.py` (+ `.out`): the counting
dispatch on d = 2/3 batteries (fires exactly on the predicted region); the q=q=2
residual instances (dispatch fails at g, fires at 2g modulo odd-harmonic handling —
demonstrated per instance); the blocking/supply ledger on residual-class adversaries.
