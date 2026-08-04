# AMM 12592 / R = 128: transportation-LP form + parity rounding (Angle 3) — CLOSED

**Status: VERIFIED-EXACT (every claim re-verified by exact integer
arithmetic; the new witness passed the hostile independent referee).
boxeph, 2026-08-03.**

Script: `04-computation/amm12592_r128_lp_rounding_boxeph.py`
Endgame widener (diagnostic): `04-computation/amm12592_r128_lp_endgame_boxeph.py`
Consolidated output: `05-knowledge/results/amm12592_r128_lp_rounding_boxeph.out`
Hunt logs: `04-computation/amm12592_r128_lp_hunt_boxeph.log` (winner, v5),
`04-computation/amm12592_r128_lp_hunt_v1_noev1_boxeph.log` (instructive failure)
**New witness (2nd independent R=128 closure):**
`04-computation/amm12592_floor_witness_R128_lp.json`
Transportation-form integer point (translated direct witness):
`04-computation/amm12592_r128_transportation_point_boxeph.json`

## 1. The parity-free claim is CONFIRMED (Phase A)

For every verified floor witness (R = 8, 16, 32, 64, 128-direct, 128-LP),
`f_{i,k} = (binom(d_i,k) - delta_{i,k})/2` is integer with
`0 <= f <= binom(d_i,k)`, and the transportation identity

```text
(**)   sum_{i<R} x^i F_i = T_R := [ (1-x^R)/(1-x) - (1-x)^{R-1} ] / 2
```

holds exactly in `Z[x]`; conversely any integer point of the capacitated
polytope is automatically parity-correct. Parity is a change of variables,
not a constraint. `T_128` is integral because every `binom(127,j)` is odd
(`127 = 2^7 - 1`, Lucas) — THM-3002 4b's dyadic clock is literally the
integrality of the target vector.

## 2. Real relaxation at R = 128: FEASIBLE, but with forced boundary faces (Phase B)

The translated direct-witness point certifies real feasibility exactly.
Saturation anatomy (14721 cells): 3127 (21.24%) on capacity; EVERY row has
a saturated cell, because the top cell `f_{i,d_i} in {0,1}`
(`Delta_i(0) = +-1`) is capacity-1: **the polytope has empty interior in
those directions in every row**, so LP margin alone can never certify
integer feasibility — the correct target is a rounding theorem.

## 3. The decisive cut: the evaluation-at-1 ballot invariant

The bottom cell is equally forced: `cap = binom(d_i,0) = 1` with odd
parity gives `delta_{i,0} = Delta_i(1) = +-1` in every row (this is F4's
sign word). Evaluating the residual recursion at `x = 1`:

```text
sigma_i(1) = sum_{j>i} Delta_j(1)   ==>   |sigma_i(1)| <= R-1-i,
sigma_i(1) == i+1 (mod 2)           (a ballot-path cone).
```

This scalar LP cut is invisible to the per-cell clamp. The first R=128
run (no cut) marched to row 125 with residual L1 = 676 but
`sigma_125(1) = 16` and 2 rows left (max 2): **every beam state was
provably dead long before the endgame** — diagnosed exactly, then fixed.
The R=64 epoch survives without the cut by luck (62 rows of slack); at
R=128 the cut is structurally necessary for this pipeline. Any future
all-R rounding argument must carry the sign-word ballot as an explicit
invariant, alongside capacity.

## 4. LP-guided rounding CLOSES R = 128 (Phase C, v5 recipe)

Per row of `p sigma_i = sigma_{i-1} - Delta_i`, `sigma_{-1} = q^{R-1}`:
Bernstein elimination is triangular, so the per-cell REAL LP optimum is
the clamp `w* = clip(res_j, +-cap)`; round to the forced parity class;
branch rounding directions/even offsets `<= 2` on the first 2 cells (cell
1 forced to leave unit constant); deterministic nearest-parity clamp
elsewhere; **sign cell k=0 branched both ways in the last 32 rows** (the
tight-cone region), deterministic clamp before that; hard-prune by the
ballot cone; beam 400 ranked `(deg, L1, |ev1|)` with sign-twin diversity
cap (max 2 states per low-part group); endgame = row R-2 offsets `<= 6`
on 3 cells + row R-1 exact Bernstein decode. Deterministic, exact.

Ladder: R=8 (0.0s), 16 (0.4s), 32 (2.9s), 64 (40.9s) all SOLVED+verified.
**R = 128: SOLVED in 351 s**, endgame decoded after 49 tries. The winner
rode the cone boundary exactly (ev1 = 13 at row 112 descending 1/row to 2
at row 125). Verified: admissibility + epoch identity + f-translation +
(**), plus the hostile independent referee (exact Fib/Lucas floor
profile, `eff_rate = 58/97 < gamma*`). The witness differs from the
direct-beam witness in 126/128 blocks — a genuinely independent second
closure of the R=128 gamma* floor epoch.

## 5. Rounding-debt anatomy: which levels resist

* **Capacity overshoot** (`sum_j max(0,|res_j|-cap_j)`, the LP-vs-target
  redistribution obstruction): astronomically large in the EARLY rows
  (~1e38, comparable to residual L1), decaying geometrically; saturation
  there is forced play, not choice. Winner total 1.15e39, last nonzero at
  row 126.
* **Parity-rounding debt** (`sum |w - w*|`): tiny throughout (winner
  total 5973, max 87 per row over 127 rows) — parity never resists.
* **Mid-run resistance** sits at the high-degree levels `j ~ d_i` where
  caps `binom(d_i, d_i - j)` are small (1, d, d^2/2, ...); mass must ride
  the shift until mid-band caps absorb it. The last-absorbed level was
  j-top at row 126.
* **The scalar level x = 1** (the ballot cone) is the one genuinely
  global constraint — see section 3.

## 6. Takeaways for the all-R statement

1. Epoch closure at the gamma* floor profile is an integer-point problem
   in a capacitated transportation polytope; its LP relaxation is
   feasible (certified 8..128) but interior-free in the forced unit
   cells; parity is free; the obstruction is rounding + the ballot cone.
2. Constructively, R=128 now has TWO independent engines: target-grid
   steering (direct beam + banded repair) and LP-clamp/parity-rounding
   with the ev1 cut — both need only O(1) branch cells per row plus a
   2-row endgame, evidence that the carry is locally absorbable once the
   sign word is steered — a concrete shape for an inductive all-R proof:
   greedy clamp + ballot-path scheduling of the sign word.
3. Golden attainment now extends to all critical n <= 255
   (epochs R = 8..128 closed at the floor profile), and per the confirmed
   D0(R) = o(R) envelope argument the all-R statement yields
   C* = 1 + gamma* = log_5(5 phi^2).
4. Search negatives remain non-evidence (THM-3029); all positives above
   are exact-arithmetic verified.
