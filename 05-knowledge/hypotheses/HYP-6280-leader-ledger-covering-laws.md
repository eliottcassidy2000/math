---
id: HYP-6280
title: Leader-ledger covering laws — staircase CONFIRMED, four-extremal climb-tightness CONFIRMED, chain-equioscillation at the covering-min REFUTED
status: MIXED (each sub-claim marked); exact computations, 50+ families, zero identity failures
source: boxeph-2026-07-12-S21
script: 04-computation/lrc14_leader_ledger_boxeph_S21.py
result: 05-knowledge/results/lrc14_leader_ledger_boxeph_S21.out
related: [THM-722, LEM-025, THM-668, HYP-6240, HYP-6250, MISTAKE-141]
---

# HYP-6280 — Leader-ledger covering laws

Companion data log to THM-722 (the conservation law ∫v_λ = 2Σ depths, chain partition,
H⁺ parity) and LEM-025 (the climb bound M ≥ v⌊q/(B+v)⌋/q at the (min,max) sum-ruler).

## CONFIRMED (exact)

1. **The deep-well staircase.** {1..12,182}: the first 14 sum-handoffs are exactly
   t = k/183, pair (1,182), depth k/183, k = 1..14; the climb is cut at the stopping time
   k* = ⌊183/13⌋ = 14 by runner 12's lander. M = 14/183 is the top of the staircase. The
   Ostrowski rung = a stopping time of the leader process.
2. **Four-extremal climb-tightness.** LEM-025 is TIGHT (bound = exact M) at: AP {1..13}
   (1/14), every ladder rung {1..12,13k} tested (k = 1,2,3,7,14 → k/(13k+1)), the deep
   well (14/183), the compressed extremal 2·{1..12}∪{13} (1/13). The known bottom of the
   covering M-spectrum lives ON the (min,max)-ruler. NOT tight at GW {1..11,13,24}
   (witness 3/14 on the (1,13) ruler; bound 1/25) — tightness is a structural signature of
   the near-dilate/single-far shapes, not universal.
3. **H⁺ parity.** H⁺ even on every family with an even element (54+ families, incl. all
   40 random primitive DC): the ι-fixed chains through 0 and 1/2 carry the parity.
4. **M = deepest handoff** and **chain mass = x_in + x_out**: exact on every family
   (internal asserts, zero failures).
5. **Efficiency gradient along the ladder.** η = ∫v_λ/(2H⁺M) and witness-chain ratio
   (x_in+x_out)/2M for k = 1,2,3,7,14:
   η = 0.785, 0.744, 0.697, 0.591, 0.534; ratio = 0.780, 0.750, 0.833, 0.929, 0.964.
   H⁺ = 58, 70, 82, 130, 220. ∫v_λ = 6.507, 7.711, 8.578, 11.698, 17.987.
   Covering (k=14) pays 3.8× the AP's handoff count; its witness chain approaches but
   does not reach equioscillation.

## REFUTED

6. **"Covering ⟹ max-chain-mass ≥ 28/183" is FALSE** — refuted by the covering-min
   extremal itself: the deep well's max chain mass is 2639/17751 = 0.148668 <
   28/183 = 2716/17751 = 0.153005. Why it fails: the two chains flanking the witness
   handoff have masses 13/183 + 14/183 = 27/183 = 0.14754 (the staircase chain entering
   the witness one step below the top) and 14/183 + 14/194 = 2639/17751 = 0.14867 (the
   exit chain, whose far end is the (12,182)-ruler depth 7/97 = 14/194) — both < 2M.
   Chain-equioscillation (mass = 2M) is an AP/compressed-extremal phenomenon (both have
   ratio exactly 1), NOT a covering-min phenomenon. Lesson: at the covering-min the extremality is the
   STAIRCASE'S TOP (a one-sided cut), not a two-sided alternation — consistent with the
   Chebyshev-equioscillation thread's finding that the deep well is a one-sided object.
   (Also: max-mass ≥ 28/183 holds trivially on all 40 random DC families, M ≥ 0.125 —
   the refutation needed the extremal; MISTAKE-127's lesson applied prophylactically.)

## LIMIT (recorded so nobody re-derives it)

7. The average-depth bound M ≥ ∫v_λ/(2H⁺) (THM-722 D2) reaches only 0.0409 at the deep
   well (M = 0.0765, factor 1.87 short): the conservation law alone cannot prove the
   crux; depth-weighted or per-ruler refinements would be needed.

## OPEN (the ledger-shaped crux)

8. **Climb-cut inverse question:** covering ⟹ no family cuts EVERY sum-ruler's climb
   below 14/183. Cutting the (v,f)-ruler at step k needs a lander within v·k/q of the
   ruler point k/q; the AP achieves all-cuts-at-1/14 with the full lattice {k/14} and is
   thereby non-covering. Quantify "cutting everywhere requires the AP lattice" = the
   lander-exclusion count (klein-S270 target (a)) with an explicit process attached.
9. Does the covering-min = min over covering-realizable (v_min, B, f) of the climb bound
   (i.e., is the crux's extremal ALWAYS climb-tight)? Deep well: yes. A counterexample
   would be a covering family with M < its climb bound's shape-minimum — none seen.
