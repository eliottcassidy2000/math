# LEM-019 — The dyadic descent lemma: at τ = (σ+1)/2 evens keep their clearance exactly and odds keep half; the doubling-chain corner collapses to a 7-runner mixed system whose tight locus is EMPTY (all-odd sets are 1/2-lonely at τ = 1/2), closing the corner boundary for o_max ≤ 5w

**Status:** PROVED (the descent lemma, the merge corollary, the 2-adic pairing wall, the
all-odd half-witness, the o_max ≤ 5w closure — elementary proofs below; descent verified
0 violations / 200k exact random checks; all-odd M = 1/2 exact over 100,919 6-sets). Honest
negatives recorded (the crude measure composition FAILS: min meas(G_{1/7}) = 0.0637 < 1/7;
full-throttle descent overflows its threshold budget).
**Source:** mac-mini-2026-07-09-S65 (cont. 10). Attacks THM-682's doubling-chain corner
(monad-S11: "the rung's content is to PROVE the cancellation or the 2-adic dispatch") — this
is the 2-adic dispatch, delivered with exact scope.

## The descent lemma (new instrument I9)

> For any real σ and τ = (σ+1)/2:
> **‖2v′·τ‖ = ‖v′·σ‖** for every integer v′ (evens descend EXACTLY), and
> **‖o·τ‖ ≥ ‖o·σ‖/2** for every odd o (odds keep at least half).
> Hence S = 2S′ ∪ O (O odd) is 1/14-lonely at τ whenever σ makes S′ 1/14-lonely and O
> 1/7-lonely — a MIXED-threshold problem one 2-adic level down.

*Proof.* 2v′τ = v′σ + v′ ≡ v′σ. For odd o: oτ = (oσ + o)/2; with x = frac(oσ) and n = ⌊oσ⌋,
oτ ≡ x/2 (n+o even) or (x+1)/2 (n+o odd), giving ‖oτ‖ = x/2 or (1−x)/2 ≥ min(x,1−x)/2. ∎

**Merge corollary.** Each ODD-ROOTED doubling (o, 2o ∈ S, o odd) merges two runners into one
(o at threshold 1/7) under one descent: the reduced distinct count is N₁ = 13 − #odd-rooted
doublings.

**The 2-adic pairing wall.** Odd-rooted doublings are disjoint pairs, so at most 6 fit in 13
runners: **N₁ ≥ 7 always** — the `13 = 2·6+1` calibration reappearing inside the 2-adic tower
(the same wall as THM-668's C3 pairing and the descent-burden parity split). Full-throttle
iterated descent overflows: a lagging odd root at depth d needs 2^d/14, so d ≤ 2; the bare
LRC(N+1) citation covers the mixed demand only for N ≤ 6 — one step short at the wall,
CALIBRATED, as everywhere in this problem.

## The corner boundary trivializes: the all-odd half-witness

> **Any set of odd integers is 1/2-lonely at τ = 1/2** (‖o/2‖ = 1/2). In particular the
> corner's maximal-collapse boundary (6 odd-rooted doublings ⟹ N₁ = 7: six odd roots at
> threshold 1/7 + one leftover) closes OUTRIGHT when the leftover is odd — all seven runners
> clear at 1/2 with margin 5/14 — no LRC citation, no composition, nothing.
(Verified: min M over 100,919 all-odd primitive 6-sets = exactly 1/2. Also immediate: no
all-odd set can be a tight dilated AP, since 2d would be even.)

> **Even leftover (w = v/2 after descent): the composition closes for o_max ≤ 5w.**
> *Proof.* The six odds keep ‖o·σ‖ ≥ 1/7 on the window σ ∈ [1/2 − δ, 1/2 + δ],
> δ = 5/(14·o_max) (Lipschitz from the 5/14 margin at 1/2). The leftover's nearest good point
> to 1/2 sits at distance exactly 1/(14w) (wσ = w/2 + wε ≡ wε). Compatible iff
> 1/(14w) ≤ 5/(14·o_max) ⟺ **o_max ≤ 5w**, a knife-edge verified exact. ∎

**The remaining sliver (w < o_max/5)** is a SEVEN-runner mixed-threshold instance — tiny,
certificate-dispatchable per instance (the reduced system's pair-sum moduli are ≤ 2·o_max; my
C0–C2 battery applies verbatim with per-runner bands). HONEST NEGATIVE: the crude measure
composition does NOT close it uniformly — min over all-odd 6-sets of meas{σ : all ‖oσ‖ ≥ 1/7}
= 0.0637 (at {3,7,9,15,23,37}) < 1/7, so the leftover's bad set CAN outweigh the odds' good
set; the closure must use the window structure (as above) or per-instance dispatch.

## Scope (honest)

This closes the ≥ 6 ODD-ROOTED-doubling corner (boundary case) and reduces the rest of the
2-adic dispatch to ≤ 7-runner mixed instances; monad's corner census counts doublings of ANY
parity (long even chains merge only at deeper levels where the threshold budget bites), so the
even-chain-heavy corner families route through descent PLUS the composed instruments, not
descent alone. The signed-cancellation arm (THM-681 option (ii)) remains monad's.

**Files:** `lrc14_dyadic_descent_macmini_S65cont10.py`,
`lrc14_allodd_halfwitness_macmini_S65cont10.py`, `lrc14_allodd_measG_macmini_S65cont10.py`
(+ consolidated `.out`). **Related:** THM-682, THM-668 (C3 pairing wall), LEM-018, THM-676,
monad THM-668-detuned (the d=1 case = descent + branch pigeonhole), LRCCoarseReduction.
