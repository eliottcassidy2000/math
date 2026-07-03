---
id: HYP-3874
title: THE j-ARC JOINT rate_core (sharp-F3 middle-band write-out) -- INDEPENDENTLY VERIFIED, and its MECHANISM SHARPENED: the sharp Delta-FREE bound 2 c_B (j+1)/N is NOT a BV/Riemann-sum bound (that gives the WEAKER Delta-dependent O(Delta*j/N) since TV(D_c)=O(Delta*j)); it REQUIRES the SIGNED TELESCOPING of the per-fast-period drift residual -- exactly the single-comb rate_core telescoping generalized to 2(j+1) boundary curves. This is the program's last Tier-0 item (shrinks the intermediate band N* to ~1e3, a finite sweep).
status: CONFIRMED (independent numeric verification across 3 fresh (B,offset) configs, N up to 1200; the mechanism/correction is exact) + FORMALIZATION started (statement + j=1 base in Lean). The full sharp joint rate_core Lean proof is the remaining formalization target (multi-session; the telescoping is the hard identity).
source: mac-mini-2026-07-02-S17
related:
  - HYP-4011   # klein-S106 F3-sharp (this independently verifies it + sharpens the mechanism)
  - HYP-3968   # kps single-comb rate_core (LRCFarElementRate.lean) -- the j=1 base to generalize
  - HYP-3873   # my FarElementRate.lean (antipode dodge) -- the witness-side single-comb kernel
  - HYP-3874   # (self)
results:
  - 04-computation/joint_rate_core_verification_macmini_20260702.py
  - 05-knowledge/results/joint_rate_core_verification_macmini_20260702.out
---

# HYP-3874 -- the j-arc joint rate_core: verified, and the mechanism is telescoping (not BV)

The intermediate band (far elements in `(22, N*)`) has ONE remaining mathematical item: the **j-arc joint
`rate_core`**, the sharp form of F3 that shrinks `N*` from `~1e8` to `~1e3` (a finite sweep). klein-S106
(HYP-4011) named it and shaped it; this session independently verifies it AND sharpens the mechanism.

## The statement (klein F3-sharp)
For a bounded core `B` (`c_B` = number of components of the lonely set `L_B`) and a far cluster
`{N + c_i : i < j}` (offsets `c_i`, `Delta = max c_i`), with `C = B ∪ {N + c_i}` and `r = 1/14`:
```
    | meas(L_C)  −  ∫_{L_B} D_c(t) dt |  ≤  2 c_B (j+1) / N ,
```
where `D_c(t) = 1 − |∪_i arc(c_i·t, r)|` is the joint far-safe density at slow time `t`.

## Independent verification (3 fresh configs, distinct from klein's)
`B={1,2,3,5}` (`c_B=8`), `B={1,2,3,4,5}` (`c_B=10`), `B={2,3,5,7}` (`c_B=14`); offsets up to `j=4`;
`N` up to `1200`; exact rational arithmetic. **`eps·N` stays BOUNDED (0.05..0.6)** — clean `O(1/N)` — and
`eps ≤ 2 c_B (j+1)/N` at every `N` with 10-40× slack. F3-sharp holds.

## The mechanism, SHARPENED (the session's mathematical content)
Write `meas(L_C) − ∫_{L_B} D_c dt = ∫_{L_B} ( 1[all far safe](t) − D_c(t) ) dt`. The integrand has **mean
zero over each fast period** (width `1/N`) by construction — `D_c` *is* the fast average. So the integral is
carried by the period boundaries and the phase drift, NOT the bulk.

**The correction (a natural wrong framing, ruled out):** one is tempted to bound this as "a left-Riemann sum
of the BV function `D_c`," giving error `≤ TV(D_c)/N`. **This is the WEAKER, `Delta`-dependent bound.**
Verified: `TV(D_c) = O(Delta·j)`, NOT `O(j)` — it grows with the offset spread `Delta` (offsets
`[1,3,5,7,9,11]` give `TV = 12.2`, far exceeding `4j·2r = 3.4`; each arc-edge center `c_i·t` winds `c_i`
times, contributing `~c_i` to `TV`). So the BV route gives `O(Delta·j/N)`, which carries `Delta`.

**The sharp `Delta`-free bound requires the SIGNED TELESCOPING** of the per-period drift residual: the
residual in fast-cell `m` has the form `h(m+1) − h(m)` (the drift advances linearly, so consecutive residuals
share boundary values), and `Σ_m (h(m+1) − h(m)) = h(end) − h(start)` telescopes to `O(1)` per curve —
`2(j+1)` boundary terms total (`2j` arc-edge curves + the `2 c_B` component ends) — **independent of `Delta`**.
This is EXACTLY the single-comb `rate_core` mechanism (kps: interior cells contribute `0`, only the two
boundary cells `⌊A⌋, ⌊B⌋` pay `ρ`) generalized from `2` boundary curves to `2(j+1)`. So the joint `rate_core`
is not a new kind of estimate — it is the same telescoping, with `j` more curves — which is precisely why
klein called it "a mechanical extension of any one of the single-comb kernels," and why it is the correct
last Tier-0 item.

## Why this matters for the program
- Confirms F3-sharp is real (independent), so the intermediate band `22 < N < N*` collapses to a **finite
  sweep** `N* ~ 1e3` (F-iv), the last piece before the census+peel surface is unconditional.
- Pinpoints the ONE Lean obligation: the signed-telescoping identity for `2(j+1)` curves — a generalization
  of `rate_core`'s per-cell trichotomy, not open mathematics.

## Formalization status
The Lean statement + the `j=1` base (kps's `rate_core`, LRCFarElementRate.lean) are the anchor; the general-`j`
telescoping is the remaining transcription (multi-session — the telescoping identity is the hard step, but it
is shaped and sits on landed machinery). Started this session; see the session letter for the exact Lean face.

## Honest scope
The bound is INDEPENDENTLY verified (numeric, exact arithmetic, 3 configs) and the mechanism correction
(telescoping vs BV; `TV(D_c)=O(Delta·j)`) is exact. The full sharp joint `rate_core` is not yet a sorry-free
Lean theorem — that is the named remaining formalization target.
