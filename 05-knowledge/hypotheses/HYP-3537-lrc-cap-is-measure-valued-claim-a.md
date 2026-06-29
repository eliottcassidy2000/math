---
id: HYP-3537
title: The LRC cap obeys a measure-valued CLAIM A -- cap(S\{s}) - cap(S) = meas(lonely(S\{s}) cap D_s) = 2 * sum_teeth mu_LRC(s,k), the exact analog of H(T)-H(T-v)=2 sum_{C} mu(C), with the s teeth of D_s as the "odd cycles through s" and the factor 2 from the complement Z_2; the chain rule cap(S)=1-sum(conditional dangers) is the measure-valued OCF
status: CONFIRMED (verified n=6,8,10,14, all deletion-contraction / factor-2 / tooth-sum / chain-rule identities exact). A structural reframe, not a new proof of LRC.
source: mac-mini-2026-06-29-S7
related:
  - THM-070   # Claim A from OCF: H(T)-H(T-v)=2 sum_{C odd cycle thru v} mu(C)
  - THM-002   # OCF (H=I(Omega,2))
  - THM-016   # the tournament alternating-sum (x=-1) identity
  - THM-582   # the two-index synthesis (Redei=odd/x=2, lonely=even/x=-1)
  - HYP-3242  # cap = Euler char of the danger-cover nerve
  - OPEN-Q-108
script: 04-computation/lrc_cap_deletion_contraction_macmini_20260629.py
result: 05-knowledge/results/lrc_cap_deletion_contraction_macmini_20260629.out
---

# HYP-3537 -- the LRC cap is a measure-valued Claim A

## The identity (verified n=6,8,10,14)
For `cap(S) = meas(lonely(S)) = meas{t : ||s t|| >= 1/n for all s in S}` and the danger comb
`D_s = {||s t|| < 1/n}` (s "teeth", arcs around `t = k/s`):
> **`cap(S\{s}) - cap(S) = meas( lonely(S\{s}) cap D_s ) = 2 * sum_{k} mu_LRC(s,k)`,**
> `mu_LRC(s,k) = meas( lonely(S\{s}) cap tooth_k cap [0, 1/2) )`.
- **Deletion-contraction:** removing a speed increases the cap by its CONDITIONAL danger
  (the danger of `s` where the rest is already safe).  Exact analog of Claim A
  `H(T) - H(T-v) = 2 sum_{C odd cycle thru v} mu(C)` (THM-070).
- **Factor 2 = the complement `Z_2`:** `lonely cap D_s` is `sigma:t->-t`-symmetric, so it is twice
  its half-circle part -- the SAME `2` as Claim A's `2 sum mu(C)`.
- **Teeth = cycles:** the `s` teeth of `D_s` (the cyclic positions `k/s` where `s` "wraps") play
  the role of the odd cycles through `s`; the interior teeth `k=1..s-1` carry the measure.
- **Chain rule (measure-valued OCF):** peeling speed-by-speed from `cap(empty)=1`,
  `cap(S) = 1 - sum_{s} meas( lonely(prefix_s) cap D_s )` -- a CONDITIONAL union bound; the floor
  `cap(S) > 0` is exactly `sum (conditional dangers) < 1` (tighter than the naive `|S|/n > 1`).

## Why it matters (and where the analogy stops)
This connects the LRC directly to the project's CORE (OCF / Claim A / Redei): the lonely measure is
built by the same deletion-contraction skeleton, with a factor `2` from the same complement `Z_2`,
and the teeth as the cyclic "through-s" pieces.  BUT the PROOF MECHANISMS diverge, exactly along the
two-index split (THM-582):
- **Tournament (x=2, odd):** `H(empty...)=1` is odd, each step adds an EVEN `2 sum mu(C)`, so `H`
  stays odd => Redei.  A PARITY argument.
- **LRC (x=-1, even):** `cap(empty)=1`, each step SUBTRACTS `2 sum mu_LRC`, and `cap` is a real
  measure with no parity; the floor is the measure INEQUALITY `sum(conditional dangers) < 1`.
So the LRC cap is "Claim A without the parity" -- the measure-valued, even-category twin.  This
CONFIRMS the two-index picture: the deletion-contraction skeleton is shared, but the tournament side
closes by oddness (Redei) and the LRC side must close by a measure bound (the conditional union
bound / the descent + SOS).

## Potential route
The chain-rule criterion `sum_s meas(lonely(prefix_s) cap D_s) < 1` is an alternative target for the
covering floor: order the peel (e.g. odd speeds first, then even, matching the 2-adic descent
THM-580) and bound the conditional dangers.  Each conditional danger is `2 * sum_teeth mu_LRC`, so
bounding the per-tooth pieces (a circular-arc / Diophantine estimate) would close the floor through
this skeleton -- a sibling of the descent route.
