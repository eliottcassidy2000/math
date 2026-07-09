---
source: opus-2026-07-09-S184
status: RESOLVED "step 2" (the theta-sum bound of the redirected forward lead) -- honestly, and as a
  synthesis rather than an independent proof. The resonance sum R (L=(6/7)^13+R, R=Sum_{t in Lambda, t!=0}
  prod h(t_i), h(m)=-sin(pi*m/7)/(pi*m)) has LEADING term = the height-3 SCHUR TRIPLES, per-triple
  coefficient c3=2*(6/7)^{k-3}*h(1)^2*h(-1)=-0.00113; verified R_lead=c3*#SchurVec captures the SIGN and
  trend of R (AP maximal at 42 vectors, |R|=0.135; sum-free 0 vectors), accounting for ~35-100%, the rest
  being higher-order relations that MODULATE per the S181 dimension effect (aligned=>constructive ratio<1,
  spread GAP=>destructive ratio>1). But the FULL bound |R|<(6/7)^13 is NOT provable from the resonance sum
  -- the absolute/signed tail diverges (opus-S172 Mertens wall). It is SUBSUMED by the fleet's already-
  PROVED moment-LP density floor D3 (THM-661: mu>=max{Sum c_i E[W^i]: Sum c_i w^i<=1_{w>0}}, degree<=4
  clears all six legs). So step 2 = "the density floor is proved by MOMENTS, not the resonance sum; the
  Schur-triple leading order EXPLAINS why it is order-3 and AP-extremal."
tags:
  - lrc14
  - theta-sum
  - schur-triples
  - moment-lp
  - density-floor
  - mertens-wall
  - honest-negative
---

# Step 2 resolves into the moment-LP — Schur leading order explains D3

**opus-2026-07-09-S184.** Owner: prove step 2 (the theta-sum bound `|R| ≤ f(E₃) < (6/7)¹³`). Working it
gives an honest resolution: the theta-sum's *leading order* is Schur-triple-governed (as S182 predicted and
LEM-015 makes extremal), but the *full* bound cannot come from the resonance sum — it is the fleet's
already-proved moment-LP density floor `D3` (THM-661). Step 2 is a **synthesis**, not a new proof.

## The leading order is the Schur triples (proved, quantified)

`L(S) = Σ_{t∈Λ} ∏ᵢ h(tᵢ)`, `Λ = {t : Σ tᵢvᵢ = 0}`, `h(m) = −sin(πm/7)/(πm)` (Fourier of `1_safe`,
`θ=1/14`), `h(0)=6/7`. The `t=0` term is `(6/7)¹³`; `R = Σ_{t≠0}`. Expand by support size `s=|supp t|`:

- **`s=1`: none.** `t = m·eᵢ` needs `m·vᵢ=0 ⟹ m=0` (`vᵢ≠0`). The first correction is `s≥2`.
- **`s=2`: the ratio relations.** `t = vⱼeᵢ − vᵢeⱼ` (always present) — but with LARGE coordinates
  `~spread`, so `|h(vⱼ)| ~ 1/(πvⱼ)` is small.
- **`s=3`: the Schur triples.** When the *speed set* is sum-closed (`v_a+v_b=v_c`), `t=e_a+e_b−e_c` has
  SMALL coordinates `±1`, so factor `h(1)=−0.138` is the LARGEST available. Per unordered triple (with
  `±t`): `c₃ = 2·(6/7)^{k−3}·h(1)²·h(−1) = −0.00113`.

Measured (`lrc14_theta_sum_leading_order`): `R_lead = c₃·#SchurVec` captures the **sign** (`R<0`, covering
below the independence value) and the **trend** — the AP has the most Schur vectors (42) and the largest
`|R|=0.135`; the sum-free `{1,3,…,25}` has 0 and `R_lead=0`. The ratio `R_lead/R` runs `0.31–1.13`: the
Schur triples are the leading term, and the remainder is **higher-order relations that modulate exactly per
the S181 dimension effect** — aligned triples (AP/near-AP) interfere constructively (`ratio<1`, `R_lead`
undershoots), spread 2-D-GAP triples interfere destructively (`ratio>1`, `R_lead` overshoots). This is the
mechanistic content of S182's `corr(E₃, deficit)=0.79` and S181's "dimension, not scalar."

So the leading-order half of step 2 is real: **`R`'s dominant term is `c₃·E₃`, and the AP uniquely
maximizes `E₃` (LEM-015), hence the leading `|R|` — the additive-triple analogue of S180's Freiman `E₂`.**

## Why the full bound can't come from the resonance sum

The naive step-2 bound `|R| ≤ Σ_{t≠0} ∏|h(tᵢ)|` **diverges** past `(6/7)¹³`: the number of relations of
each support size grows, and `Σ_s (6/7)^{k−s}|h(1)|^s N_s` is an L¹/TV sum with no sign cancellation —
exactly the **Mertens wall** (opus-S172: `TV(W) ~ spread²`, the resonant tail L¹-diverges while the true
`R` is small only by *signed* cancellation). A clean monotone `f(E₃)` bounding `|R|` does not exist (S181:
GAPs break scalar monotonicity). So the theta-sum route hits the same 2025-core signed-cancellation
difficulty it was meant to avoid.

## What actually proves it: the moment-LP (already done)

The density floor is proved by a *different* route that never forms the resonance sum: the **one-sided
moment LP** (THM-661, mac-mini-S57):
```
   μ_{1/7}(E) ≥ B_d(E) = max{ Σ_{i≤d} c_i E[W^i] : Σ_i c_i w^i ≤ 1_{w>0} on [0,6/7] },   W = Σ(g_j−1/7)_+.
```
`B_2 =` Paley–Zygmund clears `k=11,12,13`; `B_4` clears `k=8,9,10`; the closed-form `D3` clears the binding
`k=12,13` with margins `+0.157 / +0.252`. This uses only the first `≤4` moments `E[W^i]` (exact via
Farey-cell integration) and Markov–Krein positivity — **no resonance sum, no Mertens issue**. The
"order-3" of `D3` (the third moment `E[W³]`, triple gap-overlaps) is the *continuous* shadow of the same
order-3 additive structure my `E₃` counts discretely: both say **the covering floor is an order-3
phenomenon, extremized by the AP** — my Schur analysis is the resonance-side explanation of why `D3`
(degree ≥ 3), not `D2`, is where the floor becomes robust.

## Resolution (the honest statement of step 2)

> **Step 2 is subsumed, not independently proved.** The density floor `|R| < (6/7)¹³ ⇔ μ > 0` is a theorem
> (THM-661, moment-LP `D3`/`B4`, PROVED per-shape; uniform floor = compact check + decorrelation tail).
> The theta-sum route cannot prove it (Mertens wall). Its value is EXPLANATORY: the resonance sum's leading
> term is the Schur triples (`c₃·E₃`), the AP maximizes them (LEM-015), and the higher-order modulation is
> the S181 dimension effect — so the forward lead's "Schur/sum-free structure theorem" is the *understanding*
> of the density floor, while the *proof* is the moments.

## Ledger

- RESOLVED step 2: leading order of `R` = Schur triples (`c₃=−0.00113` per vector, `R_lead=c₃·#SchurVec`
  captures sign+trend, AP-maximal); full `|R|<(6/7)¹³` NOT provable from the resonance sum (Mertens wall,
  S172), SUBSUMED by the proved moment-LP density floor `D3` (THM-661).
- The Schur (order-3) leading order is the resonance-side explanation of why the floor is degree-3
  (`D3` robust, `D2` marginal) and AP-extremal; unifies LEM-015 / S182 / S181 with THM-661.
- File: `lrc14_theta_sum_leading_order_opus_S184` (+out). -> LEM-015, opus-S182/S181/S172, THM-661
  (moment-LP, PROVED), THM-515 (theta-sum), LEM-011 (arc coeffs).
