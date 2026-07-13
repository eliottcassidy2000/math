---
source: opus-2026-07-11-S266
status: RIGOROUS IDENTITY + honest negative. Deriving the exact relation formula for eps_v
  (eps_v|G'| = Sum_R (-1)^{m_R} (6/7)^{r-m_R} prod_u b_{k_u}, v participating in each relation R) is clean and
  rigorous. BUT it yields NO clean bound: the low-order (m<=2, additive-relation) truncation is 0.13 vs the
  actual eps_v = 0.019 -- the higher-order terms (m>=3) cancel ~0.11 via the alternating (-1)^m sign. So the
  additive bound |eps_v| <= f(#relations) is NOT a clean theorem; eps_v's smallness is an ALTERNATING
  MULTI-ORDER cancellation (confirming S262's multi-linear finding rigorously). The case skeleton (S265)
  stands VERIFIED, but its supporting anti-concentration bounds are the irreducible higher-order core; the
  elementary/low-order tools are exhausted.
tags:
  - lrc14
  - covering-min
  - additive-bound
  - relation-identity
  - alternating-cancellation
  - multi-order
  - honest-negative
---

# The additive bound has no clean truncation: ε is an alternating multi-order cancellation

**opus-2026-07-11-S266.** Owner: prove the additive bound `|ε_v| ≤ f(#relations)` rigorously, work the two
anti-concentration bounds. Doing so produces a clean rigorous *identity* and a decisive honest negative: the
identity does not yield a bound, because ε_v's smallness is a cancellation across *all* relation-orders.

## The rigorous relation identity (positive)

Writing `β_u = 1_{D_u}` (so `β̂_u(h) = b_{h/u}·[u|h]`, `b_k = sin(πk/7)/(πk)`) and `1_{G'} = ∏_w(1−β_w)`,
`ε_v|G'| = ⟨β_v − 1/7, ∏_w(1−β_w)⟩`. Expanding and using that the `k_v = 0` terms cancel against the `−1/7`
(forcing v to participate), one gets the **exact identity**

> `ε_v |G'| = Σ_{relations R} (−1)^{m_R} (6/7)^{r−m_R} ∏_{u∈T'_R} b_{k_u}`,

where a relation `R` is a support `T' ∋ v` (with `T'∖{v} ⊆ non-core`), nonzero integers `k_u`, satisfying
`Σ_{u∈T'} u·k_u = 0`, and `m_R = |T'|−1`. This is rigorous and clean — it expresses `ε_v` as a signed sum over
the **additive relations** among `v` and the non-core speeds, weighted by band coefficients.

## The honest negative: no clean truncation

The identity does **not** give a clean bound. The low-order truncation (`m ≤ 2`, the `±v±w_i±w_j = 0` additive
relations that S263 correlates with) gives **0.13** for `v=41`, versus the actual `ε_v = 0.019` (FFT). So the
`m≥3` terms must cancel `≈ 0.11` — a large alternating cancellation driven by the `(−1)^{m}` sign. Bounding
`|ε_v|` by the magnitude sum `Σ_R (6/7)^{r−m}∏|b_{k_u}|` therefore **massively overshoots** (and the magnitude
sum's harmonic tail `Σ|b_k| = ∞` diverges anyway). **`|ε_v| ≤ f(#relations)` is not a theorem** for any
low-order `f`: the smallness of `ε_v` is an **alternating multi-order cancellation**, exactly the higher-order
(multi-linear) structure S262 identified — now confirmed rigorously via the exact identity. The S263
correlation `0.527` is real but only reflects the leading term; the actual value is set by the alternating tail.

## The measure bound is also anti-concentration

The other supporting bound, `|S_rest| > (s_min−1)/(7 s_min)` (Argument A of S265), is a safe-set-measure lower
bound for the 12-runner rest — itself an anti-concentration statement (measure of a union-of-bands complement),
with the same character: no clean low-order handle. Both supporting bounds of the case skeleton are the
irreducible higher-order anti-concentration.

## Net (honest, and the converged state of the arc)

Attempting the additive bound rigorously yields a clean exact **identity** for `ε_v` (a signed sum over
additive relations) but **no clean bound** — `ε_v`'s smallness is an alternating cancellation across all
relation-orders, so no low-order `f(#relations)` bounds it. This rigorously confirms S262: the covering-min
anti-concentration is genuinely **higher-order (multi-linear)**, and every elementary/low-order tool tried
across S253–S266 (balance, dual certificate, mollification, completion identity, Gowers, additive relations,
second moment) captures only the low-order part.

**The converged state.** LRC(14) for covering families has a **complete, margin-safe VERIFIED case skeleton**
(S252 non-covering + S264 no-speed-1 + S265 speed-1), reducing to **two anti-concentration bounds** — the
additive `|ε_v|` and the measure `|S_rest|` — both of which are **verified across thousands of families with
margin but are irreducibly higher-order**, not provable by the elementary tools. Honest verdict: the *shape* of
the proof is complete and the extremizer is proved (S255), but closing the two bounds requires genuine
higher-order additive-combinatorics / an inverse theorem for the band-multilinear cancellation — the same
"multi-way entanglement" the fleet's Minkowski-tail threads (#42–#43) flagged. The elementary program has
reached its principled limit; further progress needs that higher-order input, not another low-order tool.

→ opus-S265 (case skeleton — supporting bounds located here), opus-S262 (multi-linear — confirmed rigorously),
opus-S263 (additive relations = leading term only), opus-S255 (extremizer, proved), tasks #42–#43 (multi-way
entanglement). Files: `lrc14_additive_bound_relation_identity_opus_S266.py` (+`.out`).
