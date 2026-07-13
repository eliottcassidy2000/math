---
source: opus-2026-07-11-S256
status: PROGRESS + DECISIVE HONEST FINDING. Attempting to prove the full beta-balance (0<beta<14/183) yields
  a clean algebraic reduction (proved via LRC(13) for 69% of hard cases, all s=1) BUT the decisive finding is
  that the beta-balance is an equioscillation HEURISTIC, NOT a rigorous lower bound on M(family): it EXCEEDS
  M(family) in real cases (e.g. {1..11,84}: balance 15/183 > M 15/184), because the true optimum sits at a
  different denominator and a third runner obstructs. So "balance>=14/183" cannot prove the covering-min. This
  re-confirms mac-mini S40 (need a DUAL certificate) and corrects the S253-S255 reliance on the balance as a
  bound. The deep-well tight case (S255 via S252) stands, independent of the balance.
tags:
  - lrc14
  - covering-min
  - beta-balance
  - heuristic-not-bound
  - dual-certificate
  - honest-negative
  - course-correction
---

# The β-balance is a heuristic, not a bound; the general covering-min needs a dual

**opus-2026-07-11-S256.** Owner: prove the full β-balance for `0 < β < 14/183`. Working it gave real algebraic
progress and then a decisive, honest negative that re-orients the whole route.

## Clean form and algebraic progress

At the core optimum `t₀ = a/q`: `M_core = m/q`, killer clearance `β = ‖182 t₀‖`, binding speed `s`. The
balance witness value is `(βs + 182 M_core)/(182+s)`, and requiring it `≥ 14/183` rearranges to the **exchange
inequality**

> `182·Δ ≥ s·ε`,  where `Δ = M_core − 14/183` (core excess), `ε = 14/183 − β` (killer deficit).

Using only `M_core ≥ 1/13` (LRC(13), so `182 M_core ≥ 14`), this reduces to **`β ≥ (14/183)(1 − 1/s)`**:
- **`s = 1`: `β ≥ 0` — always.** All s=1 hard cases clear it.
- **General: proved via LRC(13) whenever `β ≥ (14/183)(1−1/s)`** — 69% of hard cases (including all s=1).
- **Remaining zone `β < (14/183)(1−1/s)`** (deeply resonant + large s): needs `M_core` bounded above `1/13` —
  the binding-speed↔core-value coupling (verified 182/182: large `s` ⟹ large `M_core` excess).

## The decisive finding: the balance is not a lower bound

But the balance **value** is **not** a rigorous lower bound on `M(family)`. Verified — it **exceeds**
`M(family)` in real cases:

| core | balance value | `M(family)` | `q_core` | `q_fam` | balance ≤ M? |
|---|---|---|---|---|---|
| `{1..11,84}` | `15/183 = 0.08197` | `15/184 = 0.08152` | 85 | 184 | **No** |
| `{1..12}` (deep well) | `14/183` | `14/183` | 13 | 183 | yes (exact) |
| `{1..10,12,36}` | `46/465 = 0.0989` | `18/187 = 0.0963` | 40 | 187 | **No** |

The true `M(family)` is attained at a **different denominator** (`q_fam ≠ q_core`) — the operative witness is
elsewhere, and a **third core runner obstructs** the local perturbation, so the local balance **overestimates**
the achievable clearance. Consequently **"balance ≥ 14/183" does not imply "M(family) ≥ 14/183"**. The
β-balance is an **equioscillation heuristic** — exact at the deep well (where `q_fam = q_core` structure lines
up), not a general bound.

## Consequence and course-correction

The β-balance is the **wrong vehicle** for the rigorous covering-min lower bound. Sessions S253–S255 used it as
if it bounds `M(family)`; that holds **only at the extremizer** (the deep well), where it is exact. This
**re-confirms mac-mini S40** from a fresh angle: the max-min is non-convex, a **local / greedy witness has no
shortcut**, and the general covering-min bound needs a genuine **dual certificate** (Delsarte / de la
Vallée-Poussin positive trigonometric polynomial), not the local balance.

## What stands

- `M(family) ≥ 14/183` is **verified** throughout (klein S267, zero counterexamples) via family-specific
  witnesses — the covering-min is not in doubt, only the *proof route*.
- The **deep-well tight case is rigorous and independent** of the balance's heuristic nature: S255 proved it
  via S252 (`M_core = 1/13 ⟹` interval `⟹ s = 1 ⟹` equality), and there the balance *is* exact. So the
  **extremizer and its uniqueness stand.**

## Net (honest)

Working to prove the β-balance produced a clean exchange form and an LRC(13) reduction covering the majority of
hard cases — but also the important negative that **the balance does not rigorously bound `M(family)`** (it
overestimates; the true witness is elsewhere). So the general single-killer covering-min bound is **not**
reachable through the local balance; it belongs to the **dual-certificate route** (S40). The honest status of
the covering-min lower bound: **extremizer proved (deep well, S255), general bound open and correctly located
on the dual/Delsarte side** — the local-witness program (S253–S256) has reached its limit, exactly as S40
predicted.

→ mac-mini S40 (equioscillation / greedy-no-shortcut / dual certificate — confirmed here), opus-S253 (the
balance), opus-S255 (deep-well tight case via S252 — stands), klein S267 (14/183 verified), LRC(≤13). Files:
`lrc14_beta_balance_is_heuristic_not_bound_opus_S256.py` (+`.out`).
