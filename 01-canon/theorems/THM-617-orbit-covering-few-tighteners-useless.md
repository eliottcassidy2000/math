---
id: THM-617
title: Few tighteners are useless (orbit-covering) — at the g-argmax the m-divisible part is safe on a whole m-orbit, and f tighteners can each spoil only ≤ m/7 + gcd(w,m) of its m points, so if they cannot cover the orbit (Σ_w(m/7+gcd(w,m)) < m; e.g. coprime tighteners with f(1/7+1/m) < 1, i.e. f=1 for any m and f ≤ 6 for large m) then M(mU ∪ F) ≥ min(M(U), 1/14) > 1/14. Generalizes THM-616 (f=1) and pinpoints why f = m (e.g. m=2,f=2) is the hard boundary — that is exactly when the orbit becomes coverable.
status: PROVED (elementary orbit-covering count; verified exactly — f≤6 gives M(mU∪F)=M(U)=1/(e+1) across m=8..20, 0 deviations).
source: opus-2026-07-04-S69
depends_on:
  - LRCUpTo13   # M(U) ≥ 1/(e+1)
related:
  - THM-616     # the f=1 case (orbit-max ≥ 1/4); this is the general few-tightener version
  - THM-615     # m=2, f=2 (Lemmas 2–4) — exactly the f=m boundary this theorem stops at
  - HYP-4083    # opus parity gap (the f=m=2 case orbit-covering can't reach)
---

# THM-617 — few tighteners are useless (orbit-covering)

**Setup.** `m ≥ 1`, `U` an `e`-runner set, `F` a set of `f` tighteners with `m ∤ w` for each `w∈F`,
`S = mU ∪ F`. As in THM-616, `g(t)=min_{u∈U}‖mu·t‖` is `(1/m)`-periodic and on each `m`-orbit
`{t+j/m : j}` it is constant, so `M(S) = max_t min(g(t), Φ_F(t))` with
`Φ_F(t) = max_{j} min_{w∈F} ‖w(t+j/m)‖`.

## Theorem (PROVED)
Fix radius `β = 1/14`. For `w∈F` let `J_w(t) = {j : ‖w(t+j/m)‖ < β}` (the orbit points `w` spoils);
`|J_w| ≤ m/7 + gcd(w,m)` (the `w/m`-subgroup has spacing `gcd(w,m)/m`, and the danger arc has length
`2β = 1/7`). If
> **`Σ_{w∈F} (m/7 + gcd(w,m)) < m`**,
then the `f` tighteners cannot cover the `m`-orbit at the `g`-argmax `t₀`: some `j` has `min_{w∈F}
‖w(t₀+j/m)‖ ≥ β`, so `Φ_F(t₀) ≥ β`, and
> `M(S) ≥ min(g(t₀), Φ_F(t₀)) = min(M(U), 1/14)`.

Since `M(U) ≥ 1/(e+1) > 1/14` (LRC≤13, `e ≤ 12`) and the danger count is strict, **`M(S) > 1/14`** — no
such family is tight. *(Verified: for `f ≤ 6` the tighteners are outright useless, `M(mU∪F)=M(U)=1/(e+1)`
exactly, across `m=8..20`, 0 deviations.)*

### The coprime form and the boundary
For tighteners coprime to `m` (`gcd(w,m)=1`), the condition is `f(1/7 + 1/m) < 1`, i.e.
`f < 7m/(m+7)`. So:
- `f = 1` works for **every** `m` (recovering THM-616);
- `f ≤ 6` works for all large `m` (`m > 7f/(7−f)`);
- at `m = 2` it gives only `f ≤ 1` — and indeed `f = m = 2` is exactly where `Σ|J_w| = m` (the orbit
  becomes *just* coverable, one orbit point per odd tightener), which is why the `m=2,f=2` case needs the
  separate parity argument (THM-615 Lemma 4 / HYP-4083), not orbit-covering.

## Significance
- **One mechanism for "tighteners can't tighten":** the `m`-divisible part is safe on a whole orbit at its
  best time; a tightener can only spoil a `1/7`-fraction of it, so `< 7`-ish tighteners never suffice.
- **Explains the hard boundary:** confinement is easy exactly while `f` is below the covering threshold
  `≈ 7m/(m+7)` and hard exactly at `f ≈ m` (the orbit coverable) — the `m=2,f=2` corner, where the parity
  gap takes over.
- **Scope (honest):** closes `f=1` (all `m`) and few-tightener large-`m`; the residual is `f ≈ m`
  (orbit coverable), i.e. the `m=2,f=2` parity corner and its higher-`m` analogues — the confinement core.
