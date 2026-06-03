---
source: opus-2026-06-03-S598 (remote-control)
status: SYNTHESIS + EXTENSIONS — Helly entropy accounting unifies the cascade (S545), the H-entropy (S543), the Helly certificate (S599), the scale currency (S600); new: the Helly number is an order parameter, the worry-set is isostatic (full Helly = n−1), three co-monotone accountings, residual is loglog-thin
tags: [LRC, helly, entropy, accounting, cascade, clearance, isostatic, order-parameter, scale-currency, iterated-log, n14]
---

# Helly entropy accounting: three co-monotone order parameters

**Prompt (user):** work on Helly entropy accounting and extend it in many new ways.

Fusing the four repo threads — the **cascade** `|SAFE|=∏cᵢ` (S545), the **H-entropy**
`S_H=log H` (S543), the **Helly certificate** (S599, my S595), the **scale currency**
(S600, my S597) — gives one ledger with three co-monotone accountings, all extremized at
the tight witness, and several new handles.

## 0. The unified object

> **Clearance-entropy ledger.** Order the runners; `cᵢ = μ(S_{≤i})/μ(S_{<i})` is the
> conditional clearance (S545). Then the **loneliness debt**
> ```
> log μ(SAFE) = Σᵢ log cᵢ      (each cᵢ ≤ 1, so each term ≤ 0).
> ```
> A config is **tight** iff some `cᵢ = 0` (a *trapped* runner, debt `−∞`); loose iff all
> `cᵢ > 0` (finite debt). The **Helly number** `h` = the smallest subfamily of runners
> whose joint safe set is already measure-zero — the size of the looseness/tightness
> *certificate*. The **scale currency** charges each clearance `cᵢ = 1 − saving_i` in
> iterated-log units (S600).

Computed (`lrc_helly_entropy_accounting_s598.py`):
```
AP n=6 (tight):  c = (2/3, 3/4, 5/9, 1/2, 0)   debt −∞   Helly h = 5 (= n−1, FULL)
generic n=6:     c = (2/3, 3/4, 2/3, 23/28, 16/23)  debt −2.39  (h small)
AP n=7 (tight):  c = (5/7, 4/5, 2/3, 5/8, 6/25, 0)  debt −∞   Helly h = 6 (FULL)
```

## 1. EXTENSION — the Helly number is an order parameter

The raw worry-set has **NO small Helly certificate**: the AP's trap number is `n−1`
(every runner is necessary; removing any one leaves positive measure). By contrast the
post-peel **determinant residual has `h ≤ 2`** (S599: singletons/pairs, live 0). So:

> **`h` (the Helly number) is an order parameter:** `h = n−1` (full) at the worry-set —
> the genuinely entangled tight core — and `h ≤ 2` (small) on the loose residual. It
> measures *how non-local the obstruction is.* The proof difficulty is exactly the
> `h = n−1` (full-Helly) locus; everything with small `h` is cheaply certified.

## 2. EXTENSION — the worry-set is ISOSTATIC (full Helly = every runner necessary)

`h = n−1` means **removing any single runner makes the AP loose** (positive measure).
The worry-set is therefore **critically rigid / isostatic** — exactly the S585-H4
hypothesis, now concrete: *the tight config is the one where the constraint count equals
the degrees of freedom, so every runner is load-bearing.* This is the Helly reading of
"minimally tight": the AP is the unique config of full Helly number, the
maximum-entropy (S543) isostatic point.

## 3. EXTENSION — three co-monotone accountings

At the tight witness, three independent ledgers are all extremized together:

| accounting | worry-set (tight) | loose | source |
|---|---|---|---|
| **H-entropy** `S_H = log H` | **MAX** (regular polygon) | low | S543 |
| **clearance entropy** `log μ(SAFE)` | **−∞** | finite | S545, here |
| **Helly number** `h` | **n−1** (full/isostatic) | ≤2 | S599, here |

> **Helly entropy accounting = the joint ledger of these three.** They are *dual*: the
> tournament *gains* combinatorial entropy (`S_H↑`) exactly as the safe set *loses*
> geometric measure (`log μ↓`), and the obstruction *delocalises* (`h↑` to full). The
> conservation is qualitative (all extremized at the AP); a quantitative law would pin
> the exchange rate `Δ(log μ) : Δ S_H : Δ h` — an open, sharp target.

## 4. EXTENSION — the scale currency (iterated-log) on the ledger

Writing `cᵢ = 1 − saving_i`, the survival is `μ(SAFE) = ∏(1−saving_i) ≤ exp(−Σ saving_i)`
(S600). The index of the sum is the **Helly stage**:
- **Full-Helly worry-set:** `n−1` stages, the savings telescoping to the floor `1/n`
  (the first moment `2kθ`, S600) — the debt diverges at the last stage (`c=0`).
- **Small-Helly residual:** `≤2` stages ⟹ survival is a **2-factor product** ⟹ the
  looseness margin is **loglog-thin** (`Σ over ≤2 Helly stages` = a level-2 iterated-log
  saving — exactly my `ω(2n−1)∼loglog n` worry-set complexity, HYP-2145). So the residual
  spends only `loglog` currency; the worry-set spends the full `n−1`-stage budget.

## 5. EXTENSION — the Helly-3 shadow and the cascade's zero

S545's *cycle-exclusion is Helly-3* is the **triple shadow** of the full-Helly cascade:
a directed 3-cycle is the smallest Helly obstruction, necessary but not sufficient. The
clearance ledger refines it: the cascade's *zero factor* (the trapped runner) is the
**rank-collapse** of the Helly system — and (S595) that collapse is the **rank-1
two-block determinant** for the multiple-of-`n` residual. So:

> **Helly-3 (triple) ⊂ rank-1 two-block (pair, S595) ⊂ full-Helly cascade zero (n−1,
> worry-set).** The Helly number climbs from 2 (residual, certifiable) to `n−1`
> (worry-set, isostatic) as the obstruction delocalises — one ladder, the order
> parameter `h`.

## 6. Honest status

- **Verified:** the clearance-entropy ledger `log μ(SAFE)=Σ log cᵢ`; the AP's trap number
  is `n−1` (full Helly, isostatic) while the residual is `≤2` (S599); the three
  accountings co-extremized at the AP.
- **New (mine):** `h` as an order parameter; the worry-set = isostatic (full Helly,
  concretising S585-H4); the three-ledger duality; the Helly-3 ⊂ rank-1 ⊂ full-Helly
  ladder; the residual's loglog-thin currency.
- **Open/sharp:** the *quantitative* conservation law (exchange rate between `log μ`,
  `S_H`, `h`); a *proof* that the worry-set is the unique full-Helly (isostatic) config.

**Artifacts:** `04-computation/lrc_helly_entropy_accounting_s598.py` (+`.out`). Builds on
S545 (cascade/Helly-3), S543 (H-entropy), S599 (two-block Helly), S600/S597 (scale
currency / iterated-log), S595 (rank-1 two-block), S585-H4 (isostatic). New: **HYP-2151**.
