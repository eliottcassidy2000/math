---
id: THM-400
title: A speed relation is balanced (augmentation 0) ⟺ translation-invariant (observer-blind)
status: PROVED (elementary); the LRC stratification it grades is VERIFIED (k=6,8,10)
source: opus-2026-06-03-S581
related:
  - HYP-2114   # additive folds carry hardness (3-term not 4-term) -- refined here
  - HYP-2117   # binary IFS / doubling (c=2a is unbalanced)
  - THM-398    # C'
---

# THM-400 — balanced ⟺ translation-invariant

## Setup

A **speed set** `S = {v_1,…,v_k}` with the observer `v_0 = 0`. A **relation** is a
vector `c ∈ ℤ^S`, `c ≠ 0`, with `Σ_i c_i v_i = 0`. The **augmentation** (the group-ring
augmentation `ε : ℤ[S] → ℤ`) is
```
ε(c) := Σ_i c_i      (the coefficient sum, a.k.a. the augmentation index j).
```
A relation is **balanced** if `ε(c) = 0`. Translation acts by `S ↦ S + t·𝟙`
(`v_i ↦ v_i + t`).

## Theorem

> **A relation `c` holds for *every* translate `S + t·𝟙` (all `t`) ⟺ `ε(c) = 0`.**

*Proof.* `Σ_i c_i (v_i + t) = Σ_i c_i v_i + t·Σ_i c_i = 0 + t·ε(c) = t·ε(c)`. This is
`0` for all `t` iff `ε(c) = 0`. ∎

So **balanced = translation-invariant**. Balanced relations constrain only the
**inter-runner differences** `{v_i − v_j}` (they are unchanged by a global shift); they
are **observer-blind**. Unbalanced relations (`ε(c) = j ≠ 0`) are **observer-coupled**:
appending the observer `v_0 = 0` rewrites `Σ_i c_i v_i = 0` as `Σ_i c_i v_i + j·v_0 = 0`
— the relation "references the origin `j` times", and translation breaks it (the
observer is the fixed point that distinguishes a translate from the original).

## Grading of the relation lattice

`L(S) = { c ∈ ℤ^S : Σ c_i v_i = 0 }` is graded by `ε`. The **balanced sublattice**
`L_0 = ker(ε|_L)` is the translation-invariant / observer-blind part — shared by the
whole translation orbit `{S + t·𝟙}`. The quotient records the **observer coupling**.

Examples (support is *not* the invariant):
| relation | `ε` | type |
|---|---|---|
| `a+b = c` (the fold) | **+1** | unbalanced — observer-coupled |
| `c = 2a` (doubling) | **+1** | unbalanced — observer-coupled |
| `2a = b+c` (AP-triple, support **3**) | **0** | **balanced** — observer-blind |
| `a+b = c+d` (energy, support 4) | **0** | balanced — observer-blind |

## LRC relevance (the stratification — VERIFIED)

Loneliness is *distance from the observer at `0`*, so it feels only the
**observer-coupled** relations. Binning `M(S) − δ` (`δ = 1/(k+1)`) by the count of
**unbalanced** small relations `u` (a+b=c, c=2a) versus **balanced** relations
(`lrc_augmentation_stratification_s581.py`):

| k | u=0 | u=1 | u=2 | u=3 | u≥4 | (u=0) any #balanced |
|---|---|---|---|---|---|---|
| 6 | +0.088 | +0.050 | +0.031 | — | **+0.000** (tight) | stays ≥ +0.088 |
| 8 | +0.111 | +0.082 | +0.065 | +0.045 | +0.022 | b≥4 → +0.111 |
| 10 | +0.262 | +0.182 | +0.113 | +0.063 | +0.020 | b≥4 → +0.262 |

> **Hardness tracks the augmentation, not the support.** Unbalanced relations drive
> `M → δ` (tight); **any number of balanced relations leaves the margin at the
> observer-blind maximum.** This *supersedes* "3-term vs 4-term" (HYP-2114): the fold
> `a+b=c` is hard because `ε = 1` (it puts `c` at the observer when `a,b` pinch), while
> the balanced 3-term `2a=b+c` is as harmless as 4-term energy.

**Corollary (doubling is observer-coupled).** `c = 2a` has `ε = 1`, so the binary-IFS
generator `D : x↦2x` (S580) is observer-coupled — the hard branch — whereas the
"midpoint" relation `2a = b+c` is observer-blind. The 2-adic fractal generator and the
augmentation grading agree.

**Artifacts:** `04-computation/lrc_augmentation_stratification_s581.py` (+`.out`). The
theorem is elementary and proven; the stratification is verified `k=6,8,10`. See
reflection `07-reflections/lrc-augmentation-index-observer-coupling-s581.md`.
