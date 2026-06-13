---
source: opus-2026-06-03-S581 (remote-control)
status: REFRAME — the real dichotomy is the augmentation index j=ε(c), not support-3 vs support-4; balanced (j=0) ⟺ translation-invariant ⟺ observer-blind (THM-400); hardness tracks j≠0
tags: [LRC, augmentation, observer-coupling, translation-invariant, balanced, fold, THM-400, additive, n14]
---

# The augmentation index, not the support, is the dichotomy

**Prompt (user):** the real dichotomy isn't support-3 vs support-4 — it's the
augmentation index `j`: `j=0` relations are observer-blind (constrain only inter-runner
differences), `j≠0` couple to the observer. Verify the stratification, then formalize
the clean theorem (balanced ⟺ translation-invariant).

Both done. The reframe is correct and it sharpens every prior additive finding.

## 1. The clean theorem (THM-400)

For a relation `c` (`Σ c_i v_i = 0`), the augmentation `ε(c) = Σ c_i = j`. Then
`Σ c_i(v_i + t) = t·ε(c)`, so

> **`c` survives every translation `S ↦ S + t·𝟙` ⟺ `ε(c) = 0`** — *balanced ⟺
> translation-invariant.*

Balanced (`j=0`) relations see only the inter-runner **differences** `{v_i − v_j}` —
**observer-blind**, shared by the whole translation orbit. Unbalanced (`j≠0`) relations
**couple to the observer** `v_0=0`: `Σ c_i v_i = 0` is `Σ c_i v_i + j·v_0 = 0` — the
relation references the origin `j` times, and a shift breaks it.

## 2. Why this beats "3-term vs 4-term"

Support size is *not* the invariant. The two support-3 relations split by augmentation:
- `a + b = c` — `ε = 1`, **unbalanced** (the fold, observer-coupled, HARD);
- `2a = b + c` — `ε = 0`, **balanced** (the AP-triple/midpoint, observer-blind, HARMLESS).

So last session's "hardness is 3-term not 4-term" was a shadow of the real law:
**hardness is `ε≠0`, not `ε=0`.** A balanced 3-term is as harmless as 4-term energy; an
unbalanced relation of *any* support is hard.

## 3. The stratification — verified (k=6,8,10)

Binning `M − δ` by the unbalanced count `u` (a+b=c, c=2a), and — holding `u=0` — by the
balanced count:

| k | u=0 | u=1 | u=2 | u≥3/4 | u=0, many balanced |
|---|---|---|---|---|---|
| 6 | +0.088 | +0.050 | +0.031 | **+0.000** | ≥ +0.088 |
| 8 | +0.111 | +0.082 | +0.065 | +0.022 | +0.111 |
| 10 | +0.262 | +0.182 | +0.113 | +0.020 | +0.262 |

**Unbalanced relations drive `M → δ`; balanced relations never do** — `u=0` configs sit
at the observer-blind maximum no matter how many balanced relations they carry. The
augmentation, not the support, is the order parameter of LRC hardness.

## 4. The picture it unifies

- **Folds (S577) / hardness (oracle S578o):** the fold `a+b=c` is hard precisely because
  `ε=1` — it places `c` at the observer (`0`) at the `(a,b)`-pinch. *That* is the
  "observer coupling".
- **Doubling / binary IFS (S580):** `c = 2a` has `ε = 1` — the generator `D: x↦2x` is
  observer-coupled (the hard, fractal-driving branch); the AP-midpoint `2a=b+c` (`ε=0`)
  is observer-blind. The 2-adic fractal and the augmentation grading coincide.
- **Sleeve alignment (S578):** "coherent overlap that drives saturation" = the
  observer-coupled (`j≠0`) alignments; balanced overlaps are the translation-invariant
  background that never saturates.
- **The translation orbit:** `{S + t·𝟙}` all share `L_0` (balanced structure); they
  differ only in the observer-coupled part — and exactly *one* translate of the AP is
  tight (S580), the one whose `j≠0` relations land on the observer.

## 5. Abstractly

The relation lattice `L(S)` is **ε-graded**; `L_0 = ker ε|_L` is the
translation-invariant sublattice (a `GL`-/shift-invariant of the *difference set*),
and LRC hardness is a function of `L(S)/L_0` — the **observer-coupling class**. LRC
lives on the augmentation-nonzero part; the balanced part is pure Sidon/energy
background, invisible to the lonely observer. This is the precise sense in which
loneliness is a property of how the runners sit *relative to the origin*, not of their
mutual difference geometry.

## 6. Honest status

- **THM-400 (balanced ⟺ translation-invariant):** elementary, **proved**.
- **Augmentation = observer-coupling; hardness tracks `ε≠0`:** **verified** k=6,8,10
  (unbalanced → tight; balanced harmless at any count). Supersedes HYP-2114's
  3-vs-4-term as the sharp form.
- The grading `L/L_0` as *the* order parameter is a structural claim; the proof that
  "only balanced ⟹ `M ≥ c > δ`" is the augmentation-graded form of Lemma A (open).

**Artifacts:** `04-computation/lrc_augmentation_stratification_s581.py` (+`.out`),
`01-canon/theorems/THM-400-...md`. Builds on HYP-2114 (folds), HYP-2117 (doubling/IFS),
THM-398 (C′), S578 (sleeves). New: **HYP-2118**.
