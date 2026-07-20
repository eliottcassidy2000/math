---
id: THM-1455
title: "WHY ELEVEN SPECTRA: c3 = 35 + 8k where k counts the 4-subsets containing an ODD number of cyclic triples, and k is EVEN — PROVED, because each triple lies in n-3 four-subsets and n-3 is even exactly when n is ODD. The Pfaffian dichotomy is the engine: for a 4-vertex tournament |Pf| = 3 iff it has EXACTLY ONE cyclic triple (scores (0,2,2,2) or (1,1,1,3)) and |Pf| = 1 iff it has zero or two (transitive (0,1,2,3) or strong (1,1,2,2)). Hence k in {0,4,6,8,10,12,14} — seven values of c3, with k=2 empirically absent and unexplained — and the eleven spectra decompose as 1+1+1+3+3+1+1, c3 = 99 and c3 = 115 each splitting into three by c1. Paley is the unique k = 14 maximum. REFUTED: my hypothesis that the switching-invariant triple statistic t = #{triples with s_ab s_bc s_ca = +1} determines the spectrum — there are 24 values of t and t = 25 alone yields SEVEN different spectra. SCOPED OUT: the 'eleven' is a count (1+1+1+3+3+1+1), not the prime 11; no Kakeya or Mandelbrot connection was found"
status: PROVED (the congruence, the Pfaffian dichotomy, and the evenness of k) + VERIFIED-BY-EXHAUSTION (all 32768 switching classes)
author: opus-2026-07-20-S408
depends_on: [THM-1450 (the eleven-spectra census), THM-474 (switching classes = tilings)]
---

# THM-1455 — Why eleven spectra

## 1. The Pfaffian dichotomy (the engine)

For a `7×7` skew `±1` matrix, `c₃ = e₄ = Σ_{|T|=4} Pf(S_T)²`, and for a 4-subset
`Pf = s_ab s_cd − s_ac s_bd + s_ad s_bc` is a sum of three `±1` terms, so
`Pf ∈ {−3,−1,1,3}` and `Pf² ∈ {1, 9}`.

**Classification (exhaustive over all 64 four-vertex tournaments):**

| `|Pf|` | cyclic triples in the 4-subset | score sequences | count |
|---|---|---|---|
| **3** | **exactly 1** | `(0,2,2,2)`, `(1,1,1,3)` | 8 + 8 |
| **1** | **0 or 2** | `(0,1,2,3)` transitive, `(1,1,2,2)` strong | 24 + 24 |

So `Pf² = 9` **iff the 4-subset contains an ODD number of cyclic triples**. Therefore

```
c₃ = 35 + 8k ,      k = #{ 4-subsets with an ODD number of cyclic triples }
```

which immediately gives the congruence `c₃ ≡ 3 (mod 8)` — satisfied by all seven observed
values `35, 67, 83, 99, 115, 131, 147`.

## 2. `k` is even — proved, and it is an odd-`n` phenomenon

> **Claim.** For odd `n`, `k ≡ 0 (mod 2)`.

*Proof.* Let `ε(T) ∈ {0,1}` be the parity of the number of cyclic triples in the 4-subset
`T`, so `k = Σ_T ε(T)`. Each **triple** lies in exactly `n − 3` four-subsets, hence

```
Σ_{|T|=4} (#cyclic triples in T)  =  (n − 3) · (total # cyclic triples).
```

For `n = 7`, `n − 3 = 4` is even, so the left side is even. And
`k = Σ_T ε(T) ≡ Σ_T (#cyclic in T) (mod 2)`. Hence `k` is even. ∎

The hypothesis is exactly `n − 3` even, i.e. **`n` odd** — the same odd-`n` condition that
makes the switching-class ↔ even-graph bijection single-valued (THM-1430 §2). The parity of
`n` keeps deciding these questions.

**Verified:** `k ∈ {0, 4, 6, 8, 10, 12, 14}` over all 32768 switching classes. All even. ✓

**Unexplained:** `k = 2` is **absent**. Evenness permits it; the tournament structure does
not realise it. This is the one genuine gap left, and it is a sharp small question: *why can
a 7-tournament not have exactly two odd 4-subsets?*

## 3. The eleven, decomposed

| `k` | `c₃ = 35+8k` | `c₁` values | # spectra | switching classes |
|---|---|---|---|---|
| 0 | 35 | `{7}` | 1 | 720 |
| 4 | 67 | `{39}` | 1 | 1680 |
| 6 | 83 | `{23}` | 1 | 5040 |
| **8** | **99** | `{7, 71, 135}` | **3** | 7280 |
| **10** | **115** | `{55, 119, 183}` | **3** | 12768 |
| 12 | 131 | `{231}` | 1 | 5040 |
| 14 | 147 | `{343}` | 1 | **240 — Paley** |

```
11  =  1 + 1 + 1 + 3 + 3 + 1 + 1
```

Seven values of `c₃` (from the seven admissible `k`), of which two — `k = 8` and `k = 10` —
split into three by `c₁`. **Paley is the unique `k = 14` class**: every 4-subset has an odd
number of cyclic triples, the maximum, which is exactly what "doubly regular" buys. The
modal spectrum is `k = 10` (12768 classes), and within it `c₁ = 119` (10080) — the owner's
polynomial of THM-1450.

## 4. Refuted en route

I predicted that the switching-invariant **triple statistic**
`t = #{triples with s_ab·s_bc·s_ca = +1}` would determine the spectrum. It does **not**.

- `t` really is a switching invariant (switching flips 0 or 2 arcs of any triangle, so the
  product around it is preserved) — that part is right.
- But `t` takes **24** distinct values against 11 spectra, and the map is **not a refinement
  in either direction**: `t = 25` alone yields **seven** different spectra.

So the controlling invariant is the *4-subset* parity statistic `k`, not the *triple*
statistic `t`. Recorded because "the obvious switching invariant" is a natural guess and it
is wrong.

## 5. Scoping: the eleven is a count, not the prime

Per HYP-8230's rule that coincidence of index is not evidence of mechanism:

- **The prime 11.** `11 = 1+1+1+3+3+1+1` is a sum of block sizes in a `(c₃, c₁)`
  stratification. Nothing arithmetic happens at 11; no residue, order, or character mod 11
  appears anywhere in the derivation. The repo's existing *"eleven-cores"* thread (subsets of
  size 11 inside the LRC-13 speed set, codex `lrc13_scale_eleven_*`) is a different 11 —
  a cardinality of speed subsets — and is **unrelated**. Do not bridge them.
- **Kakeya.** The repo's Kakeya threads are real but live in the LRC covering geometry
  (needles = directions = speeds; "Kakeya-as-adaptive-graphic-rank"; Kakeya numbers of
  `Z_p ⋊ Z_q` multiplier groups). I found **no** connection to the spectrum census: this
  computation has no direction set, no line-in-every-direction condition, and no incidence
  bound. Searched and came up empty; saying so rather than manufacturing a link.
- **Mandelbrot.** **No connection found.** The census is a finite exact-integer count over
  `𝔽₂`-structured objects; there is no dynamical system, no parameter plane, no
  period-doubling in it. The repo's only Mandelbrot mentions sit in the
  Kaczmarz–Christoffel–Blaschke thread (HYP-3796), which concerns the unit disk and is not
  in contact with this. If there is a link it is not visible from here, and I would rather
  record a clean negative than a suggestive one.

## 6. Open

1. **Why is `k = 2` unrealisable?** The sharpest remaining question, and small.
2. **Is `k` even at every odd `n`, and what happens at even `n`?** §2's proof needs `n − 3`
   even. At even `n` the argument fails and `k` should be able to be odd — worth one run.
3. **Does the `(k, c₁)` stratification refine to a known tournament classification?** The
   splits at `k = 8, 10` into three each are unexplained.

## Verification

`04-computation/why_eleven_spectra_opus_S408.py` (full census with `t`, `c₃`, `c₁`),
`04-computation/pfaffian_parity_opus_S408.py` (the 4-vertex `|Pf|` classification and the
`k` distribution). Outputs in `05-knowledge/results/`.
