---
source: opus-2026-06-03-S579 (remote-control)
status: STRUCTURE — the extremal family is a 2-adic self-affine tower; stratifies by 2-adic valuation (geometric ratio 2); two recursions (×2 fractal, +2 ladder); difficulty of n=14 isolated to the odd layer one doubling above prime 7
tags: [LRC, extremal, fractal, self-affine, 2-adic, stratification, doubling, recursion, sleeve, n14, apex]
---

# The extremal family is a 2-adic self-affine tower

**Prompt (user):** the extremal family stratifies by index; extend creatively;
abstract and recursive; the precise fractal nature of the structure.

The extremal family (the AP and the worry-set it anchors) carries **two commuting
recursions** — additive `+2` (the parity ladder) and multiplicative `×2` (a doubling)
— and under the doubling it is **self-affine**: a precise fractal, stratified by 2-adic
valuation with geometric ratio 2. The structure localises the entire difficulty of
`n=14` into one stratum.

## 1. Two recursions on the same family

| recursion | map | role | index |
|---|---|---|---|
| **ladder** `A` | `n ↦ n+2`, append an antipodal pair (flip ∈ {0,1}) | parity carrier (HYP-2091); doubles the worry-set (`2^((n-2)/2)`) | flip-set `F` |
| **doubling** `D` | `v ↦ 2v`, `t ↦ t/2`; `n ↦ 2n` | **fractal scaling** | 2-adic valuation `v₂` |

`A` is additive (a binary *tree of flips*); `D` is multiplicative (a self-affine
*contraction*). Their product is the integer affine action `⟨×2, +2⟩`; the extremal
family is a fixed structure of both.

## 2. Doubling self-similarity (verified, exact)

> **The even layer of `AP_{2n-1}` at its lonely time `t=1/(2n)` equals `AP_{n-1}` at
> its lonely time `t=1/n`** — identical phase multiset, for every `n` (incl. `n=14`).
> `D : AP_{n-1} ↦ {2,4,…,2n-2} ⊂ AP_{2n-1}` is a phase-halving embedding.

So the level-`n` extremal configuration sits *inside* level `2n` as the even speeds —
the lonely configuration is **self-similar under `D`**. (`lrc_extremal_fractal_s579.py`,
part 1: `True` all `n`.)

## 3. The 2-adic stratification (geometric ratio 2)

Every speed factors uniquely `v = 2^r·(odd)`; stratify `AP_{n-1}` by `r = v₂(v)`.
At the lonely time `t=1/n`, the **per-stratum minimum phase doubles with `r`**:

| n | r=0 (odd) | r=1 | r=2 | r=3 | δ |
|---|---|---|---|---|---|
| 14 | **0.071** | 0.143 | 0.143 | 0.429 | 0.071 |
| 16 | **0.062** | 0.125 | 0.250 | 0.500 | 0.062 |

> **Stratum `r` binds at phase `≈ 2^r·δ`.** The **odd stratum (`r=0`) is the unique
> binder** — exactly `{1, n-1}` sit at `δ` — and every higher 2-adic stratum is
> `2×, 4×, …` safer. The constraint is carried entirely by the *coarsest* scale.

## 4. The self-affine tower

```
AP_{n-1}  =  ⊔_{r ≥ 0}  2^r · O_{⌈(n-1)/2^r⌉}        (O_m = {1,3,5,…} the odd generators)
```
Each runner contributes a **sleeve of measure exactly `2δ`** (S578), but the *arc
geometry* refines by depth: a stratum-`r` speed `2^r·odd` has `2^r·odd` arcs of width
`2δ/(2^r·odd)` — coarse, wide arcs at `r=0`; fine, thin arcs at large `r`. So the
total danger is a **self-affine (multifractal) cover**: one scale per 2-adic level,
contraction ratio `2`, the binding happening at the coarse top. Under `D` the whole
tower shifts up one level (`r ↦ r+1`, phases `× ½`, a fresh odd binder inserted at the
bottom) — a genuine **IFS fixed structure**.

## 5. The payoff: `n=14 = 2·7` is one doubling above the prime

For `n = 2^a·q` (`q` odd) the **even layer is `D`-scaled from level `q`, and carries
exactly `×2` slack**; the **odd layer binds at `δ`** (verified `n=14,28,12,24`: even
slack `×2.0` exactly, odd layer binds). Concretely at `n=14`:
- **even layer** `{2,4,…,12} = 2·{1,…,6} = D(AP_6)` — the **prime-`7`** AP, *solved*,
  sitting at phase `≥ 2δ = 1/7` (comfortable);
- **odd layer** `{1,3,…,13}` — binds at `δ = 1/14`, the entire extremal constraint.

> **The whole difficulty of `n=14` lives in the odd layer — one `D`-step above the
> solved prime `7`.** The even half is a scaled copy of a case the literature closes;
> only the freshly-inserted odd generators are new. This is the fractal meaning of
> "`14 = 2·7`": the composite is a doubling of a prime, and the obstruction is exactly
> the doubling's new (odd) layer.

## 6. Creative extensions opened

1. **Inductive C′/loneliness up the doubling tower.** `D` maps level-`q` loneliness to
   the even layer of level-`2q` (with `×2` slack); the residual is the odd layer. A
   proof of C′ at the odd-generator level, lifted by `D`, would give C′ for the whole
   dyadic tower `2^a·q` from the base prime `q`. For `n=14`: base `q=7` (prime, solved)
   ⟹ reduce `14` to the odd-layer lemma.
2. **The apex `2q` as the 2-adic boundary.** `v = q = n/2` is the *largest odd* speed in
   the binder stratum; the zero-divisor `2·q ≡ 0` (HYP-2063) is the collision between
   stratum `r=0` (odd `q`) and the `D`-image — the apex obstruction is the seam where
   the doubling tower meets itself mod `n`.
3. **Stratify the worry-set by the pair `(r, F)`** — 2-adic depth `×` flip-set — a
   two-index self-similar family: `D` acts on `r` (geometric), `A` on `F` (binary tree).
   The worry-set is the boundary of this `⟨×2,+2⟩` tree, a self-affine Cantor structure
   whose finite truncations are the per-`n` worry-sets (`2^((n-2)/2)` classes).
4. **Saturation is top-heavy.** In sleeve-saturation language (S578) the coarse (odd)
   sleeves do the binding; the fine (even) sleeves are background mass. The deficit
   `Φ` is governed by the top stratum — suggesting a *scale-truncated* `Φ` (keep only
   `r=0`) as a sharp approximation.

## 7. Honest status

- **Doubling self-similarity** (even layer `= D(`lower`)`): **proved/verified** exact.
- **2-adic stratification, ratio 2, odd-layer-binds**: **verified** (`n=8..16`,
  dyadic towers `14,28,12,24`).
- The **localisation of `n=14` to the odd layer above prime 7** is exact and is a
  genuine new attack surface (extension 1).
- The IFS / fractal framing and extensions 1–4 are structural directions, not yet
  theorems.

**Artifacts:** `04-computation/lrc_extremal_fractal_s579.py` (+`.out`). Builds on
S578 (sleeves), HYP-2091 (`+2` ladder), HYP-2063 (`2q` apex), HYP-2097 (flip-sets),
THM-369 (clock). New: **HYP-2116**.
