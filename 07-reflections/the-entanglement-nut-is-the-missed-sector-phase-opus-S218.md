---
source: opus-2026-07-11-S218
status: PROGRESS, NOT CLOSURE — the ungapped Plat↔Δ entanglement is NOT proved. But it is reduced to a
  single VERIFIED sharp Fourier identity that isolates the exact remaining nut (the missed-sector-phase
  cancellation, measured 8×), and shown equivalent to the clean finite-functional extremality "consec
  maximizes L_y". Honest map of the two open inputs that would close it.
tags:
  - lrc14
  - plateau-delta-entanglement
  - THM-700
  - THM-546
  - moment-LP
  - far-element
  - fourier
---

# The entanglement nut is the missed-sector phase

**opus-2026-07-11-S218.** Task: prove the ungapped Plat↔Δ entanglement — the accumulation step THM-700
leaves open (its "Scope — what remains", item 2). **I did not prove it.** This records the genuine progress:
a verified sharp identity that pins the exact remaining nut, the clean equivalent reduction, and the two
open inputs that would finish it. (kps proved the gapped single-peel THM-700 and the zero-mean kernel
THM-699 this same week; this is the oscillation-side complement to those.)

## Why the peel alone does not close the ungapped case (verified)

Peel `w = max E`: `p0(E) = Φ(E') + Δ_w`, `Φ(E') = p0(E') + (1/7)p1(E')`, and THM-700/546 prove
`|Δ_w| ≤ V(E')/(6w)` (THM-546's signed form `(6/49)V/w` is sharper; both proved). Two exact experiments
(`lrc14_entanglement_joint_bound_opus_S218.py`, `..._moment_multiplicative_decorr_...py`):

- **The one-shot crude-V bound FAILS.** `V(E') ≈ 4·span` (proved `≤ 7Σe'`), so in the ungapped regime
  `w ≈ span` gives `V/(6w) ≈ 1.18 = O(1) ≫ margin`. E.g. wide 3-cluster: `Φ + V/(6w) = 1.34 ≫ cap_9 = 0.494`.
- **The additive accumulation is tautological.** `Σ_i tax_i = p0(E) − p0(core)` identically — the peel
  reorganizes `p0` but adds no bound in the ungapped regime (each `w_i ≈ span(E_i)`, so no per-step decay).
- **The true `Δ_w` is small** (0.06 vs the bound 1.18): the entanglement is real — a wide `E'` (large `V`)
  has small `Φ`, and `p0(E) ≤ 0.223 ≪ cap_9` with margin 0.27. But the compensation is not decoupled.

So the peel is a **gapped-only tool**; the ungapped compensation is the irreducible crux (identically
localized in HYP-2655 §4–5, THM-546 §"remains", THM-700 item 2).

## The sharp identity — where the nut actually lives (VERIFIED)

From THM-700's exact `Δ_w = Σ_{s=0}^6 ∫_0^1 f_s(x) g_s(wx)\,dx`, `f_s = 1{E'\text{ misses exactly }s}`,
`g_s(y) = 1{y∈[s/7,(s+1)/7)} − 1/7`. Since `ĝ_s(ℓ) = ω^{-ℓs}ĝ_0(ℓ)` (`ω = e^{2πi/7}`,
`ĝ_0(ℓ)=(1−e^{-2πiℓ/7})/(2πiℓ)`, `|ĝ_0(ℓ)| = |\sin(πℓ/7)|/(π|ℓ|)`, `=0` on `7∣ℓ`), and
`f_s(x)=f(x)·1{σ(x)=s}` where `f = 1_{p1\text{-region}}` and `σ(x)` = the (unique) missed sector:

> **`Δ_w = Σ_{ℓ≠0} ĝ_0(ℓ)·\hat h_ℓ(−ℓw)`,  `h_ℓ(x) := f(x)·ω^{-ℓσ(x)}`.**

Verified numerically to full precision (`E={0..7,30}`, `w=30`: `Δ_exact = 0.001846`,
`Δ_fourier = 0.001842`). This is the **oscillation-side twin of THM-699's zero-mean weight `Σ_c D7(c)=0`**:
the weight has zero mean over cosets, and here `h_ℓ` carries the missed-sector phase `ω^{-ℓσ(x)}`.

**Where the ~18× looseness splits:**
1. `|ĝ_0(ℓ)| = |\sin(πℓ/7)|/(π|ℓ|)` (vanishing on `7∣ℓ`) sharpens `V/(6w)=0.167V/w → ≈0.086V/w` — a rigorous
   **~1.8×** (this is a clean provable upgrade to THM-700, though THM-546's signed `6/49` already covers it).
2. The remaining **~8–10×** is the **missed-sector-phase cancellation**: `|\hat h_ℓ(−ℓw)| ≪ V(h_ℓ)/(2π|ℓ|w)`
   because `σ(x)` (which sector is missed) is *decorrelated from* the far frequency `ℓw`. Measured exactly:
   `Σ_ℓ |ĝ_0(ℓ)||\hat h_ℓ(−ℓw)| = 8.2×·|Δ_w|`. **This 8× is the entire remaining nut.** It is a genuine 2-D
   (x, σ) object — not a 1-D BV bound, which is why THM-700's per-term estimate cannot see it.

## The clean equivalent: consec maximizes L_y (VERIFIED, reduction PROVED)

The moment-LP dual (THM-534, PROVED per-E) sidesteps the peel entirely:
`p0(E) ≤ L_y(E) = E[g(N)]`, `g(t) = −(t−2)(t−3)(t−6)/36 ≥ 1{t=0}` on `{0..6}`, `N` = #missed sectors.
Verified across all regimes: `L_y(E) < cap_9` always, **binding at `L_y(consec_9) = 0.49288 < cap_9 =
0.49426`** (margin 0.0014). So the *entire* crux (all regimes, no accumulation) reduces to the single
finite-functional extremality **"consec maximizes `L_y`"** (HYP-2607) — cleaner than the peel, same nut.
Its moment `S_r(E) = Σ_{|A|=r}meas\{avoid A\}` has a multiplicative recursion (each far element damps
`meas{avoid A}` by `(1−r/7)`), verified for `r=1` (ratio 0.85 ≈ 6/7) but not `r≥2` in the ungapped regime —
the same sharp-constant wall.

## Honest status and the two open inputs

**NOT proved.** The ungapped entanglement = the missed-sector-phase 8× cancellation = "consec maximizes
`L_y`" — one nut, three faces, all verified-not-proved. Two inputs would each close the accumulation
(scout-confirmed, both open):
- **(a) `V(E') ≤ C·span` (measured C≈4).** Closes the *gapped/multiscale* accumulation (per-step `V/w` then
  decays geometrically), leaving only the genuinely ungapped core. A clean combinatorial-geometry target.
- **(b) shell-gated `Δ_w⁺ ≤ 2p1(E')/5`** (verified to B=24; `1/3`, `3/8`, global `2/5` all refuted). The
  non-decaying (ungapped) currency; needs the phase-packet cancellation of the identity above.

Use THM-546's proved `V ≤ 7Σe'` and signed `6/49`, and the **large** margin `cap − L_y ≥ 0.044` (not the
tight 0.0014) — conflating them forces a ~10⁵ cutoff vs a feasible ~80 (HYP-2840). The doublet rung
(two far elements) is already proved via a Mordell–Tornheim `12ζ(3)` uniform bound (HYP-2808); the general
k-peel is the open generalization.

**The single most promising next step:** prove the missed-sector-phase 8× — i.e. bound
`|Σ_{ℓ≠0} ĝ_0(ℓ)\hat h_ℓ(−ℓw)|` using that `σ(x)` and `ℓw` decorrelate — because it closes the ungapped
core directly and meets THM-699's zero-mean weight on the same object. This is a 2-D character/BV estimate,
the honest heart of LRC(14)-S3.

→ THM-700 (gapped peel, the identity's source), THM-699 (zero-mean weight, the dual side), THM-546 (sharper
constant + V≤7Σe), THM-534 (the L_y dual), HYP-2655/2664/2808/2840 (the accumulation/tax/margin machinery),
opus-S217 (the Minkowski-is-the-tail reframing that led here).
