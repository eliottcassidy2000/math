---
source: opus-2026-07-11-S262
status: DECISIVE LOCATION. Applying the LRCFourierCompletion cancellation bound (the bilinear completion
  identity) to the residual eps_v = Sum_{h!=0} b_h ghat(-hv) shows it is BILINEAR (pairwise): it cleanly bounds
  the pairwise decorrelations |Cov(D_v,D_w)| <= 1/(3vw) (verified negligible, ~1e-4), but eps_v is ~100%
  MULTI-runner (|S|>=2). So the completion identity is one order too low; the covering-min anti-concentration
  residual is a MULTI-LINEAR cancellation bound (a core arc vs the product good-set), higher-order than the
  bilinear completion machinery. Pairwise part clean; runner 1 => S255.
tags:
  - lrc14
  - covering-min
  - anti-concentration
  - completion-identity
  - bilinear
  - multilinear
  - gowers
  - decisive
---

# The completion identity is bilinear; the residual is multi-linear

**opus-2026-07-11-S262.** Owner: apply the LRCFourierCompletion cancellation bound to `Σ_h b_h ĝ(−hv)`. Doing
so decisively locates the analytic hardness of the covering-min: it is **one order above** what the completion
identity provides.

## The completion identity and the mapping

The LRCFourierCompletion completion identity (LEM-022, tasks B.1–B.3) is

> `C_w = b²/q + (1/q) Σ_{h≠0} B̂(h)·conj(B̂(w⁻¹h))`,  `|C_w − b²/q| ≤ 5q(log q)²/P(w)`,

a **bilinear** (pairwise) bound: `C_w` is a band-vs-`w`-dilate correlation = the pairwise runner overlap
`Cov(D_v, D_w)`. Expanding `1_{G'} = ∏_w(1 − 1_{D_w}) = Σ_S (−1)^{|S|} ∏_{w∈S} 1_{D_w}`,

`ε_v |G'| = Cov(1_{D_v}, 1_{G'}) = Σ_S (−1)^{|S|} Cov(1_{D_v}, ∏_{w∈S} 1_{D_w})`.

The `|S|=1` terms are exactly the completion-identity quantities `Cov(D_v, D_w) = Σ_{k≠0} b_{vk}b_{wk}`, and for
`v` coprime to `w` the **clean** bound `|Cov(D_v, D_w)| ≤ Σ_k 1/(π²vw·k²) = 1/(3vw)` holds — from `b_h` decay
alone, no fancy machinery.

## The decisive finding: ε_v is multi-linear

Verified (FFT on `G'`, `D = 13860`): **ε_v is ~100% from `|S|≥2`, not the pairwise `|S|=1`.**

| core | v | ε_full | `|S|=1` (pairwise) | `|S|≥2` (multi) |
|---|---|---|---|---|
| {41,73} | 41 | +0.0192 | −0.0001 | +0.0193 (**101%**) |
| {29,31} | 29 | −0.0405 | +0.0001 | −0.0406 (**100%**) |
| {1,17,47,53,71,89} | 17 | +0.0388 | +0.0019 | +0.0369 (**95%**) |

The pairwise `|S|=1` term is `~0.0001` (negligible), and the core-core pairwise `Cov` (e.g.
`Cov(D_41,D_73) = −3·10⁻⁵`, bound `1/(3·41·73)`) is clean and tiny. So the completion identity closes the
**pairwise** structure (core–core *and* the `|S|=1` core–noncore term) — but that is **not where ε_v lives**.

**Why multi dominates:** the `|S|=2` term is `(6/7)^{r−2} Σ_{pairs} b·b` — there are `~C(r,2) ≈ 60` pairs
(`r=11`), each non-negligible, and `(6/7)^{r−2}/(6/7)^{r−1} = 7/6`, so `|S|=2` is `~70×` the `|S|=1` term. The
**multi-runner resonances** (frequency `hv` reached as a combination of ≥2 non-core speeds) dominate.

## Net

Applying the completion identity: it **bilinearly** bounds `Cov(D_v, D_w) ≤ 1/(3vw)` cleanly — establishing
pairwise independence of all coprime runner pairs (verified negligible) — but `ε_v` (the core arc vs the
**product** good-set) is **100% multi-linear** (`|S|≥2`). So the completion identity is **necessary but one
order too low**: the covering-min anti-concentration residual is a **multi-linear cancellation bound** — the
correlation of a core arc with combinations of ≥2 non-core runners, a higher-order (Gowers-type) object beyond
the bilinear completion machinery.

This decisively locates the analytic hardness of the covering-min, and it is *consistent* across the S253–S262
arc: every bilinear/pairwise or magnitude tool (balance, single dual certificate, far-peel, second moment,
naive Erdős–Turán, mollification-magnitude, completion identity) captures the pairwise/low-order part cleanly,
but the residual is the **multi-way entanglement of the good set with the core** — the same "multi-runner
resonance" the fleet has long flagged (the Minkowski-tail / entanglement threads, tasks #42–#43). The honest
state: extremizer proved (S255); pairwise structure clean (completion identity); the open crux is a
**multi-linear (≥3-way) cancellation** for the coprime core against the good-set product, with runner 1 folded
into S255.

→ opus-S261 (signed cancellation — shown here to be multi-linear), LRCFourierCompletion / tasks B.1–B.3 (the
bilinear identity — necessary, insufficient), the entanglement threads (tasks #42–#43 — the same multi-way
object), opus-S259/S255, s558o. Files:
`lrc14_completion_identity_is_bilinear_residual_is_multilinear_opus_S262.py` (+`.out`).
