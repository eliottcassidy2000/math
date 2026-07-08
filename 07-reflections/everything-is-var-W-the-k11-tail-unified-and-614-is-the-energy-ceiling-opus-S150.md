---
source: opus-2026-07-08-S150
status: PROVED the exact unification far<=E[W]^2 <=> Var(W)<=near (both open k=11 routes are
  one Var(W) bound); reduced brick(B) to the resonance lemma Var<=c*R2 and explained PZ-fails/
  D3-succeeds quantitatively; 614 = the tail additive-energy ceiling (why the tail is easier).
  The one remaining mile is Var(W)<=c*R2. Concurrent with klein-S179/S180, kps-S78/S80, mac-mini.
tags:
  - lrc14
  - covering-floor
  - variance
  - additive-energy
  - resonance-lemma
  - k11-tail
---

# Everything is Var(W): the k=11 tail unified, and 614 is the energy ceiling

**opus-2026-07-08-S150 (HYP-5377).** Owner: prove `far <= E[W]^2` for spread families, prove
brick (B) (`R2 <= 614 => D3 >= bar`), and look back for connections to 614. The three turn out
to be **one** question, and this note pins it down.

## 1. The unification (PROVED): both open routes bound the same Var(W)

The k=11 covering-tail (the last open density-floor leg) had two competing routes:
- **mine (S149):** on mac-mini's near/far split `E[W^2] = near + far`, `far <= E[W]^2` gives
  the sharp tail `PZ >= E[W]^2/(near + E[W]^2) ~ 0.53`;
- **kps brick (B):** `R2 <= 614 => D3 >= bar` via the resonance lemma `Var(W) <= c*R2`.

They are the **same statement**. Since `far = E[W^2] - near` identically,

> **`far <= E[W]^2`  ⟺  `Var(W) <= near`.**

(`Var(W) = E[W^2] - E[W]^2`, so `far <= E[W]^2 ⟺ E[W^2] - near <= E[W]^2 ⟺ Var <= near`.)
And brick (B) is `Var(W) <= c*R2` (via `PZ = 1/(1 + Var/E[W]^2)`). So **the entire k=11 tail is
a single bound on `Var(W)`** — the covering variance of the uncovered measure. My `far <= E[W]^2`
is the `Var <= near` sub-case (`near ~ 0.025`); brick (B) is the fuller `Var <= c*R2 <= c*614 ~
0.037` case. Verified exact on every family.

This is the clean picture the leg was missing: **there is one analytic mile, `Var(W) <= c*R2`,
and it closes everything** (far, brick B, the tail). The covariance reduction that shows it:
with `p(y) = P_x(y uncovered)`,

> `far - E[W]^2 = int int_{disjoint} Cov(y1,y2) - int int_{near} p(y1)p(y2)`,

so `far <= E[W]^2 ⟺ int_{disjoint} Cov <= int_{near} p*p`. In the decorrelated (spread) limit
`p(y) -> E[W]` and `Cov -> 0`, so **`far -> (5/7)E[W]^2`** (a `(2/7)E[W]^2` buffer; verified:
spread families sit at `far/E[W]^2 = 0.59-0.67`, disjoint-`Cov` *negative* — arcs anti-correlate
like iid). Compact/2-block families have positive disjoint covariance and `far > E[W]^2`.

## 2. Why PZ fails and D3 succeeds (kps-S78, quantified)

Everything reduces to `Var(W) <= c*R2`, but the constant `c` is **not** exact:
`Var(W)/R2 ∈ [4.9, 7.0]·10^{-5}` (mean `5.8e-5`) across 243 k=11 families — a `±20%` scatter.
This is precisely kps-S78's "coupling-tight" obstruction, now numeric:

- **Via `PZ` (degree 2):** `PZ >= bar ⟺ Var <= (1/bar - 1)E[W]^2 = 2.019 E[W]^2`. The
  *decoupled* worst case — `max(Var/R2)·614 = 0.043` variance, `min E[W] = 0.144` — gives
  `PZ >= 1/(1 + 0.043/0.144^2) = 0.324 < bar` (margin `-0.007`). PZ genuinely fails, because
  the worst `Var/R2` and the `R2 = 614` config do not co-occur, and decoupling loses the gap.
- **Via `D3` (degree 3, opus-S148):** `D3 > PZ`, and the extra `+0.0735` block margin absorbs
  the `-0.007` loss: `min_E D3 = 0.458` over `R2 <= 614` (mac-mini, margin `+0.127`). D3
  recovers what PZ drops.

So the k=11 tail closes **iff** the resonance lemma `Var(W) <= c*R2` holds with `c <= ~7e-5`,
and it must be routed through `D3` (not `PZ`) — exactly why my S148 degree-3 floor was the
enabling step. The lemma itself (the `Var(W) ≈ c*R2` proportionality, remarkably stable but
`±20%`) is the one open mile, shared with THM-656 (`Var(F) = R2*V1` for the tent) and
klein-S179 (`Var(W) ~ 5.67e-5*R2`).

## 3. 614 is the tail additive-energy ceiling (why the tail is easier)

`R2(AP_n) = (n-1)n(2n-1)/3` (the block maximizes reduced additive energy):
`R2(AP_10) = 570`, `R2(AP_11) = 770`. Brick (A) (kps THM-662): a primitive 11-set of
prim-diam `>= 16` has `R2 <= 614 = R2(AP_10) + 44`, at the block-plus-far-point `{0..9} u {16}`
(`+20` from the far point's 10 new differences, `+24 = 4T` from `T = r_7+r_8+r_9 = 3+2+1 = 6`
overlaps onto the block's short differences). So:

> **`770 - 614 = 156` is the energy cost of forcing the diameter from 10 (block) to `>= 16`.**

And this is *why the tail is easier than the compact case*. `D3` **decreases** in `R2` (more
additive energy = spikier `W` = lower floor). The compact block (`R2 = 770`) is the **global**
`D3`-minimizer (`D3 = 0.4048`); forcing diameter drops `R2` by `>= 156`, which **raises** the
`D3` floor to `>= 0.458` (`+0.053`). The tail clears with more room than the compact block —
the diameter/energy complementarity (klein THM-656) made quantitative: spreading buys margin.
`614` has no arithmetic magic; it is exactly "how much additive energy survives a forced
diameter", and its being `< 770` is the whole point. (The general-`k` pattern: the tail base is
`R2(AP_{k-1})`; compact `R2(AP_k)`; gaps `200/242/288` at `k = 11/12/13`.) A curio from the
tournament half: `614/2016` is the `p=13` log-supermodularity/FKG violation count (HYP-533) —
same `R2 <= 614`-adjacent additive object, different problem; noted, not load-bearing.

## 4. Where this leaves the leg (honest)

- **PROVED (new):** `far <= E[W]^2 ⟺ Var(W) <= near`; the two k=11-tail routes are one
  `Var(W)` bound; the decorrelated limit `far -> (5/7)E[W]^2`.
- **REDUCED:** the entire k=11 tail (and far, and brick B) to the single resonance lemma
  `Var(W) <= c*R2` (`c <= ~7e-5`), routed through `D3` (PZ is coupling-tight, `-0.007`; D3
  gives `+0.127`). This is the same mile as THM-656 / klein-S179.
- **STRUCTURAL:** `614` = the tail energy ceiling; `770 - 614 = 156` = the energy cost of
  spreading = why the tail clears with `+0.053` extra `D3`.
- **NOT proved:** the resonance lemma `Var(W) <= c*R2` itself (the `Var(W) ≈ c*R2` leading
  order is empirically stable but `±20%` — the coupling kps-S78 flagged; a full Fourier
  derivation of `c` and the error is the remaining work). k=11 closes when it does.
- Files: `lrc14_far_covariance_opus_S150.py`, `lrc14_var_resonance_opus_S150.py` (+outs).
- Builds on / cites: opus-S148 (D3, the enabling floor), opus-S149 (near/far, exact near),
  kps-S78/S80 + THM-662 (bricks, coupling-tightness), mac-mini THM-661 (near/far, unified D3),
  klein THM-656 + S179/S180 (`Var ~ R2`, factorial-moment `E[W]` floor).
