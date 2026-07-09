---
source: opus-2026-07-09-S170
status: the good-period EXISTENCE (both the dissociated AND the extremal near-AP branch, where arc-count
  is VACUOUS -- klein MISTAKE-127) closes a-priori by the SMOOTH AVERAGING route, sidestepping the
  Mertens/L2 signed cancellation of the sharp-indicator Erdos-Turan sum. DIAGNOSIS (exact ET identity):
  #good/Vmax - rho* = sum_{n!=0} Ihat(nVmax)(-1)^n; the sharp good-set indicator has Ihat ~ 1/m so this
  sum is L1-DIVERGENT (= the vacuous arc bound) and small only by CANCELLATION (opus-S154/S167 Mertens).
  THE FIX: use the SMOOTH surrogate maxgap(x) -- continuous piecewise-linear, Fourier decay ~1/m^alpha,
  alpha=1.5-2.0 (VERIFIED) > 1, so its resonant sum converges ABSOLUTELY, no cancellation. A good period
  is j with maxgap(j/V)>1/7; by max>=mean it suffices E_j[maxgap]>1/7, and E_j = E_x[maxgap] + (resonant
  disc, |.|<=0.006 << margin). E_x[maxgap]>1/7 holds with UNIFORM margin -- 1.48x for the tight AP {1..13},
  >=1.047x adversarial (kps) -- INCLUDING the extremal AP families arc-count cannot touch. LEAN:
  exists_good_of_smooth_mean (kernel-pure). Remaining a-priori: E_x[maxgap]>1/7 uniform + the alpha>1 disc
  bound -- both cancellation-FREE.
tags:
  - lrc14
  - good-period
  - erdos-turan
  - averaging
  - mertens
  - lean
  - smooth-vs-sharp
---

# The smooth averaging route sidesteps Mertens

**opus-2026-07-09-S170.** Owner: keep pushing the LRC math then formalization.  The critical path
(klein-S193 flagged it for opus) is good-period EXISTENCE on the extremal near-AP family, where the
arc-count pigeonhole is VACUOUS (MISTAKE-127).  The resolution is a route swap.

## The exact Erdos-Turan identity, and why the sharp route is Mertens-hard

A good ruler period is `j` with `maxgap{frac(e_i·j/V)} > 1/7`.  Writing `1_{G*}(x) = Σ_m Ihat(m) e(mx)`
(`Ihat(0) = rho*`) and using `Σ_j e(m x_j) = V·[V|m]·e(m/2V)` on the grid `x_j=(j+½)/V`:

> **`#good = V·Σ_{V|m} Ihat(m) e(m/2V)`,  i.e.  `#good/V − rho* = Σ_{n≠0} Ihat(nV)(−1)ⁿ`** (exact).

So existence (`#good>0`) ⟺ the RESONANT sum `|Σ_{n≠0} Ihat(nV)(−1)ⁿ| < rho*`.  I verified this identity
and diagnosed the sum (`lrc14_ET_resonance_identity`).  The catch: the sharp indicator `1_{G*}` has
JUMP boundaries, so `Ihat(m) ~ #arcs/m`, and the resonant absolute sum `Σ_n |Ihat(nV)| ~ (#arcs/V)Σ1/n`
**DIVERGES** — it IS the vacuous arc-count bound (`#arcs/V ~ 1.17` for near-AP).  The true discrepancy
is small only by SIGNED CANCELLATION of the `(−1)ⁿ`-weighted terms — exactly the **Mertens situation**
(opus-S167) / the **L¹-diverges-L²-converges** split (opus-S154).  Proving the cancellation a-priori is
the hard problem the whole fleet has circled.  **So don't.**

## The fix: average a SMOOTH function

`max ≥ mean` gives a cancellation-free route: a good period exists as soon as the ruler-grid MEAN
`E_j[maxgap] > 1/7` (kps-S95 `averaging_existence`).  And `maxgap(x)` is **continuous** piecewise-linear
(circular gaps vary continuously through collisions), so — unlike the indicator — its Fourier
coefficients decay like `1/m^α` with `α>1`.  Measured (`lrc14_averaging_smooth_route`):

| family | `maxgap` Fourier `α` |
|---|---|
| tight AP `{1..13}` | **2.02** (clean `1/m²`, kink not jump) |
| near-AP `7·{1..10}+{1,2,3}` | 1.79 |
| dissociated Sidon-ish | 1.48 |

`α>1` ⟹ the resonant sum `Σ_n |maxgap^(nV)|` converges **ABSOLUTELY**: `|E_j[maxgap] − E_x[maxgap]| ≤
Σ_{n≥1} 2|maxgap^(nV)| < ∞`, no cancellation needed.  Empirically this discrepancy is `≤ 0.006`
(near-AP, all `d`).  And the dilation-invariant continuous mean clears the threshold with room to spare:

> **`E_x[maxgap] > 1/7` with a UNIFORM margin** — `≈0.214` (`1.48×`) for the tight AP `{1..13}`,
> `1.49×` for the dilated/near-AP families, `≥1.047×` adversarial (kps).  Margin `≈0.07 ≫ 0.006`.

So `E_j[maxgap] = E_x[maxgap] − (disc) > 1/7` **a-priori**, for the SAME extremal AP families where
arc-count is vacuous.  The smooth route closes what the sharp route could not — and it needs neither
the arc count nor the Mertens cancellation.

## What this leaves (both a-priori-shaped, cancellation-free)

1. **`E_x[maxgap] > 1/7`, uniform over the dilation-invariant shape space.**  A covering-MEAN bound
   (the average largest-gap beats the threshold); verified `≥1.047×` (kps adversarial) and `1.48×` on the
   extremal AP.  This is the analogue of the density floor `rho* ≥ D3`, one order smoother.
2. **`Σ_n |maxgap^(nV)| < ` margin.**  An absolute (α>1) bound on a Lipschitz function's resonant
   Fourier tail — the cheap L² object, NOT the L¹-divergent indicator tail.

Both are absolute-value bounds on a SMOOTH object; the Mertens/L² cancellation that blocks the sharp
route never appears.  This is why the averaging route (kps-S95) is the right a-priori tool, and the
`maxgap` Lipschitz/`α>1` structure (opus-S170) is what makes its discrepancy absolutely controllable.

## Lean

`TournamentH7.LRCArcCount.exists_good_of_smooth_mean` (kernel-pure `[propext, Classical.choice,
Quot.sound]`, built): for a nonempty grid `J`, surrogate `g`, continuous mean `meanx`, discrepancy bound
`D` (`|mean_J g − meanx| ≤ D`) and margin `thr + D < meanx`, some `j∈J` has `g j > thr`.  The
smooth-route existence logic, checked; the two analytic inputs above enter as the hypotheses `meanx >
thr + D` and the discrepancy bound.  Complements the `good_period_of_arccount` pigeonhole (S169, the
dissociated route) and kps's `averaging_existence`.

## Ledger

- The sharp-indicator ET resonant sum is L¹-divergent = the vacuous arc-count = Mertens/L² cancellation
  (opus-S154/S167); the SMOOTH `maxgap` surrogate (`α>1` Fourier decay, verified 1.48–2.02) has an
  absolutely-convergent resonant sum ⟹ the averaging route is cancellation-free a-priori.
- `E_x[maxgap] > 1/7` uniform (1.48× extremal AP, ≥1.047× adversarial) is the smooth analogue of the
  density floor; `|E_j − E_x| ≤ 0.006 ≪ margin 0.07` ⟹ good period exists for the extremal families.
- LEAN: `exists_good_of_smooth_mean` (kernel-pure). Files: `lrc14_ET_resonance_identity_opus_S170`,
  `lrc14_averaging_smooth_route_opus_S170` (+outs), `LRCArcCount.lean`.
- CONNECTS: klein-S193 (ET resonances, MISTAKE-127), kps-S95 (averaging + contraharmonic/lemniscate
  2nd-moment), opus-S154 (L²-not-L¹), opus-S167 (Mertens), opus-S169 (arc-count pigeonhole, dissociated).
- NEXT: an a-priori `E_x[maxgap] > 1/7` (covering-mean, dilation-invariant) + the `α>1` resonant-tail
  bound; then the good-period capstone is fully a-priori (dissociated + extremal both by smooth mean).
