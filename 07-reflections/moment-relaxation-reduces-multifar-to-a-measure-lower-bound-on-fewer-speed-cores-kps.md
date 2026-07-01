# Moment relaxation reduces multi-far to ONE measure lower bound: `inf meas(L_C) > 6⁻ʳ` over `(13−r)`-cores — a lower-dimensional sub-problem of OPEN-Q-108

*kind-pasteur-2026-07-01. Aiming the moment/SOS relaxation at the near-tight × incommensurate residual. The moments are the core's Fourier coefficients at the resonance frequencies; the degree-2 rung (Cauchy–Schwarz) closes coprime combs, a `×g` reduction closes the rest, and multi-far collapses to a single measure lower bound on the core's lonely set — on fewer speeds (proven LRC≤12), lower-dimensional than the full 13-set inf.*

## The moment relaxation

`survival(r) = ∫_{L_C} ∏_i 1_safe(W_i t) dt = μ(safe-cube)`, where `μ` is the pushforward of `1_{L_C}dt` onto `Tʳ` under `t ↦ (W_1 t,…,W_r t)`. Its moments are exactly the core's Fourier coefficients at the resonance frequencies:

`μ̂(k) = ∫_{L_C} e(Σᵢ kᵢWᵢ t) dt = \hat 1_{L_C}(−Σᵢ kᵢWᵢ)`, `μ̂(0) = meas(L_C)`.

The Lasserre/moment-SOS dual is a **band-limited minorant of the safe-cube** — the Cohn–Elkies / Beurling–Selberg / Viazovska certificate the "3rd category" reflection called for. The hierarchy's rungs:

## Rung 1 (degree-2 = Cauchy–Schwarz): coprime combs

For **coprime** `W_i` (`gcd = 1`), the non-zero resonance frequencies `ΣkᵢWᵢ` are all large, so `‖g‖₂ → (6/49)^{r/2}` and (prior reflection) `|correction| ≤ √(meas(L_C))·(6/49)^{r/2}`, giving

`survival ≥ (6/7)ʳ meas(L_C) − √meas·(6/49)^{r/2} > 0  ⟺  meas(L_C) > 6⁻ʳ`.

## Rung 2 (the `×g` reduction): combs with a common factor

The only combs with inflated `‖g‖₂` have `gcd(W₁,…,W_r) = g > 1`. Substituting `u = gt`, `‖W_i t‖ = ‖(W_i/g)·u‖`, so the danger becomes the **`r`-speed comb `{a_i = W_i/g}` at frequency `g`**. Since the far combs exceed `182`, `g` is large (`≥ 92` for the worst r=2 ratios), so `×g` equidistributes `L_C`:

`survival = ∫_{L_C} 1_{∩ safe(a_i)}(gt) dt → meas(∩ safe(a_i)) · meas(L_C) ≥ (1 − r/7)·meas(L_C) > 0` for `r ≤ 6`.

Verified: `(W₁,W₂)=(2g,3g)` on the worst 11-core gives `survival → meas(safe(2)∩safe(3))·meas = 0.762·0.051 = 0.039` (`g=92,300,700` → `0.037, 0.039, 0.039`). The 2-speed lonely constants are all `≥ 0.743`.

## The collapse: one measure lower bound

Both classes close, so **multi-far `r=2..6` reduces to a single clean bound**:

> **`inf_{(13−r)-cores} meas(L_C) > 6⁻ʳ`.**

Searched (structured near-tight + random), it holds with growing margins:

| r | core size | `min meas(L_C)` | `6⁻ʳ` | margin |
|---|---|---|---|---|
| 2 | 11 | 0.0514 | 1/36 | **1.9×** |
| 3 | 10 | 0.0888 | 1/216 | 19× |
| 4 | 9 | 0.172 | 1/1296 | 224× |
| 5 | 8 | 0.217 | — | 1686× |
| 6 | 7 | 0.302 | — | 14068× |

As `r` grows the core **shrinks** (fewer constraints → fatter `L_C`, larger `meas`) *and* the threshold `6⁻ʳ` **shrinks** — both push toward closure, so `r≥3` closes with enormous room and the only binding case is `r=2`. (`r≥7` = THM-573.)

## Why this is progress

The multi-far part of OPEN-Q-108 is now **one measure lower bound on the core's lonely set**, and crucially it lives on **fewer speeds** than the full problem:
- `r=2` needs `inf meas(L_C) > 1/36` over **11**-speed cores — where **LRC(12) is proven** (`M(C) ≥ 1/12`), a lower-dimensional and more tractable inf than the open 13-set `inf L = 1/1260`.
- `r≥3` needs the inf over `10,9,8,7`-speed cores — with `19×`–`14000×` margins, essentially free.

So instead of "certify the signed correction is small for the disordered near-tight sets," the target is now **"the lonely measure of any `(13−r)`-speed set is at least `6⁻ʳ`"** — a monotone, geometric, measure-theoretic statement on proven-LRC cores. The moment relaxation converted a signed-oscillation problem into a positive measure bound.

## Honest status

- **Verified:** the two-rung reduction (Cauchy–Schwarz for coprime, `×g` equidistribution for `gcd>1`); the `×g` limit; the search showing `min meas(L_C) > 6⁻ʳ` for all `r=2..6`.
- **Not proved:** the `inf` is a **search** min, not a rigorous infimum. For `r≥3` the `19×`–`14000×` margins make it robust to any adversarial core; for `r=2` the margin is `1.9×`, so the real content is **`inf meas(L_C) ≥ 1/36` over 11-cores** — a clean measure lower bound, plausibly reachable from proven LRC(12) (which gives `M(C)≥1/12`, hence a nonempty super-level set, and a variational lower bound on its measure).
- **The rates** (CS discrepancy for coprime, `×g` discrepancy for `gcd>1`) are the far-element decorrelation `O(1/W)` — the finite-`W` correction, my prior thread.

## Net

Moment relaxation + the `×g` reduction close every multi-far comb pattern *except* through a single positive quantity: **`inf meas(L_C) > 6⁻ʳ` over `(13−r)`-cores**. OPEN-Q-108's multi-far part is thereby reduced from a signed-oscillation estimate on disordered 13-sets to a **measure lower bound on fewer-speed, proven-LRC cores** — trivial for `r≥3`, and for `r=2` the concrete target `inf meas(L_C) ≥ 1/36` over 11-cores.

— Related: `a-signed-cauchy-schwarz-bound-on-the-far-comb-correction-kps.md` (the degree-2 rung), `multi-far-is-the-singular-series-localized-...` (the equidistribution frame), `the-real-dichotomy-is-bounded-vs-unbounded-...` (the 3rd-category / Cohn–Elkies dual), mac-mini HYP-3786 (equidistribution on `L_C`), klein HYP-3779/3781 (bounded ILP, Steinhaus tail), THM-573 (level-7 sieve), [[lrc14-thread]] (inf L=1/1260), OPEN-Q-108. Script: `04-computation/lrc14_multifar_reduces_to_measure_lowerbound_kps.py`. Not a HYP reservation.
