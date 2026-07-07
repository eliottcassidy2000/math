# The reverse-Markov reduction: the density floor reduces to a mean max-gap `E[maxgap] > 1/7`

**opus-2026-07-07-S133.** Owner: *work the crux, think reverse-Markov reduction.* That hint cracks
the density-floor node open: the **reverse Markov inequality** reduces the floor from a hard *tail
probability* to a **mean** `E[maxgap] > 1/7`, whose binding value (the AP) is proven exactly with a
~48% margin. This is strictly cleaner than the S130 `E[U]` route and the raw `μ_{1/7}` tail.

## The reduction

Let `Y = maxgap(x)` = the largest circular gap of the cluster phases `{frac(e·x) : e∈E}`, with
`Y ∈ [0, B]` (`B ≤ 1`). The **reverse Markov inequality** (for a bounded nonnegative `Y`): for
`a < E[Y]`,
> `P(Y > a) ≥ (E[Y] − a) / (B − a)`.

*(Proof: `E[Y] = E[Y; Y≤a] + E[Y; Y>a] ≤ a·P(Y≤a) + B·P(Y>a)`, solve for `P(Y>a)`.)* With `a = 1/7`,
`B = 1`:

> **`μ_{1/7}(E) = P(maxgap > 1/7) ≥ (E[maxgap] − 1/7) / (1 − 1/7) = (7/6)(E[maxgap] − 1/7)`.**

So **density-floor positivity ⟸ `E[maxgap] > 1/7`** — a *first-moment* bound on a single order
statistic, not the tail. (Compare: S130's `μ_{1/7} ≥ E[U]` used `E[Σ(gap−1/7)_+]`, a truncated sum;
`E[maxgap]` is the cleaner mean, and the reverse-Markov constant `7/6` is explicit.)

## The AP value, proven exactly — but the AP is NOT the minimizer (corrected)

> **⚠ Correction (census, same session).** Unlike `μ_{1/7}` (which the AP uniquely *minimizes*, S130),
> `E[maxgap]` is **NOT** AP-minimized. Exhaustive census (k=8,9,10) finds shapes strictly below the AP
> (38 at k=8, e.g. an endpoint-stretched AP `{0..6,9}`); the true `k=13` minimum found is
> `E[maxgap] ≈ 0.2047` at `{0,2,3,…,10,12,17,28}`, below the AP's `0.2114`. A structural distinction:
> **the AP minimizes the max-gap TAIL (`μ_{1/7}`) but not the max-gap MEAN.** The reduction survives
> regardless — the true `inf_E E[maxgap] ≈ 0.205` is still `≫ 1/7` (margin `0.062`), giving
> `μ_{1/7} ≥ (7/6)(0.205 − 1/7) ≈ 0.072 > 0`. So the target is a **direct** `inf_E E[maxgap] > 1/7`,
> NOT AP-minimality.

The AP value is still a clean exact reference (`lrc_Emaxgap_exact_opus_S133`, rational three-gap
piecewise integration — breakpoints = Farey `m/d`, `d≤k`, order changes/wraps, plus pairwise gap
crossings; a rigorous computation, not a grid):

| k | `E[maxgap(AP_k)]` exact | ≈ | margin over `1/7` | reverse-Markov `μ_{1/7} ≥` |
|---|---|---|---|---|
| 8 | 43/140 | 0.3071 | 0.1643 | 0.1917 |
| 9 | 47/168 | 0.2798 | 0.1369 | 0.1597 |
| 10 | 19/72 | 0.2639 | 0.1210 | 0.1412 |
| 11 | 151/630 | 0.2397 | 0.0968 | 0.1130 |
| 12 | 796/3465 | 0.2297 | 0.0869 | 0.1014 |
| **13** | **93/440** | **0.2114** | **211/3080 ≈ 0.0685** | **1477/18480 ≈ 0.0799** |

So `E[maxgap(AP_13)] = 93/440 > 1/7` **with a 48% margin**, and the reverse-Markov density floor is
`μ_{1/7}(AP_13) ≥ 1477/18480 > 0` — **exact and rigorous** at the binding orbit.

## What remains (one clean step)

The floor `μ_{1/7}(E) > 0` for all clusters now needs only:

> **`inf_E E[maxgap(E)] > 1/7`** — the mean max-gap over every integer cluster exceeds `1/7`.
> (Empirically `inf ≈ 0.205`, at a stretched shape — *not* the AP; margin `≈ 0.062`.)

This is a **direct mean bound**, and it is the honest open step. AP-minimality is *not* available
(the AP is not the `E[maxgap]`-minimizer). But a mean is far more amenable to a moment/Fourier
argument than the `μ_{1/7}` tail was:

- `E[maxgap] ≥ E[Σ gap²]` (length-biased mean gap) is **too weak** — `E[Σ gap²] = 2/(k+1) = 1/7`
  exactly at `k=13` for independent points, and *less* for regular clusters. So the bound must use the
  max-gap's **excess over the length-biased mean**, which is the three-gap "a few big gaps near
  small-denominator rationals" contribution (`near x = p/q`, `q ≤ 6`, the phases cluster into `≤ q`
  groups leaving a gap `≥ 1/q ≥ 1/6 > 1/7`; the measure of such `x`, weighted by that gap size, is the
  excess). Making *that* a uniform `> 1/7 + ε` over all `E` is the crux — a mean over the modular
  clustering of `E`, i.e. squarely the additive↔multiplicative mediation of mac-mini-S15.

## Why this is the right frame (ties to mac-mini-S15)

This is the quantitative three-gap rigidity S15 called for, now with an explicit lever: the AP orbit
has `E[maxgap] = 93/440`, a *continued-fraction/Farey* quantity (the exact rationals `43/140, …,
93/440` are three-gap integrals), and reverse-Markov converts that mean into the density floor. The
`n=13`-specific number `1/7 = 2·(1/14)` is the "gap wide enough for a `1/14`-margin runner," and
`2/(k+1) = 1/7` at `k=13` is exactly the independent-points break-even — the AP beats it because it is
*more* regular in the mean but has *heavier* max-gap tails near rationals. The margin `93/440 − 1/7 =
211/3080` is the quantitative "detuning room."

## Bottom line

Reverse Markov reduces the LRC(14) density-floor crux to **`E[maxgap] > 1/7`**, a mean bound; the AP
orbit value `E[maxgap(AP_13)] = 93/440` is **proven exactly** (48% margin). The AP is *not* the
minimizer (census: true `inf_E E[maxgap] ≈ 0.205` at a stretched shape), but `inf ≫ 1/7`, so the
honest open piece is a **direct** `inf_E E[maxgap] > 1/7` — a statement about *averages of orbit
max-gaps* (not the tail, not an extremal-at-AP), squarely in classical three-gap territory and open
to a moment/Fourier attack the `μ_{1/7}` tail resisted.
