---
id: HYP-9290
title: "Frontier energy cone: stop every source of the signed ROOT transport at its first seed <= Y; for a fixed frontier height Y_0 (about 100-1000) the stopped energy E_(Y_0)(N) = sum_u c_(Y_0)(u)^2 stays below C N for every N. By the one-line frontier lemma R(N)^2 <= F_Y E_Y(N), this gives |R(N)| <= sqrt(F_(Y_0) C N) without any input from the termination regime"
status: >
  OPEN (CONJECTURE). Lemma F is PROVED (Cauchy-Schwarz on the stopped invariant).

  Evidence (FINITE-EXACT, every N < 2^23):
  - sup E_1000(N)/N = 2.08 and sup E_100(N)/N = 2.35 over 2^12 <= N < 2^23;
  - the unstopped ratio R(N)^2/N reaches 9.08 (N = 7910).

  The decomposition E_Y = E_0 + 2(MM_hi + MM_lo + MF), synchronous in depth, was computed at four N with N >= 6Y:
  - merges above the frontier remove 54-75% of E_0;
  - same-step merges at frozen seeds add at most 0.011 N;
  - arrivals at frozen seeds add at most 0.24 N for Y >= 100.

  The coherent climbs that make R(N)^2/N large (the termination regime of THM-4667) are arrivals at frozen seeds
  below about 30-100.

  Pre-registered 2026-10-17, before the scan beyond 2^23: tests P1-P6 below. ALL SIX PASSED (scan of every
  N < 2^26, invariant checked at every N). New record E_100/N = 2.62 at N = 48,251,574 (inside P2's 3.0).
  P7-P10, pre-registered for 2^26 <= N < 2^28 before that scan, ALL PASSED: over every N < 2^28,
  E_1000/N <= 1.30 and E_100/N <= 1.80 in octaves 2^26-2^27, while R^2/N reaches 7.48 (N = 133,970,036).
session: opus-2026-10-17-S42
source: 05-knowledge/results/root_frontier_penrose_20261017.md
scripts:
  - 04-computation/experiments/root_transport_frontier_energy_20261017.py (every N, several Y; checks R(N) at every N)
  - 04-computation/experiments/root_transport_frontier_profile_20261017.py (profiles in Y at fixed N)
  - 04-computation/experiments/root_transport_frontier_decomposition_20261017.py (cross-term decomposition)
related: [THM-4667, HYP-9288, HYP-9283, THM-4658, THM-4581]
---

# HYP-9290: frontier energy cone of the stopped ROOT transport

## Setting

The notation is that of THM-4667 and HYP-9288. The ROOT current is `R(N) = H(4N) - H(N)`, and `c_0` are its initial
seed coefficients. The odd transport is `s -> P(s) = seed((3s+1)/2^a)`, with sign `+` if `a = 2` and `-` if `a = 1`.
`f(u) = (-1)^(total stopping time of u)` on the ROOT basin and0 otherwise.

**Stopped transport.** Fix a height `Y >= 1`. Every source is transported until its first seed `u <= Y`, and stays
there. Seeds `<= Y` never move. Define `c_Y(u)` by summing only sources that actually hit the frontier; nonhitting sources remain in a separate moving remainder. Since initial support is finite, all hitting contributions arrive by some finite depth, without requiring every source to hit. The frontier energy is
`E_Y(N) = Σ_u c_Y(u)^2`, and `F_Y = #{seeds <= Y}` (`F_100 = 38`, `F_1000 = 375`).

**Lemma F (PROVED).** `R(N) = Σ_(u <= Y) f(u) c_Y(u)`; hence `R(N)^2 <= F_Y E_Y(N)`.

*Proof.* At finite depth, the invariant is the sum of the frozen and moving parts. Every initially ROOT-supported source eventually hits the frontier; every source that never hits has f=0. Thus the f-weighted moving remainder vanishes after the finitely many ROOT-supported arrivals, even if nonhitting sources persist. The frozen vector stabilizes after all hitting arrivals. Cauchy-Schwarz over the F_Y frozen seeds proves the inequality. This does not assert universal entry. ∎

## Statement

- **(A, frontier cone)** There are `Y_0` and `C` with `E_(Y_0)(N) <= C N` for all `N >= 1`. The data suggest
  `Y_0 = 1000`, `C = 2.1` on the tested large-source range N>=4096. The original C=2.4 suggestion at Y_0=100 is REFUTED by the incoming N=48,251,574 record E_100/N>2.62; C=3 is the revised large-source test. These are not all-N constants: at N=1, E_Y(1)=10 for every Y>=3, so the all-N formulation needs C>=10 and a separately paid finite prefix. Consequence of the stated all-N hypothesis: `|R(N)| <= sqrt(F_(Y_0) C N)`. The numerical bound `28 sqrt(N)` uses Y_0=1000,C=2.1 only on its stated large-source range.
- **(A', logarithmic)** The same with `C` replaced by `C log log (N + 16)`.
  - Heuristic: if the `F_Y` frontier coefficients fluctuate like independent Gaussian sums, `E_Y/N` deviates from its
    mean by a factor `1 + O(F_Y^(-1/2))`.
  - The running maximum over about one effectively independent value per octave then drifts by
    `O(sqrt(log log N / F_Y))`.
  - This is a model, not a theorem.
- **(B, mechanism)** At finite depth d, the exact identity is `E_frozen,d + E_moving,d = E_0 + 2(MM_hi,d + MM_lo,d + MF_d)`. The displayed terminal E_Y identity applies when the moving vector vanishes:
  - `MM_hi` is the merges of two moving currents above `Y`;
  - `MM_lo` is the same-step merges at frozen seeds;
  - `MF` is the arrival of moving current at a frozen seed already carrying current.

  Require `MM_hi,d + MM_lo,d <= 0` and `MF_d <= c N` uniformly at every finite depth d, for Y>=Y_0. Then E_Y<=E_0+2cN by nonnegativity of the moving energy and eventual stabilization of the frozen vector. This gives C=34/7+2c+o(1) at large N; an all-N constant must also cover the finite prefix. Infinite cross-term sums are not assumed to converge.

`E_0(N)/N -> 34/7 = 4.857` (the octave trace `T_0(X) = (51/7)X^2` of HYP-9288).

## Why

1. **The termination regime.** Transport after freezing never enters `E_Y`; arrival terms at frozen seeds still do. The coherent merging that makes `R(N)^2/N` large happens when a
   few dozen frontier currents finally coalesce into ROOT. In the decomposition it is an *arrival* term against frozen
   seeds below about 30-100. For example, at `N = 6904411`, `Y = 30`: `MF = +2.92 N`, which accounts for the climb
   from `E_100 = 1.37 N` to `E_30 = 4.09 N`.
2. **Above the frontier the transport is dissipative.** `MM_hi` lies between `-3.65 N` and `-2.64 N` at every computed
   `(N, Y)` with `Y <= 10^4 <= N/6`. It is concentrated in the first third of the stopped depth (the ideal regime of
   THM-4667, sibling anti-correlation).
3. **Relation to HYP-9288.** (A) is not implied by HYP-9288 and does not imply it. Its advantage: it only needs the
   transport *down to a fixed scale*. The old S41 matching test omitted ideal-only support. The [complete repair](../results/root_sparse_boundary_module_20261017i.md), R1--R5, gives O(d^2) dyadic error sites, a growing-depth window and a first-return renewal identity; it does not yet bound this stopped energy uniformly.
4. **Conditional power saving.** With `Y = N^θ`, `|R(N)| <= sqrt(F_Y E_Y)`, which is `O(N^((1+θ)/2))` whenever
   `E_(N^θ)(N) = O(N)`.

## Evidence (FINITE-EXACT; `R(N)` checked at every `N`)

| octave | `2^12` | `2^13` | `2^14` | `2^15` | `2^16` | `2^17` | `2^18` | `2^19` | `2^20` | `2^21` | `2^22` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `max E_100/N` | 2.35 | 1.68 | 1.80 | 1.65 | 1.85 | 1.22 | 1.60 | 1.48 | 1.24 | 1.46 | 2.03 |
| `max E_1000/N` | 2.08 | 1.92 | 1.75 | 1.41 | 1.42 | 1.44 | 1.38 | 1.45 | 1.50 | 1.52 | 1.31 |
| `mean E_1000/N` | 1.91 | 1.74 | 1.49 | 1.30 | 1.27 | 1.31 | 1.20 | 1.28 | 1.34 | 1.30 | 1.19 |
| `max E_30/N` | 2.40 | 2.65 | 1.27 | 2.55 | 1.66 | 2.29 | 2.52 | 2.25 | 1.16 | 1.57 | 4.64 |
| `max R^2/N` | 9.08 | 7.53 | 2.55 | 3.36 | 4.23 | 2.77 | 3.92 | 1.66 | 1.56 | 0.98 | 6.06 |

Decomposition (`/N`):

| `N` | `Y` | `E_0` | `MM_hi` | `MM_lo` | `MF` | `E_Y` | `R^2` |
|---|---|---|---|---|---|---|---|
| 7910 | 100 | 4.856 | -3.250 | -0.059 | +0.235 | 1.783 | 9.080 |
| 7910 | 1000 | 4.856 | -2.739 | -0.366 | +0.164 | 1.916 | 9.080 |
| 66463 | 1000 | 4.857 | -3.267 | -0.103 | -0.240 | 1.248 | 4.226 |
| 1048576 | 1000 | 4.857 | -3.489 | -0.036 | +0.041 | 1.374 | 1.533 |
| 6904411 | 30 | 4.857 | -3.653 | -0.031 | **+2.915** | 4.088 | 6.055 |
| 6904411 | 100 | 4.857 | -3.643 | -0.037 | +0.195 | 1.372 | 6.055 |
| 6904411 | 1000 | 4.857 | -3.547 | -0.039 | -0.101 | 1.170 | 6.055 |

(The table entries are `2·MM/N`, i.e. the actual energy contributions.)

## Pre-registered tests (written before the scan of `2^23 <= N < 2^26`)

- **P1.** `max E_1000(N)/N <= 2.08` over `2^23 <= N < 2^26` (no record beyond the `[2^12, 2^23)` maximum).
- **P2.** `max E_100(N)/N <= 3.0` over the same range (log log allowance over 2.35).
- **P3.** `max E_10000(N)/N <= 1.55` over the same range (the decline toward a plateau continues).
- **P4.** Each octave mean of `E_1000/N` lies in `[1.0, 1.5]`.
- **P5.** The `Y = 30` octave mean (3.19 at `2^22`) returns to `<= 2.0` in at least two of the three new octaves.
- **P6.** At the arg-max of `E_1000/N` and at the arg-max of `R^2/N` in each new octave, the decomposition at
  `Y = 1000` has `MM_hi <= -3.0 N` and `|MF| <= 0.5 N`.

A failure of P1 or P2 would not refute (A'), but would weaken the case for an absolute constant.

## Outcome of P1-P6 (scan of every `N < 2^26`, `root_transport_frontier_energy_big_20261017.py 25 ...`)

| octave | `2^23` | `2^24` | `2^25` | test | verdict |
|---|---|---|---|---|---|
| `max E_1000/N` | 1.337 | 1.238 | 1.449 | P1 `<= 2.08` | PASSED |
| `max E_100/N` | 2.082 | 1.586 | **2.621** | P2 `<= 3.0` | PASSED (new record, above 2.35) |
| `max E_10000/N` | 1.400 | 1.415 | 1.380 | P3 `<= 1.55` | PASSED |
| `mean E_1000/N` | 1.164 | 1.038 | 1.281 | P4 `in [1.0, 1.5]` | PASSED |
| `mean E_30/N` | 2.054 | 1.088 | 0.998 | P5 two of three `<= 2.0` | PASSED (two of three) |
| `max R^2/N` | 3.078 | 2.634 | 4.362 | | |

**P6** (decomposition at `Y = 1000`, `root_transport_frontier_decomposition_p6_20261017.out`). At the six arg-max
points `N = 9666088, 32532811, 47930131` (`E_1000`) and `N = 8429844, 32184394, 37675352` (`R^2`):
- `MM_hi/N` lies in `[-3.633, -3.575]`, all `<= -3.0`;
- `|MF|/N <= 0.245`, all `<= 0.5`.

PASSED. The `E_Y` values agree with the incremental scan to four decimals; the two codes are independent.

`E_100/N` behaves like `E_30/N` did at `2^22`: its octave-`2^25` mean is 1.69, against 1.0-1.4 elsewhere. The
coherent arrival climb reaches up to height about 100 in that octave. Both the maxima at `Y = 1000` and the record at
`Y = 100` stay inside the pre-registered bounds.

## Pre-registered P7-P10 (written before the scan of `2^26 <= N < 2^28`)

- **P7.** `max E_1000(N)/N <= 2.08`.
- **P8.** `max E_100(N)/N <= 3.0`.
- **P9.** `max E_10000(N)/N <= 1.55`.
- **P10.** Each octave mean of `E_1000/N` lies in `[1.0, 1.5]`.

## Outcome of P7-P10 (every `N < 2^28`, `root_transport_frontier_energy_big_20261017.py 27 100 1000 10000`)

| octave | `2^26` | `2^27` | test | verdict |
|---|---|---|---|---|
| `max E_1000/N` | 1.299 | 1.245 | P7 `<= 2.08` | PASSED |
| `max E_100/N` | 1.385 | 1.794 | P8 `<= 3.0` | PASSED |
| `max E_10000/N` | 1.263 | 1.268 | P9 `<= 1.55` | PASSED |
| `mean E_1000/N` | 1.214 | 1.099 | P10 `in [1.0, 1.5]` | PASSED |
| `max R^2/N` | 7.483 | 7.370 | | |

The invariant `R(N) = Σ f(u) c_Y(u)` was checked against `H(4N) - H(N)` at every one of the `2^28 - 1` values of `N`.

In octaves `2^26-2^27` the ROOT current reaches `R^2/N = 7.48`, close to the 9.08 record of `2^12`. Meanwhile the
`Y = 1000` frontier energy never exceeds `1.30 N`. All ten pre-registered tests have passed.

## Pre-registered P11-P12 (closed frontiers; written before the scan)

clock_duality's R6 (root_sparse_boundary_module_20261017i.md) certifies forward-closed frontiers:
- `V_100`: 74 seeds, largest 2051;
- `V_1000`: 627 seeds, largest 63665.

Stopping on first entry into `V` leaves no return kernel, and `R(N)^2 <= |V| E_V(N)`. Predictions for every `N` in
`[2^12, 2^26)`:
- **P11.** `max E_(V_1000)(N)/N <= 2.5`.
- **P12.** `max E_(V_100)(N)/N <= 3.5`.

### Outcome of P11-P12 (every `N < 2^26`, `root_transport_closed_frontier_energy_20261017.py 25`)

| octave | `2^12` | `2^13` | `2^14` | `2^16` | `2^18` | `2^20` | `2^22` | `2^23` | `2^24` | `2^25` | test | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `max E_(V_1000)/N` | **2.468** | 2.109 | 1.790 | 1.513 | 1.393 | 1.534 | 1.376 | 1.328 | 1.228 | 1.421 | P11 `<= 2.5` | PASSED |
| `max E_(V_100)/N` | 2.470 | 1.944 | 1.955 | 1.810 | 1.483 | 1.376 | 1.787 | 1.921 | 1.262 | **2.554** | P12 `<= 3.5` | PASSED |

The invariant `R(N) = Σ_(u ∈ V) f(u) c_V(u)` was checked at every `N`.

The renewal-free closed frontiers behave like the height frontiers:
- `E_(V_1000)` is within about 15% of `E_1000` octave by octave;
- the `2^25` octave that carried the `E_100` record carries the `E_(V_100)` maximum too.

So the frontier cone is equally supported in clock_duality's R6 setting, where stopped and unstopped transport agree
outside `V`. The Cauchy-Schwarz constant there is `sqrt(627 · 2.47) < 40`.

## Pre-registered P13-P14 (statistics of the plateau; written before the measurement)

**Model (the HEURISTIC behind (A')).** In octaves `2^24-2^27` the frontier coefficients are treated as approximately
independent centred Gaussians: `c_Y(u) ≈ sqrt(N s_u) Z_u`. Then
- `E_Y/N ≈ Σ s_u Z_u^2`, with mean `m = Σ s_u`;
- standard deviation `sqrt(2 Σ s_u^2)`;
- effective degrees of freedom `F_eff = m^2 / Σ s_u^2`.

Predictions, with `N` sampled on a regular grid:
- **P13.** The kurtosis of `c_Y(u)/sqrt(N s_u)`, pooled over seeds carrying 90% of the energy, lies in `[2, 5]` (a
  Gaussian has 3), for `Y = 100` and `Y = 1000`.
- **P14.** The model standard deviation `sqrt(2 Σ s_u^2)` is within a factor 2 of the observed standard deviation of
  `E_Y/N` over the samples, for both `Y`.

### Outcome of P13-P14 (`root_transport_frontier_statistics_20261017.py 27`; `N` sampled every `2^12` (`Y = 100`) or `2^14` (`Y = 1000`) in `[2^24, 2^28)`)

| `Y` | samples | `F_Y` | `F_eff` | mean `E_Y/N` | model sd | observed sd | max | pooled kurtosis | mean `\|corr\|` |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 61,440 | 38 | 12.2 | 1.241 | 0.502 | 0.315 | 2.617 | 2.07 | 0.42 |
| 1000 | 15,360 | 375 | 113.9 | 1.150 | 0.152 | 0.088 | 1.445 | 2.25 | 0.38 |

- **P13 PASSED.** The kurtosis is in `[2, 5]`, but at the low end: the coefficients are sub-Gaussian.
- **P14 PASSED.** The model standard deviation is 1.6-1.7 times the observed one, within the factor 2.
- **Reading.**
  - The independent-Gaussian model captures the order of magnitude. `F_eff` is about a third of `F_Y`, because the
    energy concentrates on a third of the seeds: for `Y = 100`, ROOT itself carries 20%, and seed 81 carries 12%.
  - The real plateau is about 40% *more stable* than independence predicts. The normalised coefficients are
    substantially correlated (mean `|corr|` 0.4, up to 0.95), and some have persistent means of order their RMS. The
    correlations reduce rather than inflate the fluctuations of the sum.
  - For the (A') heuristic, the maxima sit 4.4 (`Y = 100`) and 3.4 (`Y = 1000`) observed standard deviations above
    the mean over these samples.
