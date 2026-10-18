---
id: HYP-9325
title: "One-tick and pair-sign laws of first arrivals into the 627-seed frontier. In the exact pair ledger E_V = 24N/7 + (G0 - E0) + 2 Gam1 + HL (THM-4690(4)), (A) the one-tick ratio G(1)/G(0) stays in [-0.30, -0.24], (B) same-clock pairs and lag-1 pairs other than sibling sources have non-positive sums, and (C) lags >= 2 add at most 0.07N on large frontiers. Together these give E_V <= 2.6N on V_1000, the closed-frontier form of HYP-9290"
status: >
  OPEN (CONJECTURE). The consequences listed under "Consequences" are PROVED implications. Evidence is FINITE-EXACT, with
  the invariant R(N) = H(4N) - H(N) checked at every N: every N < 2^28 on V_1000 (627 seeds), every N < 2^26 on V_100
  (74 seeds) and on the height frontier H_1000 (the 375 seeds <= 1000).
  Pre-registered tests, each committed before its scan:
  - L1-L7, pair ledger, octave 2^27: ALL PASSED.
  - T1-T4, one-tick law, octave 2^26: ALL PASSED.
  - Q1-Q5, port orthogonality, octaves 2^23-2^25: Q3 passed; Q1, Q2, Q4 and Q5 FAILED. The cross-port term stays bounded
    but does not decay, so the cancellation is clock-structured rather than port-orthogonal. That reading is refuted and
    recorded below.
  Renumbered from HYP-9292, which is ceded to opus-2026-10-18-S43; the committed pre-registrations still say HYP-9292 / THM-4670.
session: mac-mini-2026-10-18-sixcolour
source: 05-knowledge/results/first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md
scripts:
  - 04-computation/experiments/first_arrival_20261018/fap.c (port, pair, depth, clock, spectrum and ledger quantities at every N)
  - 04-computation/experiments/first_arrival_20261018/fap_lean.c (memory-lean pair ledger; the 2^27 test)
  - 04-computation/experiments/first_arrival_20261018/fap_levels.c (backward ladder)
  - 04-computation/experiments/first_arrival_20261018/ledger_table.py (octave tables)
prereg:
  - 05-knowledge/results/first_arrival_port_orthogonality_prereg_20261018.json (Q1-Q5)
  - 05-knowledge/results/first_arrival_one_tick_prereg_20261018.json (T1-T4)
  - 05-knowledge/results/first_arrival_pair_ledger_prereg_20261018.json (L1-L7)
related: [THM-4690, HYP-9290, THM-4667, HYP-9288]
---

# HYP-9325: one-tick and pair-sign laws of first arrivals into the 627-seed frontier

## Setting

The notation is that of THM-4690. Sources of `R(N) = H(4N) - H(N)` are stopped at their first seed in `V`.
- `h_t(u)` is the signed mass entering `u` with literal shortcut clock `t`.
- `G(D) = Σ_u Σ_t h_t(u) h_(t+D)(u)` is the integer-clock lag Gram.
- `E_0 = Σ_x c_0(x)^2 = 34N/7 + O(log N)` is the source energy.
- `SIB_V = 5N/7 + O(|V| + log N)` is the sum over sibling source pairs `(x, 2x+1)`. Such a pair always arrives one tick
  apart with opposite signs.
- `Γ_1 = G(1) + SIB_V` collects the other lag-1 pairs.
- `HL = 2 Σ_(D>=2) G(D)` collects the longer lags.

THM-4690(4) gives the exact ledger

`E_V = E_0 - 2 SIB_V + (G(0) - E_0) + 2 Γ_1 + HL = 24N/7 + (G(0) - E_0) + 2Γ_1 + HL + O(|V| + log N)`.

## Conjecture

For the certified frontier `V_1000` (627 seeds) and the height frontier `H_1000` (375 seeds):

- **(A) One-tick band.** For every `N >= 2^13`, `-0.30 <= G(1)/G(0) <= -0.24`.
- **(B) Pair-sign laws.** For every `N >= 2^11`:
  - **(B1)** `G(0) <= E_0`: same-clock pairs of distinct sources have a non-positive sum;
  - **(B2)** `Γ_1 <= 0`: lag-1 pairs other than sibling sources have a non-positive sum.

  (B) also holds on `V_100` for `N >= 2^11`.
- **(C) Higher lags.** For every `N >= 2^12`, `HL <= 0.07 N`.

  (C) is false-looking on the small frontier `V_100`, where HL reaches `+1.29N` (octave 2^25). Small frontiers are
  therefore excluded; there the long lags carry the termination-regime climbs of HYP-9290.

## Consequences (PROVED implications)

1. **(A) + (B1) + (C) ⇒ E_V <= 2.60N + O(log N).** Indeed `E_V = (1 + 2 G(1)/G(0))·G(0) + HL`, with
   `1 + 2G(1)/G(0) <= 0.52` and `G(0) <= E_0`. So

   `E_V <= 0.52·(34/7)N + 0.07N + O(log N) <= 2.60N + O(log N)`.
2. **(B) + (C) ⇒ E_V <= (24/7 + 0.07)N + O(log N + |V|) <= 3.51N + O(log N + |V|).** This needs no ratio band.
3. **Bound on R(N).** With clock_duality's R6.4 (`R(N)^2 <= |V| E_V(N)`), statement 1 gives
   `|R(N)| <= sqrt(627 · 2.60 N) + o(sqrt N) ≈ 40.4 sqrt(N)`. This is HYP-9290's energy cone for the closed frontier,
   with `C = 2.60`, reduced to three statements about pair correlations of first arrivals with their integer clocks.

## Evidence (FINITE-EXACT)

### V_1000, every N < 2^28

Each row is one octave `2^k <= N < 2^(k+1)` (`fap.c` to 2^26, `fap_lean.c` for 2^27; the two programs agree exactly on
the octaves both cover). All values are divided by `N`, except the ratio `G1/G0`. `E0/N` lies within `10^-5` of 34/7,
and `SIB/N` within `10^-5` of 5/7, from `k = 20` on.

| k | E_V mean / max | G0 mean | G1/G0 mean [lo, hi] | (G0-E0) mean / max | Γ1 mean / max | HL mean / max | E_V/E+ |
|---:|---|---:|---|---|---|---|---:|
| 11 | 2.769 / 3.299 | 4.528 | -0.181 [-0.206, -0.138] | -0.331 / -0.218 | -0.308 / -0.221 | -0.121 / +0.068 | 0.502 |
| 12 | 2.182 / 2.468 | 4.414 | -0.231 [-0.244, -0.204] | -0.443 / -0.402 | -0.401 / -0.325 | -0.189 / -0.054 | 0.395 |
| 13 | 1.863 / 2.109 | 4.311 | -0.260 [-0.269, -0.242] | -0.546 / -0.449 | -0.443 / -0.409 | -0.208 / -0.007 | 0.351 |
| 14 | 1.535 / 1.790 | 4.184 | -0.265 [-0.272, -0.259] | -0.673 / -0.567 | -0.412 / -0.368 | -0.433 / -0.154 | 0.284 |
| 15 | 1.309 / 1.457 | 3.968 | -0.264 [-0.269, -0.256] | -0.889 / -0.724 | -0.346 / -0.321 | -0.564 / -0.450 | 0.264 |
| 16 | 1.342 / 1.513 | 3.841 | -0.276 [-0.284, -0.266] | -1.016 / -0.975 | -0.353 / -0.319 | -0.380 / -0.164 | 0.298 |
| 17 | 1.360 / 1.517 | 3.745 | -0.279 [-0.286, -0.274] | -1.112 / -1.040 | -0.334 / -0.304 | -0.295 / -0.122 | 0.296 |
| 18 | 1.258 / 1.393 | 3.652 | -0.273 [-0.278, -0.267] | -1.205 / -1.156 | -0.285 / -0.259 | -0.398 / -0.250 | 0.283 |
| 19 | 1.282 / 1.465 | 3.579 | -0.276 [-0.280, -0.271] | -1.278 / -1.230 | -0.274 / -0.253 | -0.321 / -0.112 | 0.264 |
| 20 | 1.400 / 1.534 | 3.506 | -0.277 [-0.282, -0.270] | -1.351 / -1.307 | -0.255 / -0.227 | -0.167 / -0.048 | 0.310 |
| 21 | 1.368 / 1.557 | 3.401 | -0.275 [-0.280, -0.270] | -1.456 / -1.392 | -0.221 / -0.198 | -0.163 / +0.016 | 0.324 |
| 22 | 1.229 / 1.375 | 3.340 | -0.273 [-0.279, -0.266] | -1.517 / -1.453 | -0.198 / -0.178 | -0.286 / -0.171 | 0.276 |
| 23 | 1.158 / 1.328 | 3.289 | -0.281 [-0.284, -0.276] | -1.568 / -1.512 | -0.208 / -0.192 | -0.286 / -0.150 | 0.281 |
| 24 | 1.049 / 1.228 | 3.201 | -0.284 [-0.288, -0.277] | -1.656 / -1.615 | -0.194 / -0.178 | -0.337 / -0.165 | 0.228 |
| 25 | 1.278 / 1.421 | 3.152 | -0.277 [-0.285, -0.272] | -1.705 / -1.637 | -0.159 / -0.134 | -0.128 / +0.027 | 0.243 |
| 26 | 1.163 / 1.277 | 3.086 | -0.280 [-0.285, -0.274] | -1.771 / -1.725 | -0.151 / -0.140 | -0.193 / -0.100 | 0.240 |
| 27 | 1.095 / 1.183 | 3.079 | -0.282 [-0.287, -0.277] | -1.778 / -1.758 | -0.155 / -0.137 | -0.245 / -0.134 | 0.248 |

**The cancellation budget at 2^27, per `N`.**

| Term | Value | Status |
|---|---:|---|
| Sources | +34/7 = 4.857 | exact |
| Sibling-source pairs | -10/7 = -1.429 | exact |
| Same-clock pairs | -1.778 | |
| Other lag-1 pairs (`2Γ_1`) | -0.310 | |
| Longer lags | -0.245 | |
| **Total, `E_V`** | **1.095** | |

Unsigned arrivals keep `E+/N ≈ 4.1-5.4` (octave means), about `E_0`. The signs `(-1)^clock` remove about three quarters of it:
`E_V/E+ ≈ 0.23-0.33` from 2^14 on. This is the Nyquist dip of the clock spectrum (THM-4690(3)).

The clock spectrum `S(kπ/12)/N`, `k = 0..12`, is low-pass at all six dumps (`N = 7910` to `10^8`): it stays near 4-6
for `θ <= π/3` and falls to about 1.2 at `θ = π`. At `N = 10^8` the values are
`4.85, 4.01, 4.11, 4.60, 4.09, 3.81, 3.63, 2.80, 2.04, 1.91, 1.52, 1.33, 1.18`. A pure one-tick model
`S(θ) = G(0) - 2G(1) cos θ = G(0)(1 + 0.56 cos θ)` predicts `4.78, ..., 3.06, ..., 1.35`.

### Other frontiers (every N < 2^26)

- **H_1000** (the 375 seeds <= 1000; not forward closed):
  - G1/G0 lies in `[-0.293, -0.259]` from 2^12 on;
  - (G0-E0) has octave maxima `< 0` from 2^11 on, mean -1.72 at 2^25;
  - Γ1 has octave maxima `< 0` from 2^9 on;
  - HL has maximum +0.061 (2^25) and is otherwise `<= +0.020` from 2^12 on.
- **V_100** (74 seeds):
  - G1/G0 lies in `[-0.329, -0.246]` from 2^11 on;
  - (B) holds throughout;
  - HL reaches +0.742 (2^12), +0.603 (2^23) and +1.290 (2^25), with E_V max 2.554 at 2^25.

  Small frontiers are excluded from (C).

### The port-orthogonality reading is refuted (Q1-Q5)

Before the scan of 2^23-2^25 it was predicted that the first-arrival cross term `X = E_V - E^(1)` decays on V_1000.

| Item | Prediction | Observed | Result |
|---|---|---|---|
| Q1 | V_1000 octave maxima of \|X\|/N <= 0.10 | 0.124, 0.155, 0.060 | FAILED |
| Q2 | mean of those maxima < 0.08 | 0.113 | FAILED |
| Q3 | V_1000 octave maxima of E_V/N <= 2.5 | 1.33, 1.23, 1.42 | PASSED |
| Q4 | V_100 octave maxima of \|X\|/N <= 0.35 | 0.232, 0.364, 0.435 | FAILED |
| Q5 | V_1000 octave means of X/N in [-0.06, 0.06] | -0.076, -0.085, +0.010 | FAILED |

`X/N` stays within `[-0.66, +0.11]` from 2^11 on, but it does not tend to zero. The data pointed instead to the integer
clocks, which led to T1-T4 and L1-L7.

### Pre-registered outcomes of the clock tests

- **T1-T4, octave 2^26, V_1000.**
  - T1: mean G1/G0 = -0.2803, in [-0.290, -0.265].
  - T2: every N in [-0.2845, -0.2744], inside [-0.300, -0.250].
  - T3: G0/N mean 3.086, in [2.85, 3.20].
  - T4: E_V max 1.277 <= 2.5, and cross-port mean -0.024 in [-0.15, 0.10].

  All passed.
- **L1-L7, octave 2^27, V_1000.**
  - L1: mean G1/G0 = -0.2823, per N in [-0.2873, -0.2773].
  - L2: (G0-E0)/N max -1.758 < 0, mean -1.778.
  - L3: Γ1/N max -0.137 < 0, mean -0.155.
  - L4: HL/N max -0.134 <= 0.07.
  - L5: E_V/N max 1.183, mean 1.095.
  - L6: G0/N mean 3.079.
  - L7: E_V/E+ mean 0.248.

  All passed. The theorem check also holds: `E0/N = 4.85714`, `SIB/N = 0.71428` at every N.

## Remarks

- **Weakest link.** `Γ_1/N` is shrinking: -0.44 at 2^13, -0.15 at 2^27. So (B2) may fail at larger N. Consequence 1
  does not use (B2).

  The meet-depth anatomy (`pair_series_meets.py`, `N = 2^18`, V_1000) shows the mechanism. Equal-depth meets contribute
  `-0.71, -0.31, -0.12, -0.05` at lag 1 for depths 1-4, but `+0.005..+0.018` per level at depths 6-14. As `N` grows and
  deeper levels fill in, `Γ_1` drifts up. (B2) is therefore expected to FAIL eventually.
- **Where (B1) comes from.** Same-clock pairs are mostly equal-depth "clock twins". They have the same number of
  3-divisions and the same total power of 2, hence the same scale. They contribute about -0.04..-0.06 per level at
  depths 2-8, decaying beyond. Unequal-depth meets are nearly uncorrelated (each cell at most 0.013 in absolute value).
- **Is G(0) stabilising?** The decline of `G(0)/N` slowed to -0.007 between 2^26 and 2^27, and the same-clock term moved
  only from -1.771 to -1.778.
- **Why the one-tick ratio sits near -0.28.** About 82% of `G(1)` at 2^27 is the exact sibling constant: `-5/7` of
  `G(1)/N ≈ -0.869`. The ratio is stable because `G(0)` and `Γ_1` drift in compensating directions. This is HEURISTIC;
  if `Γ_1 → 0` while `G(0)/N` keeps falling, the ratio would slowly grow in magnitude.
- **Proof routes (HEURISTIC).** (B1) asks that sources sharing an entry seed and a literal clock be anticorrelated in
  total. Same-clock pairs need paths of different depth with equal clock sums. Their `c_0` values are governed by the
  residues mod 3 along B-chains (THM-4690(4)), so a residue-class computation in the style of the `34/7` derivation is
  the natural first attempt.
