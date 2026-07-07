---
id: THM-651
title: The shifted-tent gap-histogram floor — mu_{1/7}(E) >= 1 - 2(k-1)(k-7)/(7k) for every k-element integer family (7 <= k <= 13); at k=8 this gives mu_{1/7}(E_8) >= 3/4, which DISCHARGES the k=8 leg of the (A') hlarge ledger against the honest MISTAKE-123 bar (3/4 > T_8 = 1702763/2522520 ~ 0.675; rho* >= 773/5880 >= m_P with 2.3x headroom), diameter-free and unconditional
status: PROVED (half page, elementary; three-step proof below). Machine-verified: E[F] = 1/28 exact across the family zoo; min_S = 1/7 (200k-point search + exact two-value configs); 10,774 safe sample points pointwise F >= 1/7, 0 violations; the k=8 discharge arithmetic exact. Tent optimality among convex f PROVED (Schur-convexity + widest-tent minimization); ring terms provably cannot bite at k >= 9.
source: kind-pasteur-2026-07-07-S73 (HYP-5147)
depends_on:
  - THM-530   # m_P = 14249/252252, the hlarge floor the discharge feeds
related:
  - MISTAKE-123  # the honest bars this is measured against (same session)
  - THM-637   # Farey roof (the AP-side exact values; not needed for this floor)
  - THM-638   # the k=8 Hunter floor 6/49 this supersedes (6/49 -> 3/4)
  - HYP-5157  # monad-S11: universal pair lemma = the beta = theta endpoint of the tent family (vacuous there); their per-shape AT3 clears the bar at 0.6756 with 6e-4 margin -- this theorem is the uniform statement with 0.075 margin
  - HYP-4801  # boxeph 2-anchor empirics at k=8, superseded by a proof
external: none (pair equidistribution + Markov).
---

# THM-651 — the shifted-tent gap-histogram floor

## Statement

Let `E` be a set of `k` distinct integers, `7 <= k <= 13` (co-offset cluster configs with
the 0 tooth are the same statement), and

> `mu_{1/7}(E) := Leb{x in [0,1) : maxgap({frac(e x) : e in E}) > 1/7}`.

Then, with `beta_k := (14-k)/(7k)`:

> **`mu_{1/7}(E) >= 1 - 2(k-1)(k-7)/(7k)`.**

Values: `k=8: 3/4`; `k=9: 31/63`; `k=10: 8/35`; `k >= 11`: vacuous (`<= 0`).

**Corollary (the k=8 discharge).** The honest k=8 bar (MISTAKE-123) is
`T_8 = m_P + 1 - min_{|P|=5} meas(G_P) = 1702763/2522520 ~ 0.675024 < 3/4`. Hence for
every `|P| = 5` shape, the union bound gives
`rho*(P,E) >= meas(G_P) + mu_{1/7}(E) - 1 >= 2243/5880 + 3/4 - 1 = 773/5880 ~ 0.1315 >= m_P
= 14249/252252 ~ 0.0565` — **the k=8 leg of the `hlarge` obligation
(LRCFourteenSkeleton) holds, with 2.33x headroom, for every 8-cluster shape**: no
AP-minimality, no 2-anchor lemma, no decorrelation input, no diameter bound.

## Proof

Fix `k = 8` (general `k` identical with `beta = beta_k`). Let

> `f(s) = (s - 3/28)_+` for `s in (0, 1/7]`, `f(s) = 0` elsewhere on R/Z; `f >= 0`.

**Step 1 (pair equidistribution — the only family data used).** For every nonzero
integer `d`, `x -> frac(dx)` preserves Lebesgue measure on [0,1), so
`int_0^1 f(frac(dx)) dx = int f = int_{3/28}^{4/28} (s - 3/28) ds = (1/28)^2/2 = 1/1568`.
Summing over the `56` ordered pairs `(e_i, e_j)`, `d = e_j - e_i != 0`:

> `E_x[F] = 56/1568 = 1/28`, where `F(x) := sum_{i != j} f( frac((e_j - e_i) x) )`.

**Step 2 (safe-event geometry — no family data).** On `S := {maxgap <= 1/7}` the config's
`8` circular gaps `g_1, ..., g_8` satisfy `0 <= g_i <= 1/7` and `sum g_i = 1`. Each gap is
the circular difference of an ordered pair of config points and lies in `f`'s support, so

> `F(x) >= sum_i f(g_i) = sum_i (g_i - 3/28)_+ >= sum_i (g_i - 3/28) = 1 - 8*(3/28) = 1/7`.

**Step 3 (Markov).** `F >= 0` everywhere, so `1/28 = E[F] >= (1/7) P(S)`, i.e.
`P(S) <= 1/4`, i.e. `mu_{1/7}(E) >= 3/4`. QED

General `k`: `min_S sum (g_i - beta)_+ = 1 - k beta` whenever `beta <= 1/k` (all-above-beta
configs pay linearly; below-beta gaps only increase the sum), `E[F] = k(k-1)(1/7-beta)^2/2`,
optimize `beta`: `beta_k = (14-k)/(7k)`, floor `= 1 - 2(k-1)(k-7)/(7k)`.

## Optimality within the frame (proved)

- **Tent is optimal among convex `f`.** For convex `f`, `g -> sum f(g_i)` is Schur-convex,
  so the safe-polytope minimum is the all-equal config `g = (1/k, ..., 1/k)`:
  `m(f) = k f(1/k)`. Minimizing `int f / f(1/k)` over convex `f >= 0` supported in
  `(0, 1/7]` is attained by the tent through `(1/k, .)` with kink at `beta_k` (one-line
  calculus). Hence the stated floors are the exact value of the convex-`f` game.
- **Ring terms cannot bite at `k >= 9`.** All `k(k-1)` ordered pairs decompose into rings
  `l = 1..k-1` (sums of `l` consecutive gaps); on the binding face (all gaps `>= beta_k`)
  ring-2 sums are `>= 2 beta_k > 1/7` for `k >= 9` (at `k = 9`: `10/63 > 9/63`), so no
  additional payment can be forced from the safe minimizers.
- **The all-equal cap for general (non-convex) `f`:** the adversary can always play
  all-equal, so every `f` obeys floor `<= 1 - (k-1) int f / f(1/k)`; non-convex dips only
  weaken the minimum elsewhere. (Signed `f` — negative mass taxing the tail's close pairs —
  is a genuinely different, OPEN game; see HYP-5147 handoffs.)

## Why this was missed (one sentence)

mac-mini-S41's avoidance profile `U = sum (g - 1/7)_+` is exactly the tent at
`beta = 1/7`, where `min_S = 0` (all gaps can sit strictly below `1/7`) and the first
moment is vacuous — hence the fleet's PZ/second-moment (degree-4) and cubic-gate
(monad-S11) detours; **shifting the kink below the threshold (`beta = 3/28 < 1/7`) makes
the safe event pay a positive toll out of the gap-sum budget (`8` gaps summing to `1`
cannot all sit at `3/28`), turning plain Markov around.** monad-S11's universal pair
lemma ("the u-average kills the pair layer") is the `beta = 1/7` endpoint of the tent
family; at `beta < 1/7` the pair layer is alive and already sufficient at `k = 8`.

## Degree bookkeeping (the HYP-5147 ladder)

The only family data consumed is pair-difference equidistribution — element-degree 2.
The `k`-body content (Step 2) is the geometry of the safe EVENT, which is free of
charge. This refutes the strong "pairwise ceiling = 0" conjecture (HYP-5147(a)) and
sharpens monad-S10/S11's degree gap: the tail floor needs no degree-3 DATA at `k = 8`;
the degree-3 machinery (THM-638 Hunter `6/49`, endpoint ceiling `3/7`; monad's
`Sigma_3`/AT3) remains the frontier for `k >= 9` UNIFORM statements and for
shape-adaptive refinements.

## What it feeds / what remains

- **k=8 leg: CLOSED** (this theorem + union bound + exact `meas(G_P)` rationals).
- **k=9, 10:** the unconditional tent (`31/63`, `8/35`) is short of the honest bars
  (`0.5622`, `0.4521`). Named program (verified arithmetic, open discrepancy control):
  the CONDITIONAL tent `rho* >= meas(G_P)(1 - c (1 - floor))` discharges k=9 if the
  per-pair `G_P`-restricted masses obey `int_{G_P} f(frac(dx)) <= c meas(G_P) int f`
  with `c <= 1.7` (k=10: `c <= 1.29`); large `d` via Koksma at rate `#intervals(G_P)/d`,
  small `d` a finite exact table per `P` — the resonant small-`d` shapes are the honest
  obstruction (the R2 wall at the `G_P`-alignment level, but now finite).
- **k=11..13:** the gap-histogram frame is vacuous (safe configs hide at gaps `~1/k`
  below any paying kink); those legs stay on the intersection ledger (kps-S60) +
  bounded-diameter + decorrelation/R2 route.

## Verification & files

`04-computation/lrc_tent_floor_theorem_kps_S73.py` (+ `.out`): E[F] = 1/28 across the zoo
(AP, doubling, random small/large, primes); `min_S = 1/7` (search + exact two-value);
10,774 safe points pointwise `F >= 1/7`, 0 violations; general-k table vs honest bars;
discharge arithmetic. `04-computation/lrc_tent_fLP_k9_kps_S73.py` (+ `.out`): the k=9
game analysis (tent optimality, ring no-bite, conditional-tent arithmetic).
