---
source: opus-2026-07-09-S168
status: the owner's BROKEN-CLOCK insight is the conceptual spine of the good-period dichotomy. A clock
  at rate r is "exactly right" at rate |r-1| (crossings of the exact value): the stopped clock (r=0) is
  right often; the near-correct clock (r~1) almost never -- "exactly right" is measure ZERO, being CLOSE
  means crossing it RARELY. Runners ARE clocks; relative speed d_i = v_i - w. This gives the TWO good-
  period branches as the TWO clock regimes: (near-AP = SLOW-relative = LINGERING = few long bad blocks =
  the positive-measure obstruction, closed by klein LEM-012 Dirichlet-clustering) vs (dissociated = FAST-
  relative = many SHORT bad blocks, #arcs sublinear in Vmax => closed by mac-mini arc-count). Verified:
  dissociated (longest-AP<=k-6, which REQUIRES spread>>k) has c=#arcs/spread ~0.04 median (max->0.05 at
  spread 1920) << rho*~0.96, #arcs ~ spread^0.92 -- mac-mini's route CLOSES it, confirming opus-S167.
  Records a regime error I caught before filing (small-spread k=11 is forced near-AP, not dissociated).
tags:
  - lrc14
  - good-period
  - broken-clock
  - near-resonance
  - measure-zero
  - arc-count
  - dichotomy
---

# The broken clock: the two good-period branches ARE the two clock regimes

**opus-2026-07-09-S168.** The owner's observation -- "a broken clock is right twice a day, one ticking
correctly is right 0 times" -- is not a joke; it is the conceptual spine of the good-period DICHOTOMY
(near-AP vs dissociated) that the fleet's two proof branches split along. Developed, connected, and
computationally confirmed below.

## The clock lemma (the owner's insight, made precise)

Model a clock as a point on the circle `R/1` displaying `b + r t` at true time `t` (rate `r`, offset
`b`).  It is "exactly right" when `b + r t ≡ t`, i.e. `(r-1) t + b ≡ 0 (mod 1)` -- crossings of a fixed
value by a point moving at relative speed `r-1`, happening at **rate `|r-1|`**:

- stopped clock (`r=0`): rate `|0-1| = 1` -- "right twice a day" (the reference/maximal rate);
- near-correct clock (`r ≈ 1`): rate `|r-1| ≈ 0` -- **almost never exactly right**;
- fast/slow clock (`|r-1|` large): right often.

The paradox dissolves because "exactly right" is a MEASURE-ZERO event (a single real value); a clock
set close to true time and ticking at nearly the true rate CROSSES that value rarely -- proportional to
its relative speed -- and so is exactly right *less* often than a broken one.  Being a good
approximation (`r ~ 1`, self-consistent) GUARANTEES almost never being exactly correct.

## The runners ARE clocks -- and the two branches are the two clock regimes

Every runner is a clock: runner `i` at speed `v_i`, observer at speed `w`, relative speed
`d_i = v_i - w`.  Runner `i` crosses the observer's position at rate `|d_i|` and sits in the danger arc
(width `1/7`) for duration `~ (1/7)/|d_i|` per crossing, so total blocked time `= 1/7` for EVERY runner
-- but distributed oppositely.  This is exactly the good-period proof dichotomy:

| clock regime | runners | bad set | good-period branch |
|---|---|---|---|
| **near-correct** (`d_i` small, `r~1`) | **near-AP** cluster | FEW LONG blocks -- runners LINGER | klein **LEM-012** (Dirichlet clustering) |
| **fast** (`d_i` large) | **dissociated** cluster | MANY SHORT blocks | mac-mini **arc-count** (`#arcs < rho* V`) |

So the LRC obstruction -- the persistent, positive-measure blocking -- is the LINGERING of the
near-resonances (small relative speed), which is *why* the near-AP branch is the delicate one (a
Dirichlet clustering / pigeonhole argument, LEM-012).  The measure-zero EXACT alignment
(`d_i t ≡ 0`) is negligible -- **exactly mac-mini-S61's `Corr_N` decomposition** (dominated by the
near-resonance; the exact-resonance defect `E_grid[W] - (6/7)^k ~ 0.001` is negligible).  Chasing the
exact/`r_N` resonance (opus-S167, deprioritized) is chasing when the clock is *exactly* right -- almost
never, and not the point.  The clock lemma is the conceptual root of BOTH mac-mini's near>>exact
decomposition AND the near-AP / dissociated split.

## Confirmation: dissociated closes by arc-count (and a regime error I caught)

I first "rechecked" mac-mini's arc-count with the naive pigeonhole `#g < mu*V` at spread `12-35` and
got `R = #g/(mu*spread) ~ 1-2.3` -- an apparent FAILURE.  That was a **regime error**: at spread
`12-35`, eleven points cannot avoid a long AP -- they are forced NEAR-AP (the lingering regime, high
`#arcs`), NOT dissociated.  A genuinely dissociated cluster (longest-AP `<= k-6`) REQUIRES spread `>> k`
to spread out -- the fast-clock regime.  Testing THAT regime (the clock insight tells you which one):

> k=11 dissociated, spread `30 -> 1920`: median `c = #arcs/spread` stays `~0.04` (max shrinks `0.40 ->
> 0.05`), `rho* = mu ~ 0.96 -> 0.998`, and `#arcs ~ spread^{0.92}` (SUBLINEAR).  So `c = #arcs/spread ->
> 0`, and `#arcs < rho* * V` holds with a HUGE margin (`0.04 << 0.96`).

So mac-mini's arc-count route DOES close the dissociated branch -- confirming opus-S167 (the a-priori
`rho* >= D3_inf^{(L)} >= 0.46 >> c`; the true `rho* ~ 0.96` is even more generous than the S158 bound).
The many short good intervals of a fast-clock cluster are sublinear in `V`, so plenty carry a good grid
point -- the crude worst case (all good intervals `< 1/V`) never happens for genuinely spread clusters.
The clock insight is what told me my first recheck was in the wrong regime; the small-spread `c ~ 1`
cases are the near-AP FINITE regime, closed by LEM-012 / direct enumeration, not a counterexample.

## Ledger

- The BROKEN-CLOCK insight (owner): "exactly right" is measure zero; near-correct clocks cross it
  rarely.  Runners are clocks (`d_i = v_i - w`); the two good-period branches ARE the two clock regimes
  -- near-AP = slow-relative = LINGERING = positive-measure obstruction (LEM-012); dissociated =
  fast-relative = many short blocks, `#arcs` sublinear (arc-count).  The conceptual root of mac-mini's
  near>>exact `Corr_N` decomposition and the deprioritize-`r_N` (opus-S167) conclusion.
- CONFIRMED: dissociated (longest-AP `<= k-6`, needs spread `>> k`) closes by arc-count -- `c ~ 0.04
  << rho* ~ 0.96`, `#arcs ~ spread^{0.92}` sublinear (spread `30-1920`).  opus-S167 synthesis stands.
- REGIME NOTE (caught before filing): the naive `#g < mu*V` at spread `12-35` "fails" only because
  small-spread k=11 is forced near-AP, not dissociated -- the clock insight identifies the regime.
- CONVERGENCE (mac-mini-2026-07-09-S62, independent): same two findings. mac-mini closed the
  dissociated branch as the exact dilation-aware inequality `c = #arcs/spread < D3(E)` (`c/D3` monotone
  DECREASING in spread, `<1` throughout: `0.90` at spread 80 -> `0.44` at spread 1000) + a `Vmax <= 234`
  finite check -- and wrote the SAME broken-clock reflection ("cluster = a bank of clocks ticking near
  true rate `Vmax`; near-true clocks coincide rarely").  My angle adds: (a) the mechanism `#arcs ~
  spread^{0.92}` (sublinear) that MAKES `c/D3 -> 0`; (b) the two-branches = two-clock-regimes table; (c)
  the crossing-rate `|r-1|` formalization; (d) the regime-error caution.  mac-mini's `c < D3` framing is
  cleaner than my `c < rho*` (`D3` exact + dilation-invariant, no `rho*` estimate needed) -- adopt it.
- Files: `lrc14_arccount_recheck_opus_S168`, `lrc14_arccount_growth_recheck_opus_S168` (+outs),
  `lrc14_arccount_crossover_opus_S168`.  -> mac-mini-S61/S62 (arc-count `c < D3` + clock reflection),
  klein LEM-012, opus-S158 (`D3_inf^{(L)}`), opus-S167 (interlock), the near-resonance = lingering.
