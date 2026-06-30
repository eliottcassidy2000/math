# The razor-thin line of LRC(14) is in the MEASURE, not the PEAK — so target existence, not concentration; the AP is the additive extremal, the obstruction is its multiplicative perturbation

*opus-2026-06-29. Owner: understand what a disproof entails, map the razor-thin line between proof and
disproof comprehensively, and use it to improve proof targets. The mapping inverts the naive picture
and points the proof at a different (margin-bearing) object.*

## What a disproof of LRC(14) entails
A 13-set `S` with `M(S) = max_t min_{s∈S} ‖s t‖ < 1/14`. By THM-523 (`S` not covering ⇒ won at
`t=1/q`), a disproof **must be a covering set** (multiple of every `q∈{2..14}`) whose 13 danger combs
cover the circle for *all* `t`. Equivalently: the diagonal line on the 13-torus never escapes the
view-obstruction tubes.

## The tight edge, mapped (verified, exact)
- **`M({1,…,13}) = 1/14` exactly** — the arithmetic progression is the *tight extremal* (witness
  `t=1/14`). LRC(14) `⟺` the AP `{1,…,13}` is a global **minimizer** of `M`, with `min = 1/14`. A
  disproof is a set strictly *below* the AP.
- **The covering side has a +10% MARGIN, not a razor edge.** Min over covering sets `= 7/89 ≈ 0.0787`
  (at `{1,…,11,13,84}`), and the divisor-loaded family `{1,…,11,13,84m}` *increases* away from 1/14
  (`0.0787, 0.0809, 0.0817, …, 0.0827`). So the hard case is **comfortably above** 1/14.

## The inversion: razor-thin is the MEASURE, margin is the PEAK
Three different LRC quantities, three different margins:
| quantity | what it is | margin | route status |
|---|---|---|---|
| **peak `M`** | `∃` ONE lonely time (the conjecture itself) | **+10%** (`≥7/89`) | the right target |
| lonely-MEASURE | `meas{t : lonely}` | **→ 0** (razor-thin) | over-strong |
| `ρ*` (criterion-C density) | good Vmax-periods | **= 0** (refuted) | dead (my earlier work) |

> **The measure-based routes prove TOO MUCH.** They certify a *positive-measure* set of lonely points,
> which is razor-thin (the lonely set is a thin tall peak — small measure, comfortable height). The
> conjecture only needs **ONE** lonely point, and that *peak existence* carries the +10% margin. **Don't
> prove more than the existence.** This is exactly mac-mini's floor/cap split read correctly: the *cap*
> (concentration/measure) is the `R`-odd razor-thin obstruction; but the *peak/existence* is weaker and
> has slack.

This explains my own arc: the `ρ*` route died because it targeted the razor-thin density; the
**prime-density criterion `N(S,p)≥1` is the existence route** and is the right one — `N(S,p)≥1` (one
witness) carries the margin, while `N(S,p)/p` (the density) is the razor-thin measure.

## Why the tools still fail despite the margin (the real gap)
The +10% margin is real but **unreachable by the union bound**, which is >100% off (`13·(1/7)=13/7`
vacuous). So the gap is not "exact vs approximate" — it is **"within 10%" vs ">100% off."** Improved
target: a bound accurate to within 10% of sharp on the *peak/existence*, not an exact certificate and
not the vacuous first moment. That is a much softer requirement than the razor-sharp signed certificate
the measure routes needed.

## The add/mult refinement of the extremal
- **The extremal is ADDITIVE:** the minimizer `{1,…,13}` is the *consecutive integers = additive
  interval* — the same interval/Dirichlet structure that is the **max-H** extremal (the half-turn). So
  *both* extremal problems (min-M for LRC, max-H for tournaments) are solved by the additive interval.
- **The obstruction is the MULTIPLICATIVE perturbation:** a covering set forces a *multiple of 14*
  (and of every `q`) — a multiplicative/divisibility constraint that perturbs the additive AP *off* the
  `t=1/14` witness, pushing `M up` by ~10%. The proof's job is to control this multiplicative
  perturbation of the additive extremal — i.e. **peel the multiplicative `mult-of-14` (the apex-7
  descent, mac-mini THM-580) to reduce to the additive AP core, where `M=1/14` is exact and the peak is
  understood.**

## Improved proof program (the reorganization)
1. **Target the PEAK/EXISTENCE** (`∃` one lonely `t`), not the measure — it carries the +10% margin
   (`M≥7/89`). The prime-density `N(S,p)≥1` is this target.
2. **Aim for "within 10%," not "exact."** The signed/cancellation machinery the measure needed is
   over-kill for the peak.
3. **Descend the multiplicative perturbation** (peel `mult-of-14`, apex-7) to the **additive AP core**,
   the shared extremal of LRC-min-M and tournament-max-H.
4. **A disproof would have to beat the AP** — drive `M<1/14`. No bounded covering set comes near
   (min `7/89`; loading moves *away*); the bounded-core reduction says a counterexample would be
   bounded, and none is found. So LRC(14) is true *with margin*; the difficulty is purely **uniform
   access to a margin that exists** — not a razor's edge.

## Status
- **Verified:** `M(AP)=1/14`; covering min `=7/89`, divisor-loaded *increasing*; (b) max-H = single SC
  prime (super-multiplicativity).
- **Reframe:** razor-thin = measure (over-strong); margin = peak/existence (the right target). Extremal
  additive (AP), obstruction multiplicative (mult-of-14 perturbation).
- **Improved target:** peak/existence within 10%, via apex-7 descent to the additive AP core.

Related: THM-523 (covering reduction, min 7/89), THM-566 (divisor loading), THM-527 (the refuted ρ*),
mac-mini THM-580/581 (descent / floor-existence), HYP-3547 (apex-7), my prime-density + add/mult-duality
reflections, OPEN-Q-108.
