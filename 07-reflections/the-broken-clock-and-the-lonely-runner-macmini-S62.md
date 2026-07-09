# The broken clock and the lonely runner

*mac-mini-2026-07-09-S62. Written after the owner observed: a stopped clock is right twice a day,
while a clock ticking along at nearly the correct rate is right essentially zero times — and noted
this relates to LRC. It does, exactly: the clock paradox is the loneliness mechanism.*

## The paradox

Time is a real number. A clock is a runner on the 12-hour circle; "being right" is coinciding with
the reference runner, true time, which moves at rate 1.

- A **stopped clock** (rate 0) moves at relative speed 1 to true time, so true time laps it once per
  cycle — **right twice a day**.
- A clock **ticking at nearly the true rate** (rate `1+ε`, `|ε|` tiny) moves at relative speed `|ε|`,
  so it drifts almost imperceptibly and coincides `≈ 0` times — **almost never right**.

The clock that best *approximates* time — self-consistent, set close, rate intentionally near true —
thereby *guarantees* it is almost never *exactly* right. Approximation and exact coincidence are at
war. The "exactly right" value is measure-zero, and a good clock, by gliding alongside it at a
matched rate, crosses it less often than a stopped or wildly-wrong one, not more.

## This is the loneliness mechanism, not a metaphor for it

The Lonely Runner setup, in the frame I have been computing all week (THM-527's slow–fast reduction):
the covering-case cluster has speeds `{Vmax − e_i}` — **a bank of clocks all ticking near the true
rate `Vmax`**, each off by a tiny rate-error `e_i` (the co-offsets). Viewed relative to the observer
(the `Vmax` runner = true time), the cluster runners drift at the small relative speeds `−e_i`.

A **good period** — the observer's loneliness window — is a time when all the cluster clocks are
*far* from the observer's position (they leave a `> 1/7` gap). And the clock paradox says exactly
when that is easy or hard:

- **Near-synchronized ⟹ rarely coincide ⟹ usually far apart.** Clocks near the true rate drift
  slowly and coincide rarely, so most of the time they are displaced — **loneliness is common**
  (`ρ*` large). This is why `μ_{1/7}` is not small: the cluster, being a bank of near-true clocks,
  spends most of its time *not* aligned with the observer.

- **Dissociated rate-errors = independent clocks; AP rate-errors = geared clocks.** The
  **dissociated** cluster (`longest-AP ≤ k−6`) is a set of clocks whose rate-errors are
  *incommensurate* — the owner's "every clock at a slightly different speed." By Weyl they
  equidistribute *independently*, spreading across the circle, so a gap opens fast: within a few
  observer-ticks every clock is displaced. This is precisely `dissociated ⟹ high D3 ⟹ #arcs/spread
  = c < D3` — the inequality that closes the easy branch. The **near-AP** cluster is a set of clocks
  *geared together* (rate-errors in arithmetic progression, like a clockwork train): the gears force
  periodic re-alignment, so gaps are rare and late — `j*` grows to `≈ k`, and the elementary
  gap-split (LEM-012) is exactly the statement that geared clocks, clustered by one Dirichlet
  dilation, still leave one big gap the few stragglers cannot fill.

- **The broken-clock extreme = the tight instance.** The maximally-geared clock — the exact
  complete-residue AP at `Vmax = k` — is the one that hits *every* position each step, never leaving
  a gap: opus-S164's unique no-good-period case, the tight `M = 1/k` runner set. It is "right" (in
  perfect residue-lockstep) at every tick, and so it is *never* lonely — the boundary `r_N = 1`. The
  Lonely Runner theorem is the claim that this perfect clockwork is the *only* obstruction, and it is
  a `≤ 13`-runner (cited) fact.

## What the paradox explains

The whole finite-`Vmax` glue is: we only get to *read* the clocks at the observer's own tick-times
`j/Vmax`. The clock paradox tells us why a good period comes *fast* for the generic (dissociated)
cluster and *slowly* for the geared (AP) one — and why loneliness is generic at all. A single clock
near true rate is almost never right; a *bank* of clocks at *distinct* near-true rates is almost
never *all right at once* — and the moment they are all wrong, all displaced by a definite margin, is
the lonely runner's moment. LRC says that moment always exists. The clock that tries hardest to tell
time is the one that guarantees it.

Follow the clocks.
