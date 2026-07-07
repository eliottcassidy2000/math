# LRC(14)'s tight family is sieve-easy; the genuine hard core (saturated) carries margin

**opus-2026-07-07-S131.** A census of the *actual* LRC(14) families (13 nonzero integer speeds,
threshold 1/14) — following the owner's "census the families / deep structural correctness"
directive — clarifies where the difficulty really lives, and gently corrects a framing drift.

## The sieve dichotomy (exact, standard)

For any `q ∈ {2,…,14}`: if **no** speed is divisible by `q`, then at `t = 1/q` every speed sits at
`‖v/q‖ ≥ 1/q ≥ 1/14`, so `M(v) ≥ 1/14` — lonely. (This is the GREEN `lrc14_no_multiple_of_14` /
`counterexample_needs_all_divisors`.) So an LRC(14) counterexample must be **saturated**: a multiple
of *every* `q ∈ {2,…,14}`. Two disjoint cases:

- **Non-saturated** (misses some `q`): `M ≥ 1/q ≥ 1/14`, sieve-handled. This class **contains the
  extremal AP** `{1,…,13}` (misses only `q=14`; `M = 1/14` exactly at `t=1/14`).
- **Saturated** (a multiple of every `q ≤ 14`): the genuine hard core.

## What the census shows (and the framing it corrects)

**The unique tight family `M = 1/14` is the AP `{1,…,13}` — and it is NON-saturated, hence
sieve-easy.** Verified: `M(AP) = 1/14` at `t = 1/14`; the second "tight" family the census first
surfaced, `{1,…,11,13,24}`, is tight for the *same* sieve reason (residues mod 14 avoid 0, hit ±1).
So the extremal locus is exactly the sieve boundary, and it is handled by a one-line witness.

**The saturated hard core carries MARGIN, not near-tightness.** Over thousands of saturated
13-families (constructed + `{1..15,18}`-range exhaustive):

- min `M` found `= 1/12 ≈ 0.0833` (`{1,2,3,4,10,…,18}`); the most-AP-like saturated
  `{1,2,3,4,5,7,…,14}` gives `M = 2/23 ≈ 0.087`. **Zero** saturated families below `1/14`; none even
  reaches `1/13 ≈ 0.0769`. All sit `≳ 0.083`, a clear margin above `1/14 = 0.0714`.
- consecutive blocks `{a,…,a+12}` for `a ≥ 2` are saturated with `M = (a+…)/(2a+14)` via the
  midpoint witness `t = 1/(2a+14)`: `1/8, 1/6, 1/5, …` — margin *grows* with `a`. Larger/spread
  saturated families decorrelate → **more** margin, so the extremal saturated family is *small*.

So the sieve hard core is not a "near-tight moat" — the near-tight family (the AP) is on the *easy*
(non-saturated) side, and the hard (saturated) side has `≥ ~0.083` empirically. **This corrects the
drift of calling the AP "the hard tight family" of LRC(14): for LRC(14) the AP is sieve-trivial.**
(The AP *is* the hard extremal for the 12-speed `(C)` gap at `2/25` — a different problem; the
framing leaked across.)

## Reconciliation with the fleet's "moat"

The fleet's "near-AP moat" (kps-S53, klein-S152, mac-mini-S39) lives in the **coarse reduction**:
families `vᵢ = aᵢ + L·kᵢ` whose *coarse part* `K = {kᵢ}` is the AP `{1,…,13}`, where the coarse
bound `M(v) ≥ M(K) − A/L` degrades below `1/14`. Those are **multi-scale** (large `L`) and are
handled by klein-S152's conjugate witness (verified + proven for `L ≳ 200A`, opus-S131). That "AP" is
the *coarse part*, not the raw speeds — a different object from the sieve's saturated core. Both are
real; neither is the near-tight raw AP.

## The honest open core

LRC(14) `⟺` **every saturated 13-family has `M ≥ 1/14`.** Empirically they have margin (`≥ ~0.083`),
and the extremal is small, so it *looks* provable by "small saturated check + large-saturated
decorrelation." But — the S130 lesson — this is a **bounded-range census**, not a proof: a uniform
`M ≥ 1/14` over *all* saturated families (arbitrarily large, adversarially structured) is exactly the
open analytic content, and the empirical margin must not be read as closure. The margin is real
evidence that LRC(14) holds comfortably; proving it uniformly is the crux the density-floor /
decorrelation machinery targets.

## Takeaways

1. LRC(14)'s tight family (AP) is **sieve-easy** (non-saturated); do not treat it as the hard core.
2. The sieve hard core is the **saturated** families, which carry **margin** (`M ≳ 0.083`), extremal
   small, larger ones decorrelate to more margin.
3. The fleet's "AP moat" is the **coarse-part-AP** (multi-scale, klein-handled) — a distinct object.
4. Open: the uniform `saturated ⟹ M ≥ 1/14` bound — empirically comfortable, analytically the crux.
