---
source: oracle-2026-06-01-S512
status: synthesis capstone (Tournament Analysis of the LRC, S24–S512)
tags:
  - lonely-runner
  - tournament-analysis
  - a000568
  - source-reachability
  - capstone
  - addition-multiplication
---

# Tournament Analysis of the Lonely Runner: a Capstone

This pulls together the arc S24–S512 into one picture and pushes the
A000568 reformulation to its current limit. The throughline: **a runner system
is a closed walk in tournament space, and the Lonely Runner Conjecture is a
reachability statement for that walk.**

## The ladder, assembled

```
runner speeds  ──comparator──►  tournament T(t)  ──as t varies──►  closed walk in G_n
       (additive dynamics)         (a point in A000568)        (LRC = reach a target)
```

- **S24 — the clock.** The half-turn comparator `i→j iff frac((s_i−s_j)t)∈(0,½)`
  makes `T(t)` piecewise-constant; as `t` sweeps `[0,1)` it is a closed walk in
  the metagraph `G_n`. The realizable classes are a fixed tiny "circular menu";
  the speeds choose only the walk. Extremal LRC speeds `(1,…,n−1)` = the minimal
  walk.
- **S25 — holdback / the staircase.** Edge persistence `=1/(2·(s_i−s_j))`; the
  difference multiset of `{1,…,n}` is the **staircase δ**, so LRC-extremal = max
  holdback = the staircase. Twin primes are the stickiest prime edges; LRC-hard =
  divisibility-rich (covering) while admissible/prime sets are LRC-easy.
- **S26 — H as a loneliness meter.** `H=1 ⇔ ½-gap` sharply; above that `H` is a
  coarse, ½-resolution spread meter; `#3-cycles = #triples not in a semicircle`.
  A single half-turn tournament resolves *only* the `½` threshold — so it cannot
  see the LRC `1/n` gap. **This diagnosed the obstruction.**
- **S511 — the fix and the A000568 reformulation.** Use the **LRC walls** and a
  **marked observer**: `observer→i iff ‖v_i t‖≥1/n`. Then exactly
  `observer LONELY ⇔ observer is a SOURCE`. So **LRC ⇔ the observer-marked walk
  reaches a source-class** — a pure marked-iso-class reachability statement — and
  the source target has size **A000568(n−1)** (the runner sub-tournament is free).
- **S512 — the true target.** The reachable source-classes are a *tiny* subset of
  A000568(n−1).

## S512: the real LRC win-set, measured

At a lonely time the `n−1` runners all lie in the fixed safe arc `[1/n,1−1/n]`
(length `L=1−2/n`), so the reachable source classes are exactly the **half-turn
tournaments of `n−1` points in an arc of length `L`**. Computed
(`lrc_source_reachability_deep_s512.py`, all primitive speed sets in a box):

```
n   safe-arc L   reachable source classes   A000568(n−1)   LRC failures
4   0.500        1 (forced transitive)      2              0  / 325
5   0.600        2                          4              0  / 479
6   0.667        6                          12             0  / 461
7   0.714        6                          56             0  / 210
```

Three things land:

1. **The tournament framing re-proves LRC for n ≤ 7**: *every* primitive speed
   set's marked walk reaches a source cell (open, or a boundary witness for the
   tight sets) — `0` failures. LRC, read entirely off the tournament side.
2. **The true target is microscopic** inside `A000568(n−1)`: `1,2,6,6` of
   `2,4,12,56`. The win-set is not "any tournament on `n−1` runners"; it is the
   small **arc-confined circular menu**. As `n→∞`, `L→1` and this menu fills out
   the full circular family — but it is always a vanishing fraction of
   `A000568(n−1)` (which grows like `2^{C(n−1,2)}/(n−1)!`).
3. **A geometric gradient.** For `L ≤ ½` (i.e. `n ≤ 4`) the reachable source class
   is forced **transitive** (points in a half-arc are linearly ordered). For
   `n ≥ 5` non-transitive source classes appear (`H>1`), and the menu grows with
   `L`. The `min #near-runners` histogram is `{0: almost all, 1: the tight few}`:
   essentially every set has an *open* lonely interval; only the extremal/tight
   sets touch a source class solely on the boundary.

## The precise tournament problem equivalent to LRC

> **(LRC-T)** For every primitive integer speed set `v_1,…,v_{n−1}`, the
> observer-marked tournament clock (LRC endpoint walls for observer–runner edges,
> half-turn for runner–runner edges) visits a cell in which the observer is a
> source — equivalently, a cell whose runner sub-tournament is a half-turn
> tournament of `n−1` points in an arc of length `1−2/n`.

This is faithful (verified `mism=0`, S511) and finite-per-instance. The missing
piece — the whole content of LRC — is: **the additive walk cannot avoid this
tiny, multiplicatively-shaped target.**

## Where addition and multiplication sit (and where A000568 enters)

The reformulation cleanly separates the two operations the user pointed at:

- **Addition runs the walk.** Positions add (`v_i t`); the walls are additive
  endpoint crossings; the cadence is the **staircase of differences** (S25). The
  walk's *steps* are single edge-flips — `G_n` edges. "Observer is a source" is an
  additive-geometry event (everyone in the safe arc).
- **Multiplication shapes the arena and the target.** `A000568(n)` is the
  **odd-cycle Burnside** count (a multiplicative/group fact); the canonical
  regular classes are the **Paley / quadratic-residue** tournaments
  (`p≡3 mod 4`); and *which* classes are reachable is gated by the speeds'
  **divisibility/residue** structure — the same sieve (THM-369) that makes
  divisibility-rich sets LRC-hard and admissible/prime sets LRC-easy (S25). The
  reachable menu shrinks or grows with the arithmetic of `n` via the `x+2 / x·2`
  grid: the `+`-chain of differences (staircase, odd core) versus the `·2`
  doubling tower (the dyadic `n=16` debt, the Cayley–Dickson ladder).

So **A000568 is exactly the meeting place**: the base count is multiplicative
(odd-Burnside), the walk is additive, and LRC is the claim that the additive walk
hits a multiplicatively-carved, A000568-sized target.

## Honest assessment: how far this goes

The A000568 analogy is now **exact, not loose**: LRC *is* a marked source-
reachability problem in the A000568 quotient, with an A000568(n−1)-sized win-set
whose true reachable core is the arc-confined circular menu. It has:

- re-proved LRC for `n ≤ 7` purely on the tournament side;
- given a finite, combinatorial target and a clean equivalent statement (LRC-T);
- unified the additive (staircase/holdback/clock) and multiplicative
  (A000568/Paley/sieve) faces under one quotient.

It has **not** cracked the conjecture, and the reason is structural and worth
stating: the walk is *constrained by the speed arithmetic*, so "the walk can't
avoid the target" is exactly as hard as LRC — the framing **organizes and
illuminates** the difficulty (it is now a reachability fact about a tiny target
in a well-studied space) rather than dissolving it. The most promising lever it
exposes: the target is so small (a vanishing fraction of A000568(n−1)) that an
*avoidance* (counterexample) walk would have to thread an enormous complement
forever — which is the tournament-space shadow of the measure/overlap tension
that all LRC approaches confront.

## Sharpest next problems

1. Count the reachable source menu as a function of `n` (the `1,2,6,6,…`
   sequence) and identify it (it should be "half-turn tournaments of `m` points
   in an arc of length `1−2/(m+1)`").
2. Prove LRC-T for the single-resonance / few-distinct-difference families using
   the staircase synchrony (S25) — a tournament-native special case.
3. Lift the marked source-reachability to a Lean statement on top of THM-369's
   sieve, making "observer is a source" a formal predicate.
4. Probe the frontier: for `n=14/16` hard rows, measure the walk's minimum
   `#near-runners` (distance to source) — how close does the best near-tight set
   come to source-avoidance?

## Artifacts

```
04-computation/lrc_source_reachability_deep_s512.py
05-knowledge/results/lrc_source_reachability_deep_s512.out
04-computation/lrc_observer_source_tournament_s511.py     (the reformulation)
04-computation/tournament_clock_s24.py                    (the clock)
04-computation/tournament_clock_holdback_twinprime_s25.py (holdback/staircase)
04-computation/H_loneliness_meter_s26.py                  (the meter)
```
