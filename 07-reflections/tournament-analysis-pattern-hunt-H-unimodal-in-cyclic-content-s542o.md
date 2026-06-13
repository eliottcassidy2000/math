---
source: oracle-2026-06-01-S542o
status: synthesis + computation (Tournament Analysis pattern hunt; H unimodal in cyclic content; comparator/threshold as axes)
tags:
  - tournament-analysis
  - basketball
  - lonely-runner
  - hamiltonian-path
  - H-loneliness
  - cut-cycle
  - patterns
---

# Tournament Analysis, Run Many Ways: H Is Unimodal in Cyclic Content, and the Comparator Is a New Axis

Tournament Analysis (the repo's central method, s471/s23): **pairwise data + a
switch (comparator) + a tie Hamiltonian path → a tournament-valued observable**,
studied by its fingerprints and trajectories as variables change. This session is
the empirical **pattern hunt** — running it many ways and connecting the dots to the
recent arc (S535–S541).

## The framework, sharpened

```
TA(data, comparator C, tie-path π)  =  the tournament  T  where  i→j  iff
   C(m_ij, m_ji) says so, ties resolved along π.
```

**The tie-path π is the basketball jersey order 1..5 = the tiling-model BASE
HAMILTONIAN PATH (S530/S531).** Resolving ties (`↔`, the S539 trienerment) along a
fixed Hamiltonian path is *exactly* how the tiling model fixes the base path and lets
the **tiles** be the deviations. So every TA output decomposes as

> **T = base-path backbone (π) ⊕ tiles (the metric's real signal) = cut ⊕ cycle =
> tension ⊕ flow (S537).**

The **tiles = the pairs where the data disagrees with the jersey/base order = the
cyclic content.** That is the object to watch.

## Six runs, six patterns (`tournament_analysis_pattern_hunt_s542.py`)

**P1 — basketball (discrete anchor).** Pass-count tournaments, jersey tie-break:
`star/ball-dominant → H=1` (a source/sink, transitive), `two-man-game → H=1`,
`balanced → H=15` (regular, max). So `H` (directed Hamiltonian-path count) is the
**team-balance meter**: hierarchical lineups have `H=1`, democratic ones max it —
the basketball face of `H = loneliness/spread` (S26).

**P2 — comparator zoo (the geometric/arithmetic dichotomy, S541, in TA).** For
runners, arc-distance and chord-distance (`chord = 2sin(π·arc)`, monotone) give the
**same** ranking — geometric metrics collapse. The signed cyclic **half-turn**
(`ahead`) is genuinely distinct from the **distance-rank** (`approaching`, which is
transitive). *Comparator choice = which geometry you expose.*

**P3 — the H-trajectory and the walk.** Over `t`, the half-turn tournament visits
**exactly `2·Fib(n-2) = 4` iso-classes** at `n=5` (the circular menu, S518 — confirmed
empirically), with `H` oscillating in `[1, 15]`. Fingerprint: the **AP/regular
speeds have higher mean `H` (10.0) than generic (8.0)** — arithmetic regularity
spends more time *spread*. The walk's statistics fingerprint the speed set.

**P4 — the comparator is a new axis.** Rotating the half-turn window by a phase `φ`
on **fixed** positions traces an orbit of **6** iso-classes (12 transitions) — *more*
than the `4` reachable by time-evolution. So beyond time and speeds, the **switch
functional itself is a variable**, with its own iso-class orbit.

**P5 — threshold percolation (trienerment).** As the tie-threshold `θ` grows, ties
percolate `0 → 4 (at θ=1/n) → 7 → 10 (all-tie at θ=1/2)`. `θ=0` is a pure tournament;
the **loneliness threshold `θ=1/n`** is the distinguished point where `near = tie`
(S539).

**P6 — the gem: `H` is unimodal in cyclic content.** Tabulating `H` by **tile-count**
(pairs disagreeing with the base path):

```
tiles:   0    1    3    4    5    6    7    9   10
avg H:   1    9   13.7 15   11   15  13.7  9    1
```

`H` is **symmetric and unimodal in the tile-count**, `= 1` at both extremes
(`tiles=0`: fully agrees with the base path = transitive; `tiles=10`: fully reverses
it = the opposite transitive) and **maximal in the middle** (`tiles≈4–6`: maximally
cyclic/mixed). So:

> **The loneliness meter `H` peaks at intermediate cyclic content** — the
> tournament is "loneliest/most spread/most balanced" exactly when its tiles
> (cut⊕cycle's cyclic part) are *half-on*, and is rigid (`H=1`) at both the aligned
> and anti-aligned extremes. `H` measures *how much the data has departed from its
> base-path backbone, capped at half-departure.*

## Connecting the dots

- **What the method is:** TA turns any pairwise observable into a point/trajectory in
  tournament space; the tie-path is the base Hamiltonian path; the output splits as
  base-path ⊕ tiles (cut ⊕ cycle, S537). Basketball and runners are the discrete and
  continuous faces of one method.
- **The invariant to read is `H`** (= loneliness/balance, S26), and P6 shows it is a
  **unimodal function of the cyclic content (tiles)** — a clean, general law of TA,
  not specific to runners.
- **The axes of variation:** the data (speeds/passes), **time** `t` (the walk through
  the `2·Fib(n-2)` menu, S518), the **comparator** `C`/phase `φ` (P4, a new axis), and
  the **threshold** `θ` (P5, the trienerment percolation, with `1/n` distinguished).
- **The dichotomy (S541) lives inside TA:** geometric comparators collapse to the
  circular order; arithmetic ones (differences, p-adic, Grundy) expose the channels.
- **LRC** is one TA observable: the runner half-turn tournament, asking the observer
  to reach the loneliness end of the walk (`H`/tile extreme).

## Verdict / next
- Tournament Analysis demonstrated across discrete (basketball) and continuous
  (runner) data, with the jersey = base-path tie-rule unifying them.
- New empirical law: **`H` is symmetric-unimodal in tile-count** (cyclic content),
  peaking at half-departure from the base path; rigid at both ordered extremes.
- New axis: the **comparator/phase** sweep (more iso-classes than time-evolution).
- Concrete next: (1) prove the `H`-unimodal-in-tiles law (it should follow from the
  cut⊕cycle structure + the `H=1+2^{r-1}` apex law, S531); (2) the comparator-orbit
  (P4) as a group action on iso-classes; (3) the AP-higher-mean-`H` fingerprint as a
  detector of arithmetic structure in any pairwise data.

## Artifacts
```
04-computation/tournament_analysis_pattern_hunt_s542.py
05-knowledge/results/tournament_analysis_pattern_hunt_s542.out
```
Related: s471/s23 (TA framework), s512 (TA capstone), S518 (the 2·Fib menu), S26
(H=loneliness), S530/S531 (base path/apex/recursion), S537 (cut⊕cycle/flows), S539
(trienerment/ties), S541 (geometric/arithmetic dichotomy).
