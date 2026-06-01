---
source: oracle-2026-06-01-S536o
status: synthesis + computation (dual mapping: sectors as nodes; the sector-tournament; the DFT duality with the resonance picture)
tags:
  - lonely-runner
  - sectors
  - dual-mapping
  - occupancy
  - discrepancy
  - dft-duality
  - tournament
---

# Sectors as Nodes: the Dual LRC Mapping, the Most-Restricted Sector-Tournament, and the DFT Duality

Every previous mapping put the **runners** at the nodes. This one dualizes: the
nodes are the `n` **evenly-spaced sectors** of the circle (the regular `n`-gon
cells), and the runners are moving markers that **toggle edges as they cross sector
boundaries**. The sector width is `1/n` — exactly the loneliness threshold — and
this single alignment makes the whole picture snap together.

## The set-up and the loneliness reformulation

Sectors `S_k = [k/n,(k+1)/n)`, `k = 0..n-1`. The observer at `0` sits on the
boundary of `S_0` and `S_{n-1}`. Runner `i` (speed `v_i`) occupies sector
`⌊n·frac(v_i t)⌋`; the occupancy vector is `c(t) = (c_0,…,c_{n-1})`, `Σ c_k = n-1`.

> **Observer lonely ⟺ the two sectors touching it are empty: `c_0(t) = c_{n-1}(t) = 0`**
> (all runners in `[1/n, 1-1/n]`). So **LRC@n ⟺ every speed set reaches an occupancy
> vector with the observer's two cells empty** as `t` varies.

**Pigeonhole gives the baseline.** `n-1` runners in `n` sectors ⇒ there is **always
≥ 1 empty sector**. The conjecture is the steering problem: can the empty cell be
driven to the observer, and can a *second* adjacent cell be emptied there too (the
`2/n` gap straddling `0`)? The regular `n`-gon (speeds `1..n-1`) is the tight
extremal — its empty cell sits exactly at the observer only at the times `t = k/n`.

The dynamics is literally "edges change when runners cross boundaries": runner `i`
crosses a sector boundary `n·|v_i|` times over `[0,1)`, each crossing moving one unit
of occupancy between adjacent cells — a closed walk through occupancy states (total
`Σ n|v_i|` steps, the holdback of S25 in the sector frame).

## The sector-tournament: the most restricted tournament mapping yet

Rank sectors by occupancy with a cyclic tiebreak: `a → b` iff `c_a > c_b`, or
`c_a = c_b` and `a→b` is the short way round. Edges flip exactly on boundary
crossings. Computed realizable iso-classes
(`lrc_sector_occupancy_dual_mapping_s536.py`):

```
 n   sector-tournament realizable / A000568(n)    R       occupancy-iso / all-compositions
 4            3 / 4                               0.75            13 / 20
 5            4 / 12                              0.33            38 / 70
 6           10 / 56                              0.18            136 / 252
```

This is **more restrictive than any prior tournament mapping**: at `n=6` only
`10/56 ≈ 0.18` of tournament iso-classes are realizable as sector-occupancy
tournaments, versus `0.35` for the near-graph (S535) and the unrestrictive circular
menu of the runner mapping (S518). And `LRC = membership`: every speed set reaches
the class where the observer's two sectors are joint minima (empty). Verified
`120/121, 120/121, 60/61` — the lone misses are again the **tight AP set** (lonely
only at the measure-zero `t = k/n` the grid skips), the boundary extremal.

The reason for the strong restriction is structural: the realizable occupancy vectors
are compositions of `n-1` into `n` cyclically-arranged parts arising from a
*linear flow*, and the occupancy-ranking tournament collapses many of them — the
sector node-set is a tiny, rigid stage compared to the `A000568` of all tournaments.

## The payoff: sectors are the DFT-dual of the resonance picture

The discrete Fourier transform of the occupancy vector is **exactly** the S529
exponential/character sum (verified to `~10^{-16}`):

```
ĉ_m = Σ_k c_k e^{-2πi m k/n}  =  Σ_j e^{-2πi m ⌊n·frac(v_j t)⌋ / n}      (max error 6e-16)
```

So the **sectors are the real-space dual of the Fourier / resonance / covering
picture (S529–S535)**. Two faces of the same object:

- **real space (this mapping):** occupancy of the `n` cells; loneliness = the
  observer's cells are empty;
- **frequency space (S529):** the character sums `ĝ(m)`, resonances `Σ m_i v_i = 0`,
  the inside debt and the three-channel parity law (S533).

"Empty observer cells" (a real-space discrepancy: a `2/n`-wide hole) is the DFT-dual
of "the character/covering condition fails." This unifies the sector mapping with the
entire resonance program: the discrepancy of `{v_j t}` against the equipartition
(the occupancy fluctuation) *is* the magnitude of the exponential sums, and LRC asks
that this discrepancy be one-sidedly large at the observer at some time.

## Why this is the right kind of restriction

The sector mapping retains the metric (the cells *are* the `1/n` grid), so by the
S535 principle it restricts strongly — and it does so while keeping a clean node set
(the `n`-gon) and a clean target (observer cells empty). It also reframes LRC as a
**discrepancy / empty-cell-steering** problem and exposes the exact bridge to the
Fourier side. The conjecture's content is now visibly the same on both faces: steer
the empty cell to the observer (real space) ⟺ make the exponential sums align to
clear the observer band (frequency space).

## Verdict / next
- New dual mapping (sectors as nodes); LRC = "observer's two cells empty"; pigeonhole
  guarantees ≥1 empty cell always.
- The sector-tournament is the most restricted tournament mapping computed
  (`R = 0.75, 0.33, 0.18`); LRC = membership; tight AP = boundary extremal.
- DFT duality verified exactly: sectors ⟷ S529 characters. Loneliness ⟷ character
  condition.
- Concrete next: (1) characterize the realizable sector-tournament classes (the
  `3,4,10,…` sequence) and the occupancy compositions a linear flow can realize;
  (2) phrase LRC as a one-sided discrepancy bound on the occupancy at the observer
  cell and connect to the Erdős–Turán / cascade discrepancy gap (S527/S534);
  (3) the "empty-cell walk" as an exclusion process — does the empty cell visit every
  position with multiplicity 2?

## Artifacts
```
04-computation/lrc_sector_occupancy_dual_mapping_s536.py
05-knowledge/results/lrc_sector_occupancy_dual_mapping_s536.out
```
Related: S535 (mapping spectrum / metric restriction), S529 (characters/inside debt),
S533 (three-channel parity — the frequency dual), S530 (apex = largest gap = the empty
cell), S525/S527 (covering / cascade / discrepancy).
