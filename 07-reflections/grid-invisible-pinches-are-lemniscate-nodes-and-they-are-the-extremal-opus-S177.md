---
source: opus-2026-07-09-S177
status: the grid-invisible pinches that broke mac-mini's widest-arc closure (MISTAKE-130) ARE the
  lemniscate nodes, and they ARE the LRC extremal (maxgap=1/7 <=> M=1/14 exactly). DEMONSTRATED: uniform
  grids over-merge across them (dissociated Sidon: coarse grid sees 8 arcs / maxIntG*spread=53.3, fine
  grid 22 arcs / 30.9 -- 72% over-estimate; strict and non-strict arc-counts IDENTICAL on every grid
  because the grid never samples the exact-1/7 point). SYNTHESIS: pinch = lemniscate node (r=0 on
  r^2=cos2theta, the figure-eight self-crossing) = tight loneliness (M=1/14) = measure-zero exact
  resonance (S168 broken clock) = the Eisenstein/CM certificate point (Phi_6, OPEN-Q-110). RESOLUTION:
  LRC(14) is a NON-STRICT statement (M>=1/14); the pinch belongs to the CLOSED good set and is good
  non-strictly -- so the correct a-priori route uses the closed good set (E_grid / non-strict), NOT the
  open widest-arc maxIntG (grid-unmeasurable). The lemniscate is the geometric archetype of "the node
  the coarse view merges but the exact view separates."
tags:
  - lrc14
  - grid-invisible-pinches
  - lemniscate
  - node
  - non-strict
  - complex-multiplication
  - out-of-box
---

# Grid-invisible pinches are lemniscate nodes — and they are the extremal

**opus-2026-07-09-S177.** Owner: consider how grid-invisible pinches relate to the lemniscate of
Bernoulli; get abstract.  Same day, mac-mini-S64 retracted the widest-arc dissociated closure
(MISTAKE-130) on exactly this: "the strict boundary maxgap=1/7 is hit at rationals (measure-zero pinches
invisible to every grid) ⟹ the 117k search over-merged arcs."  Here is the connection, grounded and
followed to its resolution.

## The pinch, demonstrated

The good set `G = {x : maxgap{frac(e_i x)} > 1/7}` is PINCHED wherever `maxgap(x) = 1/7` exactly: a
measure-zero point where `G` touches its complement and splits into two arcs meeting at that point.
Uniform grids cannot see it (measure zero + exact-boundary value), so they OVER-MERGE the adjacent arcs.
Measured (`lrc14_grid_invisible_pinches`):

| cluster | coarse grid (2k) | fine grid (20k+) |
|---|---|---|
| dissociated Sidon (spread 140) | **8 arcs, maxIntG·spread = 53.3** | **22 arcs, 30.9** |
| 7-structured (spread 82) | 64 arcs, 4.84 | 72 arcs, 4.41 |

The coarse grid merged 14 arcs and over-estimated the widest arc by 72%.  And the STRICT (`>1/7`) and
NON-STRICT (`≥1/7`) arc counts are IDENTICAL at every resolution — because no uniform grid ever samples
the exact `1/7` point, so it literally cannot distinguish "pinched-open" from "closed."  That is why no
grid measures `maxIntG`, and why the widest-arc pigeonhole failed.

## The pinch IS the lemniscate node

The lemniscate of Bernoulli `r² = cos 2θ` is a figure-eight with a NODE at the origin: `r = 0` at
`θ = π/4, 3π/4` — a single point where the curve self-crosses, pinching it into two loops.  The good
set's `maxgap = 1/7` points are the SAME thing: the self-crossing where `G` pinches.  The node is a
measure-zero SINGULAR point; the "coarse" (grid / real-locus) view merges the branches, and only the
"exact" (algebraic / normalized) view separates them.  Resolving the lemniscate node is a BLOW-UP
(normalization) that separates the two branches; resolving an LRC pinch is deciding strict-vs-non-strict
there — which the grid cannot do, but the exact rational structure can.

This unifies the project's measure-zero threads: the pinch = the node = the **tight loneliness point**
(`maxgap = 1/7 ⟺ M(S) = 1/14` exactly, the LRC extremal) = the **measure-zero exact resonance** (opus-S168
broken clock: the runner exactly at `1/14`) = the **cyclotomic/CM certificate** (OPEN-Q-110: the deep-well
tight set is lonely only at `14/183 = Φ₆(14)` denominator).  Every "the difficulty is at a measure-zero
point" phenomenon in this project is a pinch, and the lemniscate is its geometric archetype.

## The out-of-box arithmetic: ℤ[i] vs ℤ[ω], and the non-constructible 7

The lemniscate has complex multiplication by the Gaussian integers `ℤ[i]` (square lattice, `j = 1728`,
curve `y² = x³ − x`); its arc length is the QUARTIC elliptic integral `∫dr/√(1 − r⁴)`, and Abel's theorem
divides it into `n` equal arcs by ruler-and-compass iff `n` is a power of 2 times distinct Fermat primes
— exactly Gauss's polygon condition.  Two threads bite here:

- **The quartic.**  `√(1 − r⁴)` is quartic; the LRC's `ADDITIVE ENERGY` `E(S) = Σ_k|Ŝ(k)|⁴` (the L⁴ of
  the Fourier transform = `‖autocorrelation‖²`, THM-441/515B) that governs the singular series `L` and the
  density-floor variance is ALSO quartic.  The lemniscate's arc-length metric and the LRC's additive
  energy are the same species of 4th-power object — the natural "length" on the pinch-bearing locus.
- **`7` is not a Fermat prime.**  The lemniscate cannot be divided into 7 equal arcs (analogue of the
  non-constructible heptagon) — the arithmetic obstruction of the number `7`, which is the LRC's `θ = 1/7`
  threshold.  But the LRC's actual pinch denominators are Eisenstein (`Φ₆`, `ℤ[ω]`, `j = 0`, hexagonal),
  not Gaussian: `7 ≡ 3 (mod 4)` is INERT in `ℤ[i]` (lemniscate) yet `7 ≡ 1 (mod 3)` SPLITS in `ℤ[ω]`
  (Eisenstein) — and the LRC's `14 = 2·7` `Φ₆`-resonance lives in `ℤ[ω]`.  So the lemniscate is the
  ARCHETYPE (the node, the CM, the non-constructible-7), and the LRC's true home is its hexagonal cousin,
  the equianharmonic curve.  The invitation was the square lattice; the answer is the triangular one.

## The resolution (the LRC takeaway)

LRC(14) is a **non-strict** statement: `M(S) ≥ 1/14`, and the tight extremal (`= 1/14`) is a pinch.  So
**the pinch belongs to the CLOSED good set `{maxgap ≥ 1/7}` and is good non-strictly** (mac-mini's
surviving `LRCGoodPeriodNonStrict`, the knife-edge `M = 1/14`).  The widest-arc route failed because it
worked with the OPEN good set `{maxgap > 1/7}`, excluding the very points (pinches) that are the tight
loneliness — and the open set is grid-unmeasurable at its pinched boundary.  Concretely:

> **The correct a-priori good-period route must use the CLOSED (non-strict) good set, where the pinches
> are included; the geometric maxIntG (open) route is fundamentally grid-unmeasurable.**  This vindicates
> the `E_grid[W] ≥ 0` / non-strict routes (kps, mac-mini) over the widest-arc route, and explains
> MISTAKE-130 structurally: it measured an open-set quantity that lives only on the resolved (algebraic)
> locus, not the grid.

## Ledger

- DEMONSTRATED grid over-merging across pinches (dissociated: 8→22 arcs, maxIntG·spread 53→31; strict =
  non-strict on every grid) — grounds mac-mini MISTAKE-130.
- SYNTHESIS: pinch = lemniscate node = tight LRC point (`M=1/14`) = exact resonance (S168) = `Φ₆`/CM
  certificate (OPEN-Q-110).  Lemniscate = `ℤ[i]` archetype; LRC's true home = `ℤ[ω]` equianharmonic (`7`
  inert in `ℤ[i]`, splits in `ℤ[ω]`); additive energy `Σ|Ŝ|⁴` = the quartic arc-length species.
- TAKEAWAY: use the CLOSED (non-strict) good set; the open widest-arc `maxIntG` is grid-unmeasurable.
  Vindicates the `E_grid`/non-strict routes.  Files: `lrc14_grid_invisible_pinches_opus_S177` (+out).
- -> mac-mini-S64 (MISTAKE-130, LRCGoodPeriodNonStrict), opus-S168 (broken clock)/S169 (lemniscate
  node=collision)/S176 (two-scale), OPEN-Q-110 (Eisenstein `Φ₆`), THM-441/515B (additive energy).
