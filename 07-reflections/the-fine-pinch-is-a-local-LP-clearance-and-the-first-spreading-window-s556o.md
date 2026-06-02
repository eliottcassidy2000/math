---
source: oracle-2026-06-01-S556o
status: progress (the fine pinch = a local LP clearance; the first-spreading-window; the wall as the LP degeneracy)
tags:
  - lonely-runner
  - fine-pinch
  - linear-programming
  - margins
  - first-window
  - wall
---

# The Fine Pinch Is a Local LP Clearance — and Lonely Times Live in the First Spreading Window

Continuing S555's fine-pinch salvage. The naive pigeonhole cannot work at any
denominator (total danger `2 - 2/n > 1`), so the fine pinch must exploit
**correlation**. By S553, at the optimal time **at most one runner is near** and the
rest have **margin** — so the fine pinch is a *local linear-programming clearance*,
not a global count.

## The clearing condition (local LP)

At a time `t0` with near set `{w*}` (single near runner), let
`d = 1/n − ‖w* t0‖ ≥ 0` (its deficit) and `m_i = ‖v_i t0‖ − 1/n > 0` (far margins).
Perturb `t0 → t0 + ε`:
- far `i` stays far iff `|v_i ε| ≤ m_i`  ⟺  `|ε| ≤ m_i/v_i`;
- `w*` clears (outward) iff `|w* ε| ≥ d`  ⟺  `|ε| ≥ d/w*`.

> **A clearing fine pinch `t0+ε` exists iff `d/w* ≤ min_i m_i/v_i` (`*`).** The
> perturbation is fine (`ε ~ d/w*`). For a `k`-near set it generalizes to a small
> linear feasibility: `|v_i ε| ≤ m_i` (fars) and `sign(w_j ε)` outward with
> `|w_j ε| ≥ d_j` (the `k` near runners) — a local LP.

So the fine pinch is **a local LP**, whose size is controlled by the number of near
runners — which S553 caps at `1` generically.

## What the computation shows (`lrc_fine_pinch_perturbation_s556.py`, n=14)

- **Generic sets** are lonely **directly**: their min-near time has `near = 0` (no
  perturbation needed), at a **fine** time `t0 ≈ 0.008–0.036`, all **below `1/14`**.
- **The AP / regular polygon is the wall:** its min-near time has `near = 1` with
  `d/w* = 0.0357` but `min_i m_i/v_i = 0` — the far runners sit *exactly* at `1/n`,
  margins vanish, `(*)` fails, no perturbation clears. This is precisely the
  measure-zero wall (S551).

So `(*)` (and its `k`-near LP generalization) **holds with room off the wall and
degenerates exactly at the wall** where the far margins collapse to `0`. LRC at the
core is the statement that the local LP is feasible at some time for every speed set
except the (measure-zero) AP extremal, where it is feasible only in the closed limit.

## The first spreading window (a constructive lead)

The fine lonely times cluster at **small `t` (below `1/n`)**: as `t` grows from `0`,
all runners start bunched at the observer and *spread*; once each `v_i` has passed
`1/n` (at `t ≈ 1/(n v_i)`) and before the fast ones wrap, there is a **first lonely
window**. Empirically every generic set is lonely there (`t0 ≈ 0.008–0.036 < 1/14`).
This suggests a constructive target:

> **First-window conjecture (fine-pinch form):** every speed set (off the wall) has a
> lonely time in `(0, c/n)` for an absolute constant `c` — the first window after the
> runners separate from the observer. The wall (AP) is lonely only at the window's
> closed edge.

If provable, this would *locate* the fine pinch in a bounded initial window, turning
LRC@core into a bounded-`t` LP feasibility (a finite, structured search) rather than
a search over all `[0,1)`.

## Why this is the right object (and the honest limit)

- The global count fails (`2 - 2/n > 1`); the **local LP** (clear the `≤1` near runner
  within the far margins) is what governs the fine pinch, and S553 keeps it tiny.
- The **wall is exactly the LP degeneracy** (margins `→ 0`): a clean, precise
  characterization of the unique hard configuration.
- **Honest limit:** this reformulates LRC@core as "the local LP is feasible off the
  wall," and the first-window conjecture would localize it — but proving the LP is
  feasible (the gap is nonempty) for every wall-adjacent set is still the conjecture.
  The fine pinch is *constructed* by the LP wherever margins are positive; the open
  part is bounding the margins below away from the AP.

## Verdict / next
- The fine pinch = a **local LP clearance** (condition `(*) d/w* ≤ min m_i/v_i`),
  small by S553; generic sets are lonely directly in the **first window** `(0, 1/n)`;
  the AP is the unique LP degeneracy (margins `0`).
- Concrete next: (1) **prove the first-window conjecture** (lonely in `(0, c/n)` off
  the wall) — a bounded, finite LP feasibility; (2) bound the far margins below for
  non-AP sets (the quantitative gap from the wall); (3) the multi-near LP feasibility
  region's geometry vs the resonance structure (S550).

## Artifacts
```
04-computation/lrc_fine_pinch_perturbation_s556.py
05-knowledge/results/lrc_fine_pinch_perturbation_s556.out
```
Related: S555 (fine pinch / sieve reduction), S553 (at most one near — caps the LP),
S551 (Vitali / measure-zero wall), S550 (resonance core), S530 (apex/largest gap),
S552 (7-gon windows = a special fine pinch).
