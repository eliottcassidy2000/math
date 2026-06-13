---
source: oracle-2026-06-01-S552o
status: proof attempt (n=14 via "the 7 impossibility") — CRT 7-reduction confirmed; the 7-gon-window construction provably fails on the singleton class
tags:
  - lonely-runner
  - n14
  - crt
  - seven
  - construction
  - honest-negative
---

# "The 7 Impossibility" for LRC@14, Tested: the CRT 7-Reduction Is Right, the 7-Gon-Window Construction Provably Fails

Following the prompt to make "the 7 impossibility" the key to LRC@14. The structural
half is right and confirmed; the natural *constructive* route (post-S551) is tested
and **provably insufficient** — which sharpens exactly where the difficulty sits.

## The structure (confirmed, S524): LRC@14 = the 7 mod-7 classes

`n = 14 = 2·7`, `n* = 7` prime (S546). The 13 runners partition into the **7 mod-7
CRT classes**: six pairs `{i, i+7}` (residues `1..6`) and the **singleton `{7}`**
(the multiples of 7). **LRC@14 ⟺ all 7 classes are simultaneously safe.** Computed
class-safe measures for the AP reproduce S524 exactly: singleton `{7}` safe `6/7 ≈
0.857`; each pair-class safe `≈ (6/7)² ≈ 0.73`. It is a coupon-collector over 7
classes; the open part is the **7-way correlation** (the AP gives a measure-zero,
wall-only intersection).

## The constructive attack and why it fails

Per S551, the AP core is measure-zero and needs **construction**, not measure. The
natural construction is the **7-gon-vertex windows** `t = j/7 + s` (`j = 1..6`):

- a **pair-class** runner `v = 7u + c` (`c ≢ 0`) has `‖v·(j/7)‖ = ‖cj/7‖ ≥ 1/7 =
  2/14` (safe with margin), and stays safe for `|s| ≤ 1/(14V)` (`V = max speed`);
- a **singleton** runner `7w` has `7w·(j/7) = wj ∈ ℤ` (blocked at the vertex), with
  `‖7w·(j/7+s)‖ = ‖7ws‖` — so near the vertex the *whole* problem collapses to
  clearing **only the multiples-of-7 sub-system** in the small `s`-window.

**Computed (`lrc_n14_seven_impossibility_s552.py`):** the window construction found a
lonely time for only **5/25** sets — and **all five have `r = 0` multiples of 7**
(where the singleton class is empty and the sieve `t = 1/7` already wins). For every
set with `r ≥ 1`, the windows **failed**, even though a lonely time exists elsewhere
(`24/25` lonely somewhere).

**And this is provable, not just sampled.** The window half-width is forced to
`δ = 1/(14V)` (to keep the pair-classes safe). Within it, `7w·s` ranges over
`[0, 7wδ] = [0, w/(2V)] ⊆ [0, 1/14]` (since `w ≤ V/7`), so `‖7w s‖ = 7w s ≤ 1/14`
with **equality only at the very edge** — the windows **cannot** clear even a single
generic multiple of 7, let alone several at once.

> **The 7-gon-vertex windows are too small to clear the singleton class.** The naive
> constructive form of "the 7 impossibility" handles only the no-multiple-of-7 case
> (`= the sieve`); for `r ≥ 1` the lonely times live *away* from the vertices, where
> the singleton clears **while** the pair-classes remain safe — and that simultaneity
> is exactly the uncontrolled 7-way coupling.

## What this actually tells us (the honest verdict)

- **"7" is the right structural lens** (the CRT 7-class reduction, S524), but it does
  not, by itself, give an impossibility: neither the **measure** bound (S550/S551,
  blind to the measure-zero core) nor the **7-gon-window construction** (provably too
  small for the singleton class) closes LRC@14.
- The residual is pinned precisely: **the singleton `{multiples of 7}` sub-system,
  coupled to the six pair-classes.** A counterexample would need the multiples-of-7 to
  block every region where the pair-classes are simultaneously safe — the 7-way
  correlation, restated as a coupling between the singleton and the pairs.
- The constructive path is *not* the small windows; it must be a **wider** witness
  that clears the multiples-of-7 at a time when the six pair-classes are jointly safe.
  The pair-classes are jointly safe on a set of measure `≈ (6/7)^? `; the singleton's
  danger set has measure `2/14 = 1/7`; the real question is whether the singleton
  danger can be *contained* inside the pair-classes' joint-unsafe set for all `t`
  (covering) — which it cannot for non-AP sets (they are lonely somewhere), and only
  marginally (wall-only) for the AP.

## Verdict / next
- Confirmed the CRT 7-class reduction; tested "the 7 impossibility" via the 7-gon
  windows and showed it is **provably insufficient** (handles only `r=0` = the sieve;
  the windows are too small for the singleton class).
- The crux is the **singleton ↔ pair-class coupling** (the 7-way correlation, S524),
  untouched by both measure (S551) and the naive construction.
- Concrete next: (1) a **wider** singleton witness — clear the multiples-of-7 on the
  pair-classes' joint-safe set (not just near the vertices); (2) bound the 7-way
  correlation via the near-independence (S524 ratios `> 0.99`) with an explicit error
  term — does `(6/7)^? − error > 0`?; (3) treat the singleton sub-system under the
  fast clock `τ = 7t` (it becomes the reduced speeds `{w_k}`) and ask for a *joint*
  witness with the pair-classes — the genuine residual.

## Artifacts
```
04-computation/lrc_n14_seven_impossibility_s552.py
05-knowledge/results/lrc_n14_seven_impossibility_s552.out
```
Related: S524 (CRT 7-class reduction / coupon collector), S525/HYP-2003 (covering
wall), S546 (n=14 = 2·7, n*=7 prime), S550 (resonance-energy measure bound), S551
(Vitali: the core needs construction), THM-369 (sieve = the r=0 witness).
