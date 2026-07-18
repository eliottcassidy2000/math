# The density-route discharge (separated far element) is now kernel-pure Lean

*boxeph-2026-07-18-S107. Owner: formalize the density-route discharge. Done — `LRCDensityDischarge.lean`
adds three kernel-pure theorems (`[propext, Classical.choice, Quot.sound]`, no sorry): the density
route's Φ>0 (good-set-nonempty) mechanism as a self-contained far-element extension, with a proof
DISTINCT from the descent floor's round+kick. Built into the corpus (8477 jobs) and registered.*

## What was formalized

The density route (S96–S100) closes a far-element family by keeping the good set nonempty. Its geometric
core, now in Lean:

**`density_far_extension`** (PROVED, general `ι`):
```
(v : ι → ℤ) (vstar : ι) (t0 V : ℝ) (hV : 0 < V)
(hframe : ∀ i ≠ vstar, ∀ m, 1/13 ≤ |v_i·t0 − m|)   -- frame 1/13-lonely at t0
(hbound : ∀ i ≠ vstar, |v_i| ≤ V)                   -- frame speeds bounded by V
(hfar   : 91·V ≤ v vstar)                            -- far element separated
⊢ ∃ t, Lonely 14 v t
```
Two-step proof, exactly the density argument:
1. **Good interval.** The frame is `1/14`-lonely on the whole interval `[t0−δ, t0+δ]`, `δ = 1/(182·V)`:
   by the reverse triangle, `|v_i·t − m| ≥ |v_i·t0 − m| − |v_i|·|t−t0| ≥ 1/13 − V·δ = 1/13 − 1/182 = 1/14`.
2. **Far completion.** Since `2δ = 1/(91V) ≥ 1/(v vstar)` (from `hfar`), the interval contains a
   half-integer point `t = (k+1/2)/d` (`d = v vstar`, `k = ⌈d(t0−δ)−1/2⌉`), where the far runner sits at
   `‖d·t‖ = 1/2 ≥ 1/14`. So the whole family is lonely at `t`.

Supporting lemma **`half_integer_far`**: `1/2 ≤ |k + 1/2 − m|` for all integers `k, m` (a half-integer is
≥ 1/2 from every integer). And **`density_far_bridge`**: the complete LRC(14) rung — the LRC(≤13) citation
(`LRCUpTo13`) supplies the frame's `1/13`-loneliness, so a 13-family with 12 speeds bounded by `V` and far
element `≥ 91·V` is `1/14`-lonely, unconditionally on the citation.

## Why this is a distinct proof (and what it records)

This is the density route's actual mechanism — `Φ(E') > 0` (a positive-length good interval) plus the far
runner completing inside it — not the descent floor's `round + minimal kick` (S105). Same "separated far
element" regime (`d ≥ 91·V` here vs `d ≥ 13·V` for descent; both are "far enough"), two independent Lean
proofs. Having both formalized certifies the separated-far-element case of LRC(14) from two directions.

The Lean corpus now holds all three elementary/analytic branches of the LRC(14) dispatch:

> **`LRC(14)` ⟸ `CoveringCase`** (non-covering = sieve, S106) — PROVED
> **compact/AP-core:** `ap_core_bridge` (ρ≥13 + LRC(≤13), S105) — PROVED
> **separated far element:** `density_far_bridge` (d ≥ 91V + LRC(≤13), S107) — PROVED

Only the single open inverse theorem `INV` (= LRC(14) covering crux = Tao n=12) and the analytic `M`-split
remain unformalized; every constructive/geometric far-element and dispatch bridge is now kernel-checked.

## Net

- **Delivered:** `LRCDensityDischarge.lean`, 3 kernel-pure theorems (built, registered, no sorry/axiom).
  The density-route Φ>0 discharge for a separated far element, formalized via interval-completion.
- **Records:** the geometric density argument in the kernel — the third and last of the elementary
  far-element/dispatch bridges of LRC(14). Only `INV` remains.
- **Honest:** this formalizes the density route's *geometric core* (separated far element), not the sharp
  analytic `κ'R_G/w` Fourier bound (S100) — the elementary threshold is `d ≥ 91V`, which suffices for the
  separated regime the density route targets; the sharper constant needs the Fourier machinery (harder).

Cross-links:
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]],
[[o-R-is-false-S-is-theta-R-but-the-explicit-O-R-bound-closes-the-density-route-for-separated-far-elements-boxeph-S100]],
THM-1008/1010 (descent floor), LRC13Citation (`LRCUpTo13`).
