# The observer abstraction: the LRC is a marked ORIGIN (the deleted speed-0 runner = the empty tooth = the basepoint) escaping the harmonic blockings 1/q where the covering's runners return to it; the covering-min is the observer's Farey escape strictly between its TWO TIGHTEST blockings 1/n and 1/(n−1) — the construction = the continued fraction [0; n−1, n] (the two blockings NESTED), the drop-2 = the Farey mediant 2/(2n−1), and the n=7 transition is pure REALIZABILITY (which escape a covering set can attain); the margin ~1/n² is the observer clearing 1/n by a single Farey hair

*opus-2026-06-30. Owner: "when i typed god i meant grid" — the triangular GRID — and "work to see things
abstractly in terms of observers like this." So: develop the observer. It is the right lens, and it turns the
covering-min into one clean continued-fraction statement about a marked origin escaping its blockings.*

## The observer (the abstraction)
An **observer** is a marked point — the origin, the basepoint, the deleted element, the frame that watches.
The structure is seen RELATIVE to it; the question is about its LOCAL view (is it isolated? what's its gap?).
In the LRC the observer is the **origin** — the *deleted speed-0 runner*, the **empty tooth**, the point on
the triangular grid that was pulled out (the visual's coral ring). A runner of speed `s` **returns to the
observer** (lands on it, position `≡0`) at every `t = k/s`. So:
> **The observer is BLOCKED at `t = k/s` for each runner `s`.** A COVERING set has a multiple of every
> `q∈{2,…,n}`, so it blocks the observer at **every harmonic point `t = 1/q`** (q=2..n). The lonely time is
> the observer's **ESCAPE** — a `t` where every runner is far (the observer sits in a gap of radius `≥1/n`).

## The covering-min IS the observer's escape between its two tightest blockings (computed)
The two tightest blockings are `1/n` and `1/(n−1)` (largest `q`). **Both covering-min families live STRICTLY
between them** (verified n=4..14):
> `1/n  <  2/(2n−1)  <  n/Φ₆(n)  <  1/(n−1)` — e.g. n=14: `0.07143 < 0.07407 < 0.07650 < 0.07692`.
And the two families are the two canonical "between" points of the observer's blockings:
| family | escape | observer reading |
|---|---|---|
| **construction** (n≥7) | `n/Φ₆(n) = [0; n−1, n]` | the **continued fraction whose partial quotients ARE the two tightest blockings** `n−1, n` — *the escape is the blockings, nested* |
| **drop-2** (n≤6) | `2/(2n−1)` | the **Farey mediant** of `1/(n−1)` and `1/n` — the simplest point between |
> The covering-min is **which escape a covering set can REALIZE**: the mediant for `n≤6`, the convergent for
> `n≥7`. **The n=7 transition is a pure realizability question about the observer** — not geometry, not
> number theory, just: can the runners be arranged so the observer escapes at the mediant? (Yes for small
> `n`, no for large.) The margin `escape − 1/n = (n−1)/(nΦ₆) ~ 1/n²` is the observer **clearing `1/n` by a
> single Farey hair** — the razor-thinness is the gap between consecutive Farey escapes.

## Why this is the right lens (what the observer buys)
- **The forcing step (the one open node) is reframed:** the covering constraint blocks the observer at the
  harmonic points; the question "does the covering force the ζ₆-line?" becomes "is the observer's smallest
  REALIZABLE escape the continued fraction `[0; n−1, n]`?" The ζ₆-line is the realization of that convergent
  (its binding modulus is `Φ₆`, the convergent's denominator); the drop-2 mediant (modulus `2n−1`) is the
  only alternative, unrealizable for `n≥7`. **The whole open node is: the convergent beats the mediant for
  realizability when `n≥7`.**
- **The empty tooth = the observer = the deleted origin** — my Dirac-comb, klein's antipodal `(1,−1)`, the
  visual's coral ring, all the same marked point. The "lonely set" is exactly *where the observer is
  unblocked*.
- **The continued fraction `[0; n−1, n]` is the cleanest statement of the covering-min** — no `Φ₆`, no
  Eisenstein, no hexagon needed to STATE it: the construction's escape is just the two tightest blockings
  written as a continued fraction. (The hexagon/`Φ₆`/Kershner is the *why-optimal*; the CF is the *what*.)

## The observer across the project (the generalization)
The marked-origin lens recurs — every core object is a *pointed* structure seen from a basepoint:
| object | the observer (basepoint) | its local view |
|---|---|---|
| **LRC** | the origin (deleted speed-0 runner) | escape from harmonic blockings (this) |
| **Dirac comb Ш** | the missing tooth | the gap in the comb |
| **covering-min** | the deleted origin on the ζ₆-grid | the `2n` 0-gap (tridiagonalized) |
| **tournament** | a marked vertex | who it beats / its odd cycles (OCF) |
| **metagraph** | the transitive class `H=1` (the identity) | the H-gradient from the basepoint |
| **2-adic descent** | the origin through the levels | the renormalization zoom |
> In each, the invariant is the observer's LOCAL view (isolation, gap, escape, parity). This is the
> "pointed/based" perspective — reduced/relative invariants taken at the basepoint. The observer is the `0`
> / the identity / the frame; the structure is its orbit; the question is whether the observer is *lonely*.

## Status
- **Computed/verified (opus):** the covering-min lies strictly between the observer's tightest blockings
  `1/n, 1/(n−1)` (n=4..14); construction `= [0; n−1, n]` (the blockings nested as a CF); drop-2 `=` the
  Farey mediant `2/(2n−1)`; transition = realizability; margin `~1/n²` = one Farey hair.
- **The reframing:** the LRC is the marked origin escaping its harmonic blockings; the open forcing node is
  "the realizable escape is the convergent `[0;n−1,n]` for `n≥7`."
- **The lens:** observer = deleted origin = empty tooth = basepoint, unifying LRC / comb / covering-min /
  tournament / metagraph / descent as pointed structures read from their `0`.

Related: the-cyclic-kershner-attack (the grid/empty-tooth visual), the-covering-min-as-a-function-of-n (the
transition), the-covering-min-three-distance-gaps (the {1,n,2n}), my Dirac-comb/empty-tooth reflections;
klein HYP-3715 (ζ₆-line, antipodal), mac-mini HYP-3702 (taxonomy); OPEN-Q-108.
