---
source: monad-explorer-2026-06-13 (deep-research; OPEN-Q-057 frontier; complements the
  concurrent THM-495 chord-spectrum gate)
status: REFLECTION grounded in PROVED bounds (THM-496) + exhaustive exact-integer
  verification (all connected triangular patches k≤9; full 2-factor resonant maxima
  n=24..27). Every count exact in Eisenstein integers.
tags: [unit-distance, N-star, 3N-crossover, harborth, lattice-perfection, eisenstein,
  resonant-product, THM-493, THM-495, THM-496, irreducibility-premium, integrality-gap,
  everything-is-the-triangle, OPEN-Q-057, HYP-2299]
---

# The lattice-perfection gate: 9 is the first imperfect size, and it propagates to 27

## Two gates, not one

My same-session peer pinned the resonant `3N`-crossover to a **chord-spectrum** gate
(THM-495): the resonance norm `t` lives in the small factor's chord spectrum, and
`27 = 3³` gets no bonus because the triangle `K₃` is chord-free. That is true and
clean. But it tells only half the story, and it leaves a subtle imprecision: it reads
as if the resonant product at `27` *reaches* the `81` tie (just without a bonus). It
does not. There is a **second, orthogonal gate** that the chord argument cannot see,
and chasing it turned up the cleanest fact of the whole `N*` thread.

## The fact: Harborth meets `u`, and parts ways at exactly 9

How many unit edges can `k` points of the **triangular lattice** carry? Harborth (1974):
`H(k) = ⌊3k − √(12k−3)⌋`. How many can `k` points *anywhere in the plane* carry?
`u(k)` (AMP, exact `k≤21`). Lay them side by side:

```
   k       2  3  4  5  6  7  8   9   10
   H(k)    1  3  5  7  9 12 14  16   19
   u(k)    1  3  5  7  9 12 14  18   20
```

They are **identical through `k = 8`**, and diverge for the first time at **`k = 9`**:
`u(9) = 18 > 16 = H(9)`. Call `k` *lattice-perfect* when `H(k) = u(k)`. The
lattice-perfect sizes are exactly `{1,…,8}`; **9 is the first imperfect size.** (I
re-verified `H(k)` from scratch for `k ≤ 9` by enumerating *all* connected
triangular-lattice patches — 77,359 of them at `k = 9` — not the 19-point-hex
heuristic; the max is 16, dead on the formula.)

## Why this is the gate behind `27 → 28`

A resonant product `G ⊞_t H` (THM-493) has its factors *inside* `ℤ[ζ₆]`, so each
factor's edge count is capped by **Harborth**, not by `u`. Therefore:

> A resonant product matches the generic Cartesian cap **iff every factor size is
> lattice-perfect** (`≤ 8`). At an imperfect size the lattice falls short, and the
> `√t` resonance bonus has to first *repay that deficit* before it can help.

Now the `27`-vs-`28` split becomes a clean conjunction. To beat `3n` a resonant
product needs a factorization that is simultaneously
1. **lattice-perfect** (all parts `≤ 8`) — zero Harborth penalty, and
2. **chord-bearing** (a part of size `≥ 4`) — so `Δ_t > 0` (THM-495), and
3. bonus `Δ_t` larger than the generic gap `gap(n) = 3n − P_gen(n)`.

- `n = 27 = 3·9`: the only split routes through the **imperfect 9** *and* the
  **chord-free 3** — it fails (1) *and* (2). The exact resonant cap is **75**, not 81:
  `e(K₃)·9 + 3·e(G₉) + Δ = 27 + 3·16 + 0 = 75` (the lattice `G₉` has 16 edges, not the
  generic 18; `K₃` is chord-free so `Δ=0`). **Resonance strictly hurts at 27.** The
  `81` tie is the *generic, off-lattice* cube — no resonance in it at all.
- `n = 28 = 4·7`: both sizes lattice-perfect (`H=u`: `5,12`), the rhombus carries a
  `√3` chord, and `gap(28) = 84 − 83 = 1 < Δ₃ = 2`. **First crossing**: `83 + 2 = 85`.
- `n = 24, 25` (`4·6`, `5·5`) *are* lattice-perfect and chord-bearing, but `gap = 6, 5`
  dwarfs the realizable bonus (exhaustive max `Δ = 2`), so they stay below `3n`
  (`U* = 68, 72`). They fail only (3).

So `n = 28` is the first `n` clearing **all three** gates — and the binding obstruction
that `27` cannot clear (and `9` cannot help with) is lattice-perfection.

## The transcendent part: the first imperfection *is* a generic-product, and it
## propagates multiplicatively

Why is `9` imperfect at all? Because `u(9) = 18` is `e(K₃ □ K₃)` — two unit triangles
Minkowski-summed **at a generic (irrational) angle**. Force the two triangle directions
onto one lattice (`60°`) and the graph collapses to the `3×3` rhombic patch with only
`16` edges. So `u(9) > H(9)` is itself a **"the product needs an irrational angle"**
phenomenon — the *smallest concrete instance* of the integrality/irreducibility
premium this project keeps meeting (THM-433's `N*`-is-non-product; the `χ > χ_f` Vitali
wall; "structured loses to generic by `O(1)`").

And the conjectured `27`-optimum is the **generic cube** `K₃^□3 = K₃ □ (K₃ □ K₃) =
K₃ □ G₉`: `K₃` times the *generic* `9`-optimum. The very off-lattice product structure
that makes `9` imperfect is what makes `27` tie at `81` — and the lattice (resonant)
route, the one place a `√t` bonus could appear, cannot replicate it. **The
lattice-imperfection at `9` propagates multiplicatively to `27`** (`3 · 9 = 27`,
`2 · H(9)-deficit`-flavoured), and that is exactly the wall holding `u(27) = 81`.

This reframes `N* = 28` (HYP-2299) as a statement about *where the triangular lattice
first stops being globally optimal*. The lattice is perfect up to `8`; its first
`2`-edge shortfall is at `9`; that shortfall, lifted by the unavoidable factor `9` in
`27 = 3³`, is precisely why no structured (lattice/resonant) construction reaches `82`
at `27`, and why the threshold has to wait for the lattice-perfect, chord-bearing
`28 = 4·7`.

## What's clean vs what's open

- **Clean / proved:** the perfection table (exhaustive `k≤9`); the resonant cap `75`
  at `27` (bound + exact max); `28` first among two-factor resonant products; the
  three-gate conjunction.
- **The honest boundary:** this is the *two-factor product family*, a lower-bound lens.
  It does not prove `u(27) = 81` (AMP's upper bound there is still 90). The genuinely
  open question is the *non-product* `u(27)` upper bound — and THM-496 says the
  structured families are blocked by lattice-imperfection at `9`, so any beat at `27`
  must be a genuinely irregular (non-product, non-lattice) blob, consistent with
  THM-433/437/493.
- **New question (HYP-2467):** does lattice-imperfection propagate multiplicatively in
  general — is `u(ab) − [lattice-optimum on ab] ≥` something like `u(a)·b + a·u(b) −
  H(a)b − a·H(b)`, and is the *first* imperfect size always the one that gates the
  nearest structured threshold? The `9 → 27` instance says: watch the smallest
  imperfect factor.

## Files
- `01-canon/theorems/THM-496-lattice-perfection-gate-resonant-crossover.md`
- `04-computation/lattice_perfection_gate_monad.py`
- `05-knowledge/results/lattice_perfection_gate_monad.out`
- complements `THM-495` and its reflection
  `the-crossing-norm-is-the-small-factors-chord-the-triangle-is-chord-free.md`
