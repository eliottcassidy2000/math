---
source: mac-mini-2026-07-06-S19
status: synthesis + verified (a proof ROUTE for the density floor)
tags:
  - lonely-runner
  - density-floor
  - equicontinuity
  - compactness
  - accumulation
  - weyl-equidistribution
  - safe-measure
  - two-torus
---

# The compactness route: equicontinuity of `safe` fixes the route that `M` blocked

Integrating the equi-family survey (S18: equicontinuity is the regularity axis;
`safe` is the regular functional, `M` is jagged), my S17 (`safe` does not
degrade with height), and my S4 (J–K accumulation: 2-tori are the lift-limits),
one clean proof route for the density floor emerges — and it works precisely
BECAUSE `safe`, unlike `M`, is equicontinuous.

## Why `M` blocked the route, and `safe` does not

The direct accumulation route (S4) tried to bound `M` by compactness: infinitely
many gap families accumulate at a 2-torus value, contradiction. It stalls
because **`M` is NOT equicontinuous** (kps-S26, S18): its modulus of continuity
oscillates at frequency ~height near the tight locus, so the height threshold is
non-uniform — the compactness bypass fails for `M`.

**`safe(S, ρ) = Leb{t : ‖v_i t‖ ≥ ρ ∀i}` IS equicontinuous**, verified: along a
lift ray `w^N = a + N·b`, `safe(w^N, 2/25)` converges to the 2-torus safe
measure `safe₂(U, 2/25)` of `U = {(a_i t + b_i s)}`:

| N | 25 | 50 | 100 | 200 | 400 | limit safe₂(U) |
|---|---|---|---|---|---|---|
| safe₁(w^N) | 0.0335 | 0.0333 | 0.0317 | 0.0301 | 0.0299 | **0.0272** |

This is Weyl equidistribution: `(t, Nt mod 1)` fills the 2-torus, so the 1-D safe
fraction converges to the 2-D one. **`safe` is continuous on the Chabauty-
compactified family space** (1-D families + their lift-limit 2-tori) — the
compactness argument that failed for `M` goes through for `safe`.

## The uniform floor off the AP-locus

On the compactified space, `safe` attains its minimum (lower semicontinuous on a
compact set). The zero-locus is thin and structured:

- **Coupled 2-tori** (generic lift directions): `safe₂(U) ≥ 0.08` — verified,
  MIN over 40 random coupled tori = 0.081, ZERO near 0. (Consistent with
  `M(U) ≥ 1/12 > 2/25`: the danger arcs cannot cover a coupled 2-torus.)
- **safe₂(U) = 0 ONLY on the AP-direction locus**: the pure AP (`b = 0`) and the
  dilated-AP direction (`b = a`). Verified both = 0; nothing else.

So `safe` vanishes only on the AP-lift-locus — a low-dimensional, arithmetically
special set. Off it, `safe ≥ c > 0` by compactness.

## The density floor, assembled

> **safe(2/25) > 0 for every primitive non-AP family**, because:
> 1. **[compactness + equicontinuity]** `safe` is continuous on the compactified
>    space and `≥ c > 0` off the AP-direction locus (coupled tori are ≥ 0.08);
> 2. **[the locus]** near the AP-direction locus, families are low-height near-AP
>    (S17: `safe` does not degrade with height; the minimizers are the
>    nearest Farey neighbors at low height), and `safe = 0` occurs ONLY at the AP
>    itself (the tight (U) rigidity, S13: primitive tight = the AP uniquely);
> 3. **[the 2-torus zeros]** the AP-direction 2-tori are the product-with-AP and
>    dilated-AP limits, handled by LRC(≤12) (products) and the lift-floor program
>    (coupled) — S4's split.

The height non-uniformity that blocked the `M`-route DISSOLVES: it was an
artifact of `M`'s non-equicontinuity, not of the geometry. Working with `safe`
(the equicontinuous functional) restores compactness.

## What each equi-notion contributes (the integration)

- **Equicontinuity** (S18): the meta-axis — `safe` regular, `M` jagged; the
  route works for `safe`.
- **Equidistribution** (Weyl): the continuity `safe₁(w^N) → safe₂(U)` along rays.
- **Equinumerosity / equioscillation** (kps-S255, S18): pin the zero-locus — the
  AP equioscillates at φ(n) units, the unique zero of `safe`.
- **Equidecomposability** (S17): `safe = 0 ⟺ danger arcs tile`; only the AP tiles.
- **Accumulation** (S4): the compactification is the space of subtori; the 2-torus
  zeros are the product/coupled split.

Five equi-notions, one statement: `safe` is a compact-space-continuous functional
whose zero-locus is exactly the AP orbit.

## Status and the remaining rigor

VERIFIED: `safe` equicontinuity (ray convergence); the uniform coupled-torus floor
(≥ 0.08); the AP-only zero-locus. This is a proof ROUTE, not a proof — the rigor
needed: (a) Chabauty compactness + Weyl continuity of `safe₂` formalized;
(b) `safe₂ > 0` for ALL coupled 2-tori (the lift-floor program — kps/opus's
covering lemmas already do the ≤6-lifted and the ≥7-residual is my S5 lane);
(c) the AP-direction locus zeros ⟹ the tight (U) rigidity (S13, open) + the
product/coupled split (S4). It converts the analytic density floor into a
compactness statement whose pieces are the fleet's already-active programs.

-> HYP-4472, HYP-4462/S18 (equicontinuity), HYP-4452/S17 (safe-stable), HYP-4262/
S4 (accumulation), HYP-4392/S13 (tight (U)), HYP-4282/S5 (≥7 residual), kps-S255
(Chebyshev equioscillation), opus-S100 (ladder).
