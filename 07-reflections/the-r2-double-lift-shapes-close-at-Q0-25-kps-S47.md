# The r=2 double-lift shapes close at Q₀ = 25 — a bounded, height-uniform covering

*kps-2026-07-06-S47 — taking the r=2 (double 13-lift) shapes of the (G) residual,
mapping their coverings, and formalizing the q≤12 backbone.*

## The task

S46 laid out the proof path: the (G) residual is the AP `{1,…,12}` with `r` speeds
13-lifted; `r=0` is the AP, `r=1` (single lift) is `d=1` (mac-mini THM-633, GREEN),
and `r≥2` is the remaining program — for each lift-shape, a fixed finite covering
clears every non-AP member, height-uniform. This session does `r=2`.

## The result: Q₀ = 25 for every r=2 shape

For each of the `C(12,2) = 66` lifted pairs `(i,j)` — the AP with speeds `i,j` sent
to `i+13a, j+13b` — over all lift heights `a,b ∈ [0,25]`, **every non-AP member
clears at a modulus `q ≤ 25`**. A single fixed covering `{q ∈ [6,25]}` handles all 66
shapes and all heights (`lrc_r2_double_lift_shapes_kps_S47.out`). The covering splits
by transversality:

- **non-transversal** members (some `±`-pair missed mod 25) clear at **`q = 25`** —
  case 1, kps `LRCMod25Floor` / mac-mini THM-634 (GREEN);
- **transversal** members (blockers) clear at **`q ≤ 24`** — the small-`q` covering
  (`loose_of_no_multiple` for `q ≤ 12`; avoid-band `rational_point_margin` certs for
  `13 ≤ q ≤ 24`);
- the **AP** clears at no modulus — the tight-locus exception.

The earlier "max modulus 37" (S47 first pass) was an artifact of *excluding* `q=25`:
the non-transversal members, denied their `q=25` witness, cleared next at 37. With
`q=25` in the covering, `Q₀ = 25`. The clarifying example: `{1,2,3,5,7,8,9,10,11,12,17,19}`
(shape `(4,6)`, `a=b=1`) has `M = 2/25` *exactly* — not in the open gap — and misses
the pair `{4,21}` mod 25, so it is cleared by `LRCMod25Floor` at `q=25`.

**Height-uniformity is intrinsic:** the covering `{q ≤ 25}` is fixed, independent of
the lift heights `a,b`; clearing at `q` depends only on `{v_i mod q}`. A lift of any
size is inert at some covering modulus.

## Structure: the hard shapes lift speeds 6 and 12

Among the 66 shapes, those lifting speed `6` or `12` are the "hardest" (their
transversal members need moduli up toward the high teens/low twenties before
`q=25`); the others clear at smaller `q`. This is the divisibility signature —
`6 = 2·3` and `12 = 2²·3` are the composite small speeds, whose lifts stay
divisible-adjacent longest. It localizes the residual work to a few shapes.

## Lean: the q≤12 covering backbone (GREEN)

`LRCSmallModFloor.loose_of_no_multiple` (added this session, kernel-pure): for any
`q` with `2q ≤ 25` (i.e. `q ≤ 12`), if no speed is divisible by `q`, then
`M ≥ 1/q ≥ 2/25` at `t = 1/q` — LOOSE. This is the general `q ∈ {7,…,12}` layer of
the covering (previously only `q=12` was stated), a direct `rational_point_margin`
instance at `μ=1`. Every transversal r=2 double-lift that misses a small multiple is
now closed by one certificate.

## What remains for r=2

- The `13 ≤ q ≤ 24` avoid-band certificates (for transversal members that carry all
  of `{7,…,12}` as multiples and clear only at a mid modulus). These are
  `rational_point_margin` certs at `s = q`, `μ = ⌈2q/25⌉` — the same atom, higher
  band. Formalizing the specific `(q, μ)` for the residual members closes the
  transversal r=2 case.
- Then `r=2` is GREEN, joining `r=0` (tight-locus) and `r=1` (THM-633). The `r≥3`
  shapes follow the same template (bounded covering, height-uniform), with fewer
  viable shapes each higher `r` (more lifts ⟹ more ways to miss a pair ⟹ more
  members fall to the mod-25 branch).

## Honest scope

- The `Q₀ = 25` bound is verified over `a,b ∈ [0,25]` (heights to ~325) for all 66
  shapes; the covering is fixed and residue-only, so the bound is height-uniform by
  the S44/S46 mechanism, but a proof (every member of every shape clears at `q ≤ 25`)
  is a finite residue check per shape, not yet fully in Lean.
- The `q ≤ 12` layer is proved and formalized (`loose_of_no_multiple`); the `q=25`
  layer is `LRCMod25Floor` (GREEN); the `13 ≤ q ≤ 24` layer is the remaining certs.

## Pointers

- `lrc_r2_double_lift_shapes_kps_S47.out` (66 shapes, Q₀=25, hard shapes = speeds 6,12,
  the boundary example), `LRCSmallModFloor.lean` (`loose_of_no_multiple`, GREEN).
- kps S46 (proof path), S44 (bounded modulus, `LRCSmallModFloor`), S41
  (`LRCMod25Floor`); mac-mini THM-633 (r=1), THM-634; klein S126 (discrete ladders);
  opus S126 (finite covering).
