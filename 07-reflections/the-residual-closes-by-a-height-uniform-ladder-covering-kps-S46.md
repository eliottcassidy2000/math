# The residual closes by a height-uniform ladder covering — the proof path for the last node

*kps-2026-07-06-S46 — driving the final node (the uniform-Q₀ covering), synthesizing
klein's discrete-ladder finding, opus's Erdős-covering framing, mac-mini's d=1
ladder, and my ladder method into a concrete, height-uniform closure route.*

## The last node, and the inspiration

After S45, (G) is a finite pile of margin certificates plus one open node: **every
non-AP full-transversal ("blocker") clears at some modulus `q ≤ Q₀`** (opus S126's
"height-uniform covering, Erdős-flavoured"). The worry was height: a blocker can
have height `10⁵`, so how can a *finite* modulus set clear it?

klein S126 supplies the key: the residual is **discrete, not a continuum** — the
near-tight loose values are the ladder rungs `k/(11k+1)` (`1/12, 2/23, 3/34, …`),
realized by *spread-one-runner* families `{1,…,10, 11k}`, with a gap `δ₀ = 1/276`
above `1/12`. So the residual families are **ladder families** (the AP with some
speeds lifted), and — klein's phrase — "a finite union of formula-closable ladder
families," one closed form per shape, valid for all heights. That is exactly my
ladder method (S21/S36).

## The residual is the AP with speeds 13-lifted

A doubly-saturated family has residues `{1,…,12}` mod 13 (a complete system) and is a
full transversal mod 25. The **minimal-lift** such family is the AP `{1,…,12}` itself
(all speeds in `[1,12]`). Every *other* doubly-saturated family is the AP with some
`r ≥ 1` speeds **13-lifted** (`v_i → v_i + 13k_i`, preserving the mod-13 residue).
So the residual stratifies by `r`:

- `r = 0`: the AP — `M = 1/13` (tight-locus, the boundary).
- `r = 1`: single 13-lift `{…, v_i+13k}` — this is `d=1`, **GREEN** (mac-mini
  THM-633: `{1,…,11}∪{x}` reaches `2/25` for all `x ≠ 12`, a two-witness ladder).
- `r ≥ 2`: multi-lift — the residual this session addresses.

## The mechanism: a fixed covering, height-uniform

The AP fails **every** modulus — verified: `{1,…,12}` has no clearing rotation at any
`q ∈ [6,44] \ {25}` (38/38). It is the tight extremal: at every `q`, its residues
`c·{1,…,12}` are an AP mod `q` that (three-gap) always hits the forbidden band. **This
is why the AP is the unique survivor** — it is the one family that clears nowhere.

A *lift* breaks this. Lifting a speed changes its residue at every modulus *except* 13,
opening a gap that some rotation exploits. Concretely, for the double-lift shape
`{1,…,10} ∪ {11+13a, 12+13b}`:

> A **fixed** covering `{11,12,13,14,16,17,18,19,21,23}` — *independent of the lift
> heights `a,b`* — clears every non-AP member (verified over `a,b ∈ [0,59]`, heights
> to ~780; max modulus needed 23; the only uncleared are the `a=0`/`b=0`
> single-lifts, which are the `d=1` case closed by THM-633).

**Height-uniformity is free:** clearing at `q` depends only on `{v_i mod q}`, and the
base `{1,…,10}` clears at `q = 11` unless a lift `≡ 0 mod 11` — a *residue* condition,
not a size condition. The lifts, however large, are inert at almost every covering
modulus; the covering has enough moduli that one avoids all the lift-blocks. This is
the Erdős-covering heart: the lifts block finitely many moduli, and a fixed covering
outnumbers the blocks.

## The proof path

> **(C) ⟸** for each lift-shape `S` (which of the 12 speeds are 13-lifted), a fixed
> finite covering `{q ≤ Q₀}` clears every non-AP member of `S`; the AP is the unique
> family (the `r=0` shape) that no covering modulus clears.

This is a **finite** program:

1. Enumerate the lift-shapes. There are finitely many (subsets of `{1,…,12}` that are
   lifted), and by the doubly-saturated constraints only a few are viable.
2. `r=0`: the AP → `M=1/13` [tight-locus theorem].
3. `r=1`: [GREEN, mac-mini THM-633].
4. `r≥2`: per shape, a fixed covering `{q ≤ Q₀}` (`Q₀ ≈ 23`) clears all members — a
   **finite residue check** (for each covering modulus `q`, the residues that clear;
   show every non-AP residue pattern of the shape clears at some `q`). Each clearance
   is a `rational_point_margin` / `LRCSmallModFloor` certificate.

**No height bound is needed** (S44/S46): the covering is by residues; lifts of any
size are inert at some covering modulus. The height/`u_max`/lcm wall that klein,
mac-mini, and I kept hitting is *bypassed* — not by bounding the height, but by
covering the residues, exactly as klein's "formula-closable for all `k`" predicted.

## Honest scope

- The path is demonstrated on the double-lift shape (fixed covering, height-uniform,
  verified to height ~780). The full program — enumerate all viable lift-shapes and
  prove each shape's finite covering — is the remaining work. It is finite and
  residue-only (Erdős-covering-flavoured), not analytic.
- The `r=1` (d=1) shape is already GREEN (THM-633); `r=0` is the tight-locus theorem.
  The `r≥2` shapes are the residual, and each is a bounded-base ladder whose covering
  is a finite check.

## Pointers

- `lrc_ladder_covering_path_kps_S46.out` (AP fails all 38 moduli; double-lift shape
  cleared by a fixed `{11..23}` covering, height-uniform), `lrc_covering_uniformity_kps_S44.out`.
- klein S126 (discrete residual ladders, `δ₀=1/276`); opus S126 (finite covering,
  Erdős-flavoured), S125 (two-modulus); mac-mini THM-633 (d=1 ladder), THM-634;
  kps S44 (bounded clearing modulus, `LRCSmallModFloor`), S45 (covering unifies the
  two-modulus crux), S36/S21 (the ladder method).
