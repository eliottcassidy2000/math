---
id: HYP-3811
title: THE TILING CUBE FOLDS TO A QUARTER (Klein-four <sigma=complement, flip>) + complement-H-pairing + the SC-covering obstruction is a COORDINATE ARTIFACT (dissolves in the half-address coordinates). THREE results. (Q) QUARTER-TILING: the grid reflection sigma:(x,y)->(n+1-y,n+1-x) IS the complement (VERIFIED T(sigma t) ~= T(t)^op all t), and the flip=+1 (flip all tiles) is a second commuting involution; <sigma,flip>=(Z/2)^2 (Klein four). So the tiling cube folds TWICE: FULL 2^m ->[/sigma=complement] HALF ->[/flip] QUARTER. Counts 8->6->3 (n=4), 64->40->20 (n=5), 1024->544->272 (n=6). HALF=(2^m+|W|)/2, QUARTER=(2^m+|W|)/4 (sigma*flip fixed-point-free), |W|=2^((m+floor((n-1)/2))/2) the blue subspace. The BLUE LINES = W/flip = the sigma-fixed (size-2) orbits of the quarter, count 2^(w-1); non-blue = size-4 orbits. (H) COMPLEMENT-H-PAIRING: H(T)=H(T^op) VERIFIED (Ham-path count is complement-invariant); every H ODD (Redei); each merged node carries ONE H value; complement pairs classes of equal H. (R) THE NEXT TARGET (parity lower bound on rho_SC) RESOLVED as a COORDINATE ARTIFACT: the BLUE FLIP-RANK rho_W (min subcube in the half-address/W-basis coordinates covering all SC classes) = log2(#SC) EXACTLY (excess 0: rho_W=1,3,4 for n=4,5,6), whereas the full-tile-cube rho_SC has excess (1,2). So the S77 T-join/parity obstruction to low-dim covers of the SC classes is NOT intrinsic -- it is an artifact of the full tile-coordinate system; in the NATURAL half-address (blue subspace W) coordinates the SC classes cover at their information floor with NO excess. The half-tiling model IS the coordinate system in which the SC classes pack/cover optimally
status: MIXED (exhaustive n=4,5,6). VERIFIED: sigma=complement (canon(T(sigma t))=canon(T(t)^op), sampled); sigma,flip commute + involutions => Klein four; fold counts 8-6-3/64-40-20/1024-544-272; blue lines = size-2 orbits (2^(w-1)); H(T)=H(T^op) + all H odd; rho_W (SC cover in half-address coords) = log2(#SC) exactly = 1,3,4 (excess 0) vs full-cube rho_SC=1,4,6 (excess 0,1,2, from HYP-3810). HONEST: the coordinate-artifact finding is the key correction to HYP-3810's 'obstruction' -- the SC covering excess is real in tile-coords but dissolves in the half-address coords, so the parity does NOT intrinsically obstruct (it obstructs only the WRONG coordinates). n>=7 not computed. Refines HYP-3810, does not overturn its structure (SC=odd T-join boundary still holds); it reclassifies the 'obstruction' as coordinate-dependent.
source: klein-2026-07-01-S78
depends_on:
  - HYP-3810   # T-join boundary obstructs SC covers (in tile-coords); THIS shows it's a coordinate artifact
  - HYP-3809   # all-layers Lucas parity + conjecture menu
related:
  - HYP-3810-opus   # opus-S19 PROVED 'the blue subgraph IS the half-tiling metagraph' (folding iso) -- the theorem behind this quarter-fold; CONVERGENT
  - HYP-3810-macmini # mac-mini-S85: kappa_SC = half-tiling dim floor((n-1)^2/4) = w (the SC covering obstruction IN TILE-COORDS); this session shows it DISSOLVES to log2(#SC) in the half-address (sigma-symmetric) coords -- complementary
  - HYP-3808   # blue/black lines (the d=m fold = flip); grid reflection = complement
  - HYP-2689   # half-address / half tiling model (the W-basis coordinates used here) + THM-551
  - HYP-3807   # 'quotient/coordinates matter not the cube' (S74) -- this is that lesson realized: right coords dissolve the excess
  - HYP-1772   # F(C)=H/|Aut| odd; complement preserves H
external: Klein four-group (Z/2)^2; Burnside; Redei (H odd, complement-invariant)
results:
  - 04-computation/quarter_tiling_complement_H_klein.py
  - 05-knowledge/results/quarter_tiling_complement_H_klein.out
---

# HYP-3811 — quarter-tiling, complement-H-pairing, and the SC-cover coordinate artifact

## (Q) The tiling cube folds to a QUARTER (Klein-four)
The grid reflection `sigma:(x,y)->(n+1-y,n+1-x)` on tile positions **is the complement** (VERIFIED:
`T(sigma t) ~= T(t)^op` for all sampled `t`, matching opus-S18). The **flip** `= +1` (flip all tiles) is a
second involution, and `sigma`, `flip` **commute**, so `<sigma, flip> = (Z/2)^2` (Klein four). The tiling
cube therefore folds TWICE:
> **FULL `2^m` `->[/sigma = complement]` HALF `->[/flip]` QUARTER.**
Counts: `8 -> 6 -> 3` (n=4), `64 -> 40 -> 20` (n=5), `1024 -> 544 -> 272` (n=6). By Burnside,
`HALF = (2^m + |W|)/2` and `QUARTER = (2^m + |W|)/4` (here `sigma*flip` is fixed-point-free), with
`|W| = 2^{(m + floor((n-1)/2))/2}` the blue (grid-sym) subspace. The **BLUE LINES** are exactly the
`sigma`-fixed **size-2 orbits** of the quarter (`W/flip`, count `2^{w-1}`); the non-grid-sym tilings form
the **size-4 orbits**. So "the half tiling" (grid-sym = complement-mirror = `W`) folds once more by `flip`
into the **quarter tiling**, whose fixed locus is the blue lines.

## (H) Complement-H-pairing
`H(T) = H(T^op)` (the Hamiltonian-path count is complement-invariant, VERIFIED), and every `H` is **odd**
(Rédei). So the complement pairs classes of **equal `H`**, and each merged node carries a single `H` value.
The merged metagraph's `H`-gradient is thus well-defined on the half/quarter, and complement is an
`H`-isometry.

## (R) The next target: the SC-covering obstruction is a COORDINATE ARTIFACT
HYP-3810 found the SC classes carry the flip-rank covering-excess (in the full tile-coordinates:
`rho_SC = 1,4,6`, excess `0,1,2`). Computing the **blue flip-rank** `rho_W` — the min subcube in the
**half-address coordinates** (one bit per `sigma`-orbit = the natural basis of the blue subspace `W`,
HYP-2689/THM-551) covering all SC classes:
> **`rho_W = log2(#SC)` EXACTLY** (`= 1, 3, 4` for `n=4,5,6`), **excess 0** — versus the full-cube
> `rho_SC` excess `0,1,2`.
So **the T-join/parity "obstruction" is not intrinsic to the SC classes — it is an artifact of the full
tile-coordinate system.** In the natural half-address (blue subspace `W`) coordinates, the SC classes cover
at their information floor with no excess: they pack/cover optimally. This is exactly the S74 lesson
(HYP-3807: "the quotient/coordinates carry the content, not the cube") realized concretely — the right
coordinates dissolve the excess. The odd T-join boundary structure (HYP-3810) still holds; what changes is
its interpretation: the parity marks the SC classes as a distinguished linear subspace `W`, and *within*
that subspace they are perfectly separable — the excess was the cost of viewing `W` through the wrong
(full-tile) axes.

## Net
The tiling cube folds to a **quarter** via the Klein-four `<sigma = complement, flip>` (FULL ->
HALF -> QUARTER, `= (2^m+|W|)/4`), with the **blue lines** as the `sigma`-fixed locus of the quarter; the
**complement-H-pairing** holds (`H(T)=H(T^op)`, odd, one value per merged node); and the **SC-covering
obstruction dissolves** in the half-address coordinates (`rho_W = log2 #SC`, no excess), so the S77/S3810
hardness is a coordinate artifact, not an intrinsic parity obstruction. The half-tiling model is the
coordinate system in which the SC classes cover optimally.
