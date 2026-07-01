---
id: HYP-3810
title: THE T-JOIN BOUNDARY PARITY OBSTRUCTS LOW-DIMENSIONAL COVERS OF THE SC CLASSES -- the SC classes carry the flip-rank covering-excess, and it is their odd T-join-boundary structure (odd blue-fiber in the linear blue subspace) that makes them the bottleneck. Pursuing the HYP-3809 conjecture menu. The blue (grid-sym, d=m) lines form a graph on merged nodes with every SC node ODD-degree => SC nodes = the boundary T of a T-JOIN (|T|=#SC must be EVEN, the handshake). VERIFIED (n=4,5,6): (A) the grid-sym (blue) tilings form a LINEAR SUBSPACE W (sigma-fixed; sigma is a linear tile-coord involution), dim=(m+floor((n-1)/2))/2 = 2,4,6; each SC class has an ODD blue-fiber |C∩W|, every NS class has 0 (W = disjoint union of the SC classes' odd blue-fibers). (B) FLIP-RANK by class-subset: rho_all=2,4,7; rho_SC=1,4,6 (own log2-bounds 1,3,4 => EXCESS 0,1,2); rho_NS=0,1,6 (bounds 0,1,5 => excess 0,0,1). At n=5 rho_SC=4=rho_all (covering the SC classes ALONE is as hard as covering everything) with excess 1 while NS has excess 0; at n=6 SC has the LARGEST excess (2 vs NS 1). So the covering-excess (flip-rank above the information floor, HYP-3803/3804) is CONCENTRATED ON THE SC CLASSES. (C) So YES: the T-join boundary parity obstructs low-dim covers of the SC classes -- they are the odd-clustered boundary of the blue T-join in the linear subspace W, they cannot be packed into a rainbow subcube at their information floor (rho_SC > ceil(log2 #SC)), and they carry the covering-excess. The parity forces #SC even (handshake); the quantitative excess aligns with the high-|Aut| SC classes (opus-S15) -- the odd T-join boundary is the STRUCTURAL marker of the covering bottleneck
status: MIXED (exhaustive n=4,5,6 + structural interpretation). VERIFIED exhaustive: blue W linear (sample-checked closed under XOR + contains 0), dim (m+fix)/2 = 2,4,6; SC blue-fiber odd, NS zero; #SC even; rho_all=2,4,7, rho_SC=1,4,6 (excess 0,1,2), rho_NS=0,1,6 (excess 0,0,1). So SC carries the largest covering-excess (the low-dim-cover bottleneck). HONEST: the T-join parity is a PROVED structural feature (odd blue-fiber = odd boundary, #SC even = |T| even); that it OBSTRUCTS is EMPIRICAL (SC has the excess) + interpretive -- the parity per se yields the handshake, not a tight excess lower bound (the quantitative excess is the symmetry-folding/|Aut| mechanism of opus-S15, concentrated on SC). n>=7 not computed. A confirmed structure + a located bottleneck, not a closed lower-bound proof.
source: klein-2026-07-01-S77
depends_on:
  - HYP-3809   # the all-layers Lucas parity law + conjecture menu (this pursues the T-join/SC-cover items)
  - HYP-3803   # flip-rank / covering-excess (the covering side)
related:
  - HYP-3804   # packing/covering asymmetry (rainbow vs flip-rank); SC = the hard-to-pack side
  - HYP-3805   # opus-S15: flip-rank excess = symmetry-folding (high-|Aut| = Paley/SC classes) -- the quantitative mechanism
  - HYP-2250   # metagraph blue-black boundary parity (prior)
  - HYP-1772   # bucket parity (SC odd / NS even) -- the T-join boundary is the odd side
  - HYP-3808   # d=m blue/black lines (the T-join itself)
external: T-joins (odd-boundary edge sets, |T| even); linear codes / affine covers; Lucas
results:
  - 04-computation/tjoin_sc_cover_obstruction_klein.py
  - 05-knowledge/results/tjoin_sc_cover_obstruction_klein.out
---

# HYP-3810 — the T-join boundary parity obstructs low-dim covers of the SC classes

## The T-join boundary (confirmed structure)
The blue (grid-symmetric, `d=m`) lines form a graph on merged nodes in which **every SC node has odd
blue-degree**, so the SC nodes are the **boundary `T` of a T-join** (and `|T| = #SC` must be even — the
handshake, verified `#SC = 2,8,12` for `n=4,5,6`). Underneath: the grid-symmetric (blue) tilings are the
fixed set of `sigma` (a linear involution on tile-coordinates), so they form a **linear subspace `W`**
(verified closed under XOR, contains `0`), of dimension `(m + floor((n-1)/2))/2 = 2,4,6`. And:
> `W` = the disjoint union of the SC classes' blue-fibers, each `|C ∩ W|` **ODD**; every NS class has
> `|C ∩ W| = 0`.
So the SC classes are exactly the **odd-sized clusters** of the linear blue subspace — the odd T-join
boundary made concrete.

## The SC classes carry the flip-rank covering-excess (the bottleneck)
Flip-rank restricted to a class-subset (min subcube whose completions cover all classes in the subset),
with each subset's own information floor `ceil(log2 |subset|)`:
| `n` | `rho_all` (exc) | `rho_SC` (exc) | `rho_NS` (exc) |
|---|---|---|---|
| 4 | 2 (0) | 1 (0) | 0 (0) |
| 5 | 4 (0) | **4 (1)** | 1 (0) |
| 6 | 7 (1) | **6 (2)** | 6 (1) |
- At `n=5`: **`rho_SC = 4 = rho_all`** — covering the SC classes *alone* is as hard as covering everything —
  and the SC subset has excess `1` while NS has excess `0`. The entire covering-excess is on the SC side.
- At `n=6`: the SC subset has the **largest excess (`2`)**, more than NS (`1`); `rho_all = 7` needs both.
So the covering-excess (the flip-rank rising above the information floor, HYP-3803/3804) is **concentrated on
the SC classes** — precisely the odd T-join boundary.

## The answer: yes, the T-join parity obstructs low-dim SC covers
The SC classes are the **odd-clustered boundary** of the blue T-join in the linear subspace `W`. They cannot
be packed into a rainbow subcube at their information floor (`rho_SC > ceil(log2 #SC)` for `n>=5`), and they
carry the covering-excess. The **odd blue-fiber parity** is the structural fingerprint of this bottleneck:
an SC class, being an odd cluster of `W`, resists being isolated by a low-dimensional affine cover. The
parity yields the handshake `#SC` even; the *size* of the excess is the symmetry-folding / high-`|Aut|`
mechanism (opus-S15: the Paley/SC classes have fewest labeled reps, hardest to catch in a thin subcube).
The two threads meet: **the odd T-join boundary and the high-symmetry covering bottleneck are the same SC
classes.**

## New conjectures (from this)
- **C-SC-excess**: the covering-excess `rho_all - ceil(log2 A000568(n))` equals the SC-subset excess for
  small `n` (SC is the sole/dominant obstruction); test `n=7` (opus-S15: excess from the `|Aut|=21` SC
  heptagon).
- **C-parity-lower-bound**: is there a genuine PARITY (not just information) lower bound `rho_SC >=
  log2(#SC) + f(parity structure of W)`? The odd-cluster partition of the linear `W` may force excess via a
  Fourier/character argument on `GF(2)^w`.
- **C-blue-subspace-rank**: the SC classes are separated within `W` (dim `w=(m+fix)/2`); the min # of
  GF(2)-linear coordinates on `W` separating the SC classes is a "blue flip-rank" worth computing.

## Net
The SC classes are the odd boundary `T` of the blue T-join, sitting as odd-sized clusters in the linear
blue subspace `W`; they carry the flip-rank covering-excess (`+1` at `n=5` where `rho_SC = rho_all`, `+2` at
`n=6`), so the T-join boundary parity is the structural marker of the low-dimensional-cover obstruction. The
parity gives the handshake (`#SC` even); the excess magnitude is the high-`|Aut|` symmetry-folding — and it
lands on the same SC classes. Confirmed structure + located bottleneck; a tight parity lower bound on
`rho_SC` is the open next step.
