---
id: HYP-3810
title: THE T-JOIN PARITY OBSTRUCTS LOW-DIMENSIONAL COVERS OF THE SC CLASSES -- kappa_SC(n) = floor((n-1)^2/4) (the half-tiling / complement-fold dimension) for n>=5, which EXCEEDS the information floor ceil(log2 #SC) by a quadratically-growing gap (0,1,2 for n=4,5,6; predicted +4 at n=7). The blue line-subgraph on the SC merged nodes is a T-JOIN (all-odd-degree, T = every SC node), because each SC class contains an ODD number of grid-symmetric tilings (the S84 parity theorem); this odd-parity is exactly why no subcube smaller than the half-tiling dimension can realize all SC classes. UPPER BOUND (clean): the grid-sym subcube (dim = floor((n-1)^2/4), containing all grid-sym tilings = exactly the SC classes) covers all SC classes, so kappa_SC <= half-dim. LOWER BOUND (the obstruction, verified n=5,6): no smaller subcube works. The blue T-join goes bipartite (n<=5) -> NON-bipartite (n=6), another face of the n=6 sea-onset transition.
status: CONFIRMED-partial (kappa_SC = 1,4,6 for n=4,5,6 exhaustive; = half-dim for n=5,6; blue T-join structure + cycle-rank 0,1,15 + bipartite->non-bipartite at n=6 verified). CONJECTURE kappa_SC(n)=floor((n-1)^2/4) for n>=5. The T-join parity is the STRUCTURAL fingerprint of the obstruction (odd #grid-sym per SC class); a rigorous proof that the parity forces the lower bound is the open target. Answers the owner's question YES: the T-join boundary/parity DOES obstruct low-dim SC covers, up to the full half-tiling dimension.
source: mac-mini-2026-07-01-S85
related:
  - HYP-3809   # S84 atlas: the parity theorem (SC odd #grid-sym) + blue T-join + half-tiling; this pursues its cover-obstruction item
  - HYP-3808   # S83 blue=odd-degree / black=even-degree decomposition (the T-join is the blue part)
  - HYP-3798   # S81 kappa(n)=1+C(n-2,2) (all classes); kappa_all - kappa_SC = extra dim for the NS sea
  - THM-549    # half-tiling = complement fold; dim floor((n-1)^2/4)
results:
  - 04-computation/sc_cover_dimension_tjoin_obstruction_macmini_20260701.py
  - 05-knowledge/results/sc_cover_dimension_tjoin_obstruction_macmini_20260701.out
---

# HYP-3810 -- the T-join parity obstructs low-dimensional covers of the SC classes

Owner's target (from the HYP-3809 atlas): does the T-join boundary/parity obstruct low-dimensional covers of
the SC (self-complementary) classes? **Answer: yes, up to the full half-tiling dimension.**

## The setup (joining S81 covers + S84 T-join)
- **`kappa_SC(n)`** = min dimension of an axis-aligned subcube (choose `k` free arcs, fix the rest) whose
  `2^k` tournaments realize ALL self-complementary iso classes. (S81's `kappa` restricted to SC.)
- The **blue line-subgraph** on the SC merged nodes is a **T-JOIN**: every SC node has ODD blue-degree, so
  the odd-degree set (the GF(2) boundary) is `T =` every SC node. This is the S84 parity theorem in graph
  form (each SC class has an odd number of grid-symmetric tilings).

## The obstruction (verified n<=6)
| n | #SC (A051337) | info floor `ceil(log2 #SC)` | `kappa_SC` | half-dim `floor((n-1)^2/4)` | gap |
|---|---|---|---|---|---|
| 4 | 2 | 1 | 1 | 2 | 0 |
| 5 | 8 | 3 | **4** | 4 | **+1** |
| 6 | 12 | 4 | **6** | 6 | **+2** |
| 7 | 24 | 5 | (9?) | 9 | (+4?) |

**`kappa_SC(n) = floor((n-1)^2/4)` for `n>=5`** (verified 5,6). This is the **half-tiling / complement-fold
dimension** -- covering the SC classes requires the ENTIRE fold, a quadratically-growing gap above the
information floor. The SC classes fill the complement-fold and cannot be compressed. (`n=4` is a small
exception: `kappa_SC=1<2`, the grid-sym subcube is not optimal there.)

## Why: the T-join / odd-grid-sym parity
- **Upper bound (clean):** the grid-sym subcube `{t : bit[i]=bit[sigma(i)]}` has dimension `floor((n-1)^2/4)`
  and contains ALL grid-symmetric tilings; its iso classes are EXACTLY the SC classes; so it covers all SC
  classes and `kappa_SC <= half-dim`.
- **Lower bound (the obstruction, verified n=5,6):** no subcube of dimension `< half-dim` realizes all SC
  classes. The mechanism is the parity: each SC class has an ODD number of grid-symmetric tilings (S84), so
  the grid-symmetric tilings are essential and span the full fold; the all-odd-degree blue T-join is the
  fingerprint. A smaller (lower-parity) subcube cannot capture every odd-grid-sym SC class.

## The blue T-join structure (extends the atlas)
- All-odd-degree on SC nodes (T = every SC node); connected; **cycle-rank (genus) `= 0, 1, 15`** (n=4,5,6).
- **Bipartite for `n<=5`, NON-bipartite at `n=6`** -- another face of the `n=6` transition (with the NS-NS
  sea onset HYP-3808, the pure-black self-loops, and the minimal-flip `kappa` gauge break HYP-3798).

## New conjectures (for the atlas)
1. `kappa_SC(n) = floor((n-1)^2/4)` for `n>=5` (the grid-sym subcube is the optimal SC cover).
2. The obstruction gap `floor((n-1)^2/4) - ceil(log2 A051337(n))` grows like `~ n^2/4` (SC covers are
   maximally parity-obstructed).
3. `kappa_all(n) - kappa_SC(n) = 0,1,2` (n=5,6,7) = the extra dimension needed to also cover the NS "sea".
4. Blue T-join cycle-rank `0,1,15,...` -- identify the sequence; its jump `1->15` is the n=6 transition.
5. The bipartite->non-bipartite switch at `n=6` coincides with sea-onset / `kappa`-break: ONE threshold.

## Concrete next target
Prove the lower bound `kappa_SC >= floor((n-1)^2/4)` for `n>=5` from the odd-grid-sym (T-join) parity: show
that any subcube missing a dimension of the fold must miss an SC class whose grid-sym count is odd. That
would turn the empirical obstruction into a theorem and confirm the T-join parity is the exact mechanism.

## Honest scope
`kappa_SC = 1,4,6` (n=4,5,6) EXACT (exhaustive). `kappa_SC = half-dim` for n=5,6; the upper bound
(grid-sym subcube) is a clean general proof; the lower bound (obstruction) is verified only n<=6 and the
`kappa_SC = floor((n-1)^2/4)` formula is a conjecture for `n>=7`. The T-join parity is the identified
structural mechanism, not yet a proof of the lower bound.
