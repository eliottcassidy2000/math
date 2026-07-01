---
id: HYP-3813
title: SIX BRIDGES integrating the recent merged-metagraph structure (S81-S86) with the LRC work (S60-S80) -- both live on the staircase triangle folded by a Z_2, with a Gram/spectral fine-structure layer. GROUNDED BRIDGE (2): the tournament spectral-twin separator (the adjacency matrix A = the tournament relation) and the LRC floor's danger relation composed with itself (D D^T, HYP-3571) are the SAME move -- the Gram / 2nd-moment spectrum of a 0/1 relation; the elementary counts are cospectral (tournament H/|Aut|/c3/scores; LRC M/trace) while the fine structure is spectral. VERIFIED: the LRC danger-Gram (D D^T) of construction/AP/GW all have identical trace 1.857 (=|S|*2r, the 1st moment) but DISTINCT eigenvalue spectra (top eig 0.5006/0.4992/0.4803) = the LRC analog of the S86 invariant-cospectral spectral twins; and the S76 near/far split is visible IN the Gram (far element off-diag ~ (2r)^2 independent, near-core correlated).
status: SYNTHESIS (Bridge 2 grounded computationally; the other 5 are structural correspondences, one already canon THM-589). A cross-domain integration proposing testable bridges, not a theorem. The Gram-spectrum bridge is the sharpest concrete anchor; the Z_2-fold, parity, apex-7, and covering-obstruction bridges are structural.
source: mac-mini-2026-07-01-S87
related:
  - HYP-3811   # S86 spectral twins (the tournament side of Bridge 2)
  - HYP-3571   # LRC floor: the danger relation D D^T (the LRC side of Bridge 2)
  - HYP-3810   # S85 T-join covering obstruction (Bridge 5)
  - HYP-3789   # S76 moment relaxation / near-far split (Bridge 2 mechanism)
  - THM-589    # W(n) = metagraph H-variance <-> LRC floor CV(N_R)^2 (Bridge 6, already canon)
  - HYP-3771   # S65 (2,3,n) angle-defect, apex-7 (Bridge 4)
results:
  - 04-computation/metagraph_lrc_integration_macmini_20260701.py
  - 05-knowledge/results/metagraph_lrc_integration_macmini_20260701.out
---

# HYP-3813 -- six bridges: the merged metagraph and the LRC are two folds of one triangle

Owner: find novel intriguing ways to integrate the recent merged-metagraph structure with LRC. The unifying
picture: **both the tournament metagraph and the LRC live on the staircase triangle `delta_{n-2}`, each
folded by a `Z_2`, with a 2nd-moment / Gram spectral layer carrying the fine structure.**

## Bridge 1 -- THE Z_2 FOLD (complement <-> iota)
- Tournament: the merge is by **complement** `sigma` (reverse all arcs = transpose = the `y=x` staircase
  mirror, THM-549). Self-complementary (SC) classes = `sigma`-fixed; the half-tiling is the fold's
  fundamental domain.
- LRC: the **iota** involution `t -> 1-t`. The lonely set is `iota`-symmetric; the covering-min binding pair
  is `{t* , 1-t*} = {14/183, 169/183}` (the 2 atoms, S76).
- **Same `Z_2` on the same triangle.** SC classes (fold-fixed tournaments) `<->` the `iota`-symmetric
  extremal lonely times. Covering the fold-fixed family is obstructed on both sides (Bridge 5).

## Bridge 2 -- THE GRAM / SPECTRAL SEPARATOR (GROUNDED)
The S86 tournament **spectral twins** are invariant-cospectral (identical `H, |Aut|, c3, scores`) but
separated by the **adjacency spectrum** of `A` (the tournament relation). The LRC **floor** (HYP-3571) is
"the danger relation `D` composed with itself, `D D^T`" -- the pairwise-overlap Gram (S76). VERIFIED: the
LRC danger-Gram of `construction / AP / GW` all have the SAME trace `1.857 = |S|*2r` (the 1st moment) but
DISTINCT spectra (top eigenvalue `0.5006 / 0.4992 / 0.4803`). So the LRC sets are **"trace-cospectral,
spectrally distinct"** -- the exact analog of the tournament spectral twins. And the S76 near/far split is
visible IN the Gram (far element off-diagonal `0.0201 ~ (2r)^2` independent; near-core `0.0299` correlated).
**One move: the Gram / 2nd-moment spectrum of a 0/1 relation carries the structure the elementary counts
miss** -- tournament `A` (adjacency), LRC `D D^T` (danger). This is the sharpest concrete integration.

## Bridge 3 -- THE EVEN/ODD `Z_2`-GRADING
- Tournament: SC odd tiling count / NS even (S84 parity theorem); the **black lines form an even-degree
  (Eulerian) graph, the blue lines an odd graph** (S83); RPGFD "even graph" = tournaments (S82).
- LRC: the `iota`-**even** (Eisenstein / `E_2` bulk) vs `iota`-**odd** (cusp form `f_14`; the margin = the
  Dedekind sum, an `iota`-odd sawtooth object) split.
- The even/odd grading under the fold is shared; the tournament black-even-graph IS an even graph in the
  same sense the LRC parity dual chases.

## Bridge 4 -- APEX-7 and the `n=6` EUCLIDEAN BOUNDARY
- Tournament: forbidden `H in {7, 21}`; `G_5` = icosahedron (spherical); `n=6` spherical->hyperbolic
  transition (5 faces, S86); odd holes / diameter jump at `n=7`.
- LRC: `LRC14 = 2*7`; apex-7 = Klein quartic `PSL(2,7)`; the `(2,3,n)` angle defect makes `n=6` EUCLIDEAN
  (hexagonal `A_2` = the covering-min `Phi_6`), `n=7` hyperbolic (genus jump, HYP-3771/S65).
- **`7` is the apex and `n=6` the flat/Euclidean boundary on BOTH sides.** The metagraph goes hyperbolic at
  `n=6`, one step before the LRC apex at `n=7`; the covering-min itself lives at the `n=6` Euclidean node.

## Bridge 5 -- COVERING OBSTRUCTION (T-join parity <-> covering-min rigidity)
- Tournament: covering the SC classes is **parity-obstructed** to the half-tiling dimension
  `kappa_SC = floor((n-1)^2/4)` by the all-odd-degree blue **T-join** (S85).
- LRC: the **covering-min** is rigidity-obstructed (no covering set beats `n/Phi_6`; the construction is an
  isolated deep well, S69/S77).
- Both: **covering a self-symmetric extremal family is parity/rigidity-obstructed far above the naive
  bound.** The half-tiling dimension `floor((n-1)^2/4)` is the tournament analog of `Phi_6`'s rigidity.

## Bridge 6 -- 2nd MOMENT (already canon, THM-589)
`W(n)` = the metagraph `H`-variance = the simplicial-Redei succession count, and `Var(H) = (n!/4^{n-1})
(W(n)-n!)`; the LRC analog is `CV(N_R)^2` (the floor's 2nd moment). The S86 spectral twins refine the `H`
structure; the LRC danger-Gram refines the floor -- **the same 2nd-moment layer, the same `Z_2`-folded
triangle.**

## New conjectures / next targets
1. **Spectral-twin LRC test:** do covering sets that are M-cospectral AND covering-profile-identical but
   danger-Gram-spectrally distinct exist (true LRC spectral twins)? The construction vs a tuned near-miss.
2. **Fold-fixed correspondence:** do the SC (fold-fixed) tournament classes map to the `iota`-symmetric LRC
   witnesses under a shared staircase coordinate?
3. **Euclidean-`n=6` unification:** is the metagraph's `n=6` spherical->hyperbolic transition the same event
   as the LRC `(2,3,n)` Euclidean node (both `n=6`), and the apex-`7` the same on both?
4. **Half-tiling `<-> Phi_6`:** relate `floor((n-1)^2/4)` (SC-cover obstruction) to `Phi_6 = n^2-n+1`
   (covering-min) -- both quadratic rigidity dimensions of the fold.

## Honest scope
Bridge 2 is GROUNDED (the LRC danger-Gram is trace-cospectral / spectrally-distinct, the exact analog of the
tournament spectral twins). Bridge 6 is canon (THM-589). Bridges 1,3,4,5 are structural correspondences on
the shared `Z_2`-folded staircase -- intriguing and testable, not yet theorems. A cross-domain synthesis
opening concrete tests, not a proof.
