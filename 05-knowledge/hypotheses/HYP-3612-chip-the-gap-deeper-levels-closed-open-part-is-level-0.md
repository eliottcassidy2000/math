---
id: HYP-3612
title: CHIPPING THE GAP -- THM-580's per-level CS bound rho_j >= 1-CV(N2_Oj).sqrt((1-m')/m') is POSITIVE for every deeper level j>=1 (verified CS floor 0.65-0.87 at levels 1,2 over standard + binding coverings) and FAILS only at the TOP level j=0; the reason is sharp: at j>=1 the cores O_j are small (CV(N2)<=0.5) AND the descended sets are small (m'>=0.69 >= the union bound 1-|S'|/7>0 once |S'|<=6), but at j=0 the top core O_0 (the full odd part, lonely-poor) has CV(N2)~1.4 (blows up); so the descent's CS machinery CLOSES the deeper levels, localizing the open part to level 0 = the original existence question rho_0>0 (HYP-3599); the least-eigenvalue certificate (HYP-3604/mac-mini-3606, bounded m-independent g in [4cos^2(3pi/7),2]) is the robust discrete object for the top, degenerating only at the exact cusp O_0=Z_7 (g=0, the tight extremal)
status: VERIFIED (regime split over consec{1..13}, tightest{1..12,182}, skip-12, binding {1..13}\7; CS floor by level). RIGOROUS pieces: the union bound m(lonely P)>=1-|P|/7 (each danger comb measure 1/7) makes m' bounded below for |S'|<=6; THM-580's CS inequality. CONDITIONAL: uniform small-core CV(N2) bound for full rigor. Refines/localizes HYP-3599. Converges with mac-mini-S37.
source: klein-2026-06-29-S20
depends_on:
  - HYP-3599   # the bridge splits by level; level 0 is the existence question
  - THM-580    # the per-level CS bound rho_j >= 1 - CV(N2).sqrt((1-m')/m')
related:
  - HYP-3604   # klein: the least-eigenvalue certificate = signless Laplacian of the odd cycle
  - HYP-3606   # mac-mini-S37: same certificate = non-bipartiteness (convergent)
  - HYP-3554   # the 14-sheet CV is unbounded (the whole-covering object the descent avoids)
  - THM-576    # the odd-set caps
results:
  - 04-computation/chip_gap_regime_split_klein.py
  - 05-knowledge/results/chip_gap_regime_split_klein.out
---

# HYP-3612 — chipping the gap: the deeper levels close, the open part is level 0

## The regime split (verified)

THM-580's per-level Cauchy-Schwarz bound is `rho_j >= 1 - CV(N2_{O_j}) . sqrt((1-m')/m')`,
`m' = m(lonely S^{(j+1)})`. Computing the CS FLOOR (the RHS) by descent level over standard + binding
coverings:
```
covering            lvl0      lvl1     lvl2
consec{1..13}       -0.52     0.74     0.87
tightest{1..12,182} -0.31     0.65     0.87
skip12{1..11,13,84} -0.59     0.74     0.84
binding {1..13}\7   -0.03     0.74     0.87
```
**The CS bound is POSITIVE at every deeper level `j >= 1` (`>= 0.65`) and negative only at the TOP level
`j = 0`.** So THM-580's standard certificate already closes the deeper levels; the open part is level 0.

## Why level 0 fails and the deeper levels close (the sharp reason)

`rho_j >= c` needs BOTH `CV(N2_{O_j})` small AND `m'` bounded below.
- **Deeper levels (`j >= 1`):** the descent halves sizes, so the cores `O_j` are small (`<= 3-4` residues),
  giving modest `CV(N2_{O_j}) <= 0.5`; and the descended sets are small (`|S^{(j+1)}| <= 6`), so the UNION
  BOUND `m(lonely P) >= 1 - |P|/7 > 0` (each danger comb has measure `1/7`) makes `m' >= 0.43-0.86`. Both
  good => CS floor `> 0`.
- **Top level (`j = 0`):** the core `O_0` is the FULL odd part (`6-7` residues, lonely-poor), so
  `CV(N2_{O_0}) ~ 1.4` BLOWS UP (`m(lonely O_0)` small => `CV^2 ~ 1/m -> infinity`). The measure CS bound
  cannot hold at the top regardless of `m'`. This is exactly HYP-3599's level-0 obstruction, now with its
  cause pinned: the large top core's variance.

The descent's value is precisely this: it never takes the CV of the WHOLE covering (HYP-3554's unbounded
14-sheet object); it uses per-level 2-sheet CVs of SMALL cores, which are tame -- isolating the single
large-CV difficulty to the top level.

## The least-eigenvalue certificate is the top-level tool

Where the measure CS bound dies (top, large core / `m' -> 0`), the least-eigenvalue certificate
(HYP-3604; mac-mini-S37/HYP-3606) is the robust object: `g(O) = lambda_min(Gram(O)) in [4cos^2(3pi/7), 2]`
is BOUNDED and `m`-INDEPENDENT (it never blows up), and positive because the apex cycle is odd. It
certifies the DISCRETE / existence content at the top (the odd cycle is present, non-degenerate) -- per
mac-mini-S37 and HYP-3599, the side the proof needs once the measure was retired (klein-S16, inf=0). It
degenerates ONLY at the exact cusp `O_0 = Z_7` (`g = 0`), which is the tight extremal (the disproof
boundary), not a generic open case.

## Net (what is chipped)

- **Closed (deeper levels `j >= 1`):** THM-580's CS bound is positive (verified `>= 0.65`); rigorous given
  the union bound on the small descended sets + a uniform small-core `CV(N2)` bound (finite to establish).
- **Open (level 0):** `rho_0 = m_S/(m_{O_0} m_{S'}) > 0` is the original existence question; the measure
  CS bound provably fails (top core CV blows up). The least-eigenvalue certificate covers its DISCRETE
  side; the measure->existence passage at the top (HYP-3599) is the residue.

So the gap is now localized to a SINGLE level (the top) with a named cause (the large top core's variance),
and the least-eigenvalue certificate is positioned exactly there. Converges with mac-mini-S37 (HYP-3606):
the discrete floor is discharged (odd cycle, non-degenerate); the top-level measure->existence passage
remains.
