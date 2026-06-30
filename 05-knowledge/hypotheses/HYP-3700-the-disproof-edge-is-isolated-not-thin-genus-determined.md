---
id: HYP-3700
title: PUSHING THE DISPROOF MODE + the EXACT NATURE of the proof/disproof edge -- the disproof boundary (apex gap=0) is a SINGLE POINT (the full Z_p core = the complete mod-p covering), and its ISOLATION from all proper cores (the min positive gap = how close a disproof can get) is GENUS-DETERMINED: the DOUBLET 4sin^2(pi/2p) for genus<=1 (p=3,5,7: 1.0, 0.382, 0.198), but a LARGER near-vanishing-character-sum core for genus>=2 (p=11,13: {1,2,4,6,8,9}-type, 0.0078, 0.0049). So LRC14 (p=7, genus 1) is the LAST case where the discrete edge is ISOLATED (by the doublet 0.198, a spectral gap with NO cores in (0,0.198)) -- its 'razor-thinness' is purely a MEASURE artifact (the per-level product ->0 over deep descents, klein-S16), NOT a discrete approach to the boundary; for genus>=2 the edge genuinely goes razor-thin (proper cores approach 0). Pushing the disproof: a disproof must land EXACTLY on the full Z_7 (it cannot slide toward it -- isolated jump 0.198), where the lonely measure is 0 but EXISTENCE (the phi(n) units) gives a lonely point -- so the edge is closed by COUNTING, not measure or forbidden-H; the forbidden-H 7 mirrors the special-ness
status: COMPUTED + verified (edge=full Z_p only; min-positive-gap by p: 1.0,0.382,0.198,0.0078,0.0049; p=11 minimizer = {1,2,4,6,8,9} size 6, NOT a doublet). Rediscovers + sharpens klein-S11's genus story (genus>=2 => larger cores bind below). The edge-isolation = genus characterization is the new framing; the disproof is closed by existence (established, klein-S16/S38), not by a new mechanism.
source: mac-mini-2026-06-30-S41
related:
  - HYP-3617  # the forbidden-H disproof mode (pushed here)
  - HYP-3587  # klein-S11: genus = the obstruction; genus>=2 => larger cores bind below (0.0078,0.0049)
  - HYP-3548  # the two razor-thin lines (refined here: the gap-edge is ISOLATED for LRC14, not thin)
  - HYP-3615  # the small-measure regime: existence (the units) closes the edge, measure vanishes
  - THM-590   # the apex gap law (doublet min for p=7; gap=0 only at full Z_p)
results:
  - 04-computation/forbidden_H_certificates_macmini_20260630.py  # (related; the H-spectrum certificates)
---

# HYP-3700 -- the disproof edge is isolated (not thin) for LRC14; razor-thinness is genus-determined

The owner asked to push the disproof mode while being precise about the razor-thin edge between proof and
disproof and its patterns. The precision overturns the intuition: for LRC14 the edge is NOT razor-thin.

## The edge is a single point
The disproof boundary is `gap(O)=0`. By the irreducibility of `Phi_p` (THM-590), the apex gap
`min_{k!=0}|sum_{x in O} w^{kx}|^2` is `0` ONLY at the full core `O=Z_p` (the complete mod-`p` covering) --
a SINGLE point, measure 0. Everywhere else (every proper non-empty core) the gap is positive. So the
proof/disproof boundary is one isolated configuration: the complete covering.

## The isolation is GENUS-determined (the key pattern)
How close can a PROPER core (a genuine partial covering) get to the edge? The min positive gap, by apex
prime `p` (and `genus(X_0(2p))`):
| `p` | `genus` | min positive gap (edge isolation) | binding core |
|---|---|---|---|
| 3 | 0 | `1.000` | doublet |
| 5 | 0 | `0.382` | doublet |
| **7** | **1** | **`0.198 = 4cos^2(3pi/7)`** | **doublet** |
| 11 | 2 | `0.0078` | `{1,2,4,6,8,9}` (size 6, NOT a doublet) |
| 13 | 2 | `0.0049` | larger core |
- For **genus <= 1** (`p<=7`): the closest proper core is the DOUBLET, `4sin^2(pi/2p)`. The edge is
  ISOLATED -- there is a spectral GAP, NO proper core in `(0, doublet)`.
- For **genus >= 2** (`p>=11`): a LARGER core with a NEAR-VANISHING character sum (e.g. `{1,2,4,6,8,9}` mod
  11) approaches `0` (`0.0078`, `0.0049`) -- the doublet is no longer the minimum. The edge becomes
  RAZOR-THIN.
This rediscovers and sharpens klein-S11 (HYP-3587: genus>=2 => larger cores bind below): the GENUS literally
measures how thin the disproof edge is.

## LRC14 is the LAST isolated-edge case
`p=7` is `genus 1`: the discrete edge is ISOLATED by the doublet `0.198`. A disproof cannot continuously
slide toward the boundary -- between any partial covering (gap `>= 0.198`) and the complete covering
(gap `0`) there is a JUMP of `0.198`, with nothing in between. So **for LRC14 the proof/disproof edge is NOT
razor-thin at the discrete (gap) level -- it is a clean, isolated point.** The apparent "razor-thinness" is
entirely a MEASURE artifact: the TOTAL floor is the product of per-level gaps, and over a deep 2-adic
descent that product `-> 0` (klein-S16, `inf R'=0`) -- a soft, accumulative thinness, not a discrete
approach to the disproof boundary. The two razor-thin lines (HYP-3548) refined: the gap-`M` line is far
(`+10%`), the measure line is thin (`inf=0`), but the DISCRETE gap-edge is ISOLATED (`0.198`).

## Pushing the disproof mode
A would-be LRC14 disproof must reach the edge `gap=0`, i.e. land EXACTLY on the full `Z_7` core (the
complete mod-`7` covering) at a binding descent level -- it cannot approach it (isolated by `0.198`). But
the full `Z_7` is the apex CUSP: the lonely MEASURE there is `0`, yet EXISTENCE survives -- the `phi(n)`
unit touch-points are a lonely set of measure 0 but nonempty (HYP-3615, klein-S16). So a disproof at the
edge fails by COUNTING (the units are there), not by measure. The forbidden-H value `7 = Phi_3(2) = |Fano|`
(HYP-3616) MIRRORS this -- `7` is the "pure" value a complete structure cannot realize -- but the actual
closer of the edge is the existence/counting of the units, not a forbidden-H certificate. So the disproof
mode, pushed to the edge, lands on the same wall the rest of the program found: the edge is a measure-0
complete-covering point, closed by the existence of the unit witnesses.

## The patterns (summary)
1. **Edge = a single point**: the full `Z_p` core (complete covering), measure 0. (`Phi_p` irreducible =>
   no other gap-0 core.)
2. **Isolation = genus-determined**: doublet `4sin^2(pi/2p)` for genus `<=1`; a near-vanishing-sum larger
   core `-> 0` for genus `>=2`.
3. **LRC14 (p=7, genus 1) = the frontier**: the LAST case with an ISOLATED edge (`0.198`); razor-thinness
   begins at genus `>=2` (`p>=11`).
4. **Two kinds of thinness**: DISCRETE (the gap-edge -- isolated for LRC14, thin for genus>=2) vs MEASURE
   (the product `->0` -- thin for all, a soft artifact). LRC14 is thin in measure, isolated in gap.
5. **The edge is closed by EXISTENCE** (the units at the cusp), not by measure (vanishes) or by a
   forbidden-H certificate (mirrors but doesn't close).
6. **Primality + oddness of the apex** are the edge-keepers: primality (`Phi_p` irreducible) makes the edge
   a single point; oddness (`C_p` non-bipartite) makes the per-level gap positive.

## What it buys
The exact nature of the proof/disproof edge: a single measure-0 point (the complete covering), with
genus-determined isolation. LRC14's edge is ISOLATED (the doublet `0.198`), not razor-thin -- the thinness
is a measure artifact; genuine razor-thinness is a genus`>=2` (`p>=11`) phenomenon, so LRC14 is the last
tractable frontier. Pushing the disproof mode to the edge confirms it is closed by the existence of the
unit witnesses (the small-measure/counting regime), with the forbidden-H `7` as its counting-side mirror.
