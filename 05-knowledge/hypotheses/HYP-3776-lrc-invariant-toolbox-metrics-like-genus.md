---
id: HYP-3776
title: A DIVERSE METRICS TOOLBOX for LRC(2p) -- beyond the genus, which niche invariants detect the apex frontier, and along which axis. Computed a broad panel across apex primes p=3,5,7,11,13,17,19,23 and sorted every metric by WHAT it detects. TWO orthogonal axes emerge. AXIS 1 (the FRONTIER, spherical->hyperbolic at p=7, all = the genus jump 0->1): genus g(X0(2p))=1+psi/12-nu2/4-nu3/3-cusps/2 (=dim S2^new = #cusp-form obstructions f14); the (2,3,p) TRIANGLE-GROUP order 4/(1/2+1/3+1/p-1), FINITE (24,120) for p<=5, INFINITE (hyperbolic) p>=7; the GONALITY of X0(2p), 1 (rational P^1) for p<=5, 2 (elliptic/hyperelliptic) p>=7; the HURWITZ THRESHOLD -- (2,3,7) is the minimal hyperbolic triangle group and Klein quartic genus 3 has Aut=PSL(2,7)=168=84(g-1), so p=7 is THE Hurwitz prime. AXIS 2 (ORTHOGONAL: QR / p mod 4, the Paley-graph-vs-tournament dichotomy on 2p-1): nu2(2p)=1+(-1/p) =2 if p=1mod4, 0 if p=3mod4; h(-p)=class number of Q(sqrt-p), nonzero for p=3mod4 (h(-23)=3 the first >1); 2p-1 = a PALEY GRAPH (Ramanujan) if 1mod4, a rotational TOURNAMENT (Redei/OCF) if 3mod4 -- the covering-min circulant's type. THE UNIQUE-ELLIPTIC METRIC: p=7 is the ONLY apex prime where X0(2p) has genus EXACTLY 1, so X0(14) IS an elliptic curve (14a): conductor 14, RANK 0, TORSION = Z/6 -- and |torsion|=6=phi(14)=#(Z/14)^* = the lonely set (the S66 triple-6 recurs in the modular curve's torsion). NEW FRAMING: genus = dim(un-regularizable residual) (S67); the frontier-detectors are its geometric shadows, the QR-detectors an orthogonal covering-min axis. SCALE metrics (monotone, no transition): psi, phi(2p), J2, hyperbolic covolume psi*pi/3, apex gap 2-2cos(pi/p) ~ pi^2/p^2.
status: TOOLBOX / synthesis. The panel is exact/computed (genus, triangle-group, gonality, nu2, h(-p), psi, covolume, apex gap for p=3..23); the 14a facts (rank 0, torsion Z/6, conductor 14) are standard (LMFDB). The organization into two orthogonal axes (frontier vs QR) + the new framings (genus=residual dim; torsion=phi(14); p=7=Hurwitz/unique-elliptic) are the contribution. NOT a proof; a metrics catalog for the team, extending S65/S66/S67.
source: mac-mini-2026-06-30-S68
related:
  - HYP-3775   # klein-S59: the modular-invariant ZOO + newform metric (CONVERGENT, same task) -- primary companion
  - HYP-3771   # S65 the (2,3,n) angle-defect spine; genus tracks the geometry
  - HYP-3772   # S66 the triple 6 (phi(14)=6=lonely set); dim>=phi(n)
  - HYP-3774   # S67 genus = dim(un-regularizable residual); zeta-regularization
  - HYP-3586   # X0(14) cusps=Klein V4, apex cusp d=7, genus obstruction = f14
  - HYP-3715   # covering-min = zeta6-line; the circulant / Paley side
results:
  - 04-computation/lrc_invariant_toolbox_macmini_20260630.py
  - 05-knowledge/results/lrc_invariant_toolbox_macmini_20260630.out
---

# HYP-3776 -- a diverse metrics toolbox for LRC(2p)

The owner asked for "other formulas for metrics like genus -- search all around, diverse niche possibilities."
I computed a broad panel across apex primes `p = 3,5,7,11,13,17,19,23` and sorted every metric by **what it
detects**. The upshot: there are **two orthogonal axes**, plus a unique-elliptic coincidence at `p=7`.

## The panel

| p | genus | (2,3,p) grp | gonality | nu2 | h(-p) | psi | covol/pi | phi(2p) | apex gap | 2p-1 |
|---|-------|-------------|----------|-----|-------|-----|----------|---------|----------|------|
| 3 | 0 | 24 | 1 | 0 | 1 | 12 | 4 | 2 | 1.000 | 5 (Paley) |
| 5 | 0 | 120 | 1 | 2 | 0 | 18 | 6 | 4 | 0.382 | 9 (comp) |
| **7** | **1** | **INF** | **2** | 0 | 1 | 24 | 8 | 6 | 0.198 | 13 (Paley) |
| 11 | 2 | INF | 2 | 0 | 1 | 36 | 12 | 10 | 0.081 | 21 |
| 13 | 2 | INF | 2 | 2 | 0 | 42 | 14 | 12 | 0.058 | 25 |
| 23 | 5 | INF | ? | 0 | 3 | 72 | 24 | 22 | 0.019 | 45 |

## Axis 1 -- FRONTIER-DETECTORS (spherical -> hyperbolic at p=7 = the genus jump)
All of these fire at `p=7`, the first hyperbolic apex (S65):
- **genus** `g(X0(2p)) = 1 + psi/12 - nu2/4 - nu3/3 - cusps/2` `= dim S_2^new` `=` the number of cusp-form
  obstructions (`f_14`). `0,0,1,2,2,...` -- jumps `0->1` at 7.
- **`(2,3,p)` triangle-group order** `= 4/(1/2+1/3+1/p-1)`: finite (`24, 120`) for `p<=5`, **infinite**
  (hyperbolic) for `p>=7`.
- **gonality** of `X0(2p)`: `1` (rational `P^1`) for `p<=5`, `2` (elliptic/hyperelliptic) for `p>=7`.
- **Hurwitz threshold**: `(2,3,7)` is the minimal hyperbolic triangle group; the Klein quartic (genus 3) has
  `Aut = PSL(2,7) = 168 = 84(3-1)`, attaining the Hurwitz bound `84(g-1)`. `p=7` is **the Hurwitz prime**.

## Axis 2 -- QR / `p mod 4` DETECTORS (orthogonal: the covering-min circulant type)
A **different** structure, invisible to the genus:
- **`nu2(2p) = 1 + (-1/p)`**: `2` if `p = 1 (mod 4)`, `0` if `p = 3 (mod 4)` (`p=7 -> 0`).
- **`h(-p)`** (class number of `Q(sqrt(-p))`): nonzero for `p = 3 (mod 4)`; `h(-23) = 3` is the first `> 1`.
- **`2p-1`**: a **Paley GRAPH** (Ramanujan) if `= 1 (mod 4)`, a rotational **tournament** (Redei/OCF) if
  `3 (mod 4)` -- the type of the covering-min circulant (`2p-1 = 13` at `p=7` is a Paley graph, HYP-3715).

`LRC14` sits at `p=7`: **hyperbolic** on axis 1 (genus 1), **`3 mod 4`** on axis 2 (`nu2=0`, `2p-1=13` Paley).

## The unique-elliptic metric (the sharp coincidence)
`p=7` is the **only** apex prime where `X0(2p)` has genus **exactly 1** -- so `X0(14)` **is an elliptic curve**,
`14a` (`p<=5` give `P^1`, `p>=11` give genus `>=2`). And `14a` has conductor 14, **rank 0**, **torsion `Z/6`** --
so `|torsion| = 6 = phi(14) = #(Z/14)^* =` the LRC14 **lonely set** (the units). The **triple-6** of S66 recurs
as the torsion order of the modular curve. `p=7` is the "Goldilocks" apex: genus not `0` (would be rational,
trivial), not `>=2` (would be higher-genus, generic) -- exactly `1`, an elliptic curve, rank `0`, tractable.

## New framing: genus = dimension of the un-regularizable residual
Tying to S67: `genus = dim S_2^new =` the dimension of the **cusp-form part**, which is exactly the
**un-regularizable** piece of the floor (the Eisenstein/`E_2` bulk regularizes to `zeta`-values; the cusp form
does not). So the "metric like genus" *is* the genus, correctly read as the **residual dimension**:
`0` (fully `zeta`-regularizable, easy) for `p=3,5`; `1` (one obstruction, first hard) at `p=7`. The other
frontier-detectors (triangle group, gonality, Hurwitz) are its **geometric shadows**; the QR-detectors
(`nu2`, class number, Paley/tournament) are an **orthogonal axis** measuring the covering-min circulant's type,
not the residual.

## Honest scope
The panel is exact/computed; the `14a` facts (rank 0, torsion `Z/6`, conductor 14) are standard (LMFDB). The
contribution is the **organization into two orthogonal axes** (frontier vs QR) and three framings: genus `=`
residual dimension (S67), the unique-elliptic torsion `= phi(14)` (S66 triple-6), and `p=7 =` the Hurwitz /
unique-elliptic apex. This is a metrics toolbox for the team -- diverse invariants and what each one sees --
not a new proof step.
