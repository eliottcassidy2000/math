---
id: HYP-3773
title: RIEMANN-HURWITZ DESCENT extends Gauss-Bonnet -- the repo's 2-adic parity descent (THM-580, peel n=2p -> odd apex p) IS the degeneracy branched cover X_0(2p) -> X_0(p) (degree 3), and Riemann-Hurwitz chi(2p)=3 chi(p)-R (= Gauss-Bonnet for covers) removes the level-2 curvature, LOWERING the genus to the g=0 (spherical, obstruction-free) apex. The removed curvature R (ramification) IS the iota-odd RESIDUAL: for n=14, X_0(14)->X_0(7) drops genus 1->0 with R=6, and the removed genus-1 = the cusp form f_14 = elliptic curve 14a (S56/HYP-3768). Verified genus(X_0(N)) via the standard formula (g=0,0,1,2,2 for X_0(6,10,14,22,26)) and R = 4-6g(p)+2g(2p) = 4,4,6,2,8 for p=3,5,7,11,13. This unifies THM-580 (analytic measure descent), HYP-3768/S57 (curvature=genus, the residual), mac-mini HYP-3771 (the (2,3,p) apex spine), and opus HYP-3770 (the O(log) reciprocity descent) as ONE descent seen four ways: measure / genus-curvature / triangle-group / continued-fraction. KERSHNER anchor: the covering-min = flat hexagonal A_2 covering is the g=1 MIDDLE of the descent tower (spherical/Platonic below, hyperbolic above); n=14 is the flat frontier
status: MIXED (rigorous RH computation + synthesis). PROVED/standard: genus(X_0(N)) formula (verified vs known g at N=7,10,11,13,14,15,22,23,26); Riemann-Hurwitz for the degree-3 degeneracy cover X_0(2p)->X_0(p): chi(2p)=3chi(p)-R, R computed 4,4,6,2,8. VERIFIED: the genus descent 2p->p (n=14: 1->0), R=4-6g(p)+2g(2p). SYNTHESIS (not proof): identifying R with the LRC iota-odd residual / f_14; the THM-580 (measure) <-> Riemann-Hurwitz (curvature) descent bridge; the covering-min = Kershner flat middle. The exact per-cusp/elliptic ramification breakdown and the rho_j <-> R correspondence are untested next levers.
source: klein-2026-06-30-S58
depends_on:
  - THM-580    # the 2-adic parity descent (the analytic side of this cover)
  - HYP-3772   # curvature = genus (S57); this descends it
related:
  - HYP-3768   # covering-min = E2/Eisenstein bulk; residual = genus cusp form (the R here)
  - HYP-3771   # mac-mini: the (2,3,p) apex-prime spine (the descent walks it)
  - HYP-3770   # opus: the O(log) reciprocity descent (the CF side of the descent)
  - HYP-3587   # genus X_0(2p) = the obstruction dimension
  - HYP-3706   # Kershner hexagonal covering (the flat g=1 middle)
  - HYP-2108   # apex descent 127/127 Z_7 cores (the g=0 base)
results:
  - 04-computation/riemann_hurwitz_descent_klein.py
  - 05-knowledge/results/riemann_hurwitz_descent_klein.out
---

# HYP-3773 — Riemann-Hurwitz descent extends Gauss-Bonnet

## The move
Last session (HYP-3772) put the LRC family on the curvature/genus trichotomy and identified the coverage
crux's obstruction with `chi` (Gauss-Bonnet). The repo's central descent is THM-580's **2-adic parity
descent**: peel `n = 2p -> p` (the odd apex), multiplicatively on the lonely measure. This session:
**that peel IS a branched cover, and Riemann-Hurwitz is its Gauss-Bonnet.**

The degeneracy map `X_0(2p) -> X_0(p)` (from `Gamma_0(2p) subset Gamma_0(p)`) has degree
`psi(2p)/psi(p) = 3`. **Riemann-Hurwitz** (Gauss-Bonnet for branched covers):
```
    chi(X_0(2p)) = 3 * chi(X_0(p)) - R,     R = total ramification (removed curvature).
```
Equivalently `R = 4 - 6 g(p) + 2 g(2p)`.

## The computation (genus + ramification, verified)
```
  p    g(p)  g(2p)   chi(p) chi(2p)   R
  3     0     0        2      2        4
  5     0     0        2      2        4
  7     0     1        2      0        6     <== n=14: genus 1 -> 0, curvature R=6 removed
 11     1     2        0     -2        2
 13     0     2        2     -2        8
```
Genus computed from the standard formula `g = 1 + psi/12 - nu2/4 - nu3/3 - nu_inf/2` (sanity-checked
against known genera: `X_0(11)=1, X_0(14)=1, X_0(22)=X_0(23)=X_0(26)=2`).

## What the descent does: removes the residual curvature
Following the cover **down** (`2p -> p`) is the descent; it drops the genus from `g(2p)` to `g(p)`,
removing `g(2p)-g(p)` units of curvature, produced by the ramification `R` via Riemann-Hurwitz. For
**`n = 14`**: `X_0(14)` (genus 1, the LRC target) is a degree-3 branched cover of `X_0(7)` (genus 0, the
apex) with `R = 6`; the descent removes exactly the **genus-1** piece -- the `iota`-odd cusp form
`f_14 = ` elliptic curve `14a` (S56/HYP-3768). The descent **bottoms at the genus-0 (spherical) apex
`X_0(7)`**, which has no cusp-form obstruction (the base case; the `127/127` `Z_7`-cores of HYP-2108).
So:
> the hard residual (the genus-1 obstruction) = the ramification `R` the 2-adic descent concentrates;
> descending to the odd apex removes it, landing in the obstruction-free genus-0 base.

## One descent, four faces
- **Measure** (THM-580): `meas(lonely S_{2p}) = rho * meas(lonely S_p)` -- the analytic peel.
- **Genus/curvature** (here): `chi(2p) = 3 chi(p) - R` -- the geometric peel (Riemann-Hurwitz).
- **Triangle group** (mac-mini HYP-3771): the peel walks the `(2,3,p)` angle-defect spine; `p=7` is
  where it turns hyperbolic (`1/2+1/3+1/7-1 < 0`) = where `X_0(2p)` turns genus 1 = the frontier.
- **Continued fraction** (opus HYP-3770): the covering-min margin is an `O(log)` reciprocity/Euclidean
  descent (`1/M = [0; n-1, rung]`).
The `rho_j` (measure decorrelation) `<->` `R` (ramification) correspondence is the natural next bridge.

## Kershner anchor
The covering-min `= n/Phi_6` is the **flat hexagonal `A_2` covering** (Kershner's thinnest plane
covering, HYP-3706) -- the `genus-1`, curvature-zero **middle** of the descent tower: spherical/Platonic
below (`n=6,10`, `g=0`), hyperbolic above (`n=22,26`, `g>=2`). The 2-adic descent moves `n=14` off the
flat middle down to the spherical base; Kershner is the geometry AT the frontier.

## Honest scope
Riemann-Hurwitz and the genus formula are exact/standard (rigorous). The identification of the
ramification `R` with the LRC `iota`-odd residual `f_14`, and the THM-580(measure)<->RH(curvature)
descent bridge, are SYNTHESIS -- they say the two descents have the same shape and the same fixed point
(the `g=0` apex), not (yet) that they are the same map. The per-cusp/elliptic breakdown of `R` and the
`rho_j <-> R` correspondence are the next levers.

## Net
The repo's 2-adic parity descent is the degeneracy branched cover `X_0(2p) -> X_0(p)`, and Riemann-Hurwitz
(`chi(2p)=3chi(p)-R`) is its Gauss-Bonnet: the descent **removes the level-2 curvature `R`**, lowering
the genus to the `g=0` obstruction-free apex, and the removed curvature is the `iota`-odd residual (`n=14`:
`g 1->0`, `R=6`, `= f_14 =` curve `14a`). Gauss-Bonnet extended by descent = Riemann-Hurwitz; the LRC-14
residual is the ramification of the level-2 cover, and the descent's base case is the spherical apex.
