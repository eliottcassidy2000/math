---
id: HYP-3722
title: THE ONE REMAINING NODE + the realizability node in OBSERVER terms -- the covering-min witness t* (the lonely runner's escape time, = M since speed-1 binds) is a continued fraction on the Stern-Brocot path from 1/(n-1): the MEDIANT 2/(2n-1)=[0;n-1,2] for n<=6 (projective/drop-2 regime, while PG(2,n-1) exists), the CONVERGENT n/Phi_6(n)=[0;n-1,n] for n>=7 (hexagonal/Eisenstein). WHY the convergent not the mediant for n>=7 (observer terms): the projective/mediant escape is BLOCKED at the first PG-plane failure PG(2,6) (n=7, opus), so the lonely runner must DESCEND to the deepest (Diophantine-BEST, Lagrange) escape -- the full convergent -- where the other runners spread into the optimal hexagonal AP (three-distance gaps {1,n,2n}). THE SPECTRAL NODE: the convergent = the SL2 CF-recurrence matrix M=[[Phi_6,n],[n-1,1]] (det 1, unimodular); the CF recurrence q_k=a_k q_{k-1}+q_{k-2} IS tridiagonal/Jacobi, so the convergent = the Pade/Gauss-quadrature spectral edge = the A2-'tridiagonalized'-Catalan(semicircle) Jacobi edge (klein); the CERTIFICATE that the construction is optimal = LAGRANGE (convergents are the best rational approximations) = MARKOV-STIELTJES (the Jacobi convergent is the optimal spectral/quadrature bound) -- CLASSICAL. So the one remaining node REDUCES to the LRC->CF (1D->2D) reduction (the covering constraint forces the convergent escape past n=7); the spectral certificate is then classical
status: COMPUTED + verified (CF expansions [0;n-1,2] vs [0;n-1,n]; the SL2 matrix M det 1; the Stern-Brocot path). Converges with klein-S28 (HYP-3717: CF=[0;n-1,n], three-gap covering-min, torus lift = 2-row hexagonal patch) + opus (three-distance, PG transition). The observer framing + the spectral-certificate classification are the distinctive pieces; the certificate is classical, the open node is the LRC->CF/1D->2D reduction (= klein-S27/S28).
source: mac-mini-2026-07-01-S46
related:
  - HYP-3704  # mac-mini S45: the three routes (lift/three-distance/spectral); this resolves the spectral node as classical CF/Jacobi
  - HYP-3717  # klein-S28: three-gap covering-min, CF=[0;n-1,n], torus lift = 2-row patch, open 1D<->2D bridge
  - HYP-3716  # klein-S27: 2D side = theorem (Kershner); the gap = LRC->2D
  - HYP-3701  # the covering-min n-dependent transition (mediant n<=6 / convergent n>=7)
  - HYP-3564  # relations-not-things / the observer (relational) viewpoint
results:
  - 04-computation/observer_convergent_vs_mediant_macmini_20260701.py
  - 05-knowledge/results/observer_convergent_vs_mediant_macmini_20260701.out
---

# HYP-3722 -- the observer's escape (convergent, not mediant) and the spectral node

The owner asked to work the one remaining node (the spectral certificate / A2-Jacobi edge) AND to push the
realizability node in OBSERVER terms -- why the convergent, not the mediant, is the observer's escape for
n>=7. They are the same object, seen two ways.

## The observer's escape: a continued fraction
The covering-min witness `t*` is the lonely runner's escape time. Take the observer to be runner 1 (speed
1); it binds (`||1*t*|| = t* = M`), so `t* = M`. As a continued fraction on the Stern-Brocot path from
`1/(n-1)`:
- **n <= 6: the MEDIANT** `2/(2n-1) = mediant(1/n, 1/(n-1)) = [0; n-1, 2]` -- the shallow escape (partial
  quotient 2), the projective/drop-2 regime.
- **n >= 7: the CONVERGENT** `n/Phi_6(n) = [0; n-1, n]` -- the full convergent (partial quotient n), the
  deepest descent from `1/(n-1)`, the hexagonal/Eisenstein escape.
The Stern-Brocot path `1/(n-1) -> [0;n-1,2] -> [0;n-1,3] -> ... -> [0;n-1,n]` climbs the mediants to the
convergent endpoint.

## Why the convergent, not the mediant, for n >= 7 (observer terms)
The lonely runner escapes by finding a time `t*` that keeps every OTHER runner far from the origin -- a
simultaneous Diophantine approximation. The MEDIANT `[0;n-1,2]` is a shallow, "first-compromise"
approximation; it spreads the other runners into a COARSE arithmetic progression that suffices only while
the covering is loose. That looseness is exactly the **projective-plane regime**: for `n<=6` (`q=n-1<=5`,
all prime powers) `PG(2,n-1)` exists, the drop-2/difference-set covering realizes the mediant escape, and it
is the tightest. At **n=7** the first projective plane FAILS (`PG(2,6)`, Bruck-Ryser/Euler; opus), the
mediant/difference-set escape is BLOCKED, and the observer must DESCEND to the **full convergent** -- the
Diophantine-BEST approximation (Lagrange) -- where the other runners spread into the OPTIMAL hexagonal AP
(three-distance gaps `{1, n, 2n}`, the `omega = mult-by-n` rotation). So: *the convergent is the observer's
escape for n>=7 because the shallow mediant escape only exists while the projective plane does; past `q=6`
the runner falls through to the deepest (best-approximation) hexagonal/Eisenstein convergent.* The
convergent is "the observer's escape" in the precise sense of Lagrange: it is the best rational a runner can
use to avoid all the others.

## The spectral node = the convergent = a Jacobi/Stieltjes edge (the certificate is classical)
The convergent `[0;n-1,n]` is computed by the CF recurrence `q_k = a_k q_{k-1} + q_{k-2}`, whose matrix is
the **SL2/unimodular** `M = [[n,1],[1,0]][[n-1,1],[1,0]] = [[Phi_6, n],[n-1, 1]]`, `det M = Phi_6 - n(n-1) =
1`; the convergent `n/Phi_6` is its slope. A continued-fraction recurrence IS a tridiagonal (Jacobi)
object: the convergents are the **Pade / Gauss-quadrature** approximants of the associated Jacobi operator.
This is exactly klein's spectral lead: the A2 lattice "tridiagonalized" is the Jacobi operator with
**Catalan moments** (the semicircle), and the `zeta_6`-line covering radius is its **spectral edge** -- which
is the convergent `n/Phi_6`. The CERTIFICATE that the construction is the cyclic-Kershner optimum is then
**classical**: convergents are the **best rational approximations** (Lagrange), equivalently the Jacobi
convergent is the **optimal spectral/quadrature bound** (Markov-Stieltjes). The "Eisenstein-symmetric
Fourier-positive certificate" of HYP-3704 IS this Markov-Stieltjes/Lagrange extremality of the
`omega`-symmetric Jacobi convergent.

## The one remaining node, clarified
So the spectral node is NOT a hard new inequality -- it is **classical** (Lagrange / Markov-Stieltjes:
convergents/quadrature-edges are optimal). What remains genuinely open is the **LRC -> CF reduction**: that
the COVERING CONSTRAINT forces the covering-min's escape to be the convergent (not a shallower mediant or
some other rational) for `n>=7`. That is exactly klein-S27/S28's "`1D <-> 2D` metric bridge" / "LRC -> 2D
reduction." Once it holds, the convergent's optimality (and hence `M = n/Phi_6 >= 1/n`, trivially) follows
from the classical theory. The observer framing makes the open node concrete: *show that, past the
projective-plane failure at `n=7`, the lonely runner's tightest escape is forced to be the hexagonal
convergent.*

## What it buys
Unifies the two halves: the observer's escape (convergent for `n>=7`, mediant for `n<=6`) and the spectral
node are the SAME continued-fraction convergent, which is a Jacobi/Stieltjes spectral edge. It RESOLVES the
spectral certificate as classical (Lagrange/Markov-Stieltjes), pinning the single open node to the
`LRC -> CF (1D->2D)` reduction; and it gives the realizability/transition (`mediant -> convergent` at `n=7 =
PG(2,6)`) an observer meaning -- the lonely runner falls from the shallow projective escape to the deep
hexagonal convergent exactly when the projective plane fails.
