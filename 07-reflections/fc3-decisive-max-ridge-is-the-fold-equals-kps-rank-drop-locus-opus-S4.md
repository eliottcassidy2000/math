---
source: opus-2026-07-31-S4 (the decisive step: the FC(3) residual lives on the fold = kps's rank-drop locus)
status: >
  DECISIVE REDUCTION. Proved: the max-modulus ridge of a C_3-covariant tilde-phi lies EXACTLY on the fold
  {J=0} (J = |tilde-phi_zeta|^2 - |tilde-phi_zetabar|^2), which IS kps's rank-drop discriminant locus
  {det J=0}. Algebra: at a critical point of |tilde-phi|^2 with tilde-phi != 0, nabla|tilde-phi|^2=0 gives
  phibar phi_zeta = -phi conj(phi_zetabar), so |phi_zeta|=|phi_zetabar|, i.e. J=0 (verified). Consequences:
  (i) GENERIC (Morse |tilde-phi|^2, isolated max) => saddle cap => nonzero leaks => NO counterexample, ALL
  degrees. (ii) The ONLY residual is a NON-Morse tilde-phi whose |tilde-phi| maxes on a critical CURVE Gamma
  ON THE FOLD {J=0}, with tilde-phi|_Gamma winding UNIFORMLY around |z|=M and CONSTANT transverse curvature
  (forced: C_3 makes the ridge density period-2pi/3 => only 3k Fourier modes; the leak kills all 3k modes =>
  the density is constant). So the FC(3) residual sits precisely on kps's rank-drop locus, where kps's
  integer-leak DISCRIMINANT and ideal-emptiness (Groebner=(1) at D=2, transversal rank=P through D=5) already
  operate. The remaining lemma is a single free-boundary/quadrature statement on that locus: no polynomial has
  a max-modulus critical curve mapping to a circle with constant transverse curvature. Geometric (mine) and
  algebraic (kps) caps now attack the SAME locus.
tags: [factorial-conjecture, FC3, fold, rank-drop, jacobian, max-ridge, non-morse, saddle-cap, kps-synthesis, quadrature-domain, decisive]
related: [fc3-saddle-cap-closes-isolated-max-all-degrees-synthesis-with-kps, kps-S161, fc3-cap-the-decisive-step-null-quadrature-and-the-edge-obstruction]
---

# The decisive step: the FC(3) residual IS the fold = kps's rank-drop locus

## 1. The max-modulus ridge lies on the fold {J=0} (rigorous)

Let `tilde-phi(zeta, zetabar)` be a `C_3`-covariant polynomial on the equilateral triangle `T`, and let
`J = |tilde-phi_zeta|^2 - |tilde-phi_zetabar|^2` be its real Jacobian (`{J=0}` = the fold). At a critical point
of `|tilde-phi|^2` where `tilde-phi != 0`:

```
   partial_zeta |tilde-phi|^2 = phibar * phi_zeta + phi * conj(phi_zetabar) = 0
   => phibar * phi_zeta = - phi * conj(phi_zetabar)
   => |phi| |phi_zeta| = |phi| |phi_zetabar|   =>   |phi_zeta| = |phi_zetabar|   =>   J = 0.
```

So **every max-modulus critical point (and hence any max-modulus ridge curve) lies on the fold `{J=0}`.**
Since the max `M = max_T |tilde-phi| > 0`, `tilde-phi != 0` there automatically. Verified numerically: for
`phi = z^2 + z̄/2`, the critical points with `phi != 0` satisfy `|phi_zeta| = |phi_zetabar|` (`J=0`); the ones
with `J != 0` are the ZEROS of `phi` (minima of `|phi|^2`, irrelevant). **`{J=0}` is exactly kps's rank-drop
discriminant locus `{det J=0}`** -- so my geometric residual and kps's algebraic locus are the SAME set.

## 2. The dichotomy, sharpened

- **GENERIC (Morse `|tilde-phi|^2`):** critical points isolated, so the max is at isolated points. The saddle
  cap (`fc3-saddle-cap-...`) gives `int_T tilde-phi^{3k} ~ (3/k) M^{3k} sum_j p_j e^{3ik theta_j} != 0` -- so no
  counterexample, at EVERY degree.
- **NON-Morse residual:** `|tilde-phi|^2` has a critical CURVE `Gamma` at its max, forced onto `{J=0}`. This is
  a degenerate, codimension-`>=1` condition (a curve of critical points requires `Re(F_zeta)` and `Im(F_zeta)`
  to share a real-curve factor).

## 3. The residual density must be CONSTANT (C_3 + leak)

On the critical ridge `Gamma`, `tilde-phi|_Gamma = M e^{i theta}` winding `n` times. The Laplace transverse to
`Gamma` gives leak `int_T tilde-phi^{3k} ~ M^{3k} oint_Gamma e^{3ik theta} tilde-rho(theta) dtheta`, where
`tilde-rho ~ 1/sqrt(transverse curvature kappa(theta))`. The `C_3` action rotates `theta -> theta - 2pi/3` and
preserves Lebesgue, so `tilde-rho` has PERIOD `2pi/3`, i.e. **only `3j` Fourier modes**. The leak vanishing
(`all k`) forces the `3k` modes to vanish for `k>=1`; combined, `tilde-rho == const`, i.e. the transverse
curvature `kappa` is CONSTANT along `Gamma`, and `tilde-phi` maps `Gamma` to the circle `|z|=M` with a
constant-density uniform winding (the disk boundary). So the residual is exactly the fold-maps-to-circle case.

## 4. The synthesis: geometric and algebraic caps hit the same locus

```
   saddle cap (opus)   : Morse => isolated max => nonzero leak. Kills the GENERIC case, ALL degrees.
   max-ridge = {J=0}    : the residual is a critical curve on the FOLD = kps's rank-drop discriminant locus.
   kps leak-ideal       : Groebner=(1) exact at D=2, transversal rank=P through D=5, ON that locus.
   edge obstruction     : the boundary (edge) version of the ridge can't be a circle.
```

The FC(3) counterexample, if any, is forced to be a NON-Morse `tilde-phi` whose max-modulus critical curve is
a fold mapping to a circle with constant transverse curvature -- a single quadrature-domain / free-boundary
configuration, sitting on kps's rank-drop locus, closed exactly at `D<=2` and transversally at `D<=5`.

## 5. The one remaining lemma (clean, sharp)

> **No non-constant polynomial `tilde-phi` has `|tilde-phi|` attaining its maximum on a critical curve `Gamma`
> that `tilde-phi` maps onto a circle `|z|=M` with constant transverse curvature.**

Equivalently: the fold `tilde-phi({J=0})` cannot be a circle carrying a Fourier-flat (`3k`-null) uniform
pushforward. This is a Schwarz-function / free-boundary rigidity on the fold. Two independent tools now bear on
it -- kps's integer-leak discriminant of `{det J=0}` (algebraic) and the constant-curvature/uniform-winding
obstruction (geometric). Proving it upgrades **bare FC(3) TRUE while GM(3) FALSE** (the `{e} < S_3 < SO(3)`
symmetry trichotomy) to a full theorem. Handed to kps as the joint endgame; the fold `{det J=0}` is the shared
object.
