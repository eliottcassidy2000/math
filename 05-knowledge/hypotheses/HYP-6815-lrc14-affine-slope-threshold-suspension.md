---
id: HYP-6815
title: LRC14 affine-slope threshold suspension
status: EXACT REPRESENTATION CLAIMED; preservation audit in progress
source: codex-2026-07-14-S2
related:
  - HYP-6780
  - HYP-6755
  - HYP-3106
  - HYP-3072
  - HYP-3025
  - THM-755
---

# HYP-6815: LRC14 Affine-Slope Threshold Suspension

This namespace is reserved for an exact reparameterization of affine LRC14
families and for the information-preservation audit it suggests.

Let `P=(p_i)` and `R=(r_i)` be integer 13-vectors and let

`V(c)=(c p_i+r_i)`,

on the positive, distinct, primitive integer scales `c` under study.  On the
two-torus set

`Phi_{P,R}(u,t)=min_i ||p_i u+r_i t||`.

For the slope-`c` closed geodesic `L_c={(ct,t):t in R/Z}`, direct substitution
gives the exact identity

`M(V(c)) = max_t Phi_{P,R}(ct,t)`.

Thus the mixed four-coordinate incidence object

`X_{P,R}={(u,t,c,lambda): u=ct, Phi_{P,R}(u,t)>=lambda}`

in `(R/Z)^2 x N x [0,1/2]` has the LRC14 assertion as nonemptiness of its
integer-slope fiber at `lambda=1/14`.  It is a stratified continuous/discrete
object, not a smooth four-manifold.  Runner owners label its strip strata.

When `R=0`, `Phi` is independent of `t`; the suspension is cylindrical and
recovers HYP-6780's dilation law.  Nonzero `R` is transverse shear/holonomy,
not a small perturbation that may be discarded.

Still missing before promotion:

1. an exact executable audit on named affine families;
2. a ledger of which quotients preserve fiber nonemptiness, exact clearance,
   peel/deletion functoriality, and endpoint ownership;
3. a dual Fourier statement identifying the slope resonance line
   `m*c+n=0` and separating finite trigonometric identities from conditional
   indicator expansions;
4. an explicit comparison with closed arc-Cech, CRT/p-adic, interaction-code,
   observer-cut, and certificate-sheaf carriers;
5. Tournament Analysis on carrier choices, with runners and raw arcs included
   as challenged vertex sets rather than assumed ones.
