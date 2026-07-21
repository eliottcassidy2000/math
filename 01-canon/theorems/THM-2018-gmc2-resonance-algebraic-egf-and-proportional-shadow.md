---
id: THM-2018
title: "GMC(2) resonance closures from an algebraic EGF and a proportional radial shadow"
status: >
  RESERVED by codex-2026-07-21-NC2-followup. Two candidate proved parts are under
  independent coefficient and analytic review before canonization. Part A concerns
  p=q=1, deg b<=1 and deg(ac)<=1; Part B concerns the arbitrary-degree hypersurface
  s*a*c=kappa*b^2. The exact identities are known; the remaining work at reservation
  time is hostile review, exact-script verification, and integration into HYP-8766.
source: codex-2026-07-21-NC2-followup
depends_on:
  - THM-1570
  - THM-2014
  - THM-2017
related:
  - HYP-8766
  - HYP-8769
---

# THM-2018 — reserved resonance closures

This file reserves the next theorem identifier while two exact mechanisms are checked.

1. For
   `P=Z*a(s)+(b0+b1*s)+Zbar*c(s)` with
   `a(s)c(s)=delta+alpha*s`, the candidate exponential generating function is
   algebraic-exponential. If all moments vanish, it forces an algebraic function and
   its exponential both to be algebraic; the compact-curve pole/essential-singularity
   argument should force that exponent to be constant and hence all two-sided data to
   degenerate.
2. On the arbitrary-degree central hypersurface
   `h=s*a*c=kappa*b^2`, every channel factors as
   `S_m(kappa)*b^m`. EMP controls `L(b^m)`. **Audit correction:** the elementary
   Legendre recurrence only forbids consecutive zeros; that does not imply the
   eventual radial vanishing needed by EMP. THM-2021 supplies the correct bridge:
   finite recurrence of a fixed nonzero Legendre root, with a self-contained
   nonexceptional-phase theorem and a clearly labelled use of the announced 2026
   Mangoubi--Kadets--Weller Weiser finiteness result.

No use should be made of this reserved stub until the status and full proofs are
installed.
