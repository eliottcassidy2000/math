---
id: THM-5292
title: "Anatomy of THM-1300's Keller map F on C^3. (A) F = Psi o tau with Psi(s,w,c) = (s^2 + sw - cs^3, 4s + 2w - 3cs^2, c), det J Psi = 2w, and tau(x,y,z) = ((1+xy)/x, 1/x, x l), l = 2 - 3xy - x^2 z, det J tau = -1/w: the zero of the fold family cancels the pole of the inversion. (B) Psi is the root cover of the cubic c S^3 - 2S^2 + bS - 2a; every Psi_r = (sw + r - s r_s, 2w - 2r_s, c) has Jacobian 2w. (C) Level surfaces of F_3 = x l: cylinders carrying the planar map Phi_c of Jacobian 2w minus its critical line; the plane x = 0 carries an automorphism. (D) JC_2 is equivalent to injectivity of spatial Keller maps on plane components of coordinate level sets. (E) The planar ansatz psi_r o (gamma x^3 y + phi(x), 1/x) is never polynomial"
status: >
  PROVED (exact symbolic computation, scripts/jacobian_levels.py): (A), (B), (C), (E); (D) is an elementary
  equivalence. Re-verified here: det JF = -2 and F(-2,1/3,-2) = F(-4,1/3,0) = F(6,-1/6,1/4) = (0,1/3,-24).
  The map and its non-injectivity are THM-1300 (owner-supplied; verified there); the fibre cubic is THM-4621; the
  tangent-sweep families are THM-4625; the torus reduction is THM-4727; the planar no-mate theorem for the
  factorized coordinate is THM-2045. Added here: the factorisation, the level-surface statement, (D), (E).
  NOT CLAIMED: anything about JC_2, which remains OPEN; (E) excludes one ansatz only.
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-tapes
date: 2026-11-15
source: 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/JACOBIAN.md
scripts:
  - 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/scripts/jacobian_levels.py
related:
  - THM-1300, THM-4621, THM-4625, THM-4727, THM-2045, THM-2722, HYP-9832, HYP-9343
---

# THM-5292 - the dimension-3 Keller map as a fold family after an inversion

Source: [the note](../../05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/JACOBIAN.md).

`F = (u^3 z + y^2 u(4 + 3xy), y + 3x u^2 z + 3x y^2 (4 + 3xy), 2x - 3x^2 y - x^3 z)`, `u = 1 + xy`,
`l = 2 - 3xy - x^2 z`, `F_3 = x l`.

**(A)** `F = Psi o tau`, `Psi(s, w, c) = (s^2 + sw - c s^3, 4s + 2w - 3c s^2, c)`, `tau = (u/x, 1/x, x l)`;
`det J Psi = 2w`, `det J tau = -x`; `tau` is an isomorphism `{x != 0} -> {w != 0}` with inverse
`(1/w, s - w, 5w^2 - 3sw - c w^3)`.

**(B)** On the image of `Psi`: `c s^3 - 2s^2 + b s - 2a = 0`, `2w = b - 4s + 3c s^2`. For every polynomial
`r(s, c)` the map `(sw + r - s r_s, 2w - 2 r_s, c)` has Jacobian `2w`; `Psi` is `r = (c s^3 - 2s^2)/2`.

**(C)** For `c != 0`, `{F_3 = c}` is `{(s, w) : w != 0}` and `(F_1, F_2)` restricts to
`Phi_c = (s^2 + sw - c s^3, 4s + 2w - 3c s^2)`. `{F_3 = 0} = {x = 0} + {l = 0}` (disjoint); on `x = 0`,
`F = (z + 4y^2, y, 0)`; on `l = 0`, `F = ((v+1)(v+2)/x^2, 2(2v+3)/x, 0)`, `v = xy`.

**(D)** JC_2 holds iff every Keller map of `C^3` is injective on every component, isomorphic to `C^2`, of a level
set of one of its coordinates.

**(E)** Let `r in C[s]`, `psi_r = (sw + r - s r', 2w - 2r')`, `tau_2 = (gamma x^3 y + phi(x), 1/x)` with `phi` a
Laurent polynomial (the shape forced by `det J tau_2 = gamma/w`). Then `psi_r o tau_2` is not a polynomial map:
the second component forces `deg r = 2`, `phi = 1/(rho x) + polynomial`, and then the first component is
`1/(2 rho x^2) + polynomial`.
