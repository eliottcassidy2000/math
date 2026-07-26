---
id: HYP-9027
title: "The odd-exponent Jelonek law for 2-jet fiber quartics (P1-prime)"
status: >
  OPEN (three-map exact evidence base + d=3 control; the d=3-style
  readings are REFUTED: mu need not be a unit, and odd disc
  components are not exclusively Jelonek). Blocking hypothesis G1:
  no field-degree-4 2-jet Keller map is known, so the law is filed
  for general dominant 2-jet maps, not Keller maps.
source: kind-pasteur-2026-07-26-S132
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
  - THM-1310-conic-pair-fibers-and-design-equations
script: 04-computation/jacobian_twojet_disc_jelonek_maps_kps_S132.py
output: 05-knowledge/results/jacobian_twojet_disc_jelonek_maps_kps_S132.out
controls_script: 04-computation/jacobian_twojet_disc_jelonek_controls_kps_S132.py
controls_output: 05-knowledge/results/jacobian_twojet_disc_jelonek_controls_kps_S132.out
script_sha256: d1175b2c0e6140ebd09859d47480bf7373bc82c06143a4866ea811f108a31ba4
output_sha256: 615a94e4dc723a1fabc1cebcaf383ef67a0b58d81acc1f6f374a27d74b552be5
controls_script_sha256: fd778bd076e5cdba9ac38b1fec7cbcedc6cd5047a16ae48b23adee4d743ae218
controls_output_sha256: b4facd9293a31765f220b6903d5e74a1d9fef4e5acb8764f65417991c6b75854
---

# HYP-9027 -- P1-prime: Jelonek components enter the fiber discriminant to odd order

Corrected successor to THM-2446 (P1), with THM-2455 supplying the
proved classical layer. Let `F = A z^2 + B z + C : C^3 -> C^3` be
dominant with generic fiber field degree 4; `N(t)` the primitive
graph polynomial of a separating coordinate, leading coefficient
`ell`; `q` its monic depressed form with invariants `I, J`
(so `27 disc q = 4 I^3 - J^2`, THM-2455); `L_4` the reduced Jelonek
polynomial.

## The conjectured law

1. Every irreducible component of `{L_4 = 0}` divides `ell` -- the
   Jelonek factor is leading-coefficient-driven (fiber degree drop /
   escape to infinity), not a ramification divisor.
2. `v_Lambda(disc q)` is **odd** for every component `Lambda` of
   the Jelonek set, with
   `N = (2*4-2) deg_Lambda(ell) - v_Lambda(disc N)` the exact odd
   bookkeeping; equivalently
   `L_4^N * disc q = mu * S^2` with `N` odd, `S` polynomial, and
   `mu` = the product of odd-order components of the honest branch
   divisor.
3. `mu` need NOT be a unit (refuting the naive d=3 transfer), and
   the square class of `disc q` computes the `A_4` Galois test but
   does NOT locate the Jelonek set (THM-2455's resolvent equality
   kills every zero-locus version of "resolvent conductor").

## Evidence (all exact, hard-checked)

| map | fiber poly | disc | Jelonek | k | mu |
|---|---|---|---|---|---|
| d=3 control (THM-1310) | cubic, lc = L | -4 Q^2 L | {L=0} | 1 | unit |
| (x^2-z^2, xz, y) proper | monic | -16 v^2 (u^2+4v^2)^2 | empty | 0 | unit |
| I: (x^2-z^2+y, xz, y) proper | monic | -16 v^2 ((u-w)^2+4v^2)^2 | empty | 0 | unit |
| II: (x, xz^2+y, yz+y^2) | lc = u^2 | u^3 Delta, Delta irred | {u=0} | 3 | Delta (S4 Galois) |
| III: (x, xz^2+y, yz^2) | lc = u | 16 u w (v^2-4uw)^2 | {u=0} | 1 | 16 w (D4 Galois) |

Master identities: `u^9 disc(monic II) = Delta`;
`u^5 disc(monic III) = 16 w (v^2 - 4uw)^2`; d=3:
`disc(monic) = -4 Q^2 / L^3`. The `N` values `9, 5, 3` are all odd.
Monodromy realizations certified in the class: `V_4` (proper),
`D_4` (non-proper), `S_4` (non-proper). Coordinate-artifact squares
are real (the witness's `v^2` covers two distinct unramified
points) and must be certified away per map (controls file, exact
unramified-fiber checks).

## Open obligations (from the program's validity gate)

- G1 (blocking): exhibit a field-degree-4 2-jet Keller map or prove
  none exists; both known 2-jet Keller examples are tame with field
  degree 1 and the THM-1310 wild map has degree 3.
- G5: prove the odd-exponent mechanism from lc-degeneration parity.
- G6 (decoupling search, floors reported): a proper 2-jet map with
  Galois not contained in A_4 would decouple "square disc" from
  "proper"; a non-proper map with even Jelonek exponent would
  refute clause 2. Neither was searched for yet.
- The invariant-form prize: a family identity expressing
  `mu L_4^N` through `4 I^3 - J^2` along a modulus path that
  genuinely crosses the cuspidal edge `{I = J = 0}` -- THM-2455
  proves the Chebyshev quadrisection path is disqualified.
