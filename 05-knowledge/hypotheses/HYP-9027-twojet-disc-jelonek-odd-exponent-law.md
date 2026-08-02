---
id: HYP-9027
title: "Odd-exponent Jelonek law: general form refuted, Keller successor open"
status: >
  REFUTED IN THE STATED GENERAL-DOMINANT SCOPE + KELLER-RESTRICTED SUCCESSOR
  OPEN.  THM-3059 gives an explicit dominant field-degree-four two-jet map
  with generic S4, exact Jelonek plane, C3 infinity inertia, and even cleared
  exponent 8.  The exact survivor is parity = inertia sign after reciprocal
  reversal.  No field-degree-four two-jet Keller map is known; whether every
  Jelonek component of such a Keller map must have odd inertia remains open.
source: kind-pasteur-2026-07-26-S132
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
  - THM-1310-conic-pair-fibers-and-design-equations
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
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

> **SUPERSEDED SCOPE.**  THM-3059 refutes the general-dominant odd-exponent
> clause below.  Retain this file for the evidence and correction lineage.
> The live successor asks whether the Keller condition forces odd infinity
> inertia at every Jelonek component; it does not assume that conclusion.

Corrected successor to THM-2446 (P1), with THM-2455 supplying the
proved classical layer. Let `F = A z^2 + B z + C : C^3 -> C^3` be
dominant with generic fiber field degree 4; `N(t)` the primitive
graph polynomial of a separating coordinate, leading coefficient
`ell`; `q` its monic depressed form with invariants `I, J`
(so `27 disc q = 4 I^3 - J^2`, THM-2455); `L_4` the reduced Jelonek
polynomial.

## The historical conjectured law

1. Every irreducible component of `{L_4 = 0}` divides `ell` -- the
   Jelonek factor is leading-coefficient-driven (fiber degree drop /
   escape to infinity), not a ramification divisor.
2. **REFUTED in this scope.** `v_Lambda(disc q)` was conjectured to be odd
   for every component `Lambda` of
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
| IV (THM-3059): (x, xz^2+y, xyz^2+z) | lc = u^2 | u^4 H, H=-27 mod u | {u=0} | 4 | H (S4 Galois) |

Master identities: `u^9 disc(monic II) = Delta`;
`u^5 disc(monic III) = 16 w (v^2 - 4uw)^2`; d=3:
`disc(monic) = -4 Q^2 / L^3`. The original `N` values `9, 5, 3` are all odd,
but THM-3059's value is `12-4=8`, so that sample pattern is not uniform.
Monodromy realizations certified in the class: `V_4` (proper),
`D_4` (non-proper), `S_4` (non-proper). Coordinate-artifact squares
are real (the witness's `v^2` covers two distinct unramified
points) and must be certified away per map (controls file, exact
unramified-fiber checks).

## Open obligations (from the program's validity gate)

- G1 (blocking): exhibit a field-degree-4 2-jet Keller map or prove
  none exists. STATUS UPDATE (THM-2465, S134): no witness found;
  G1 pinned by a proved package (total degree >= 3 forced; universal
  z-rationality; conic-cap emptied in four QQ-exact boxes; degrees
  {1,3} realized on point/line caps, so unconditional degree-4
  exclusion there implies the order-{1,3} degree-4 case -- a
  hardness floor); any witness must satisfy the (N1)-(N5) resolvent
  dossier.  Purity and the sign character force an odd Jelonek component
  somewhere on an S4 branch; they do **not** force every Jelonek component
  to be odd.  Decisive pending: the (1,2,2)+fiber4 GB box.
- G5 is repaired by THM-3059: after reciprocal reversal, the exact formula is
  `N=6v(ell)-(d_sigma+2i)`, hence `N mod 2` is the sign of infinity inertia.
  The remaining Keller task is geometric: exclude C3 infinity inertia at a
  Jelonek component, or prove that it forces a forbidden finite critical
  divisor.  THM-3059's family shows that all Jacobian valuations can remain
  zero; the missing local sidecar is the branchwise unit residue of the
  primitive-element Jacobian cofactor.
- G6 (decoupling search): a proper 2-jet map with
  Galois not contained in A_4 would decouple "square disc" from
  "proper".  The other requested decoupling now exists: THM-3059 is a
  nonproper generic-S4 map with even Jelonek exponent, so clause 2 is
  refuted for dominant maps.
- The invariant-form prize: a family identity expressing
  `mu L_4^N` through `4 I^3 - J^2` along a modulus path that
  genuinely crosses the cuspidal edge `{I = J = 0}` -- THM-2455
  proves the Chebyshev quadrisection path is disqualified.
