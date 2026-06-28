---
id: HYP-3311
title: LRC14 CRT/Galois sidecar audit for unit residue skeleton and 2-adic covering magnitude
status: EVIDENCE / exact CRT and Galois-sidecar audit; proof-target refinement, not an LRC14 proof
source: mac-mini-2026-06-28-S84 + codex-2026-06-28 extension
tangent: T1361
technique: LTI-361
tournament_technique: LTT-261
script: 04-computation/lrc_census_crt_factorization_macmini_S84.py
result: 05-knowledge/results/lrc_census_crt_factorization_macmini_S84.out
reflection: 07-reflections/the-census-factors-via-crt-7-adic-residue-c3-skeleton-times-2-adic-magnitude-doubling-hinge.md
related:
  - HYP-3400
  - HYP-3310
  - HYP-3301
  - HYP-3265
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3256
  - HYP-3255
  - HYP-3254
  - HYP-3253
  - HYP-3250
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3246
  - HYP-3243
  - HYP-3300
  - HYP-2909
  - HYP-3087
  - THM-523
  - OPEN-Q-108
---

# HYP-3311: LRC14 CRT/Galois Sidecar Audit

## Claim

This is an exact sidecar audit for the HYP-3310 C6 residue-magnitude frame.
The current LRC14 census picture has a two-prime decomposition:

```text
{1,...,13} = U union 2U union {7}
U  = (Z/14Z)^* = {1,3,5,9,11,13}
2U = {2,4,6,8,10,12}
```

The binding runners are the units.  The covering runners are the even shadow
`2U` plus the ramified singleton `7`.  The map `u -> 2u mod 14` is a bijection
from binding units to even covering classes, while `7=(1 mod 2, 0 mod 7)` is
the apex: odd, but ramified at the 7-adic coordinate, so it belongs to the
covering layer.

This is not a proof of LRC14.  It is a better proof packet.  It says the
residue skeleton, field sidecars, and height/magnitude layer must be kept
separate until a finite chamber, sheaf-exactness, or covering-flex theorem
discharges them.  In the HYP-3301 language, the first forgotten coordinate
among `C3`, `Q(sqrt(-7))`, and height/flex should become an explicit
obstruction cocycle or boundary-transfer kernel, not an invisible scalar
quotient.

HYP-3400 supplies the global conservation language for the same rule:
the C3 unit skeleton, the quadratic `Q(sqrt(-7))` sidecar, and the 2-adic
height/flex ledger are three shadow-charge reservoirs.  A quotient is legal
only if it preserves the charge, transfers it to a dual shadow, or emits a
named debt.

## Exact Galois Readout

Inside `Q(zeta_7)`,

```text
Gal(Q(zeta_7)/Q) = (Z/7Z)^* = C6 = C2 x C3.
```

The C3 quotient is the real-cubic binding-pair skeleton:

```text
(1,13) -> (3,11) -> (5,9) -> (1,13)
```

This is the S83/HYP-3259 orbit and matches the rank-3 unit-nullspace readout
of HYP-3257.  It is the right home for the HYP-2909 single binding-pair seed:
prove one orbit point, then move equivariantly through the C3 action.

The quadratic `Q(sqrt(-7))` character is not the same coordinate.  It is the
Legendre character mod `7`, fixed by the quadratic subfield of the cyclotomic
field.  The executable audit verifies:

```text
each_binding_pair_has_one_QR_and_one_NQR = True
```

Thus every C3 antipodal binding pair contains one quadratic-residue and one
quadratic-nonresidue element.  The C3 binding quotient and the
`Q(sqrt(-7))` character are transverse sidecars.  Neither can legally replace
the other unless the downstream LRC predicate is proved constant on the
forgotten coordinate.

## LRC Transfer

The proof should now split into two coupled rigidity claims.

1. Unit-contact rigidity:

   Rows retaining all six unit contacts should force the C3 unit skeleton,
   with HYP-2909 supplying the seed and C3 propagation supplying the other
   two antipodal binding pairs.  This links HYP-3257's unit nullspace and
   HYP-3265's contact graph.

2. Covering-flex rigidity:

   The nonunit layer `2U+{7}` must be controlled by a height/magnitude ledger.
   The hinge `12 -> 24` is a 2-adic doubling lift: it contains the old
   12-kill lattice and adds interleaved 24-kills.  Same-residue height moves
   such as the earlier `2 -> 16` decoy show why residue equality is not proof
   currency by itself.

The desired theorem is therefore:

```text
unit-contact skeleton survives
  => C3/HYP-2909 binding rigidity
unit contacts are killed or leave the global maximum
  => HYP-3265/HYP-3300/HYP-3301 off-unit covering chamber,
     first-obstruction exactness, strict open witness, Phi14d dilation witness,
     finite trap discharge, or named residual debt
covering layer survives as a boundary equality
  => only AP/Goddyn-Wong integer tight lifts, including the 12->24 hinge
```

This is a labelled-packet theorem target, not a scalar extremality claim.

## Tournament Analysis

The associated script uses proof carriers, not runners, as tournament
vertices.

```text
pairwise observable:
  preserved LRC predicate
  binding detection
  covering detection
  field-lattice payload
  rigidity target
  proof readiness

binary gauge:
  A -> B iff the weighted proof-carrier score of A is larger than B;
  ties use the lexicographic Hamiltonian path.
```

Exact fingerprint from the stored output:

```text
vertices = 9
score_hist = {29:1,42:1,50:1,54:1,58:1,60:1,61:1,62:1,73:1}
directed_3cycles = 0
hamiltonian_path_count = 1
selected_path =
  crt_packet_U_plus_2U_plus_ramified7
  -> c6_subfield_lattice_Qzeta7
  -> c3_real_cubic_binding_pair_orbit
  -> apex7_ramified_covering_switch
  -> two_adic_unit_to_even_shadow_2U
  -> quadratic_Qsqrt_minus7_transverse_character
  -> hinge_12_to_24_height_lift
  -> raw_mod14_residue_table
  -> slogan_binding_units_covering_evens
```

The ranking demotes the slogan.  The slogan remembers a true pattern, but it
forgets the transverse quadratic character, the same-residue height ledger,
and the finite chamber exits.

## Assumption Challenge

Do not assume the tournament vertices are runners or arcs.  This packet
explicitly considered alternate vertex sets:

```text
CRT classes
unit pairs
quadratic QR/NQR sidecars
subfields of Q(zeta_7)
covering flex directions
contact-graph cases
Morse chamber exits
proof obligations
```

Preserved LRC predicate: whether a speed or carrier belongs to the binding
unit skeleton, the even covering shadow, the ramified apex, or a proof exit
that can still certify `M(S) >= 1/14`.

Destroyed information: actual safe-set geometry, endpoint owners, off-unit
component topology, height/lift data, and finite trap status unless those
sidecars are retained.

Challenged assumption: the partition "units bind, evens cover" is enough.
It is not.  The C3 quotient, the `Q(sqrt(-7))` character, and the 2-adic
height/flex ledger must remain visible until HYP-3265/HYP-3300/HYP-3301-style
contact, chamber, or sheaf-exactness logic discharges them.

## Next Proof Obligations

A. Formalize unit-contact rigidity: tight rows retaining all six unit contacts
force the C3 unit skeleton.

B. Use the transverse `Q(sqrt(-7))` character as a sidecar, not a replacement
for the C3 quotient.

C. Prove the covering layer `2U+{7}` has a one-dimensional height flex and
that its only integer tight lift is the AP/Goddyn-Wong `12 -> 24` hinge.

D. Combine with HYP-3265 and HYP-3300: killed unit contacts must route to
off-unit covering, Morse chamber witnesses, strict open cells, finite
Toeplitz/Green/root-motion discharge, `Phi14d` dilation equality, or named
residual debt.

E. Instantiate HYP-3301 on this three-coordinate packet: prove that forgetting
the C3 skeleton, the quadratic sidecar, or the height/flex ledger produces an
exact cocycle, a zeta_7 holonomy repair, a positive boundary-moment image, a
forbidden AP/Goddyn-Wong kernel, or named K33/H7 debt.

F. Price the same forgetting through HYP-3400: each proposed scalar shadow of
this packet must preserve, transfer, or explicitly debt the C3, quadratic, and
height/flex charges.
