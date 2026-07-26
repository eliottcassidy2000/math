---
id: THM-2455
title: "The d=4 swallowtail scaffold: master invariant identity, resolvent equality, and endpoint corrections"
status: >
  PROVED (exact classical layer for the 2-jet d=4 program: the
  depressed-quartic discriminant satisfies 27 Delta_4 = 4 I^3 - J^2
  with I = p^2+12r, J = 2p^3+27q^2-72pr; all three standard
  resolvent cubics have disc EQUAL to Delta_4 identically, unit 1,
  and the square class is convention-invariant; the swallowtail
  {Delta_4=0} stratifies with A2 cuspidal edge exactly {I=J=0}
  ideal-theoretically, A1A1 curve {q=0, p^2=4r}, A3 = origin only;
  singular locus certified both directions) + REFUTED (the
  "A_3 at both ends" reading of THM-2446 (P1): the Chebyshev
  quadrisection pencil disc_T(8T^4-8T^2+1-m) = 2^17 (1-m)(1+m)^2
  provably NEVER meets the cuspidal edge or the A3 vertex; endpoint
  types are A1 at m=+1 and Z/2-symmetric A1A1 at m=-1 through the
  inner double cover; even the d=3 trisection endpoints are A1
  fiber collisions) + VERIFIED-EXACT (frame identities with signs
  pinned: on the fiber z (A x B) = A x (P-C) and
  z^2 (A x B) = B x (C-P); g2 weights w_l = (A x B)_l). The theorem
  does not decide any Jelonek factorization law (HYP-9027) and does
  not exhibit a field-degree-4 2-jet Keller map.
source: kind-pasteur-2026-07-26-S132
depends_on:
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
related:
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - THM-1310-conic-pair-fibers-and-design-equations
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/jacobian_d4_swallowtail_scaffold_kps_S132.py
output: 05-knowledge/results/jacobian_d4_swallowtail_scaffold_kps_S132.out
script_sha256: ec546ed691855050ed68f701f0868104e851272d8ce0bd96993c4a1a4cc5693b
output_sha256: 305ce12fa8745b2a72950c104ec92cf65d8b056d4770f9ddf456528e95b327fd
hash_basis: working-tree bytes (LF)
---

# THM-2455 -- the quartic's master identity, and two retired slogans

**PROVED + REFUTED + VERIFIED-EXACT** as itemized in the status.

This is the classical layer THM-2446's prediction (P1) needs, plus
the corrections it forces. All identities exact (sympy, hard
assertions, 46 checks + jet supplements).

## 1. The master invariant identity (PROVED)

For the depressed quartic `f(z) = z^4 + p z^2 + q z + r`, with the
weight-(4,6) invariants `I = p^2 + 12 r`, `J = 2p^3 + 27q^2 - 72pr`:

```text
27 * Delta_4 = 4 I^3 - J^2,                                     (1)
Delta_4 = 16p^4 r - 4p^3 q^2 - 128 p^2 r^2 + 144 p q^2 r
          - 27 q^4 + 256 r^3.
```

This is the d=4 occupant of the slot THM-1335's
`108 a^2 L = (12a - b^2)^3 + E^2` filled at d=3: discriminant =
(cube of an invariant) minus (square of an invariant), up to the
constant. The depressed resolvent is `w^3 - (I/3) w - J/27`.

## 2. Resolvent equality is exact (PROVED)

All three standard resolvent cubics (`g1` with roots `-(z_i+z_j)^2`,
`g2` with roots `(z_i+z_j)^2`, `R2` with roots `z_i z_j + z_k z_l`)
satisfy `disc = Delta_4` **identically, unit 1, no square
correction**; likewise for the general monic quartic and its cubic
resolvent over `Z[A,B,C,D]`. Mechanism: under `e_1 = 0`,
`(z_1+z_2)^2 - (z_1+z_3)^2 = (z_1-z_4)(z_2-z_3)`, so resolvent root
differences are a permutation of the products `(z_i-z_j)(z_k-z_l)`.
Renormalizations move disc by even powers, so the square class is
convention-invariant. Consequence for any "resolvent conductor"
argument: **the resolvent's discriminant zero-locus can never
separate anything from the quartic's** -- only odd-order valuations
can, and those are addressed in HYP-9027.

## 3. Swallowtail stratification (PROVED)

`Sing{Delta_4 = 0}` = (A2 cuspidal edge) union (A1A1 curve),
set-theoretically certified both directions; the edge is
parametrized `(p,q,r) = (-6s^2, 8s^3, -3s^4)` and equals `{I = J = 0}`
**as an ideal**; the A1A1 curve is `{q = 0, p^2 = 4r}`; the A3
stratum is the origin alone (tangent cone `256 r^3`). Transversal
types certified by exact jets (cusp `v^2 ~ w^3` on the edge, normal
crossings on the curve).

## 4. Two slogans retired (REFUTED)

The quadrisection pencil: `disc_T(8T^4 - 8T^2 + 1 - m) =
2^17 (1-m)(1+m)^2`. The path `(p,q,r) = (-1, 0, (1-m)/8)` has
`8p^3 + 27q^2 = -8` identically, so it **never meets the cuspidal
edge** `{I = J = 0}` and never meets the A3 vertex (`p = -1`).
Endpoints: `m = +1` is a single A1 collision; `m = -1` is
`2(2T^2-1)^2`, a Z/2-symmetric A1A1 pair factoring through the
inner double cover `u = T^2` (whose own disc is `32(1+m)`).

- Retired: "A_3 at both ends" for the d=4 modulus path. False; no
  member of the pencil is more degenerate than A1A1. Any hunting
  family for the master identity must be chosen to actually cross
  the cuspidal edge.
- Corrected reading of d=3: the trisection endpoints
  (`disc_T(T^3 - 3T - 2m) = -108(m^2-1)` in the `T = 2 cos`
  normalization) are A1 fiber collisions; THM-1335's "A_2 both
  ends" survives only as target/caustic language, not as fiber
  root-collision type.

## 5. Frame identities (VERIFIED-EXACT)

For `F = A z^2 + B z + C`, `E := A z^2 + B z + (C - P)`: on the
fiber, `z (A x B) = A x (P - C)` and `z^2 (A x B) = B x (C - P)`
(sign-flipped variants fail); `g1 = [A,B,C] - (A x B).P` vanishes on
the fiber and is affine-linear in `P`; the quadratic consequences
carry weights `w_l = (A x B)_l`, with the full bilinear family
`a_l a_m = u_l b_m` (`u = A x B`, `a = A x (P-C)`, `b = B x (C-P)`).

## 6. Scope

Nothing here proves a Jelonek factorization law (that is HYP-9027's
conjecture with its three-map evidence), and no field-degree-4
2-jet Keller map is exhibited -- the known 2-jet Keller examples
have field degree 1, so the Keller-restricted (P1) currently has an
empty verified domain (blocking hypothesis G1 of the program).

## 7. Reproduction

```bash
python 04-computation/jacobian_d4_swallowtail_scaffold_kps_S132.py
```

Exact symbolic; prints `FAILED CHECKS: NONE`.
