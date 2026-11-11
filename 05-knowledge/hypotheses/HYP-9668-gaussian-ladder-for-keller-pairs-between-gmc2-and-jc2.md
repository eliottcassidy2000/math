---
id: HYP-9668
title: "A Gaussian ladder for Keller pairs between GMC(2) and JC(2). For a Keller pair (f, g) with Stein kernel M_F (THM-5138): (A) K-GMC: M_F is a Mathieu-Zhao subspace of C[x,y]. (B) CV, cone vanishing: for every h, [h U^k V^(m+k)] = [h U^(m+k) V^k] = 0 for m large and 0 <= k <= m. (C) SN, one-sided Stein nilpotence: for every h, d_u^m h and d_v^m h lie in M_F for m large. Each is implied by JC(2) (A => B => C are proved, and JC(2) => A by GMC(2)); each is a theorem for automorphisms. (D) Detection: SN at a generic centre forces the pair to be an automorphism; then A, B, C and JC(2) are equivalent. (E) Rank-one potentials: if H^2(C^2, d - dW ^) is one-dimensional, the kernel of its functional is a Mathieu-Zhao subspace."
status: >
  OPEN (all five parts). PRE-AUDIT DRAFT.
  (A), (B), (C): no Keller pair other than an automorphism is known, so there is no direct evidence on C^2. They
  hold for automorphisms (THM-2022 and THM-5138 (A3), (D)). The implications JC(2) => (A) => (B) => (C) are PROVED
  (THM-5138 (D)).
  (D): PROVED only in a two-sided form: finite chaos (all classes of d_u^m d_v^n h vanish for m + n large) forces an
  automorphism (THM-5138 (C)). FINITE-EXACT evidence for the one-sided form on other surfaces: on C^* x C with
  (u, v) = (x^2, y/x) SN fails in the direction d_u (proved by residues) and holds in the direction d_v, along which
  no sheet is missing; on the slice F_1 = 1 of the map of THM-1300 it fails in both directions (classes of w V^m,
  w U^m nonzero for m <= 5). A heuristic derivation is in section 7 of the source note; it is not a proof.
  (E): GMC(2) is the case W = xy. For W = x + x^2 y and W = x^2 y the space is one-dimensional (FINITE-EXACT,
  THM-5138 controls); nothing is known about the Mathieu property there.
  None of (A)-(C) is known to imply JC(2), and (A)-(C) together with (D) are equivalent to JC(2), so this is a
  reformulation with intermediate rungs, not a reduction.
session: kind-pasteur-2026-11-11-S24
source: 05-knowledge/results/kps_s24_stein_modules_keller_20261111/README.md
related:
  - 01-canon/theorems/THM-5138-stein-modules-of-keller-pairs-dimension-count-and-finite-chaos-criterion.md
  - 01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md
---

Notation as in THM-5138: `U = f - a`, `V = g - b`, dual fields `d_u, d_v`, `M_F = Im(d_u - V) + Im(d_v - U)`,
classes `[h]` in `S_F = C[x,y] / M_F`.

## The ladder

    GMC(2)  =  K-GMC for automorphisms                              (THM-2022, a theorem)
    JC(2)  <=>  finite chaos  <=>  dim S_F = 1 for every Keller pair  (THM-5138 (B), (C))
    JC(2)   =>  (A) K-GMC   =>  (B) CV   =>  (C) SN                   (THM-5138 (D))
    (D) detection:  (C) => JC(2)                                      (open)

(A) is the statement of GMC(2) with the Gaussian functional replaced by all Stein functionals of an arbitrary
Keller pair. A proof of GMC(2) that only used integration by parts in the coordinates, and not the monomial basis,
would prove (A). The proof in canon (charges, lowest balanced face, Frobenius) uses the grading by charge, which a
Keller pair does not have unless the rotation field `U d_u - V d_v` is locally finite.

(B) is what the Mathieu property gives uniformly in a one-parameter family of null polynomials `V (1 + c UV)`.
(C) is the case `k = 0`.

## Why detection is plausible, and what would refute it

A counterexample of generic degree `d` has at least `2d - 1` Stein functionals (THM-5138 (B)), of which at most
`d` can be local Gaussians at the preimages of the centre; the others come with the missing sheets. In the two
examples on other surfaces these extra functionals are exactly what breaks (C), and the mechanism (moments of a
function with a singularity at finite distance along the non-properness curve) does not involve the degree of the
map. A refutation of (D) would be an etale map from a smooth affine surface, non-proper along a curve met by both
families of coordinate lines, for which SN holds at a generic centre.

## What a proof of (D) would need

For the direction `d_u`: the module `C[x,y] / T_u C[x,y]` is holonomic on the `v`-line with generic rank
`sum_i delta_i deg_u S_i >= d - 1`; SN says that all of its de Rham cohomology is carried by a quotient supported at
`v = b`. An Euler-Poincare count with the irregularity at infinity of the exponential factors
`e^(-(s(v) - a)(v - b))` attached to the branches `u = s(v)` of the non-properness curve should exclude this for
generic `a`. The formal structure at infinity of this module has not been established here.

## (E) Rank-one potentials

`dim H^2(C^2, d - dW ^) = 1` iff the generic fibre of `W` has `b_1 = 1`. Besides `xy` and its images under
automorphisms this happens for `x^p y^q` and for Broughton's polynomial `x + x^2 y`, which has no critical point. For
`x^2 y` the functional is `x^k y^l -> l! [k = 2l + 1]`, a weighted Gaussian. The question is whether the Mathieu
property of GMC(2) is a property of the Gaussian or of every one-dimensional twisted cohomology.
