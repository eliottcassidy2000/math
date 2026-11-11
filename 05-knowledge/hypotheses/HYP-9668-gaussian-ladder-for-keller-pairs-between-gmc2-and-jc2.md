---
id: HYP-9668
title: "Gaussian conditions on Keller pairs: what is left of the ladder between GMC(2) and JC(2). For a Keller pair (f, g) with Stein kernel M_z (THM-5138): K-GMC (M_z is a Mathieu-Zhao subspace), CV (cone vanishing) and SN (one-sided Stein nilpotence) follow from the automorphism property at every centre. As first filed, parts (A)-(C) proposed them for all Keller pairs as conjectures between GMC(2) and JC(2), (D) proposed that SN detects automorphisms, (E) proposed the Mathieu property for all rank-one potentials. Outcome: (A)-(C) are each equivalent to JC(2), so they are not intermediate; (D) is PROVED at centres in general position, FALSE on surfaces at special centres, and OPEN on C^2 at a prescribed centre; (E) is REFUTED (W = x^2 y). Open parts: (D') prescribed-centre detection on C^2; (E') the Mathieu property for rank-one potentials whose class of 1 is not zero; (Q) which hypothesis on a potential q makes Im(d - dq) a Mathieu-Zhao subspace."
status: >
  (A), (B), (C) [K-GMC, CV, SN for all Keller pairs and all centres]: each is EQUIVALENT TO JC(2) (THM-5138 (E),
  PROVED; the two-direction form of this collapse, conditional on the dimension count, was pushed first by
  Session269, CONCURRENT). They hold for automorphisms (THM-2022). K-GMC is the instance q = (f - a)(g - b) of the
  Image Conjecture of Zhao (J. Algebra 324 (2010), Conjecture 1.3); it was not a new conjecture.
  (D) [SN forces an automorphism]: PROVED at every centre whose coordinates are algebraically independent over the
  field of coefficients, in one direction and from first hits of x and y only (THM-5138 (E)); PROVED at every
  centre when SN is replaced by finite chaos (THM-5138 (C)). FALSE for etale maps of surfaces at special centres
  (PROVED: for the inclusion of C^2 minus {1 + xy = 0} at the centre (0, 0), SN and CV hold, the Stein space has
  dimension 2 and the class of 1 is not zero; for the double cover (x^2, y/x) at centres (0, b) the Stein space is
  zero).
  (D') [the same at a prescribed centre, on C^2]: OPEN. It is implied by JC(2), so it cannot be tested on C^2
  short of a counterexample; it is a statement about what a counterexample would look like.
  (E) [kernel of the functional is Mathieu-Zhao whenever H^2(C^2, d - dW ^) is one-dimensional]: REFUTED as filed
  (first audit): for W = x^2 y, 1 = 2 (d/dy - W_y)(y) - (d/dx - W_x)(x) lies in the kernel, which has codimension
  one. The same for x^p y^q, p, q >= 1, p != q.
  (E') [the same when the class of 1 is not zero, e.g. W = x + x^2 y]: OPEN; no evidence either way.
  (Q) OPEN question, no conjecture offered.
  First audit: PASS WITH CORRECTIONS (MISTAKE-764); second audit of the statements quoted from THM-5138 (D), (E)
  and of the special-centre example: PENDING at this commit.
session: kind-pasteur-2026-11-11-S24
source: 05-knowledge/results/kps_s24_stein_modules_keller_20261111/README.md
related:
  - 01-canon/theorems/THM-5138-stein-modules-of-keller-pairs-dimension-count-and-finite-chaos-criterion.md
  - 01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md
  - 05-knowledge/results/session269_centres_ports_switches_20261111/jacobian/proof.md
---

Notation as in THM-5138: `U = f - a`, `V = g - b`, dual fields `d_u, d_v`, `M_z = Im(d_u - V) + Im(d_v - U)`,
classes `[h]_z` in `S_z = C[x,y] / M_z`; `E_a = f - d_v`, `E_b = g - d_u`.

## What became of the ladder

    GMC(2)  =  K-GMC for the pair (x, y), and by transport for automorphisms      (THM-2022)
    at every centre:      automorphism  =>  K-GMC  =>  CV  =>  SN  =>  first hits  (THM-5138 (F))
    at a centre in general position: each of these  <=>  automorphism             (THM-5138 (E))
    for all Keller pairs: each of these  <=>  JC(2)                               (THM-5138 (E))

So there is no statement in this family strictly between GMC(2) and JC(2). The reason is THM-5138 (D): the
Stein kernels `M_z` are the fibres of one torsion-free module, a polynomial with coefficients in a field `k_0` lies
in `M_z` for no centre that is algebraically independent over `k_0`, and so "`d_u^m x` lies in `M_z`" means
"`d_u^m x = 0`" there: at such centres the Gaussian conditions are the classical criterion (a Keller pair with a
locally nilpotent dual field is an automorphism). The Gaussian side has content of its own only at special centres.

Bounded variants do not help either: "the Stein rank `dim S_z` of every Keller pair is at most `N`" is equivalent to
JC(2) for every `N`, because the generic degree is multiplicative under composition and the rank is at least
`2d - 1` (THM-5138 (B)); and "K-GMC for Keller pairs of generic degree at most `D`" restates "no Keller pair of
generic degree `2 <= d <= D`".

## (D') Detection at a prescribed centre on `C^2` (OPEN)

Let `(f, g)` be a Keller pair and `z` any centre, for instance one with coordinates in the field of coefficients.
Does SN at `z` (or CV, or K-GMC) force an automorphism?

- True with finite chaos in place of SN, at every centre (THM-5138 (C)).
- True for SN in one direction at centres in general position (THM-5138 (E)).
- False on surfaces. For `X = C^2 \ {1 + xy = 0}`, `F` the inclusion, `z = (0, 0)`: the Stein kernel is graded by
  the charge `deg x - deg y`, contains every element of nonzero charge, has codimension 2 in charge 0 and does not
  contain 1; so SN and CV hold and `F` is not an isomorphism (source note, Proposition 6.3). The centre is special:
  the missed curve `uv = -1` belongs to the pencil of hyperbolas `uv = t`. For a Keller pair on `C^2` every
  component of the non-properness curve leaves the plane at `(infinity, infinity)` only (THM-4122 (4)), so no
  component belongs to such a pencil; whether that is enough is the question.
- On the slice of the map of THM-1300 at the excluded centre (0, 0) (first audit, FINITE-EXACT): the classes of
  `w V^m` vanish for `m = 4, 6, 8` and those of `x V^m` for `m = 1, 3, 5, 7, 9`, by a parity symmetry of the
  potential, while the classes `h U^m` do not vanish; at the rational centre (2, 3) no vanishing is observed
  (`m <= 9` and `m <= 12`).

What would settle it on `C^2`: a description of the set of special centres of a Keller pair (the centres at which
some `d_u^m x` or `d_u^m y`, `m >= 1`, lies in `M_z` without being zero). By THM-5138 (D) it is a countable union of
proper algebraic subsets defined over the field of coefficients; detection at every centre would follow if for a
non-automorphism these subsets, for `m` in an infinite set, never have a common point.

## (E) REFUTED, and (E') (OPEN)

`dim H^2(C^2, d - dW ^) = 1` iff the generic fibre of `W` has `b_1 = 1`. Besides `xy` and its images under
automorphisms this happens for `x^p y^q` with `gcd(p, q) = 1` and for Broughton's polynomial `x + x^2 y`, which
has no critical point. For `x^2 y` the functional is `x^k y^l -> l! [k = 2l + 1]` and 1 lies in its kernel, so the
kernel is not a Mathieu-Zhao subspace: part (E) as filed is false. The two operators `d/dx - 2xy`, `d/dy - x^2`
commute and have order one with constant leading coefficients, so this is also an instance in which the statement
of Zhao's Conjecture 1.3, read for an arbitrary potential, fails over `C`; this is not claimed as new (one search
found no record; it needs its own sourced note).

(E') Let `W` be a polynomial on `C^2` with `H^2(C^2, d - dW ^)` one-dimensional and the class of 1 not zero. Is
`Im(d/dx - W_x) + Im(d/dy - W_y)` a Mathieu-Zhao subspace? Test case: `W = x + x^2 y` (the functional takes the
value -1 on 1). GMC(2) is the case `W = xy`.

(Q) For which potentials `q` is `Im(d - dq ^)` a Mathieu-Zhao subspace? Necessary: 1 is not in the image, unless
the image is everything. Known cases: `q = xy` (GMC(2)); `q = (f - a)(g - b)` for automorphisms; one variable
(van den Essen, Wright, Zhao). For `q = (f - a)(g - b)` with `(f, g)` a Keller pair the answer is JC(2).

## History of this file

The first version (commit d701fad1cd, pre-audit) presented (A)-(C) as new conjectures between GMC(2) and JC(2),
claimed that SN fails on the slice, and filed (E) without checking whether 1 lies in the kernel. The first audit
corrected all three (MISTAKE-764); the collapse of (A)-(C) was found the same day by Session269 and, in the form
quoted above, in this session after the audit.
