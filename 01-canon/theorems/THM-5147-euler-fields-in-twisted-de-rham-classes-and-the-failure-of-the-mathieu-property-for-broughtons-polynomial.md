---
id: THM-5147
title: "Euler fields in twisted de Rham classes, and the failure of the Mathieu property for Broughton's polynomial. (A) For a polynomial W on C^2 and M_W = Im(d/dx - W_x) + Im(d/dy - W_y): if a polynomial vector field xi has xi(W) = W and constant divergence delta, then W h = (xi + delta) h modulo M_W, so [W^m h] = (e + delta)(e + delta + 1)...(e + delta + m - 1) [h] for an eigen-polynomial xi h = e h; if eta(W) = 0 and div eta = 0 then eta(h) lies in M_W. (B) For W = x^(l+1) y + r(x), l >= 1, r(0) = 0 (r = x, l = 1: Broughton's polynomial, generic fibre C^*, no critical point): M_W is the kernel of an explicit linear form, has codimension one, contains W^m for m >= l, and never contains x^l W^n = n! x^l; so M_W is not a Mathieu-Zhao subspace, and 1 lies in M_W iff the coefficient of x^l in e^(-r) is zero. This refutes HYP-9668 (E'). (C) Obstruction: an Euler field with divergence a non-positive integer, together with an eigen-polynomial of nonzero class whose eigenvalue plus delta is not a non-positive integer, excludes the Mathieu property. (D) For an etale map (u, v) of a surface and the Stein kernel at a centre z: [(UV)^m h]_z = [(E + 1)...(E + m) h]_z with E = (U d_u + V d_v)/2."
status: >
  (A) PROVED (two lines; the description of the image as 'divergence minus derivative along W' is standard, no
  novelty claimed for the lemma). (B) PROVED by hand: for r = x by the author (residue form, Euler field
  x d/dx - l y d/dy of divergence 1 - l); for general r, and the elementary proof of codimension one, by the
  independent audit (its U1, U2; re-checked by hand by the author); machine checks for l = 1, 2 with the engine of
  THM-5138 and exact certificates over Q by the audit (23 vector fields, five potentials). (C) PROVED from (A).
  (D) PROVED from (A) in the frame du ^ dv; checked on a Keller automorphism pair, on the slice map of
  THM-5138 (G), and by the audit with explicit functionals (190 identities).
  Consequence for HYP-9668: (E') REFUTED; question (Q) has two necessary conditions (the 1-property, the Euler
  obstruction). Not settled: rank-one potentials outside the monomials x^p y^q and the family x^(l+1) y + r(x).
  Literature: (B) is an instance of the statement of Zhao's Conjecture 1.3 (arXiv:0902.0210: images of commuting
  operators h(z) + sum c_i d/dz_i), which is already REFUTED in characteristic zero in canon (THM-2801 refutes SIC,
  the special case of the form restated in arXiv:1008.3962); (B) is one more failure, in the smallest setting (two
  operators in two variables over C, quotient of dimension one), with 1 outside the image. One search and three
  papers read through summaries by the audit found no record of it (the nearest: arXiv:1004.0521, Example 2.4, a
  derivation, by the 1-property). A search, not a clearance; no priority is claimed.
  AUDIT: independent, own code, PASS WITH CORRECTIONS
  (05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/audit/AUDIT.md; all applied; MISTAKE-768).
  JC(2) remains OPEN; nothing here bears on it beyond the identity (D).
session: kind-pasteur-2026-11-11-S25
source: 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/scripts/euler_field.py
related:
  - 01-canon/theorems/THM-5138-stein-modules-of-keller-pairs-dimension-count-and-finite-chaos-criterion.md
  - 01-canon/theorems/THM-2801-sharp-special-image-boundary-and-beta-shift-witness.md
  - 01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md
  - 05-knowledge/hypotheses/HYP-9668-gaussian-ladder-for-keller-pairs-between-gmc2-and-jc2.md
---

Notation. `B = C[x, y]`, `W` in `B`, `M_W = Im(d/dx - W_x) + Im(d/dy - W_y)`, classes `[h]` in `B / M_W`.
Proofs: part A of the source note.

## (A) The lemma

For every polynomial vector field `v`, `div(v) - v(W)` lies in `M_W`, and `M_W` consists of these elements. So for
a vector field `xi` and `h` in `B`, `xi(h) + (div xi - xi(W)) h` lies in `M_W`.

- If `xi(W) = W` and `div xi = delta` is constant: `[W h] = [(xi + delta) h]` and
  `[W^m h] = [(xi + delta)(xi + delta + 1)...(xi + delta + m - 1) h]`.
- If `eta(W) = 0` and `div eta = 0`: `eta(B)` lies in `M_W`.

## (B) `W = x^(l+1) y + r(x)`, `l >= 1`, `r(0) = 0`

1. `B = C x^l + M_W` (elementary reduction; audit U1).
2. `lambda_r(x^k y^j) = j! [x^((l+1)(j+1)-1-k)] e^(-r(x))` (0 for a negative index) vanishes on `M_W`
   (`lambda_r(h) = Res_(x=0) int_0^infinity h e^(-W) dy dx`), and `lambda_r(x^l) = 1`. So `M_W = ker lambda_r`, of
   codimension one, and `1` lies in `M_W` iff `[x^l] e^(-r) = 0`.
3. `lambda_r(x^k W^m) = [x^(l-k)] (m! e^(-r) e_m(r))`, `e_m` the truncated exponential series, and
   `m! e^(-s) e_m(s) = m! + O(s^(m+1))`. So `W^m` lies in `M_W` for every `m >= l`, and `x^l W^n = n! x^l` modulo
   `M_W` for every `n >= 0`.
4. `P = W^l` has all its powers in `M_W`, and `x^l P^m` is never in `M_W`: `M_W` is not a Mathieu-Zhao subspace.
5. For `r = x` (`W = x + x^(l+1) y`; `l = 1` is Broughton's polynomial): `lambda(x^k y^j) = (-1)^n j!/n!`,
   `n = (l+1)(j+1) - 1 - k`, `lambda(1) = (-1)^l/l!`, so `1` is not in `M_W`; and (3) follows from (A) with the
   Euler field `xi = x d/dx - l y d/dy`, `div xi = 1 - l`: `[W^m] = (1-l)(2-l)...(m-l) [1]`,
   `[W^n x^l] = n! [x^l]`. For `m < l` the class `[W^m]` is not zero (for `l = 2`: `W + 1` lies in `M_W`).
6. `W = x (x^2 y + x + 1)` has no vector field with `xi(W) = W` and constant divergence (audit U3); it is covered
   by (1)-(4).

## (C) The Euler obstruction

If `xi(W) = W`, `div xi = delta` is an integer `<= 0`, and `xi h = e h` with `[h] != 0` and `e + delta` not an
integer `<= 0`, then `W^(1 - delta)` has all its powers in `M_W` while `h W^((1 - delta) m)` never lies in `M_W`.
For `W = xy` (and for `W = (f - a)(g - b)` of a Keller pair, in the frame `df ^ dg`) the Euler field has
`delta = 1`.

## (D) Radial identity for Stein modules

For an etale map `(u, v) : X -> C^2` of a smooth affine surface, a centre `z = (a, b)`, `U = u - a`, `V = v - b`,
dual fields `d_u, d_v` and `E = (U d_u + V d_v)/2`: for all `h` in `O(X)` and `m >= 1`,
`[(UV)^m h]_z = [(E + 1)(E + 2)...(E + m) h]_z`, and `[(U d_u - V d_v) h]_z = 0`. In particular
`[(UV - c)^m]_z = m! e_m(-c) [1]_z` with `e_m` the truncated exponential series; where `[1]_z != 0` no `UV - c` is a
null element, and where `[1]_z = 0` (the double cover of THM-5138 at the centres `(0, b)`) every `UV - c` is.
