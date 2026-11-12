---
id: THM-5147
title: "Euler fields in twisted de Rham classes, and the failure of the Mathieu property for Broughton's polynomial. (A) For a polynomial W on C^2 and M_W = Im(d/dx - W_x) + Im(d/dy - W_y): if a polynomial vector field xi has xi(W) = W and constant divergence delta, then W h = (xi + delta) h modulo M_W, so [W^m h] = (e + delta)(e + delta + 1)...(e + delta + m - 1) [h] for an eigen-polynomial xi h = e h; if eta(W) = 0 and div eta = 0 then eta(h) lies in M_W. (B) For W = x + x^(l+1) y, l >= 1 (l = 1: Broughton's polynomial, generic fibre C^*, no critical point): M_W does not contain 1, W^m lies in M_W for m >= l, and x^l W^n = n! x^l modulo M_W is never in M_W; so M_W is not a Mathieu-Zhao subspace. This refutes HYP-9668 (E'). (C) Obstruction: an Euler field with divergence a non-positive integer, together with an eigen-polynomial of nonzero class whose eigenvalue plus delta is not a non-positive integer, excludes the Mathieu property. (D) For an etale map (u, v) of a surface and the Stein kernel at a centre z: [(UV)^m h]_z = [(E + 1)...(E + m) h]_z with E = (U d_u + V d_v)/2."
status: >
  (A) PROVED (two lines; the description of the image as 'divergence minus derivative along W' is standard, no
  novelty claimed for the lemma). (B) PROVED by hand (the linear form is a residue of an exponential integral; the
  Euler field is x d/dx - l y d/dy, divergence 1 - l); machine check for l = 1, 2 with the engine of THM-5138
  (exact modulo a prime, membership relative to the truncation). (C) PROVED from (A). (D) PROVED from (A) in the
  frame du ^ dv; checked on a Keller automorphism pair and on the slice map of THM-5138 (G).
  Consequence for HYP-9668: (E') REFUTED; question (Q) has two necessary conditions (the 1-property, the Euler
  obstruction). Not settled: rank-one potentials without an Euler field.
  The general Image Conjecture of Zhao is already REFUTED in characteristic zero (THM-2801); (B) is one more
  failure in the smallest setting (two operators in two variables over C, quotient of dimension one), this time
  with 1 outside the image. No literature search was made for (B) beyond the sources read by the audits of
  THM-5138; no priority is claimed.
  AUDIT: pending at this commit (PRE-AUDIT DRAFT).
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

## (B) `W = x + x^(l+1) y`

1. `lambda(x^k y^j) = (-1)^n j!/n!` for `n = (l+1)(j+1) - 1 - k >= 0`, and 0 otherwise, vanishes on `M_W`
   (`lambda(h) = Res_(x=0) int_0^infinity h e^(-W) dy dx`); `lambda(1) = (-1)^l/l!`, `lambda(x^l) = 1`.
2. With `xi = x d/dx - l y d/dy`: `xi(W) = W`, `div xi = 1 - l`. So `[W^m] = (1-l)(2-l)...(m-l) [1] = 0` for
   `m >= l`, and `[W^n x^l] = n! [x^l]`.
3. `P = W^l` has all its powers in `M_W`, and `x^l P^m` is never in `M_W`: `M_W` is not a Mathieu-Zhao subspace,
   although `1` is not in `M_W`. For `l = 1` the quotient `B / M_W` is one-dimensional.

## (C) The Euler obstruction

If `xi(W) = W`, `div xi = delta` is an integer `<= 0`, and `xi h = e h` with `[h] != 0` and `e + delta` not an
integer `<= 0`, then `W^(1 - delta)` has all its powers in `M_W` while `h W^((1 - delta) m)` never lies in `M_W`.
For `W = xy` (and for `W = (f - a)(g - b)` of a Keller pair, in the frame `df ^ dg`) the Euler field has
`delta = 1`.

## (D) Radial identity for Stein modules

For an etale map `(u, v) : X -> C^2` of a smooth affine surface, a centre `z = (a, b)`, `U = u - a`, `V = v - b`,
dual fields `d_u, d_v` and `E = (U d_u + V d_v)/2`: for all `h` in `O(X)` and `m >= 1`,
`[(UV)^m h]_z = [(E + 1)(E + 2)...(E + m) h]_z`, and `[(U d_u - V d_v) h]_z = 0`. In particular
`[(UV - c)^m]_z = m! e_m(-c) [1]_z` with `e_m` the truncated exponential series.
