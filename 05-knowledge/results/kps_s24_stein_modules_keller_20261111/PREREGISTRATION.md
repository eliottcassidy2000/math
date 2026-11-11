# Pre-registration - Stein modules of Keller pairs (kind-pasteur-2026-11-11-S24)

Owner directive (2026-11-11): merge prior work on the planar Jacobian conjecture with the proof of the two-dimensional
Gaussian moments conjecture, and define conjectures between GMC(2) and JC(2) that could form a chain excluding
high-degree counterexamples resting on an undiscovered topological anomaly. This file is written and pushed BEFORE
the runs S2-S5 below. It was written AFTER the prior-art sweep returned (lesson of MISTAKE-762).

## What exists (read before writing this)

- GMC(2) = two real Gaussians = one pair `(Z, conj Z)`, `E[Z^a W^b] = a! [a = b]`: THM-2022 (PROVED, Lean endpoint
  `GMC2Main.gmc2`), through NC2 (a null polynomial has charges of one sign) and one-variable DvdK (THM-2067,
  THM-2101, THM-2111). GMC(n) is false for n >= 3 (THM-1480, THM-1490, THM-1500). No arrow from GMC(2) to JC(2) is
  in canon (MISTAKE-237); THM-1830 lists one only as a superseded programme.
- Dual fields of a Keller map, the kernel ladder, cokernel = relative de Rham cohomology, "(f, g) automorphism iff
  D_f surjective iff the generic fibre of f is a line": THM-4620, THM-4626 (crediting Gurjar-Masuda-Miyanishi,
  Cerveau, Stein, Bavula's twisted module). The C^3 Keller map with a cubic fibre: THM-1300, THM-1310, THM-1335,
  THM-4621. Jelonek set: THM-2633 (inertia generates the monodromy), THM-4122 (components are polynomial curves).
- Session267 (same owner prompt, pushed 2026-11-11 12:02-12:20): NC2 consequences, the row repair operator as
  Gaussian rotation, a rational-cone criterion, a shear-word compiler for inverse coordinates, a wreath law for the
  fixed-sheet gate, three open bridges A, B, C. It has no twisted cohomology, no Stein functional of a Keller pair,
  no dimension count. The repo has no occurrence of "twisted de Rham", "Dimca", "bifurcation set" in canon.

## Working claims (derived by hand before this file; the runs below test parts of them)

Notation. `(f, g)` a Keller pair on `C^2`, `{f, g} = 1`; dual fields `d_u = {., g}`, `d_v = {f, .}` (`d_u f = 1`,
`d_u g = 0`, `d_v g = 1`, `d_v f = 0`; THM-4620). Centre `(a, b)`, `U = f - a`, `V = g - b`, `W = U V`.

**W1 (Stein module).** `M_F := (d_u - V) C[x,y] + (d_v - U) C[x,y]` is the image of `d - dW ^` on 1-forms, and
`S_F := C[x,y] / M_F = H^2(C^2, d - dW ^)`. Linear forms on `S_F` are the functionals that satisfy the Gaussian
integration-by-parts identities in the coordinates `(f, g)` ("Stein functionals"). Identities:
`Im D_W` lies in `M_F` (rotation invariance); `U^k V^l = k! [k = l]` in `S_F`; `h V^m = d_u^m h` and
`h U^m = d_v^m h` in `S_F` (Cameron-Martin); more generally `d_u^m d_v^n h = h H_(n,m)(U, V)` with the complex
Hermite polynomials. For an automorphism `S_F` is one-dimensional, its functional is `h -> E[h o F^(-1)]`, and
GMC(2) (THM-2022) says exactly that `M_F` is a Mathieu-Zhao subspace.

**W2 (dimension; uses the theorem of Dimca-Saito 1993 that `H^k(C^n, d - dW ^)` is the reduced cohomology of the
generic fibre of `W` in degree `k - 1`).** For a Keller pair of generic degree `d` with non-properness curve `S`
(components `S_i`, `delta_i` sheets missing over a generic point of `S_i`) and a generic centre,

    dim S_F = b_1(generic fibre of W) = c + sum_i delta_i (deg_u S_i + deg_v S_i) >= 2d - 1,

`c` the number of components of the generic fibre. So `dim S_F = 1` iff `F` is an automorphism:
JC(2) is the statement that the Stein functional of every Keller pair is unique.

**W3 (Hermite finiteness).** Say that `F` has finite chaos at the centre `(a, b)` if for every polynomial `h` the
classes of `d_u^m d_v^n h` in `S_F` vanish for `m + n` large. Claim: a Keller pair has finite chaos at a generic
centre iff it is an automorphism. Sketch: `h -> sum [d_u^m d_v^n h] s^m t^n / (m! n!)` is then a homomorphism of
Weyl-algebra modules from `C[x,y]` (Bavula's twisted module) onto `S_F (x) C[s,t]`, a sum of `dim S_F` copies of
the simple module `C[u,v]`; comparing ranks over `C(f, g)` gives `dim S_F <= d`, against W2.

**W4 (ladder).** For Keller pairs: JC(2) <=> finite chaos <=> Stein uniqueness => "the Stein kernel `M_F` is a
Mathieu-Zhao subspace" (K-GMC) => one-sided Stein nilpotence (`d_u^m h` and `d_v^m h` lie in `M_F` for `m` large).
GMC(2) is K-GMC for automorphisms. Whether the two weaker rungs imply JC(2) is the open part.

**W5 (toy, by hand).** On `C^* x C` with `(u, v) = (x^2, y/x)` (an etale double cover of `C^* x C`, not proper over
`C^2`): the Stein space of `(u - a)(v - b)` is two-dimensional, spanned by residue functionals at `x = +-sqrt(a)`,
and `x^(-2) (v - b)^m` is not in the Stein kernel for any `m`: one-sided Stein nilpotence and the Mathieu property
fail.

## Controls already run (not predictions)

`scripts/stein.py` (exact linear algebra modulo a prime): `dim H^2(C^2, d - dW ^) = 1, 4, 1, 1, 1` for
`W = xy`, `x^3 + y^3`, `x + x^2 y` (no critical point), `x^2 y` (non-isolated critical line), and
`W = f g` for the automorphism `(x + y^2, y + (x + y^2)^3)`. All agree with the theorem of Dimca-Saito.
`scripts/slice_pred.py`: the slice `F_1 = c` of the map of THM-1300 is `X_c = C^2 minus {1 + xy = 0}` with
`(u, v) = (y + 3c x w, x w + x w^2 - c x^3 w^3)`, `w = 1/(1 + xy)`, Jacobian `-2 w^3`; its non-properness curve is
`27 c^2 v^2 - 18 c u v + 16 c + u^3 v - u^2 = 0` (bidegree (3, 2)), which meets a generic hyperbola
`(u - a)(v - b) = t` in 5 points. Nothing about the Stein space of the slice or of the toy has been computed.

## Predictions

- **S2 (toy).** `dim H^2(C^* x C, d - dW ^) = 2` for `W = (x^2 - a)(y/x - b)`, `(a, b) = (2, 3)`; and the class of
  `x^(-2) (y/x - 3)^m` is nonzero for every `1 <= m <= 6`.
- **S3 (slice of THM-1300, c = 1, centre (2, 3)).** `dim H^2(X_1, d - dW ^) = 11`
  (`= b_1(Y) - 1 + b_2(X)` with `b_1(Y) = 1 + 2 * 5`: a connected fibre, two sheets missing over each of the five
  points, and the meridian of the removed hyperbola restricting nontrivially to the fibre).
- **S4 (one-sided Stein nilpotence fails on the slice).** For `h = w` the classes of `h (v - 3)^m` and of
  `h (u - 2)^m` are nonzero for every `1 <= m <= 5`.
- **S5 (Gaussian subalgebra on the slice and the toy).** `(u - a)^k (v - b)^l = k! [k = l]` in the Stein space for
  all `k + l <= 4` (an identity of W1; the run checks the implementation), and the class of 1 is nonzero.

## Reading

S3 confirmed: the dimension count of W2 survives on a surface with nontrivial topology, where the correction terms
`H^1(X)`, `H^2(X)` enter. S3 refuted: the count or the comparison theorem is being misapplied; W2 must be re-derived
before anything else is used. S4 confirmed: in the only planar non-proper etale map in canon, Gaussian moments see
the missing sheets (evidence for the detection half of the ladder). S4 refuted (classes vanish from some m on): the
one-sided rung does not detect non-properness, and only the two-sided statement W3 can.

## Addendum (written after S2-S5 were run and before the run S6; pushed first)

Outcomes of S2-S5 are in the outcomes section below (all four confirmed). A further test of the dimension count
W2 on the slice, with LINEAR potentials (fibres are preimages of lines; `S_1` has bidegree (3, 2), total degree 4,
two sheets missing over a generic point; `dim = b_1(fibre) - 1 + 1` as in S3, `b_1 = 1 - chi`,
`chi = 3 (1 - n) + n = 3 - 2n` for a line meeting `S_1` in `n` points):

- **S6.** `dim H^2(X_1, d - dW ^)` is 2 for `W = u` (`n = 2`), 4 for `W = v` (`n = 3`), and 6 for `W = u + v`
  (`n = 4`).

Reading: confirmed, the count is right for a second family of potentials with different intersection numbers;
refuted, the agreement in S3 was partly accidental and W2 needs re-derivation.

---

## Outcomes (added after the runs; nothing above was changed)

| prediction | outcome |
|---|---|
| S2 (toy on `C^* x C`: dimension 2; `x^(-2) V^m` nonzero for `m <= 6`) | **CONFIRMED** (`runs/s2_toy.txt`): `c = 2` at `N' = 18, 24, 30`; six classes nonzero |
| S3 (slice of THM-1300, centre (2, 3): dimension 11) | **CONFIRMED** (`runs/s3_s5_slice.txt`): `c(N) = 11` for `N <= 36` at `N' = 40, 46` |
| S4 (`w V^m`, `w U^m` nonzero for `m <= 5`) | **CONFIRMED** at `N' = 46` (all elements of weight at most 34, where the truncation is exact) |
| S5 (`U^k V^l = k! [k = l]`, `[1] != 0`) | **CONFIRMED**: 15 of 15 on the toy and on the slice; `[1]` nonzero on both |
| S6 (linear potentials on the slice: 2, 4, 6) | **CONFIRMED** (`runs/s6_linear.txt`) |

Notes.

1. S2 and S5 were derivable by hand (W5 and W1); they test the implementation on the rings `C[x, 1/x, y]` and
   `C[x, y, w]/(w(1+xy) - 1)`. The genuine predictions are S3, S4 and S6.
2. A first slice run at `N' = 22` printed "NONZERO" for elements of weight 28 and 34, which a truncation at 22
   cannot decide; the stored run uses `N' = 40, 46`, where `c(34) = c(36) = 11` shows that the kernel is complete in
   those weights.
3. Not predicted, found in the runs: the centre (0, 0) gives dimension 9 on the slice (an excluded centre of W2:
   `v = 0` is a horizontal asymptote of the non-properness curve); every Hermite class of `w` and of `x` computed
   on the slice is nonzero.
