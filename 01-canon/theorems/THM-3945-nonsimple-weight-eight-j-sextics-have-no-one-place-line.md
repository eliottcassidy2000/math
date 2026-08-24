---
id: THM-3945
title: "Non-simple weight-eight J-sextics have no one-place line"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Every stereographic plane sextic obtained from the four-cuspidal trigonal
  curve in Sigma_2 has no line whose pullback to the normalization is
  supported at one point.  This holds uniformly for centers over a finite or
  infinite base point.  As a cited corollary, it excludes a one-place line
  for both non-simple weight-eight configurations J_{2,0}+4A_2 and
  J_{2,3}+3A_2.  It does not exclude other torus sextics or arbitrary
  rational discriminant branches.
source: jc-cohn3709 / Degtyarev stereographic four-cuspidal trigonal model, 2026-08-24
depends_on: []
related:
  - THM-3879-rational-torus-sextic-c3-one-place-tradeoff
  - THM-3882-rational-dual-one-place-wronskian-projection-criterion
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3943-rational-weight-eight-four-torus-sextics-have-no-one-place-line
script: 04-computation/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.py
output: 05-knowledge/results/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.out
script_sha256: b0d94071071d912550601789151943bf566d8e608c00033f403621d62f9bf818
output_sha256: d57a30e69516710812b3edf046481e4e54ccd99f8d103bbf0f4fa6f4dcb8841e
semantic_sha256: 42244fe4548295ba73bab808a633d595b75bf4e18d8f9cc85ddca26f5b356c31
hash_basis: raw LF bytes
---

# THM-3945 -- non-simple weight-eight J-sextics have no one-place line

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**
Work over `C`.  Let `Bbar` be the four-cuspidal trigonal curve in the
Hirzebruch surface `Sigma_2` given, on the affine `O(2)` chart, by

```text
f(x,y)=4y^3-(24x^3+3)y+(8x^6+20x^3-1)=0.             (1)
```

Contract the exceptional section to obtain the quadric cone

```text
Q={Z1^2=Z0 Z2} in P3,       (x,y) |-> [1:x:x^2:y].   (2)
```

For any nonvertex point `o in Q` not on `Bbar`, let `C_o` be the plane
curve obtained by projecting `Q` stereographically from `o`.  Then no
projective line `L` has pullback to the normalization of `C_o` supported at
one point.  Equivalently, `C_o` has no choice of line at infinity making its
affine normalization `A1`.

Degtyarev's classification and stereographic construction are used as a
**CITED** input: every irreducible plane sextic with singularities

```text
J_{2,0}+4A_2              or              J_{2,3}+3A_2              (3)
```

is obtained in this way, with the center respectively generic or sharing a
fiber with a cusp of `Bbar`.  See Section 3.6 of
[Degtyarev, *Irreducible plane sextics with large fundamental groups*](https://arxiv.org/abs/0712.2290)
and the trigonal-model construction in
[Degtyarev, *Oka's conjecture on irreducible plane sextics. II*](https://arxiv.org/abs/math/0702546).
Consequently neither family in `(3)` has a one-place line.

The classification and singularity labels in `(3)` are literature input.
The uniform one-place obstruction below is repo-derived and exact.

## 1. The normalization is a sparse binary-sextic map

The cited rational parametrization of `(1)` is

```text
x(t)=3t/(t^3+2),
y(t)=-(t^6-20t^3-8)/(2(t^3+2)^2).                    (4)
```

Substitution proves `f(x(t),y(t))=0`.  It is birational: on a dense open
set the inverse is

```text
t=(4xy-6x)/(4x^3-2y-1).                              (5)
```

Put `D=T^3+2S^3`.  After clearing the common denominator in `(4)`, the
normalization map to `Q` is

```text
Z0=D^2,
Z1=3 T S^2 D,
Z2=9 T^2 S^4,
Z3=-(T^6-20T^3S^3-8S^6)/2.                           (6)
```

These are basepoint-free binary sextics and satisfy `Z1^2=Z0 Z2`.
Dehomogenizing by `S=1`, their supports are

```text
Z0=t^6+4t^3+4,
Z1=3t^4+6t,
Z2=9t^2,
Z3=-t^6/2+10t^3+4.                                   (7)
```

The missing `t^5` coefficient is the decisive coordinate.  A binary sextic
supported at one finite point has the form `lambda(t-r)^6`; its `t^5`
coefficient is `-6 lambda r`.  Since `lambda!=0`, any such point would have
to be `r=0`.  Support only at `t=infinity` instead means that the
dehomogenized sextic is a nonzero constant.

## 2. A finite-base projection center gives two incompatible values of b

Every point of `Q` with `Z0!=0` has a unique expression

```text
o=[1:a:a^2:b].                                         (8)
```

A basis for the hyperplanes through `o`, and hence for the coordinates of
the projection from `o`, is

```text
U=Z1-aZ0,        V=Z2-a^2Z0,        W=Z3-bZ0.          (9)
```

The pullback of an arbitrary target line is

```text
F=A U+B V+C W.                                         (10)
```

Suppose first that `F=lambda(t-r)^6`.  The absent `t^5` coefficient gives
`r=0`.  The `t^4` and `t^2` coefficients of `(10)` are `3A` and `9B`, so
`A=B=0`.  Since `(10)` is nonzero, `C!=0`.  Its `t^3` and constant
coefficients now give

```text
10-4b=0,                 4-4b=0,                      (11)
```

which demand `b=5/2` and `b=1` simultaneously.

If the sole support point is infinity, `F` is constant.  The same `t^4`
and `t^2` equations give `A=B=0`, while the `t^6` and `t^3` coefficients
give

```text
-1/2-b=0,                10-4b=0,                     (12)
```

which demand `b=-1/2` and `b=5/2`.  Thus no center in `(8)` works.

## 3. Centers over the infinite base point also fail uniformly

If `Z0=0` on `Q`, then `Z1=0`.  Away from the vertex one has `Z2!=0`, so
every remaining center is uniquely

```text
o=[0:0:1:b].                                           (13)
```

Projection coordinates and an arbitrary line pullback are now

```text
U=Z0,       V=Z1,       W=Z3-bZ2,
F=A U+B V+C W.                                         (14)
```

Again `F` has no `t^5` term.  If it is a sixth power at a finite point, that
point is `t=0`.  Vanishing of the `t^4`, `t^3`, and constant coefficients
gives

```text
B=0,        4A+10C=0,        4A+4C=0.                 (15)
```

Hence `C=A=B=0`, contrary to `F!=0`.  If `F` is supported only at infinity,
its `t^6`, `t^4`, and `t^3` coefficients instead give

```text
A-C/2=0,       B=0,       4A+10C=0,                   (16)
```

again forcing `A=B=C=0`.  The parameter `b` never enters either
contradiction.  Thus the second center chart is closed, including every
fiber-collision specialization used for `J_{2,3}+3A_2`.

## 4. Why this is exactly the one-place obstruction

The center is not on `Bbar`, so the three projected coordinate sextics have
no common zero on its normalization `P1`.  Stereographic projection is the
birational elementary transformation in the cited construction, and its
image has degree six.  Therefore a target line pulls back to an effective
divisor of degree six.  If its set-theoretic support were one point `P`, the
divisor would be `6P`, and its binary form would be a sixth power of the
linear form cutting out `P`.  Sections 2 and 3 exclude all such forms.

This closes a precise many-character counterexample route.  The exceptional
families in `(3)` have four torus structures, but the extra cubic characters
do not buy a one-place affine boundary.  The theorem does **not** claim that
all rational torus sextics fail, nor that arbitrary rational one-place branch
curves cannot occur in a Keller discriminant.  It excludes exactly the two
non-simple weight-eight stereographic families and, more strongly, every
admissible center in their common trigonal construction.

## Reproduction

```bash
python3 04-computation/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.py
python3 -O 04-computation/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.py
cmp <(python3 04-computation/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.py) \
    05-knowledge/results/jc2_nonsimple_weight_eight_j_sextic_one_place_thm3945.out
```

The companion freezes `(4)`--`(7)`, both center charts, all load-bearing
coefficients, and the two full-rank contradictions.  It does not substitute
for the cited classification in `(3)`.
