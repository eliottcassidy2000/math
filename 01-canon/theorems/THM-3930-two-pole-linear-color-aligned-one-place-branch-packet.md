---
id: THM-3930
title: "Two-pole linear-color cubic has an aligned one-place line-plus-quintic branch packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Over an
  algebraically closed characteristic-zero field, an explicit unit-ideal
  binary cubic defines a normal nonmonogenic S3 cubic order whose squarefree
  discriminant is a line times an irreducible rational quintic. Both branch
  components have affine normalization A1 and the same projective infinity
  point. The quintic normalization has at most three addresses over every
  affine target point, with a sharp two-address finite root-at-infinity
  fibre. Thus the packet passes the nonmonogenic, one-place alignment, and
  cubic affine-address gates that killed the preceding constructions.
  Subsequent audited THM-3931 proves scalar units, class group Z, principal
  quintic ramification and a nonprimitive vertical class, and excludes every
  affine-plane/Keller open. Thus the packet is a sharp near miss, not a
  Jacobian counterexample.
source: root / post-THM-3927 finite repeated-root pole-divisor lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_degree6_one_place, 2026-08-23). The
  audit independently reconstructed generic irreducibility, the exact
  centered two-pole calculation, normalization birationality and finiteness,
  the squarefree/maximal-order normality step, common infinity, and the two
  incompatible divisor congruences in the Laurent-unit equation. It made
  explicit that Section 2 is complete only in the centered t(infinity)=0
  gauge and that m=e+1 is forced by cancellation against a zero at infinity,
  not by an assumed generic leading term. The visible two-address collision
  also exposes the stricter source-unibranchness question subsequently closed
  by THM-3931 and not used in the proof here. The assertion-free 35-gate
  companion LF-normalizes exactly to the frozen raw-LF output in normal and
  optimized mode; raw and semantic hashes and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
related:
  - THM-3927-unit-ideal-rational-sextic-affine-address-cap-two-place-boundary
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3931-degree-two-pole-cubic-principal-ramification-no-atlas
  - THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry
script: 04-computation/jc2_two_pole_linear_color_aligned_branch_thm3930.py
output: 05-knowledge/results/jc2_two_pole_linear_color_aligned_branch_thm3930.out
script_sha256: bb71b72a362663bf4694355b3ab7a6ef40911e5f61712578083b6213aac871fe
output_sha256: 487bed58c0db0a73bbfecf31f249f6ce27e705e96b00403c419a97e496e09c2f
semantic_sha256: bc4fefd7ee1f28b2af68456da3c5d9129c7453ac128a5706b304ecb81ee2a7fa
hash_basis: raw LF bytes
---

# THM-3930 -- finite root poles unlock the branch geometry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Choose `rho in k` with

```text
rho^2=-3.                                                        (1)
```

Put `R=k[A,C]` and consider the binary cubic

```text
Phi=a U^3+C U^2V+c UV^2+d V^3,                                (2)

a=(A-rho)^2,
c=2rho A+2/3,
d=A/2-rho/18.                                                  (3)
```

Its Delone--Faddeev algebra `S` is a normal finite-flat cubic domain, is
globally nonmonogenic, and has generic Galois group `S3`. Its discriminant is
squarefree with exactly two irreducible components: a vertical line and a
rational quintic. Both have affine normalization `A1`, their projective
closures have the same unique point at infinity, and the quintic has at most
three normalization addresses over every affine target point.

This is a positive completion packet, not a counterexample. In particular,
this theorem does not assert that `Spec S` contains an open `A2`, or even that
the natural ramification deletion has the required unit and class groups.

## 1. Unit coefficient ideal and generic cubic field

The four coefficients in `(2)` generate the unit ideal. Indeed `d` vanishes
only at `A=rho/9`, while

```text
a(rho/9)=(-8rho/9)^2 !=0.                                      (4)
```

Over `k(A)`, the dehomogenization

```text
f(T)=aT^3+CT^2+cT+d                                            (5)
```

is irreducible. It is linear in `C`; a factor independent of `C` would have
to divide both `T^2` and `aT^3+cT+d`. Their gcd is one because `d!=0` in
`k(A)`. Gauss's lemma therefore proves irreducibility over `k(A,C)`.

The unit ideal is only a local-generation passport. It does not imply that
the index form represents a unit; Section 5 proves that it does not.

## 2. Degree-two repeated-root incidence

Let `u` be a normalization parameter and define

```text
A(u)=u^3+rho u^2-u,
t(u)=(u+rho/3)/(u^2-1),                                       (6)

C(u)=-(3/2)(u^2-1)(u+rho/3)
             (u^2+(8rho/3)u-11/3).                           (7)
```

The rational function `t` has degree two. Its two finite poles are `u=1,-1`;
they are not cancelled because `rho/3` is neither `1` nor `-1`. Exact
substitution gives

```text
a(A(u))t(u)^3-c(A(u))t(u)-2d(A(u))=0,                        (8)
3a(A(u))t(u)^2+2C(u)t(u)+c(A(u))=0.                          (9)
```

Thus `t(u)` is the double root of `(5)` along the image of `(6)-(7)`.
The finite poles are load-bearing: they are precisely what escapes
THM-3929's root-regular integrality argument.

The construction comes from the centered trace-zero two-simple-pole
calculation. This calculation is complete only in the root gauge where the
value at normalization infinity is zero. After normalizing the two finite
simple root poles to `u=+/-1`, put

```text
A=u^3+p u^2+q u+r,                 t=(u+kappa)/(u^2-1).     (10)
```

Within this centered `t(infinity)=0` gauge, vanishing of the quadratic
coefficient in the minimal polynomial of `t` forces, away from the degree-one
cancellations `kappa=+/-1`, either

```text
p=0, q=-3,                         or
p=-1/kappa, q=-1.                                             (11)
```

In the first family the pulled-back color `C` has a pole at `u=-kappa`
unless `kappa=+/-1`, recovering the THM-3927 two-place debt. In the second,
polynomiality is equivalent to

```text
3kappa^2+1=0.                                                 (12)
```

Taking `kappa=rho/3` and translating `r` to zero gives `(3),(6),(7)`.
No classification of arbitrary degree-two root maps in a fixed linear-color
presentation is asserted outside this centered gauge.

## 3. Exact discriminant and normalization

Define

```text
J=243A^5-4455rho A^4-648rho A^3C-27738A^3-4104A^2C
  +18506rho A^2+432AC^2+2376rho AC+10839A
  +72C^3-48rho C^2+648C-627rho.                              (13)
```

The exact discriminant factorization is

```text
Disc(Phi)=-(9A-rho)J/324.                                    (14)
```

The parametrization `(6)-(7)` lies on `J=0`. Homogenizing it to degree five
gives a basepoint-free map `P1 -> P2`: finite parameter values lie in the
affine chart, and at parameter infinity the leading `C` coefficient is
nonzero. The map is birational. Indeed `[k(u):k(A)]=3`, while the rational
map `t(u)` has degree two, so `t` cannot lie in `k(A)`; the prime-degree
extension gives

```text
k(A,t)=k(u).                                                    (15)
```

At a generic discriminant point the unique double root is rational in the
cubic coefficients, hence belongs to `k(A,C)`. Equations `(6)-(9)` then give

```text
k(A,C)=k(A,t)=k(u).                                            (16)
```

The image is therefore an irreducible degree-five curve and `(13)` is its
equation up to the displayed scalar. The other factor in `(14)` is the line

```text
L: A=rho/9.                                                     (17)
```

It is not a component of `J`, since `A(u)` is nonconstant. Hence `(14)` is
squarefree.

The free rank-three algebra `S` is `S2`. Away from `(14)` it is etale. At
each of the two discriminant DVRs the exponent is one; the order/maximal-
order discriminant formula has an even index correction, so the local index
is zero. Thus every height-one localization is normal, and `S` is `R1` and
therefore normal. Generic irreducibility makes it a domain. The nonsquare
discriminant gives generic Galois group `S3`.

## 4. The exact one-place and address packet

Both affine branch normalizations are `A1`. For the quintic this is immediate
from the polynomial birational parametrization `(6)-(7)`. Its sole
projective normalization point at infinity maps to

```text
[A:C:Z]=[0:1:0],                                                (18)
```

because `deg A=3` and `deg C=5`. The projective closure of the vertical line
`(17)` meets `Z=0` at the same point `[0:1:0]`. Thus the two components pass
the common-infinity condition forced on a planar Keller nonproperness set.

For any affine target point, every quintic normalization address lies in one
fibre of the cubic polynomial `A(u)`. Therefore

```text
# addresses <=3.                                                 (19)
```

The finite root poles exhibit a sharp visible collision:

```text
u=1,-1  |-> (A,C)=(rho,0).                                    (20)
```

These are two distinct quintic addresses, and the double root moves to
`[1:0]` in the binary-root line. The line component is smooth and has one
address over every point. Hence neither branch violates THM-3920's cubic
affine-address cap. This numerical cap does not by itself certify that the
corresponding ramification boundary is unibranch in the source: distinct
normalization addresses may still coalesce there. That stricter successor
audit is deliberately not part of the present packet.

## 5. The index form represents no scalar unit

It remains to prove global nonmonogenicity rather than infer it from the
unit coefficient ideal. Put `K=k(A)` and let `z` be a root of `(5)`. Solving
for the target color gives

```text
C=-(az^3+cz+d)/z^2.                                           (21)
```

As a map from the `z`-line to the `C`-line, `(21)` has exactly two poles:
order two at `z=0` and order one at `z=infinity`. Consequently the integral
closure of `K[C]` in `K(z)` is

```text
K[z,z^(-1)],                 units=K* z^Z.                    (22)
```

Suppose, for contradiction, that `x,y in k[A,C]` and

```text
Phi(x,y)=mu in k*.                                             (23)
```

Since

```text
Phi(x,y)=a Norm_(K(z)/K(C))(x-yz),                             (24)
```

the integral element `x-yz` has unit norm over `K[C]`. It is therefore a
unit in `(22)`, so

```text
x-yz=lambda z^n,                 lambda in K*.                 (25)
```

If `x=0` or `y=0`, `(23)` would respectively make `d` or `a` a cube in
`K*` modulo a scalar. This is impossible because `d` has a simple zero and
`a` has a zero of order two.

Assume both are nonzero and put

```text
m=deg_C x,                         e=deg_C y.                  (26)
```

At `z=0`, the two terms in `x-yz` have valuations `-2m` and
`1-2e`. They have opposite parity and cannot cancel. At `z=infinity`, their
valuations are `-m` and `-e-1`. If `m<=e` or `m>=e+2`, direct comparison at
the two places is incompatible with the opposite valuations of a monomial.
Therefore `m=e+1`; at `z=0` the `x` term dominates and fixes `n=-2m`. At
infinity the two terms then have the same pole order `m`, whereas the right
side of `(25)` has a zero of order `2m`. Thus equality itself forces their
leading Laurent tails to cancel. In particular the degree relation is

```text
m=e+1,                         n=-2m.                           (27)
```

Taking norms in `(25)` and using `Norm(z)=-d/a` now turns `(23)` into

```text
mu=lambda^3 a^(2m+1)d^(-2m).                                  (28)
```

At the zero `A=rho` of `a`, whose order in `a` is two, equation `(28)`
requires

```text
2(2m+1)=0 mod 3,                 hence m=1 mod 3.              (29)
```

At the distinct simple zero `A=rho/9` of `d`, it requires

```text
-2m=0 mod 3,                     hence m=0 mod 3.              (30)
```

This contradiction proves that `(2)` represents no scalar unit. By the
Delone--Faddeev index-form criterion, `S` is globally nonmonogenic.

## 6. What remains

The packet simultaneously has

```text
unit coefficient ideal;             normal nonmonogenic S3 cubic;
squarefree branch;                   two A1 branch normalizations;
one common infinity point;           affine address count at most three.
                                                                    (31)
```

Subsequent independently audited THM-3931 completes that intrinsic
calculation: `S*=k*`, `Cl(S)=Zq`, the vertical ramification class is `2q`,
and the quintic ramification divisor `E` is principal. The two visible
addresses coalesce at one point of `E`, making `E` non-unibranch;
equivalently, deleting `E` makes its nonconstant defining function a unit.
Hence this packet admits no `A2`/Keller open. Euler characteristic is not
needed for its exclusion, and no broader claim about `JC(2)` is made.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_two_pole_linear_color_aligned_branch_thm3930.py
python3 -O 04-computation/jc2_two_pole_linear_color_aligned_branch_thm3930.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
