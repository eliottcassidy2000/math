---
id: THM-4430
title: "Laurent root-rotation rigidity and marked-packet genus"
status: >
  PROVED by rational-map rigidity, explicit monodromy and Riemann-Hurwitz,
  with independent exact permutation/formal-series controls. Exact rational
  root rotations are only monomial dilations; primitive width-three marked
  packet curves have explicit genus at least three. General width-three
  first-return bounds remain OPEN.
source: overnight-hexagon-sep05; Milnor section 5.1 is a cited Riemann-Hurwitz input
depends_on:
  - THM-2111-effective-compound-root-bound-for-one-variable-constant-terms
  - THM-4417-width-two-laurent-first-return-parabolic-critical-bound
proof: 05-knowledge/results/nc2_root_rotation_obstruction_overnight_hexagon_sep05.md
script: 04-computation/nc2_root_rotation_obstruction_overnight_hexagon_sep05.py
output: 05-knowledge/results/nc2_root_rotation_obstruction_overnight_hexagon_sep05.out
script_sha256: f078906ab85b5d8162b28d684137d3e4be3248b43d6bb9a8ac0fce1838607a49
output_sha256: 354c876c7aceaae0d07eb0a601f0a46cf97ddf7970c691ca529c6acc1925a661
hash_basis: raw LF bytes
---

# THM-4430 -- Laurent root-rotation rigidity and marked-packet genus

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof, primary source pin and 373-gate replay](../../05-knowledge/results/nc2_root_rotation_obstruction_overnight_hexagon_sep05.md)
are part of this theorem. This identifies a method boundary and its faithful
replacement object; it does not improve the general Laurent first-return bound.

Let R(0)!=0, deg R=M+N, M,N>=1, h(z)=z^M/R(z). Choose a local Mth
root and put w=z/R(z)^(1/M). The local root rotation is
 iota=w^(-1)(zeta*w), zeta a primitive Mth root. A power of order q>1
extends rationally if and only if it is multiplication by a qth root xi,
if and only if R(xi*z)=R(z), equivalently R(z)=S(z^q). In particular
q divides M and N. A finite-order rational map has degree one; preservation
of the two-point zero fibre forces it to be a dilation. This excludes neither
rational approximations as in THM-4417 nor higher-dimensional constructions.

For the primitive binomial f=z^(-M)+z^N, gcd(M,N)=1 and d=M+N>=3,
mark one small root z. The trace and norm of the other M-1 small roots
separately have degree binom(d-1,M-1) over C(z), as does their joint field.
The cover h has d distinct simple branch values and full S_d monodromy;
the marked stabilizer is S_(d-1). Transpositions exclude both subset-sum
and subset-product collisions, including the complementary-product boundary.

For M=3, the faithful remaining pair satisfies

```text
y^2-S*y+P divides (y^3*R(z)-z^3*R(y))/(y-z),
S/z -> -1, P/z^2 -> 1.
```

These branch conditions are essential. For every N>=1 coprime to three,
the smooth compact marked-pair curve of f=z^-3+z^N has genus

```text
g_N=[N(N+1)(N+3)-2-floor(N^2/2)]/2 >=3.
```

Over the t-line its degree is (N+3)(N+2)(N+1)/2. Each of the N+3
transposition branches contributes (N+1)(3N+2)/2 ramification, and the
zero-fibre (3)(N) permutation has 4+floor(N^2/2) cycles. Riemann-Hurwitz
gives the formula and forbids degree>1 holomorphic self-maps of this curve.
For the first interior-width hostile z^-3+z^4, the marked packet degree is
15 and genus65; its trace first deviates at order15, its norm at order23,
and its first constant-term return is7. Normalization therefore does not
turn this correspondence into the rational dynamics used at width two.
No external priority or new Lean verification is claimed.
