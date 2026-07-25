---
id: THM-2181
title: "Exact-square-prefix compression and monic depressed quartic closure"
status: >
  PROVED. In a positive-weight Rees filtration, suppose a Keller component
  is exactly f=K^2+B, with top(K)=h and with the first face A of B not
  divisible by rad(h). After reducing a mate to odd top exponent s, the full
  resonant train has a first unavoidable rational tooth
  A^((s+1)/2)/h. Every later seed reaches its first pole strictly after this
  tooth, even when K^2 has arbitrarily many earlier faces. Hence only the
  terminal s=1 equality survives. Consequently every planar Keller pair
  whose component is already a monic depressed quartic with polynomial
  coefficients is an automorphism. This does not by itself regularize the
  nonmonic approximate root of THM-2180.
source: codex-2026-07-24-JC-square-prefix-restoration
depends_on:
  - THM-2071
  - THM-2102
  - THM-2113
  - THM-2127
  - THM-2129
related:
  - THM-2071
  - THM-2102
  - THM-2129
  - THM-2158
  - THM-2180
---

# THM-2181 -- compressing an exact square prefix

Work in `C[x,y]`, with

```text
{F,G}=F_x G_y-F_y G_x.
```

The point is narrow. A quartic square completion begins

```text
P=H_0^2+(qz+delta),
```

and lower faces of `H_0^2` may occur before the first genuine remainder
face. They form one exact approximate-root prefix, not independent
obstructions. The following Rees argument compresses that whole prefix
before locating the first pole.

## 1. Exact-square-prefix theorem

Fix positive integral weights `w=(w_1,w_2)` and put

```text
W=w_1+w_2.                                            (1)
```

Let `K,B in C[x,y]` and suppose

```text
f=K^2+B.                                               (2)
```

Write `h=in_w(K)`. Assume that `h` is nonconstant, power-free, and
`w`-homogeneous of degree `d`. Encode the lower faces as

```text
K(tau)=sum_(eta>=0) K_eta tau^eta,
deg_w(K_eta)=d-eta,                 K_0=h,

B(tau)=sum_(delta>=rho) B_delta tau^delta,
deg_w(B_delta)=2d-delta,            B_rho=A!=0.       (3)
```

Assume

```text
rho>0,          a=deg_w(A)=2d-rho>0,
rad(h) does not divide A.                             (4)
```

There is no restriction on later faces of either series.

Suppose `{f,g}=1`. Subtract polynomial target shears until the mate is
reduced and has top face

```text
in_w(g)=c h^s,             c!=0,              s odd. (5)
```

Put

```text
D=(s+2)d-W,                 j_0=(s+1)/2.              (6)
```

Here `D>0`: the top faces `h^2` and `c h^s` commute, so the nonzero
constant bracket cannot occur at the top weight, while `D` is exactly its
Rees defect.

Then

```text
j_0 rho>=D,                                             (7)
d+j_0 a<=W.                                           (8)
```

Condition (4) also forces `d+a>=W`. Hence the only possible reduced
exponent is

```text
s=1,                 d+a=W,                 {h,A}!=0. (9)
```

This is a terminal reduction. A separate geometric theorem is still needed
to close the resulting mate; the quartic application below lands on
THM-2071's quadratic stratum.

## 2. The full-seed proof

Put

```text
F(tau)=K(tau)^2+B(tau),
G(tau)=sum_(t>=0) g_t tau^t,
deg_w(g_t)=sd-t.                                      (10)
```

The weighted Keller equation is

```text
{F(tau),G(tau)}=tau^D.                                (11)
```

Work in the completed rational function field
`C(x,y)[[tau]]`. For every integer `ell>=0`, take the constant-one branch

```text
S_ell(tau)
 =K(tau)^(s-ell)
  (1+B(tau)/K(tau)^2)^((s-ell)/2)
 =F(tau)^((s-ell)/2).                                (12)
```

It has leading face `h^(s-ell)` and commutes with `F`. For every integer
`T<D`, the centralizer induction of THM-2127 gives scalars `lambda_ell`,
with `lambda_0=c`, such that

```text
G(tau)=sum_(ell d<=T)
 lambda_ell tau^(ell d) S_ell(tau)
 modulo tau^(T+1).                                   (13)
```

Indeed, after earlier seeds are subtracted, a first residual coefficient
below defect `D` centralizes `h^2`, hence `h`. The weighted rational
centralizer theorem makes it zero unless its defect is `ell d`; in that
case it is a scalar multiple of `h^(s-ell)` and is removed by (12).
Impossible negative polynomial seeds simply have coefficient zero. Thus
(13) is an induction in the stated completed ring and uses no sparsity.

Assume for contradiction that

```text
T_0=j_0 rho<D.                                        (14)
```

Expand the first seed around the **full** approximate root:

```text
S_0=sum_(j>=0) binom(s/2,j) K^(s-2j)B^j.             (15)
```

For `j<j_0`, the exponent of `K` is positive, so the term is polynomial.
For `j>j_0`, its defect exceeds `T_0`. At `j=j_0`, equality of defect
forces every selected face of `B` to be `A`, and `s-2j_0=-1`. The unique
nonpolynomial contribution at defect `T_0` is therefore

```text
c binom(s/2,j_0) A^j_0/h.                            (16)
```

The coefficient is nonzero because `s` is odd.

No later seed cancels (16). If `s-ell` is even and nonnegative, (12) is a
polynomial power of `F`. If `s-ell` is odd, its first negative `K`-power
requires `(s-ell+1)/2` selections from `B`. Including the seed defect, the
first possible pole is

```text
ell d+((s-ell+1)/2)rho
 =T_0+ell(d-rho/2)
 =T_0+ell a/2
 >T_0                                                (17)
```

for `ell>0`. Moreover `T_0<(s+1)d`, since `rho<2d`, so every seed relevant
at `T_0` has `ell<=s`; the parity split above covers them all.

Thus polynomiality of the coefficient of `G` at `T_0` would force

```text
h divides A^j_0.                                      (18)
```

This contradicts (4). Hence (7) holds. Using
`rho=2d-a` and `s=2j_0-1`, it is equivalent to (8).

Finally, if `d+a<W`, then `{h,A}` has negative weighted degree and vanishes.
If `d+a=W` but the bracket still vanished, the graded centralizer theorem
would make `A` a positive scalar power of `h`. Either case would force every
irreducible factor of `h` to divide `A`, again contradicting (4). Therefore

```text
d+a>=W.                                               (19)
```

If `j_0>1`, then `d+j_0a>d+a>=W`, contradicting (8). Thus `j_0=1`,
which gives (9). This proves the abstract theorem.

## 3. Monic depressed quartics

Let `(P,Q)` be a planar polynomial Keller pair and suppose that, in
polynomial source coordinates `(x,z)`, one component is already

```text
P=z^4+p(x)z^2+q(x)z+r(x),             p,q,r in C[x]. (20)
```

Then `(P,Q)` is a polynomial automorphism.

If all three coefficients are constant, `P_x=0` and the Jacobian bracket
cannot be a nonzero constant. Otherwise choose the first upper Newton edge
from `z^4`, with slope

```text
max(deg(p)/2,deg(q)/3,deg(r)/4).                      (21)
```

Take positive integral weights clearing its denominator. The top face
contains `z^4` and another term. If it is power-free, THM-2102 closes the
pair. If it is a proper power, THM-2129's quartic face classification makes
it a square:

```text
in_w(P)=h^2,
h=z^2+p_top/2,                 deg_w(h)=d=2w_z.       (22)
```

Complete the square before inspecting lower faces:

```text
H_0=z^2+p/2,
delta=r-p^2/4,
P=H_0^2+B,                     B=qz+delta.            (23)
```

Every face of `B` lies strictly below `h^2`. If `B` is constant, an output
translation gives `{P,Q}=2H_0{H_0,Q}`, which cannot equal one. Otherwise
let `A` be its first nonconstant weighted face.

The one-variable top coefficient is a monomial, so

```text
h=z^2+c x^e
```

is square-free: a common factor with `partial h/partial z=2z` would also
divide `c x^e`, impossible. Thus `rad(h)=h`. Since `deg_z(A)<=1`, one has

```text
rad(h) does not divide A.                             (24)
```

Reduce the mate by subtracting scalar powers of `P` whenever its top face
is an even power of `h`. The descent reaches a positive odd exponent:
below terminal weight no bracket term has weight zero, while at equality a
putative identity `{h^2,in_w(Q)}=1` would be divisible by `h`.

Sections 1--2 force `s=1`. The reduced mate has weighted degree
`d=2w_z` and contains `z^2`; no term of `z`-degree above two can have that
weight. Hence the reduced mate has source-fibre degree exactly two.
THM-2071 closes every Keller pair with a quadratic member. Restoring the
subtracted target shears preserves invertibility.

### Terminal control

The equality in (9) is real. With weights `deg(x)=2`, `deg(z)=1`, put

```text
h=x+z^2,              P=h^2+z,              Q=-h.    (25)
```

Then `{P,Q}=1`, `d=2`, `a=1`, and `d+a=W=3`. The pair is triangular in
coordinates `(h,z)`. The theorem must land on quadratic closure rather than
declare every square face impossible.

## 4. Scope

The quartic application begins only after monic depression has polynomial
coefficients. THM-2180 proves `V|beta` for the nonmonic twice-odd branch,
but its remaining constant coefficient

```text
[4gamma-(beta/V)^2]/(8V)
```

may still have a pole. This theorem does not automatically turn that
nonmonic approximate root into (20), and it proves neither the remaining
congruence nor general JC(2).

The faithful carrier is the ordered Rees defect filtration, seed parity,
exact square root, factor support, and terminal weight. A tournament on
faces would erase the delay `ell a/2` in (17), which is the obstruction.
QED.
