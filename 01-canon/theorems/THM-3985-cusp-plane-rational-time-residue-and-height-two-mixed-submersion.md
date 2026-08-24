---
id: THM-3985
title: "Cusp-plane rational-time residue gate and height-two mixed submersion"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. On the
  height-two completion, every nonconstant first coordinate F(p,y) is
  excluded from a polynomial Keller pair. If its restriction to the cusp
  (p,y)=(s^2,s^3) is nonconstant, a transverse cusp-address residue proves
  the stronger exact statement J(F,k(x,t)) intersect k(F)=0. If the
  restriction is constant, F-c is divisible by H=p^3-y^2=t*p^2 and dF
  vanishes on t=0, which excludes regular mates but is not overstated as a
  rational-mate obstruction. The surprising positive panel
  A=alpha*p+gamma*y^m is a source submersion exactly at height two with
  m>=2. Its generic fibre is rational with 4m+1 source places, yet its time
  form has 3m nonzero simple residues. Adding a nonconstant h(x^2t), taking
  m=1, or raising the height creates an explicit affine critical point.
source: jc-zero-debt-lift + root / post-THM-3983 cusp-plane lane, 2026-08-24
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3983-coordinate-boundary-constancy-and-rational-place-budget
  - THM-3984-boundary-generator-coupling-criticality-and-holomorphic-time-form
script: 04-computation/jc2_cusp_plane_rational_time_residue_thm3985.py
output: 05-knowledge/results/jc2_cusp_plane_rational_time_residue_thm3985.out
script_sha256: da7f6b39a0826861ee1ca94b7c1d151f58e84868d7b1ed9ff4ae257f2dfa146a
output_sha256: 13cbb7b951a7abcec3a91c7eb4cb8db30fab4e6c5e8fc4b17772e4d9f929b2ad
semantic_sha256: 1db4bab794803c2ca49fba6f94904b9083c8a91609fd729da44c9a915b369fed
hash_basis: raw LF bytes
---

# THM-3985 -- the cusp plane has no Keller first coordinate

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.** Work over
an algebraically closed field `k` of characteristic zero. On the height-two
surface of THM-3973 put

```text
u=x^2t,                 z=1+u,
p=zt=t+x^2t^2,          y=xzt^2,                           (1)
H=p^3-y^2=tp^2.                                            (2)
```

Then the cusp-plane chart is birational and its affine modification is exact:

```text
k(x,t)=k(p,y),
x=yp/H,        u=y^2/H,        z=p^3/H,        t=H/p^2,    (3)
B_2=k[p,y,yp/H,y^2/H].                                    (4)
```

For every nonconstant `F in k[p,y]`, define its cusp restriction

```text
phi(s)=F(s^2,s^3) in k[s].                                (5)
```

There is a sharp dichotomy.

```text
phi nonconstant:
    J_(x,t)(F,k(x,t)) intersect k(F)={0};                 (6)

phi constant:
    F-phi(0) is in (H), dF vanishes on V(t), and no
    Q in k[x,t] (hence no Q in B_2) has J(F,Q) in k*.      (7)
```

The second row is deliberately only a **regular-mate** obstruction. A
rational function can carry a pole along a critical divisor, so `(7)` is
not claimed to imply the rational-image statement `(6)`. Together,
`(6)--(7)` prove that no nonconstant element of the entire cusp-plane
subring `k[p,y]` can be the first entry of a polynomial Keller pair.

The strongest positive hostile inside this closed subring is

```text
A_m=alpha*p+gamma*y^m,       alpha*gamma!=0,       m>=2. (8)
```

It has no affine critical point. Its geometric generic source fibre is
rational, but it has exactly

```text
r_U=4m+1                                                   (9)
```

places at infinity. The height-two completion adds the `m` boundary points
and therefore retains `3m+1` places. More sharply, `(6)` gives

```text
J(A_m,k(x,t)) intersect k(A_m)={0}.                       (10)
```

Thus `(8)` is a genuine all-degree family of polynomial submersions with
rational generic fibre, but not even a rational constant-Jacobian mate.

The positive stratum is exact. For `h in k[u]`,

```text
h(u)+alpha*p+gamma*y^m       is a submersion
iff h is constant                                             (11)
```

when `m>=2`. In the analogous height-`n` generators, `n>=2`, the unshifted
row `alpha*p+gamma*y^m` is a submersion exactly when

```text
n=2 and m>=2.                                               (12)
```

The endpoints `m=1` and `n>=3` have explicit critical points below.

## 1. The exact affine-modification chart

The identities

```text
v=y/p=xt,                  p=t+v^2                         (13)
```

give the inverse formulas `(3)`: indeed

```text
t=p-v^2=(p^3-y^2)/p^2=H/p^2,
x=v/t=yp/H.                                               (14)
```

They also give `u=x^2t=y^2/H` and `z=1+u=p^3/H`. Hence the
right side of `(4)` contains `x,u,p,y`, while the reverse containment follows
from `(1)--(2)`. This proves both the field equality and the ring equality;
no normalization or saturation is being guessed.

The Jacobian density of the birational chart is

```text
J_(x,t)(p,y)=(y^2-p^3)/p=-H/p.                            (15)
```

Equivalently, the source volume is

```text
eta=dx wedge dt=-(p/H) dp wedge dy.                       (16)
```

The exceptional curve `H=0` is the ordinary cusp

```text
(p,y)=(s^2,s^3),                                         (17)
```

and the kernel of `k[p,y]->k[s]` in `(17)` is exactly `(H)`.

## 2. The complete cusp-plane dichotomy

### 2.1 Constant cusp restriction

Suppose `phi(s)=c` is constant. The kernel statement after `(17)` gives

```text
F-c=H G(p,y)                                             (18)
```

for some polynomial `G`. On the source, `H=tp^2`, while `p=y=0` on
`t=0`. Thus both `H` and `dH` vanish along the whole line `V(t)`, and

```text
dF=G dH+H dG=0                  on V(t).                 (19)
```

For every `Q in k[x,t]`, the polynomial `J(F,Q)` therefore vanishes on
that line and cannot be a nonzero scalar. Since `B_2 subset k[x,t]`, the
same holds for every regular completion mate. This proves `(7)`.

As warned above, `(19)` alone says nothing about a rational `Q` allowed to
have a pole on `t=0`; no such overclaim is used.

### 2.2 Nonconstant cusp restriction

Now suppose `phi` is nonconstant. Let `q` be transcendental over `k` and
work over an algebraic closure of `K=k(q)`. Every root `s_i` of

```text
phi(s)-q=0                                                (20)
```

is nonzero and simple. The first assertion follows because `q!=phi(0)`;
the second is generic separability in characteristic zero. At the
corresponding cusp point `(p,y)=(s_i^2,s_i^3)`, direct differentiation gives

```text
J_(p,y)(F,H)
 =-2s_i^3 F_p-3s_i^4 F_y
 =-s_i^2 phi'(s_i) !=0.                                  (21)
```

Thus the generic fibre `F=q` meets the cusp transversely at every one of
these points. This argument is local on the normalization and does not
require the whole geometric generic fibre to be integral.

Using `(16)` and `(21)`, a time form on the normalized generic fibre is

```text
omega_F=-p/[H J_(p,y)(F,H)] dH.                          (22)
```

Its residue at the address `s_i` is

```text
Res_(s_i)(omega_F)=-p/J_(p,y)(F,H)
                  =1/phi'(s_i) !=0.                      (23)
```

The displayed final sign follows from `p=s_i^2` and `(21)`; only
nonvanishing matters. If `Q in k(x,t)` and

```text
J(F,Q)=b(F),                 b(F) in k(F),                (24)
```

then on the generic fibre `dQ=b(q)omega_F`. A rational differential of the
form `dQ` has zero residue at every normalization place, whereas `(23)` is
nonzero whenever `b!=0`. Hence `b=0`. This proves `(6)` without a generic-
fibre irreducibility assumption.

## 3. The mixed powers are rational submersions with many ends

Fix `(8)` and let `q` again denote the generic value. In the cusp-plane
coordinates,

```text
p_q(y)=(q-gamma*y^m)/alpha,
D_q(y)=y^2-p_q(y)^3.                                     (25)
```

Equations `(3)` become

```text
x=-y p_q/D_q,                 t=-D_q/p_q^2.              (26)
```

Thus the generic fibre has function field `K(y)`. More precisely, its
source-plane part is exactly

```text
P1_y minus (V(p_q) union V(D_q) union {infinity}).       (27)
```

Indeed, when `p_qD_q!=0`, formulas `(26)` give a source point and invert
the map. On the source, `D=y^2-p^3=-tp^2`; if `pD=0`, then `p=y=0`, which
cannot lie over the transcendental value `q`.

The polynomial `p_q` has `m` simple roots: a common zero of `p_q,p_q'`
would have `y=0` and then `q=0`. It is coprime to `D_q`, since a common
root would again have `p_q=y=0`. The degree of `D_q` is `3m`. It is
generically squarefree. To see the last assertion without a
discriminant census, a repeated root of `D_q` can be written

```text
p_q=s^2,                    y=s^3,                       (28)
```

with `s!=0`. The derivative equation is exactly

```text
2+(3gamma*m/alpha)s^(3m-2)=0.                            (29)
```

This confines `s` to a finite constant set and then makes

```text
q=alpha*s^2+gamma*s^(3m)                                (30)
```

constant, a contradiction. Thus `(27)` removes `m+3m+1=4m+1` distinct
geometric points, proving `(9)`. On `X_2`, precisely the `m` roots of
`p_q=0` are restored as the transverse points of the boundary
`D=V(x,z)=A1_y`. No root of `D_q` is restored: if `p!=0` and
`H=0`, write `s=y/p`, so `p=s^2` and `y=s^3`. The two determinantal
relations

```text
zy=xp^2,                    p(z-1)=xy
```

would give both `z=xs` and `z-1=xs`, a contradiction. Since `D` has
affine coordinate `y`, the point `y=infinity` is not restored either.
Thus the completion retains exactly `3m+1` places.

The exact specialized Hamiltonian row is

```text
delta(y)=J(A_m,y)=alpha*D_q/p_q,                          (31)
omega_m=p_q/(alpha*D_q) dy.                              (32)
```

Every root `rho` of `D_q` is simple and has `p_q(rho)!=0`, so `(32)` has
the nonzero residue

```text
p_q(rho)/(alpha*D_q'(rho)).                              (33)
```

This is the `3m`-address specialization of `(22)--(23)` and independently
reproves `(10)`.

## 4. Exact submersion and endpoint classification

For height two, `(31)` is also the shortest criticality test. In source
coordinates it reads

```text
J(A_m,y)=-alpha*t*p.                                    (34)
```

If `dA_m=0`, then `(34)` forces `t=0` or `p=0`. On `t=0`, one has
`p=y=0` and `(A_m)_t=alpha`. On the other component of `p=0`, one has
`z=0`, `t!=0`, and `y=0`; because `m>=2`,

```text
dA_m=alpha dp=alpha*t dz !=0.                            (35)
```

Thus `(8)` is a submersion.

The exponent endpoint is sharp. For `m=1`, the point

```text
x=gamma/alpha,                 t=-alpha^2/gamma^2        (36)
```

has `z=p=y=0` and

```text
d(alpha p+gamma y)=(alpha+gamma*x*t)dp=0.                (37)
```

It is an affine critical point.

The height endpoint is equally sharp. At height `n>=3`, put

```text
u=x^n t,       r=u(u+1),       f=u^2(u+1),
p=x^-n r,      y=x^(-n-1)f.                              (38)
```

At an off-color critical point, the `x` equation and then the `u` equation
reduce to

```text
x^[m(n+1)-n]
 =-m(n+1)gamma*f(u)^m/[n alpha r(u)],                    (39)

[(2-n)u+(1-n)]/(n+1)=0.                                 (40)
```

Choose

```text
u_0=-(n-1)/(n-2).                                       (41)
```

It is neither color, the right side of `(39)` is nonzero, and algebraic
closedness supplies a nonzero `x`. This gives an affine critical point for
every `m>=1`. Together with `(36)`, this proves `(12)`.

## 5. A nonconstant weight-zero shift always restores criticality

Return to height two and put

```text
A_(m,h)=h(u)+alpha*p+gamma*y^m,          g=h'(u),         (42)
N=3m-2,                                 K_0=-3m gamma/(2alpha).
```

At a point with `x u(u+1)!=0`, the two critical equations reduce exactly to

```text
x^N=K_0 u^(2m-1)(u+1)^(m-1),
x^2=alpha/(3g(u)).                                       (43)
```

For constant `h`, the original `u` equation after the first row reduces to
`-alpha/(3x^2)=0`, so it cannot hold; Section 4 proves that the map is a
submersion. Suppose `g` is not identically zero, and set

```text
e=gcd(2,N),
B(u)=K_0 u^(2m-1)(u+1)^(m-1).                           (44)
```

The two power equations in `(43)` have a common nonzero solution `x` exactly
when

```text
[alpha/(3g(u))]^(N/e)=B(u)^(2/e).                        (45)
```

After clearing the denominator, `(45)` is the polynomial equation

```text
alpha^(N/e)-3^(N/e)g(u)^(N/e)B(u)^(2/e)=0.               (46)
```

Its second term is nonconstant, while its values at `u=0,-1` are the
nonzero scalar `alpha^(N/e)`. It therefore has a root away from both colors;
no such root can be a zero of `g`. For even `N`, condition `(45)` is exactly
`x^N=(x^2)^(N/2)`; for odd `N`, the two square roots have opposite
`N`-th powers, so the coprime-power relation selects one of them. Hence a
valid nonzero `x` exists and `(43)` gives an affine critical point.

This proves the if-and-only-if statement `(11)`. It also isolates the
mechanism: the coefficient ratio `2/3` in the height-two Euler rows kills
the constant shift's critical equation, while every nonconstant weight-zero
jet supplies exactly the missing compatibility root.

## 6. Scope, connections, and reproduction

The theorem closes first coordinates in the subring `k[p,y]`, not arbitrary
elements of `B_2`. Statement `(6)` is stronger than criticality and survives
arbitrary rational companions; Statement `(7)` is intentionally restricted
to regular companions. The family `(8)` is a positive submersion/rational-
fibre control, not a Keller pair.

The place count `(9)` is far from THM-3983's extremal `r=m+1` boundary:
the cusp plane pays three additional punctures per boundary degree. The
nonconstant `h(u)` closure is complementary to THM-3984, whose cells carry
a leading `x`; neither theorem subsumes the other.

The companion verifies `(1)--(4)`, the cusp Jacobian and residue rows,
generic degree/squarefreeness controls, all endpoint critical formulas, and
both parity branches of the compatibility equation. Finite `m,n` loops are
hostile controls only; the theorem proves every degree symbolically.

Reproduce with

```bash
python3 04-computation/jc2_cusp_plane_rational_time_residue_thm3985.py
python3 -O 04-computation/jc2_cusp_plane_rational_time_residue_thm3985.py
sha256sum 04-computation/jc2_cusp_plane_rational_time_residue_thm3985.py \
  05-knowledge/results/jc2_cusp_plane_rational_time_residue_thm3985.out
python3 agents/check_docs.py
```

Independent hostile audit is requested especially for the typed constant-
versus-nonconstant cusp split, the residue sign and normalization passage,
the `4m+1` puncture exhaustion, and the even/odd common-power reconstruction.
