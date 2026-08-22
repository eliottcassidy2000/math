---
id: THM-3348
title: "Linear-z generic puncture response and one-root valuation"
status: >
  PROVED + FINITE-EXACT HOSTILE-AUDITED.  For P=f(x)+g(x)z over a
  characteristic-zero field, the Hamiltonian response cokernel on the
  generic P-fibre is the algebraic de Rham H^1 of the affine line punctured
  at the roots of g, hence has rank deg(rad(g)).  Explicit logarithmic
  classes give a free K[P]-lattice of that rank.  The unit class sees only
  -dx/g and can miss puncture channels.  In the one-root case its exact
  annihilator is computed without gradient-unimodularity by a valuation
  formula.  This completes the generic response object behind THM-3326; it
  does not add a new planar Jacobian chart.
source: root/lrc-math-2026-08-12
depends_on:
  - THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion
related:
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
  - THM-2063-one-fiber-linear-planar-keller-pairs
script: 04-computation/jc_linear_z_generic_puncture_response_thm3348.py
output: 05-knowledge/results/jc_linear_z_generic_puncture_response_thm3348.out
script_sha256: 105bed99614bde959db73c2d3b7efb3927e94b9eb234db74717f73886520f4f6
output_sha256: 5137ba9a2d33adb4e33f2d9823ef5f8e7e471230587b0e690fa0e622acd4a4de
hash_basis: LF-normalized bytes
---

# THM-3348 -- linear-`z` generic puncture response and one-root valuation

**PROVED + FINITE-EXACT HOSTILE-AUDITED.**

## 1. Generic puncture theorem

Let `K` be a field of characteristic zero and put

```text
R=K[x,z],
P=f(x)+g(x)z,                 0 != g in K[x],
D_P(q)=P_x q_z-P_z q_x,
C_P=R/D_P(R).
```

Regard `C_P` as a `K[P]`-module.  Let `s=rad(g)`, `r=deg(s)`, and
`F=K(P)`.  Then

```text
C_P tensor_(K[P]) F
  ~= H^1_dR(Spec F[x,g^(-1)]/F)
  ~= F^r.                                                  (1)
```

More precisely, the polynomial classes

```text
eta_j=[g(x)x^j/s(x)] in C_P,                0<=j<r,       (2)
```

generate a free `K[P]`-submodule of rank `r`, and its quotient is
elementwise `K[P]`-torsion:

```text
0 -> K[P]^r -> C_P -> T -> 0.                              (3)
```

No finite generation of `C_P`, uniform annihilator of `T`, or splitting of
`(3)` is asserted.  If `g` is constant, the stronger integral statement
`C_P=0` holds.

The unit class `theta=[1]` maps under `(1)` to `[-dx/g]`.  For nonconstant
`g`, its generic image vanishes exactly when

```text
g=lambda*(x-a)^d,                         d>=2.            (4)
```

Thus the scalar observer `theta` need not see the whole response.  Under the
gradient-unimodularity hypothesis of THM-3326, every root of `g` is repeated:

- with one root, the generic response still has rank one, but `theta`
  vanishes there and its polynomial obstruction is vertical torsion;
- with at least two roots, `theta` has a nonzero residue vector, so
  `Ann_(K[P])(theta)=0` already after generic localization.

When `deg(g)>=2`, the finite residues of `dx/g` sum to zero.  The unit class
therefore occupies one conserved vector in the `r`-channel puncture space
and may have zero entries even when the full response has rank `r`.

## 2. Generic-fibre identification

Write `A=K[P]` and let `t` denote its abstract generator.  A nonconstant
polynomial in `K[x,z]` is transcendental over `K`, so `A~=K[t]`.  Localization
at `A\{0}` gives

```text
F tensor_A R ~= F[x,z]/(t-f(x)-g(x)z).                    (5)
```

The polynomials `g(x)` and `t-f(x)` are coprime in `F[x]`.  Their resultant
in `x` is a nonzero polynomial in the transcendental `t`; over a splitting
field it is, up to a nonzero scalar, the product of the nonzero factors
`t-f(alpha)`, with multiplicity.  Choose `u,v in F[x]` with

```text
u g+v(t-f)=1.
```

In `(5)`, `t-f=gz`, so `g(u+vz)=1`.  Hence `g` is a unit and

```text
F tensor_A R ~= F[x,g^(-1)],              z=(t-f)/g.      (6)
```

For a function written in coordinates `(t,x)`, the chain rule at fixed `t`
gives

```text
D_P=-g(x) partial_x.                                      (7)
```

The operator is `A`-linear because `D_P(P)=0`.  Localization is flat, and
tensoring the cokernel presentation therefore gives

```text
C_P tensor_A F
 ~= F[x,g^(-1)]/(-g partial_x)F[x,g^(-1)].                (8)
```

The map

```text
[h] |-> [-h dx/g]                                         (9)
```

identifies `(8)` with

```text
Omega^1_(F[x,g^(-1)]/F)/d(F[x,g^(-1)]).                  (10)
```

Indeed, if `h=-g q_x=D_P(q)`, then `-h dx/g=dq`; the converse is identical.
Inverting `g` is the same as inverting its squarefree part `s`, proving the
first isomorphism in `(1)`.

## 3. Rank and the integral logarithmic lattice

Let `E/F` split `s`.  Faithfully flat base change commutes with the derivative
cokernel.  Over `E`, partial fractions integrate every polynomial term and
every pole of order at least two.  Only simple-pole residues remain, and they
are independent because a rational derivative has zero residue at every
pole.  Thus the residue map identifies the base-changed quotient with `E^r`.

Already over `F`, the forms

```text
x^j dx/s,                              0<=j<r,            (11)
```

form a basis.  After writing `s=prod_i(x-alpha_i)`, their residue matrix is

```text
(alpha_i^j/s'(alpha_i))_(i,j),                              (12)
```

a diagonal rescaling of a Vandermonde matrix.  It is invertible over `E`, so
the `F`-linear map defined by `(11)` is an isomorphism by faithful flatness.
This proves the dimension statement in `(1)`.

Under `(9)`, `eta_j` maps to the negative of the `j`th form in `(11)`.  If
`sum_j a_j(P)eta_j=0` in `C_P`, localization and basis independence give
`a_j(t)=0` in `F` for every `j`; injectivity of `A->F` gives `a_j=0` in `A`.
Hence the map in `(3)` is injective.  Its localization is an isomorphism, so
the localization of its cokernel is zero.  This is exactly the elementwise
torsion assertion.

For constant `g`, one has integrally

```text
R=K[P,x],                              D_P=-g partial_x.
```

Characteristic zero permits coefficientwise polynomial integration in `x`,
so `D_P` is surjective and `C_P=0`.

## 4. What the unit observer sees

The class `theta=[1]` maps to `-dx/g`.  Hermite reduction says this
differential is exact precisely when its residues at all finite poles vanish.
Factor over an algebraic closure as

```text
g=lambda product_i(x-a_i)^(m_i),          d=deg(g).        (13)
```

If a rational primitive exists, its local pole order at `a_i` is exactly
`m_i-1`; its reduced denominator therefore has degree `d-r`.  After
subtracting its value at infinity, integration of the leading term of `1/g`
shows that its order of vanishing at infinity is exactly `d-1`.  A nonzero
proper rational function with denominator degree `d-r` has infinity order at
most `d-r`.  Thus `r<=1`.  A unique simple root has nonzero residue, whereas
for a unique repeated root

```text
d/dx[-1/(lambda(d-1)(x-a)^(d-1))]=1/g.                   (14)
```

This proves `(4)`.  The root in `(4)` descends to `K` from the coefficient of
`x^(d-1)`.

Gradient-unimodularity gives the THM-3326 specialization directly.  At a root
`alpha` of `g`,

```text
P_z=0,                         P_x=f'(alpha)+g'(alpha)z.
```

Avoiding a common zero for every `z` forces `g'(alpha)=0` and
`f'(alpha)!=0`.  Finally, when `d>=2`, `dx/g` has zero residue at infinity;
the residue theorem makes the sum of its finite residues zero.

## 5. One-root valuation refinement

The exact annihilator of `theta` can be computed without
gradient-unimodularity.  Suppose

```text
g=lambda*(x-a)^d,                       d>=2,
b=f(a),
e=ord_(x=a)(f(x)-b),                    e=infinity if f=b,
nu=min(e,d),
N=ceil((d-1)/nu).                                         (15)
```

Then

```text
Ann_(K[P])(theta)=((P-b)^N).                              (16)
```

Put `u=x-a` and invert `u`.  Then

```text
R[u^(-1)]=K[P,u,u^(-1)],
D_P=-lambda*u^d partial_u,
ker(D_P)=K[P],
Q_0=u^(1-d)/(lambda(d-1)),               D_P(Q_0)=1.      (17)
```

Every localized primitive of `H(P)` is `H(P)Q_0+J(P)`.  Since `J(P)` is
already polynomial, it cannot cancel a negative `u`-valuation.  Write

```text
H(T)=(T-b)^k L(T),                       L(b)!=0.
```

In the polynomial ring with coefficient variable `z`,

```text
ord_u(P-b)=nu,                ord_u(L(P))=0,
ord_u(H(P)Q_0)=k*nu-(d-1).                               (18)
```

Thus a polynomial primitive exists exactly for `k*nu>=d-1`, proving `(16)`.
In particular `theta` remains nonzero.  Under gradient-unimodularity,
`f'(a)!=0`, so `nu=1` and `N=d-1`, exactly recovering THM-3326's unit-class
annihilator.  Its canonical class `mu`, longer annihilator `((P-b)^d)`, and
marked jet bridge are additional integral extension data erased by generic
localization.

## 6. Hostile controls and scope

For `f=x` and repeated roots, the explicit Bezout row

```text
(1-g'z)P_x+((g')^2/g)z^2 P_z=1                           (19)
```

shows that the following puncture channels coexist with a smooth gradient.

1. For `g=x^2`, `theta` has residue vector `(0)`, while
   `eta_0=x` maps to `-dx/x` and is free.
2. For `g=x^2(x-1)^2`, `theta` has residues `(-2,2)` and the positive
   logarithmic basis has residue determinant `-1`.
3. For `g=x^2(x-1)^2(x-2)^2`, `theta` has residues
   `(-3/4,0,3/4)`, while the full basis determinant is `-1/2`.
   Thus the unit observer misses the middle puncture.

The exact companion checks six split residue matrices, including simple-root
and mixed-multiplicity boundaries, and 51 sharp instances of `(16)` through
`d=10`.  These are controls; the uniform theorem is the proof above.

This theorem does not close a new `JC(2)` chart: THM-2063 already handles the
linear-source-variable Keller locus.  It instead completes the generic
response object behind THM-3326.  Generic de Rham localization preserves
puncture residues and discards vertical extension torsion; THM-3326's
one-root jet packet is exactly the sidecar restoring that lost data.  Nor
does nonzero generic `H^1` itself obstruct a mate: the mate equation is
specifically `D_P(Q)=1`, so the relevant observer is `theta`, not an arbitrary
response class.
