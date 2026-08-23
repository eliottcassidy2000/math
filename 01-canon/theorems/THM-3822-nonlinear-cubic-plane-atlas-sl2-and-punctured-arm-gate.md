---
id: THM-3822
title: "Nonlinear cubic plane atlases have an SL2 and punctured-arm gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every dominant
  etale plane atlas of the THM-3811
  nonlinear cubic surface determines a polynomial SL2 matrix whose first
  row satisfies an explicit genus-three square-discriminant condition.
  The divisor h=0 is intrinsically G_m, so every component of its etale
  plane pullback has a nonconstant unit and at least two places at infinity.
  Consequently h cannot be a coordinate, and no atlas with h=xy-1 exists,
  for any polynomial companion k.  This closes the standard elementary
  big cell, every hyperbolic-unit first row, and the first Cohn row.  It
  does not close nonstandard multi-ended arms and does not construct a
  Jacobian counterexample.
source: root / nonlinear-cubic source-atlas and Cohn-word lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).
  The audit independently reconstructed all cubic laws and the different
  from the intrinsic equations without saturation, checked the
  D-quadratic-to-polynomial-square implication including k_0!=0 and the UFD
  scalar step, reconstructed B/(h)=K[k_0,k_0^-1], and verified that etale
  base change gives a genuinely nonconstant unit on every plane-arm
  component.  It then rederived the squarefree genus-three fibre,
  properness/Riemann--Hurwitz obstruction, the Laurent proof of
  K[x,y] intersect K(xy)=K[xy], and every edge of the Pell R=0 reduced-ratio
  argument.  No proof or scope repair was found.  The 37-gate exact
  companion verifies the
  cubic-order reductions, SL2 determinant, two intrinsic lift laws,
  reconstruction of all three cubic relations and the different, the
  intrinsic D-quadratic, both discriminant identities, the Pell completion
  and its square-zero binary cubic, the G_m arm quotient, standard-cell and
  hyperbolic specializations, squarefreeness, the genus-three boundary, and
  the first Cohn row.  Normal and optimized replay agree with the frozen
  output and raw hashes.  A separate 11,707-gate finite sidecar rejects all
  5,832 named shifted/inverse triple-Cohn words at each of two exact fibres;
  that bounded scout is not used in the all-degree proof.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
related:
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
script: 04-computation/jc2_nonlinear_cubic_atlas_sl2_gate_thm3822.py
output: 05-knowledge/results/jc2_nonlinear_cubic_atlas_sl2_gate_thm3822.out
script_sha256: 92604dd1388491e656644c4a072c1259c9c9677588b224522cad7d9f11d60b2a
output_sha256: 5db8d4ae4b8a6d2251f248aa341cc7d1e7cbbd5b2515d0bb70bbe0697a163f34
semantic_sha256: a83ceddd81bf567a58c457b82978bec96f49297823c4b2374083fd96f3efc6f7
hash_basis: raw LF bytes
---

# THM-3822 -- the source-plane problem is an SL2 completion problem

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.  Let `U=Spec B` be
the affine etale surface of
THM-3811, with

```text
B=S[A/D,omega/D],
D=C omega-3A theta-14A^2,
h=A/D,                 k_0=omega/D,              m=3theta+14A.    (1)
```

The theorem gives necessary conditions for a dominant etale morphism

```text
psi:A2_(x,y) -> U.                                                (2)
```

It is an obstruction theorem, not an existence theorem for `(2)`.

## 1. The intrinsic SL2 row

The Bezout identity of THM-3811 becomes

```text
C k_0-mh=1.                                                       (3)
```

Thus the five intrinsic functions `(h,k_0,m,C,D)` carry the matrix

```text
M = [ k_0  h ] in SL_2(B).                                       (4)
    [  m   C ]
```

The cubic multiplication laws imply two further identities:

```text
D(7h^2+3k_0^2)=1+2Ck_0,                                         (5)
hD(9h+14k_0)=k_0m+3hC^2.                                        (6)
```

These laws retain more information than the determinant-one row.  In the
other direction, over any integral `K`-algebra, a quintuple satisfying
`(3)--(6)` reconstructs

```text
A=hD,          omega=k_0D,          theta=(m-14hD)/3,             (7)
```

and direct reduction gives all three cubic relations of THM-3811 and
reconstructs its different `D`.  This is an exact reconstruction on the
dense `D!=0` graph.  No sufficiency assertion about possible parasitic
components over `D=0` is needed below: every actual map `(2)` already comes
from `B` and hence satisfies `(3)--(6)`.

Since `U -> A2_(A,C)` is etale, `(2)` is etale exactly when its composite
has unit Jacobian.  In polynomial coordinates this is

```text
Jac_(x,y)(hD,C) in K*.                                           (8)
```

So the counterexample search has become a constrained polynomial
`SL_2`-completion problem: find `(h,k_0;m,C)` and `D` satisfying
`(3),(5),(6),(8)` together with the required codimension-one coverage.

## 2. Eliminating the second row produces a genus-three sidecar

Put

```text
Q=3k_0^2+7h^2,
B_0=k_0^5-7h^2k_0^3-3h^2k_0^2-6h^3k_0^2-7h^4.                  (9)
```

Eliminating `(m,C)` from `(3),(5),(6)` gives the intrinsic quadratic

```text
h^2Q^2D^2+2B_0D+h^2-2k_0^3=0.                                 (10)
```

Its discriminant is

```text
Disc_D(10)=4k_0^2 H(h,k_0),                                    (11)
```

where

```text
H(h,k)=
 84h^7+36h^6k^2+196h^6k+84h^5k^3+36h^5k^2
 +49h^4k^4+112h^4k^3-12h^3k^5-14h^2k^6+12h^2k^5+k^8.          (12)
```

The same polynomial controls the missing second row.  Once one solution
`(m,C)` of `(3)` is fixed, every determinant-one completion is

```text
(m,C) |-> (m+k_0s,C+hs).                                        (13)
```

Substitution in `(6)` gives a quadratic in `s` whose discriminant modulo
`(3)` is exactly

```text
9H(h,k_0).                                                       (14)
```

For an actual dominant plane atlas, `k_0` is not the zero polynomial:
dominance makes `B -> K[x,y]` injective.  Evaluating `(10)` at its root `D`
therefore shows that `H(h,k_0)` is a square in `K(x,y)`.  Factoriality of
`K[x,y]`, followed by algebraic closedness of `K`, upgrades this to

```text
H(h,k_0)=w^2                    for some w in K[x,y].             (15)
```

This is only a necessary gate; a square does not by itself construct the
second row or prove `(8)`.

There is a useful completed-square form:

```text
H=L^2+4h^4 R,                                                    (16)
L=k^4-7h^2k^2+6h^2k-6h^3k,
R=49k^3+(27h-9)k^2+49h^2k+21h^3.                               (17)
```

Even the tempting boundary `R=0` cannot carry a nonconstant polynomial
row.  If `h!=0`, write `k/h=p/r` in lowest terms in the UFD `K[x,y]`.
Then `R=0` forces

```text
h=9p^2r/F(p,r),                  k=9p^3/F(p,r),                  (18)
F(p,r)=49p^3+27p^2r+49pr^2+21r^3.                              (19)
```

Polynomiality makes `F` divide both `p^2r` and `p^3`.  But
`F(0,r)=21r^3` and `(p,r)=1`, so `F` is a unit.  The binary cubic `F` has

```text
disc F(p,1)=-27046348 != 0.                                     (20)
```

Over `K` it is a product of three distinct linear forms.  Their product
being a unit makes each a unit; two independent factors then make `p,r`,
and hence `h,k`, constant.  If `h=0`, `(17)` reads
`R=k^2(49k-9)`, again giving only constant solutions.  Thus annihilating
the Pell correction cannot yield a dominant atlas.

## 3. The companion divisor is a punctured arm

Set `h=0` in the actual ring `B`.  Equations `(3),(5),(6)` give

```text
C=k_0^-1,                 D=k_0^-2,                 m=0.          (21)
```

Conversely these values satisfy the specialized cubic laws, so

```text
B/(h) = K[k_0,k_0^-1],             V_U(h) isomorphic to G_m.      (22)
```

The ideal `(h)` is prime in the regular surface `U`, hence `h` is a local
parameter along this smooth divisor.  If `(2)` is etale, base change of
`(22)` gives

```text
V_A2(psi^*h) -> G_m                                                   (23)
```

etale on every component.  Each component therefore has the pulled-back
coordinate `k_0` as a nonconstant unit.  Equivalently, every projective
completion of a component has at least two points at infinity: a
nonconstant rational function with neither zero nor pole on the affine
curve must have both a zero and a pole on its boundary.

This simple geometric fact is a strong gate:

```text
psi^*h cannot be a plane coordinate, and no component of V(psi^*h)
can be A1.                                                        (24)
```

It is the missing sidecar in approaches that record only the `SL_2` row.

## 4. Every standard hyperbolic arm is impossible

Assume, after a polynomial change of source coordinates, that

```text
h=xy-1.                                                          (25)
```

Let `tau=xy`, regard `K[x,y]` inside `K(tau)(x)` via `y=tau/x`, and put

```text
P_tau(Z)=H(tau-1,Z).                                             (26)
```

This is monic of degree eight in `Z`.  At `tau=0` it specializes to

```text
P(Z)=Z^8-14Z^6+24Z^5+49Z^4+28Z^3+196Z-84.                       (27)
```

Exact Euclidean reduction gives `gcd(P,P')=1`; equivalently its resultant
is nonzero (and is `100 mod 101`).  Therefore `P_tau` is squarefree over
`K(tau)`, and

```text
W^2=P_tau(Z)                                                     (28)
```

has smooth projective model of genus three.

The necessary square `(15)` gives a `K(tau)(x)`-point
`(Z,W)=(k_0,w)` of `(28)`.  If `k_0` were nonconstant in `x`, it would
define a nonconstant rational map from `P1_K(tau)` to this genus-three
curve.  Properness extends the map across its missing points, and
Riemann--Hurwitz forbids it.  Hence

```text
k_0 in K(tau).                                                    (29)
```

The elementary intersection

```text
K[x,y] intersect K(xy)=K[xy]                                    (30)
```

then shows that `k_0` is constant on the irreducible hyperbola `xy=1`.
This contradicts the nonconstant-unit requirement `(23)`.  Consequently:

```text
There is no dominant etale plane atlas of U with h=xy-1, for any
polynomial k_0 and any determinant-one completion.               (31)
```

The same proof applies whenever a source automorphism carries the arm to
`xy=1`.

## 5. Closed cells and the remaining constructive target

Several formerly separate failures are now instances of the same gate.

* The elementary first row `(h,k_0)=(x,1+xy)` already violates `(24)`.
  Independently, setting `y=0` in `(12)` gives a degree-seven polynomial
  in `x`, of leading coefficient `84`, so it cannot satisfy `(15)`.
* For the hyperbolic unit rows

  ```text
  h=xy-1,      k_0=cx^n,      m=sum_(j=0)^(n-1)(xy)^j,
  C=c^-1 y^n,                 c!=0,
  ```

  a harmless constant rescaling makes `(3)` exact.  They are all excluded
  by `(31)`.  Directly, the boundary square would be `P(cx^n)`.  Since
  `P(0)=-84` and `P` is squarefree, the chain rule makes `P(cx^n)`
  squarefree for every `n>=1`, another all-degree contradiction.
* The first Cohn row

  ```text
  (k_0,h;m,C)=(4t^2,2xt-1;1+2xt,x^2)                             (32)
  ```

  has determinant one but lies on the hyperbolic arm `(31)`.  The direct
  specialization `x=0` is the nonsquare `P(4t^2)`.
* The separate companion
  `04-computation/jc2_nonlinear_cubic_atlas_shifted_triple_cohn_scout_thm3822.py`
  exhausts the `5,832` ordered words
  `C(x+a,y+b)^epsilon C(x+c,y+d)^eta C(x+e,y+f)^zeta`, with every shift in
  `{-1,0,1}` and every exponent in `{+1,-1}`.  Over `GF(1009)`, both fibres
  `y=2,3` independently give an odd-degree or derivative-gcd nonsquare
  certificate for every word.  This is **FINITE-EXACT** evidence only: it is
  not an all-word Cohn closure and is not a dependency of `(31)`.

What remains is now sharper.  A plane atlas, if one exists, needs all of:

```text
(i)   a genuinely multi-ended, non-coordinate plane curve h=0;
(ii)  a nonconstant unit k_0 on every component of that curve;
(iii) the global square H(h,k_0)=w^2;
(iv)  a compatible polynomial second row and D solving (3),(5),(6);
(v)   the Keller condition Jac(hD,C) in K*;
(vi)  the codimension-one coverage required by THM-3811.         (33)
```

The highest-value constructive lane is therefore not another bounded
coefficient search in `(A,C)`.  It is the classification or construction
of **nonstandard punctured plane arms** equipped with a unit and a rational
map to the genus-three fibration `(28)`, followed by the second-row and
Jacobian lifts.  Interacting non-elementary Cohn words are a concrete source
of such arms, but the standard and single-hyperbola grammars are closed.

No planar Jacobian counterexample is claimed.  **QED.**
