---
id: THM-3081
title: "Terminal toric residue-parameter Mobius rigidity and autonomous decoder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the normalized
  residue-degree-one target-coordinate-line scope inherited from THM-3074,
  assume C(M_0,R_0)=C(x,y) and that C(u) maps isomorphically onto the
  residue field.  At
  at the primitive terminal stage of THM-3080 choose A g+B e=1 and set
  S=R^A Z^B, Theta=Z^g/R^e.  Then w(S)=1 and the residue theta of Theta
  generates the full residue field: C(u)=C(theta).  Consequently theta is a
  Mobius function of u, and every terminal initial coefficient has the form
  rho^w H(theta) with H rational.  The target, two-form prefactor, and exact
  Keller coefficient reduce to one autonomous rational decoder for theta;
  Mobius differentiation further forces its target/prefactor quotient to be
  a monomial times the square of one linear form.
  At every earlier strict stage, retaining the nonprimitive relation
  m^g/r^e=c^d and a Bezout complement is a degree-d torus isogeny with exact
  kernel mu_d; retaining the primitive c makes it unimodular.  Thus the gcd
  descent of THM-3080 is also an exact root-torsor discharge ledger.  This
  identifies, but does not supply, the pointed/common-atom lift needed to
  synchronize conjugate A4 quotient charts.  No polynomial globalization,
  arbitrary-Jelonek straightening, full C3, A4/S4, or JC(2) exclusion is
  asserted.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-3080-c3-finite-toric-key-tower-depth-partition-and-gcd-descent
related:
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
  - THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law
  - THM-3077-pointed-norm-relative-line-lift-and-relation-carry-obstruction
script: 04-computation/jc_terminal_toric_residue_decoder_thm3081.py
output: 05-knowledge/results/jc_terminal_toric_residue_decoder_thm3081.out
script_sha256: 46d63117fe376e13d534707453824db2e8578e52aba7106ee9e30981431bee46
output_sha256: 6e25871afd87df7c0057871a8ab4aacf740455454757d665b45f7343f94bcd09
hash_basis: LF-normalized bytes
---

# THM-3081 -- the terminal key recovers the residue coordinate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and new object

[THM-3080](THM-3080-c3-finite-toric-key-tower-depth-partition-and-gcd-descent.md)
turns every coordinate-line `C3` cancellation branch of THM-3074 into a
finite Laurent toric tower.  At its terminal stage `N`, write

```text
Omega=U dlog(M) wedge dlog(R),
M=c(1+m(u)s^e+O(s^(e+1))),
R=r(u)s^g(1+O(s)),                                    (1)

g,e>=1,                    gcd(g,e)=1.                  (2)
```

The terminal wedge coefficient is nonzero.  THM-3080 uses this to show that
equal-weight Laurent initials cannot cancel and then proves `(2)` from the
normalized divisorial value group.

This theorem identifies the object left after that primitivity argument.
The terminal value-zero coefficient ratio is not an arbitrary rational
function of the residue parameter `u`: it is itself a coordinate on the
residue line.  Moreover the target and Jacobian coefficient become rational
functions of this one coordinate and satisfy an exact autonomous ODE.

The same exponent calculation also exposes what every earlier strict stage
forgets if it records only a nonprimitive scalar relation: an exact `mu_d`
root torsor.  This is the precise contact with THM-3077's relation/carry
formalism.

The residue parameter throughout is the actual coordinate-line parameter
`u=P` from THM-3074.  The abstract terminal exponent algebra makes sense
over other coefficient fields, but the identification `P=u`, the residue
field equality `C(u)`, and the Mobius conclusion below do **not** straighten
an arbitrary Jelonek component.

More literally, the inherited hypotheses are

```text
K=C(x,y),            w(K*)=Z,            C(M_0,R_0)=K,
C(u) -> kappa(w) is an isomorphism.                         (0)
```

Every strict toric change is unimodular, so the terminal pair still
generates `K`.  Both full-field generation and the displayed residue-field
isomorphism are used in Section 3.

## 2. A terminal value-one/residue chart

Choose integers `A,B` with

```text
A g+B e=1.                                               (3)
```

Put `Z=M/c-1` and define

```text
S=R^A Z^B,
Theta=Z^g/R^e.                                          (4)
```

The exponent matrix in the ordered coordinates `(R,Z)` is

```text
[  A   B ]
[ -e   g ],                                             (5)
```

whose determinant is one by `(3)`.  Hence

```text
C(M,R)=C(Z,R)=C(S,Theta),                               (6)

w(S)=1,                    w(Theta)=0.                  (7)
```

The leading coefficients are

```text
rho=r^A m^B,
theta=m^g/r^e in C(u)*.                                 (8)
```

The terminal coefficient in THM-3080 is nonzero, and because `gcd(g,e)=1`
it says

```text
theta'/theta=g m'/m-e r'/r !=0.                         (9)
```

Thus `theta` is a nonconstant rational function of `u`.

Inverting `(5)` at the level of leading coefficients gives

```text
r=rho^g theta^(-B),
m=rho^e theta^A.                                        (10)
```

Therefore a Laurent monomial `R^j Z^n` of value

```text
w=jg+ne                                                  (11)
```

has leading coefficient

```text
r^j m^n=rho^w theta^(-Bj+An).                           (12)
```

The residue coordinate is independent of the Bezout choice.  Replacing
`(A,B)` by `(A+ke,B-kg)` changes

```text
S to S Theta^(-k),              rho to rho theta^(-k),
```

but leaves `Theta` and `theta` fixed.  The later powers of `rho` and rational
decoders transform covariantly, so the Mobius residue coordinate is the
canonical part of the chart.

## 3. Universal terminal initial decoder

Let `F` be a nonzero Laurent polynomial in `M,R`.  Expand each Laurent
coefficient of `M` at `M=c`, equivalently at `Z=0`.  THM-3080's terminal
noncancellation shows that the least-weight coefficient is nonzero.  By
`(11)--(12)`, all monomials on that weight line share `rho^w`, and their
remaining exponents form a Laurent polynomial in `theta`.  Thus

```text
lead_w(F)=rho^(w(F)) H_F(theta),
H_F in C[T,T^(-1)] nonzero.                             (13)
```

Every element of `C(M,R)*` is a quotient of two such Laurent polynomials.
Taking the quotient of `(13)` proves the field version

```text
lead_w(F)=rho^(w(F)) H_F(theta),
H_F in C(T)*                                             (14)
```

for every nonzero `F in C(M,R)`.

In particular, a value-zero field element has residue in `C(theta)`.  But
`M,R` generate `C(x,y)`, whose residue field is `C(u)` by the coordinate-line
hypothesis.  Since `theta` itself lies in `C(u)`, `(14)` gives both inclusions:

```text
C(u)=C(theta).                                          (15)
```

This is stronger than merely saying that `theta` is nonconstant.

## 4. Mobius rigidity of the residue parameter

In the inherited coordinate-line setup, apply `(14)` to the target coordinate
`P=u`, which has value zero.  There is
a rational function `H_P in C(T)` such that

```text
u=H_P(theta(u)).                                        (16)
```

The nonconstant rational maps

```text
theta:P^1_u -> P^1_theta,
H_P:P^1_theta -> P^1_u                                 (17)
```

compose to the identity.  Degrees of rational maps multiply, so

```text
1=deg(H_P o theta)=deg(H_P)deg(theta).                  (18)
```

Consequently both maps have degree one.  Equivalently,

```text
theta(u)=(a u+b)/(c u+d),             ad-bc!=0,          (19)
```

and `H_P` is its inverse Mobius transformation.  No finite cover of the
residue line remains hidden in the terminal coefficient.

The Mobius conclusion is rational, not necessarily affine.  Earlier inverse
toric substitutions can introduce denominators, so it would be wrong to
claim in general that `H_P` is a Laurent polynomial.  The exact controls in
Section 7 happen to land in the affine or reciprocal-affine subcases.

## 5. The target and prefactor decoders

Return briefly to tame ramification index `E`, as in THM-3080.  At the
terminal stage let

```text
w(U)=sigma,                   e=E-sigma.                 (20)
```

Apply `(14)` to the target `Q=t`, whose leading coefficient is `tau`, and to
the prefactor `U`.  There are nonzero rational functions `K,L in C(T)` with

```text
tau=rho^E K(theta),
lead(U)=rho^(E-e)L(theta).                               (21)
```

The terminal coefficient equation is

```text
lead(U)[g m'-e m r'/r]=E kappa^(-1)tau.                 (22)
```

Since `gcd(g,e)=1`, equations `(8)--(10)` give

```text
g m'-e m r'/r=m theta'/theta,
m=rho^e theta^A.                                        (23)
```

Substitution of `(21)--(23)` into `(22)` cancels the entire value-one scale
`rho^E` and leaves the autonomous decoder

```text
theta'
 =E kappa^(-1) theta^(1-A) K(theta)/L(theta).           (24)
```

Combining `(24)` with `u=H_P(theta)` yields the rational identity

```text
H_P'(T) E kappa^(-1) T^(1-A) K(T)/L(T)=1.               (25)
```

The Mobius form `(19)` makes `(25)` completely explicit.  With
`Delta=ad-bc`, direct differentiation and substitution of the inverse map
give

```text
theta'=(a-c theta)^2/Delta,                              (25a)

K(T)/L(T)
 =kappa/(E Delta) T^(A-1)(a-cT)^2.                      (25b)
```

Thus away from the toric points `T=0,infinity`, the terminal
target-to-prefactor ratio has at most one zero and it has even multiplicity
two; it has no other finite zero or pole.  The square marks the point
`u=infinity` in the Mobius coordinate.  This resembles the square sidecars
in the repo's cubic resolvent/discriminant anatomy, but no identification
with a physical quartic resolvent is asserted here.

Thus the terminal stage separates into a value-one scale `rho` and one
Mobius residue coordinate `theta`; the Keller equation lives entirely on
the latter after the predictable powers cancel.

Equation `(25)` is a decoder, not yet a contradiction.  The remaining
global problem is to determine which rational functions `K,L` can actually
arise from simultaneous polynomial target coordinates after all inverse key
substitutions.

## 6. Earlier strict stages are exact root-torsor lifts

At an earlier strict stage of THM-3080, write

```text
g=d alpha,                    e=d beta,
gcd(alpha,beta)=1,                                         (26)
```

and choose `gamma,delta` with

```text
alpha delta+beta gamma=1.                                (27)
```

The differential cancellation proves that the **primitive** ratio

```text
c=m^alpha/r^beta in C*.                                  (28)
```

Consider instead retaining only the nonprimitive scalar and the Bezout
complement

```text
C=m^g/r^e=c^d,
psi=m^gamma r^delta.                                    (29)
```

The exponent matrix of the torus map

```text
(m,r) |-> (C,psi)                                       (30)
```

is

```text
[ g       -e ]
[ gamma delta ],                                        (31)
```

with determinant `d`.  In the primitive coordinates `(c,psi)`, whose
exponent determinant is one by `(27)`, map `(30)` is simply

```text
(c,psi) |-> (c^d,psi).                                  (32)
```

Hence it is a degree-`d` torus isogeny with exact kernel

```text
mu_d x {1}.                                             (33)
```

Retaining `c` rather than only `C` chooses a lift, equivalently a
trivialization of this root torsor over the packet, and makes the next key

```text
Z^alpha-c R^beta                                        (34)
```

well typed.  This is the rank-two signed-character version of
[THM-3077](THM-3077-pointed-norm-relative-line-lift-and-relation-carry-obstruction.md):
a nonprimitive norm plus a complement loses finite diagonal torsion, while a
pointed primitive root restores it.

This choice is not a canonical global section of the torsor; changing the
primitive root acts by the kernel `mu_d` in `(33)`.

Along the key tower the torsor order is

```text
d=gcd(g_i,e_i)=g_(i+1).                                 (35)
```

The gcd descent is therefore also a **root-torsor discharge ledger**.  The
terminal equality `gcd(g_N,e_N)=1` says that no finite root ambiguity remains
in the last residue coordinate.

## 7. Three exact terminal decoders

The controls from THM-3074/3080 realize both Mobius orientations.

### One-stage equality packet

For

```text
x=s^(-1),                  y=s^(-1)+3u s^4,              (36)
```

the terminal values are `(g,e)=(1,5)`.  Take `(A,B)=(1,0)`.  Then

```text
rho=1,                    theta=3u,
u=theta/3,                theta'=3.                     (37)
```

### Two-stage strict packet

At the terminal stage of THM-3080's `7=4+3` packet,

```text
(g,e)=(2,3),              (A,B)=(-1,1),
r=u,                       m=-3u.                        (38)
```

Consequently

```text
rho=-3,                   theta=9/u,
u=9/theta,                theta'=-theta^2/9.            (39)
```

Here the residue coordinate is reciprocal-affine.

### Three-stage packet

At the terminal stage of the `11=4+4+3` packet,

```text
(g,e)=(4,3),              (A,B)=(1,-1),
r=u,                       m=3u.                         (40)
```

Thus

```text
rho=1/3,                  theta=81u,
u=theta/81,               theta'=81.                    (41)
```

The exact companion verifies `(24)` in both nontrivial multistage controls,
including the target and prefactor normalization constants.

## 8. Three-view synchronization boundary

The strict-stage torsor `(33)` sharpens the relation to the incoming
three-view tomography results without manufacturing a physical map.
[THM-3072](THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md)
and
[THM-3076](THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law.md)
reconstruct additive flag signals from all three conjugate quotient views.
That linear reconstruction does not select a section of a multiplicative
`mu_d` root torsor.

If three conjugate local charts retain only the descended scalars `C_j=c_j^d`
and their complements, then each chart carries the lift torsor `(33)`.
Synchronizing the primitive keys requires either

```text
- a common physical atom/point identifying the three roots, or
- a based continuation whose cyclic holonomy in mu_d is trivial.          (42)
```

This is a typing necessity, not an impossibility theorem.  No construction
here identifies the toric charts with THM-3072's three `C2` quotients, and no
claim is made that their torsor classes are nontrivial.  The proved message
is narrower: quotient tomography of scalar/additive data cannot by itself
perform the primitive-root choice used in `(34)`.

For `C3`, a strict stage with `3|d` contains a `mu_3` sub-torsor.  THM-3080's
terminal primitivity forces the local tower eventually to break that factor,
but synchronizing the break across conjugate charts still needs the sidecar
in `(42)`.

## 9. Exact evidence and boundary

The companion, replayed byte-identically under normal and optimized Python
during an independent line audit, checks:

- every primitive terminal Bezout chart and leading exponent decoder for
  coprime `1<=g,e<=24` and a hostile exponent box;
- rational-map degree conventions and the general Mobius inverse identity;
- the affine and reciprocal-affine residue decoders in `(37),(39),(41)`;
- the exact autonomous ODE `(24)` on both multistage packets; and
- all displayed `rho`, `theta`, target, and prefactor constants.

Reproduce with

```bash
python3 04-computation/jc_terminal_toric_residue_decoder_thm3081.py
python3 -O 04-computation/jc_terminal_toric_residue_decoder_thm3081.py
```

The theorem identifies the terminal residue coordinate and the exact finite
root sidecar discarded at strict stages.  It does not make the Laurent keys
polynomial, globalize the local tower, straighten an arbitrary Jelonek curve,
construct a common `A4` atom or holonomy, or exclude `C3`, `A4/S4`, `G1`,
`JC(2)`, or `DC(2)`.
