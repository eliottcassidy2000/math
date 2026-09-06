# Exceptional triple: formal collision compensation and a second-order hostile

**PROVED in the stated formal-local compiler scope; exact symbolic audit
and independent proof/second-order substitution audit PASS.** Research date:
2026-09-05.

The first-order tradeoff below was independently recovered in this session
from THM-4404, while THM-4411 was an empty reservation. A fresh fetch at
`71caf0315` brought in its concurrent proof. It is now an **inherited** result,
not a new first-order claim of this session. The additions here are the
invertible six-variable formal lifting system and the exact second-order
failure of a collision-preserving straight pencil.

## Inheritance and exact scope

The closest proved mechanism is
[THM-4411 / first-order-collision-transgression-seminormal-tradeoff](../../01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md),
with the rank-88 image and missing line from
[THM-4404 / exceptional-quartic-descended-two-form-seminormal-cokernel](../../01-canon/theorems/THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel.md).
The canonical hostile is `h=x`: it fills the missing transgression coordinate
but separates the triple at first order. The corrected near miss is treating
a first-order direction as an integrated moving graph. The least-used
sidecar is the formal nonlinear compatibility equation after the first-order
repair system is solved.

Work over the characteristic-zero exceptional field `K` of
[THM-4381 / exceptional-quartic-seminormalization-and-conductor-fibre-classification](../../01-canon/theorems/THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification.md).
For a polynomial `q`, retain the Russell compiler

```text
D=1+x^2 q, B=(D-1)(D+2)^2, C=xD(D+2), E=q(D+3).
C^2 E=B(B+4).
```

The exceptional polynomial has the explicit form

```text
P=x^2(x^2-1)^2,
Q1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,
p(alpha)=(520/9)alpha^2-(1688/81)alpha-5717/729,
Q=Q1+P*(a+b*x+c*x^2),
a=-259/36+p(alpha)+4alpha, b=-9alpha, c=-p(alpha),
F(alpha)=72783360alpha^4-77822208alpha^3-28419741alpha^2
         +7849770alpha-1276420=0.
```

The auxiliary letters `a,b,c` in that formula are polynomial coefficients,
not surface coordinates. The symbolic verifier first allows them to vary
freely, then specializes to the exceptional field.

## Formal compensation theorem

Fix any finite-dimensional polynomial subspace `V` containing `Q` and `x`;
for example `V=K[x]_(<=8)`. Split the vector space of polynomial perturbations
as `K*x direct_sum W`. Given arbitrary formal parameters for `W`, vanishing
at the origin, there exist unique formal series

```text
z_-(s)=-1+O(s), z_0(s)=O(s), z_+(s)=1+O(s),
C_*(s)=O(s), E_*(s)=-3+O(s), c_x(s)=O(s)
```

such that the polynomial family

```text
Q_s=Q+H_s+c_x(s)*x,             H_s in s W[[s]],
```

sends all three labelled sections `z_-,z_0,z_+` to one target point. The same
statement holds with finitely many formal base parameters. Thus the formal
locus of nearby coefficient vectors retaining this labelled triple is a
smooth hypersurface, and the coefficient of `x` is a transverse coordinate.
This concerns polynomial dependence on `x` with formal coefficient series;
it does not assert polynomial termination in the base parameter.

The same recursion also accepts any **prescribed** forcing `H_s in s V[[s]]`,
even with an `x` component: it uniquely determines the additional correction
`c_x`. The total `x` coefficient then includes the prescribed contribution.
Using a complement `W` is needed only to parametrize the hypersurface without
redundant free coordinates. The fixed-pencil controls below use this more
general prescribed-forcing formulation.

**Proof.** Near `(B,C,E)=(0,0,-3)`, the surface relation has derivative `-4`
with respect to `B`. Hence `(C,E)` are formal local coordinates, and equality
of those two coordinates implies equality of the complete target point near
this point. Impose six equations

```text
C(z_i,Q_s(z_i))=C_*, E(z_i,Q_s(z_i))=E_*, i=-,0,+.
```

Their Jacobian in the ordered unknowns
`(z_-,z_0,z_+,C_*,E_*,c_x)` at the undeformed triple is

```text
J = [ 3  0  0 -1  0 -2 ]
    [-9  0  0  0 -1  2 ]
    [ 0  3  0 -1  0  0 ]
    [ 0  4  0  0 -1  0 ]
    [ 0  0  3 -1  0 -2 ]
    [ 0  0  9  0 -1 -2 ].
```

Its determinant is `-288`, a nonzero element of `K`. A direct formal proof
avoids any convergence claim: at every homogeneous base degree the unknown
coefficient vector occurs through this same `J`, with all other terms
determined at lower degrees. Invert `J` to obtain the unique next coefficient.
This recursively proves existence and uniqueness at every degree. The
remaining `W` coefficients are free, proving the hypersurface assertion.

The nonzero constant terms of `z_i-z_j` keep the sections distinct. The
nonzero pairwise tangent determinants keep their tangents transverse
formally. No assertion about additional singularities outside the three
retained branches follows merely from this local argument.

## Exact tangent and curvature of the constraint

At `x=(-1,0,1)`, write the two surface-coordinate rows as

```text
t_i=(C',E')=(3,-9),(3,4),(3,9),
n_i=(partial_q C,partial_q E)=(2,-2),(0,4),(-2,-2).
```

Their determinants are all `12`. Deleting the last column of `J` gives the
rank-five endpoint/common-target repair matrix `A`. Its one compatibility
condition for a normal polynomial `h` is

```text
L(h)=5h(-1)-18h(0)+13h(1)=0.
```

This is the THM-4411 first-order criterion. For a prescribed perturbation
`H_s=s*h`, the first compensator coefficient is `c_x'(0)=-L(h)/8`.
For collision-preserving `h`, it is zero. In particular `h=x` has obstruction
`8`; a degree-two surviving direction is `h=4x^2-9x`.

The important new hostile is that the latter direction does **not** preserve
the triple along its uncorrected straight pencil. Write the repaired family

```text
Q_s=Q+s*(4x^2-9x)+s^2*c_2*x+O(s^3).
```

The first endpoint and common target coefficients are

```text
(z_-',z_0',z_+',C_*',E_*')=(-14/3,4,2/3,12,16).
```

At second order, the unique needed compensator is

```text
c_2=-(3072*a-7424*b+8256*c+102511)/144
   =(2695680*alpha^2-1684224*alpha-1089575)/1296.
```

This element is nonzero: its numerator has degree two and is nonzero, while
the minimal polynomial defining `K` has degree four. Therefore the straight
pencil `Q+s*(4x^2-9x)` cannot maintain the triple modulo `s^3`, even though it
does maintain it modulo `s^2` with moved endpoints. The six-variable formal
theorem repairs this failure by adding the displayed `s^2*c_2*x` and then
one uniquely determined coefficient at every subsequent order.

The general second-order formula is checked by substituting into the six
equations. If `eta_i` is the first endpoint displacement, its known
second-order residual is

```text
(1/2) gamma''(i)*eta_i^2
 + (h*partial_q gamma)'(i)*eta_i
 + (1/2) partial_q^2 gamma(i)*h(i)^2,
```

where `gamma=(C,E)` and primes denote total `x` differentiation along `Q`.
Multiplication by `-J^(-1)` gives the second-order unknowns. The companion
performs these polynomial identities symbolically for arbitrary `a,b,c`,
before specializing. This is an exact coefficient calculation, not a
numerical branch fit.

**Positive boundary.** The constant direction `h=1` has first displacements
`(-2/3,0,2/3,0,4)`. Its general second-order compensator is
`-(36a+16b+36c+259)/9`; this vanishes for the exceptional coefficients.
Its third-order compensator also vanishes there. These are finite jet
controls only; an uncorrected all-order constant pencil is not proved here.

## Connection to the seminormal and graph frontiers

The source is a polynomial compiler deformation with labelled endpoint
sections. Its first-order quotient sends branch velocities to wedge
contractions. It preserves the linear collision compatibility equation and
forgets the higher-degree coefficients of the nonlinear equation. The
missing sidecar is exactly the formal `x`-coefficient compensator above.
The cheapest decisive test is the quadratic direction: it is in the tangent
kernel but has the nonzero quadratic correction displayed above.

The same first-order functional controls the fixed-fibre quotient
`E_x=K[x]/dS`: by THM-4404/4411,

```text
Lambda(tau_h(dC wedge dE))=(2/3)*L(h).
```

Consequently collision-preserving first-order motion cannot remove the
seminormal transgression obstruction at that fibre. Formal integrability
does not reverse this statement: the compensator arranges compatibility
rather than constructing a missing descended primitive.

This fills the labelled-endpoint **existence** requirement locally at the
retained triple, but not all of
[THM-4067 / seminormal-period-kernel-and-figure-eight-completeness-obstruction](../../01-canon/theorems/THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction.md).
The graph's edge-primitive module, comparison with the actual target algebra,
global coefficient restrictions, and algebraic/polynomial termination remain
separate. It also does not remove the fixed-subring compensator in
[THM-4408 / rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall](../../01-canon/theorems/THM-4408-rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall.md).
No Keller pair, chart entry, `JC(2)`, or `DC(2)` result follows.

The general viewpoint of parametrization deformations with sections is
classical; see Greuel--Lossen--Shustin, *Introduction to Singularities and
Deformations*, Chapter II, section 2
([primary book](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/greuel.pdf)).
The proof here supplies its own formal recursion and the compiler-specific
matrix and hostile; no external priority claim is made.

## Reproduction

```text
python3 -B 04-computation/synthesis_20260905_transgression.py
python3 -B -O 04-computation/synthesis_20260905_transgression.py
```

The companion checks the surface relation, retained values, exact tangent
and normal rows, rank-five compatibility, determinant `-288`, positive and
hostile normal directions, the symbolic second-order repair, exceptional
specialization, and the constant third-order control. Every test is an
explicit runtime check, active under optimization. Its frozen transcript is
`synthesis_20260905_transgression.out` alongside this report.

The independent companion
`04-computation/synthesis_20260905_wildcard_transgression_audit.py` imports
neither the primary script nor its symbolic engine. It substitutes truncated
series with `Fraction` arithmetic directly into the compiler modulo `s^3`,
on an affine basis and mixed rational controls for `a,b,c`. It checks both
second-order formulas in 70 exact gates, and independently audits the formal
recursion, distinct sections, and the limits on the graph interpretation.
