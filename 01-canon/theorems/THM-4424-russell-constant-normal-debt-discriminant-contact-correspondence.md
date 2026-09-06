---
id: THM-4424
title: "Russell constant-normal debt and discriminant-contact correspondence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT TWO-ENGINE AUDIT, in the formal-local
  labelled-triple and inherited retained-debt scope stated below. For the
  Russell compiler's constant pencil Q+s, the successive D2, D4, D6, and
  order-eight kappa residuals are, up to explicit nonzero rational factors,
  exactly the s^2,s^3,s^4,s^5 coefficients of the unique transverse x
  compensator that retains the triple. The exceptional quartic therefore
  cuts out fourth-order tangency, and its four conjugate fibres have exact
  contact order five because the next coefficient is -9*kappa/20 !=0.
  This proves no all-order debt formula, polynomial termination, global
  discriminant statement, target lift, Keller pair, JC(2), or DC(2).
source: root / JC2 and arXiv continuation session, 2026-09-05
audit: >
  PASS. A hand-written Fraction implementation works directly in the
  four-dimensional exceptional quotient and solves the six formal equations
  through order twelve, with 159 dynamic checks. An independent SymPy engine
  derives the compiler Jacobian, follows the nested symbolic coefficient
  strata, and recovers all four debt/contact ratios, tangent-relation
  periods, adjoint first-exit identities, and the same fifth coefficient in
  78 checks. Normal, optimized, and fixed-hash-seed replays
  byte-match their frozen outputs. Hashes are recorded below.
depends_on:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
  - THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification
  - THM-4411-first-order-collision-transgression-seminormal-tradeoff
related:
  - THM-4397-rank-two-poisson-counterexample-symplectic-gauge-equivalence
  - THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel
  - THM-4408-rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall
script: 04-computation/jc2_exceptional_quartic_constant_normal_high_order_scout_s616.py
output: 05-knowledge/results/jc2_exceptional_quartic_constant_normal_high_order_scout_s616.out
independent_script: 04-computation/jc2_constant_normal_debt_contact_independent_audit_s616.py
independent_output: 05-knowledge/results/jc2_constant_normal_debt_contact_independent_audit_s616.out
script_sha256: f674c16cab7bcc37d29dd3b6a2d933017cf03828d82af8408e3efd8c81862440
output_sha256: 65369a2cbd8246cb1630145561de6b1d00f60f757d71387f6d9cea62d93342c9
independent_script_sha256: 96edf5d7252bc401cbc45d3795cac1b6e959a2b2ccf01b6f5672c93d46e95d20
independent_output_sha256: 6f3e3ad831ac8e02cefbd1bdd1104a24d1d8e80159f2d53efcb84e49d3cbb4cd
hash_basis: raw LF bytes
---

# THM-4424 -- Russell constant-normal debt and discriminant-contact correspondence

**PROVED + VERIFIED-EXACT + INDEPENDENT TWO-ENGINE AUDIT, FOR THE FORMAL-LOCAL
LABELLED TRIPLE AND THE FOUR INHERITED RETAINED DEBTS BELOW. THIS IS NOT AN
ALL-ORDER IDENTITY, A GLOBAL DISCRIMINANT CLASSIFICATION, POLYNOMIAL
TERMINATION, A TARGET LIFT, A KELLER PAIR, `JC(2)`, OR `DC(2)`.**

The obstruction ladder that produced the exceptional quartic has a geometric
meaning that was previously hidden. Its equations are the successive Taylor
coefficients of one transverse graph: the local hypersurface on which the
Russell restriction curve keeps its labelled ordinary triple point.

## 1. The local collision hypersurface

Work in characteristic zero with

```text
P=x^2(x^2-1)^2,
Q1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,
Q_(a,b,c)=Q1+P(a+bx+cx^2).                              (1)
```

For any polynomial `q`, use the Russell surface compiler

```text
D=1+x^2q,
B=(D-1)(D+2)^2,
C=xD(D+2),
E=q(D+3),                  C^2E=B(B+4).                 (2)
```

Every polynomial in `(1)` maps the three labelled normalization points
`x=-1,0,1` to `(B,C,E)=(0,0,-3)`. Near that target, `(C,E)` are regular
surface coordinates because the derivative of `C^2E-B(B+4)` with respect to
`B` is `-4`.

Consider the constant normal pencil, with a transverse `x` coefficient left
free:

```text
q_s=Q+s+chi_Q(s)x,             chi_Q(s) in s k[[s]].    (3)
```

Ask for formal sections

```text
z_-(s)=-1+O(s), z_0(s)=O(s), z_+(s)=1+O(s)             (4)
```

whose three `(C,E)` images equal one formal target `(C_*(s),E_*(s))`. The
Jacobian of these six equations in

```text
(z_-,z_0,z_+,C_*,E_*,chi_Q)
```

at `s=0` is

```text
J = [ 3  0  0 -1  0 -2 ]
    [-9  0  0  0 -1  2 ]
    [ 0  3  0 -1  0  0 ]
    [ 0  4  0  0 -1  0 ]
    [ 0  0  3 -1  0 -2 ]
    [ 0  0  9  0 -1 -2 ],          det J=-288.          (5)
```

Thus `(3)--(4)` have a unique formal solution. In any finite coefficient
space containing `Q,1,x`, the nearby labelled-triple locus is a smooth formal
hypersurface, with the `x` coefficient transverse. Write

```text
chi_Q(s)=sum_(n>=1) chi_n(Q)s^n.                        (6)
```

THM-4411 gives `chi_1(Q)=0`: the constant normal satisfies the first-order
collision-period equation. The new content is the exact higher curvature of
this hypersurface.

## 2. Orienting the inherited debt residuals

The older retained calculations successively solved the pullback-Jacobian
equations for the quadratic fold `q=Q(x)+t^2`, with `s=t^2`. At each even
stage let `delta_(2m)` mean the actual remaining scalar

```text
delta_(2m)=Lambda(provisional J_(2m))                  (7)
```

after the lower coefficients have been killed. This orientation matters
because the historical symbols did not all put the residual on the same side
of their displayed identity. Concretely,

```text
delta_2 = D_2,
delta_4 = D_4,
delta_6 = -D_6,
delta_8 = -kappa,                                      (8)
```

where `D_2,D_4,D_6,kappa` are exactly the quantities of
THM-3641, THM-3677, THM-3683, and THM-4046. No new interpretation of those
theorems is being inserted into their proofs; `(7)--(8)` only align signs.

## 3. Four exact debt/contact identities

The identities hold successively on the stratum where all preceding
residuals vanish:

```text
boxed:
chi_2=(9/8) delta_2,
chi_3=(9/12)delta_4,
chi_4=(9/16)delta_6,
chi_5=(9/20)delta_8.                                   (9)
```

Equivalently, for `m=1,2,3,4`, the four verified levels have the uniform
finite pattern

```text
chi_(m+1)=9 delta_(2m)/(4(m+1)).                       (10)
```

Equation `(10)` is **not** asserted beyond these four computed levels.

### 3.0 The common conormal mechanism

The rational factors in `(9)` are not inferred from four matching values.
They come from the same moving tangent relation. Along the compensated
collision let

```text
t_i(s)=d(C,E)/dx at z_i(s),
u(s)=(det(t_0,t_+),det(t_+,t_-),det(t_-,t_0)).          (10a)
```

Then `sum_i u_i t_i=0`. At the base triple,

```text
u(0)=(15,-54,39),
lambda_0=u(0)/54=(5/18,-1,13/18).                     (10b)
```

Normalize formally by `lambda_s=u(s)/(-u_0(s))`, where `u_0` denotes the
middle component, and let

```text
Pi_Q(s)=lambda_-(s)+lambda_0(s)+lambda_+(s).            (10c)
```

This is the moving tangent-relation period of the constant source density.
On the four successive zero strata the exact expansions give

```text
[s^m]Pi_Q(s)=-delta_(2m),                 m=1,2,3,4,   (10d)
```

and all lower coefficients vanish. Differentiate the collision equation
along `(3)`. Its normal velocity is `1+chi_Q'(s)x`. Applying the
collision--transgression relation of THM-4411, then taking the first live
coefficient, gives

```text
(m+1)chi_(m+1)=-(9/4)[s^m]Pi_Q=(9/4)delta_(2m).        (10e)
```

The independent companion constructs the moving tangents directly from the
compiler, checks their raw base relation, verifies `(10d)`, and verifies
`(10e)` separately at all four levels. Thus `(9)` has a common conormal
explanation. What remains open is an all-order theorem identifying every
future target-lift debt with the corresponding moving period, especially
outside this compiler.

Here are the literal coefficient calculations.

### 3.1 Second coefficient: the zero-second-debt plane

For general `(a,b,c)` in `(1)`, exact recursion with `(5)` gives

```text
chi_2=-(36a+16b+36c+259)/9.                            (11)
```

The THM-3641 principal debt is

```text
delta_2=D_2=-8(36a+16b+36c+259)/81,                   (12)
```

so `(11)` is the first identity in `(9)`. The complete degree-eight
zero-second-debt plane is parametrized as

```text
a=-259/36+p+4r,       b=-9r,       c=-p.               (13)
```

On `(13)`, `chi_2=0` identically.

### 3.2 Third coefficient: the zero-fourth-debt parabola

Restrict the already solved lower formal coefficients to `(13)`. The next
recursion gives

```text
chi_3=16(729p-42120r^2+15192r+5717)/2187.              (14)
```

THM-3677's fourth debt is

```text
delta_4=D_4=64(729p-42120r^2+15192r+5717)/6561.        (15)
```

Thus `chi_3=(3/4)D_4`, and `chi_3=0` cuts out exactly the old parabola

```text
p=p(r)=(520/9)r^2-(1688/81)r-5717/729.                (16)
```

### 3.3 Fourth coefficient: the exceptional quartic

On `(16)`, the next coefficient is

```text
chi_4=16F(r)/177147,                                   (17)

F(r)=72783360r^4-77822208r^3-28419741r^2
     +7849770r-1276420.                                (18)
```

THM-3683 has

```text
D_6=-256F(r)/1594323,          delta_6=-D_6,           (19)
```

so `(17)` is `chi_4=(9/16)delta_6`. Hence the four roots of the irreducible
quartic `(18)` are not merely the retained sixth-debt survivors: they are
exactly the points of the parabola where the constant line has at least one
further order of contact with the labelled-triple hypersurface.

### 3.4 Fifth coefficient: the order-eight exit

Let

```text
K=Q[alpha]/(F(alpha)).                                  (20)
```

At the exceptional point, exact reduction in the basis
`(1,alpha,alpha^2,alpha^3)` gives

```text
chi_1=chi_2=chi_3=chi_4=0,

chi_5=
 259188338368/3^17
 -(46584993664/3^14)alpha
 +(23019960448/3^12)alpha^2
 +(9180348416/3^11)alpha^3                             (21)

 =64(4049817787-19653044202alpha+87403912326alpha^2
       +104569906176alpha^3)/3^17.
```

THM-4046's order-eight constant response is

```text
kappa=
 -5183766767360/3^19
 +(931699873280/3^16)alpha
 -(460399208960/3^14)alpha^2
 -(183606968320/3^13)alpha^3.                          (22)
```

Coordinatewise comparison of `(21)--(22)` gives

```text
boxed: chi_5=-(9/20)kappa=(9/20)delta_8.               (23)
```

The degree-three representative in `(21)` is nonzero modulo the degree-four
minimal polynomial. Therefore it is nonzero at every embedding of `K`.

## 4. Exact fifth-order contact

Equations `(6)` and `(21)` say that the straight constant line `Q+s`, with
no `x` correction, meets the formal labelled-triple hypersurface with exact
intersection multiplicity five. More explicitly, its three labelled
sections and common target can be solved with equal images modulo `s^5`, but
they cannot be lifted to equal images modulo `s^6`. The latter statement is
not a choice-of-gauge artifact: deleting the `x` column from `(5)` leaves
rank five, and the unique augmented solution at order five has the nonzero
last coordinate `(21)`.

The full compensated solution nevertheless exists uniquely to every formal
order by `(5)`. The primary companion carries it through order twelve and
finds nonzero coefficients at every order from five through twelve. Those
extra coefficients are finite positive controls, not a claim that the
series is algebraic, convergent, or polynomial.

## 5. Why this is useful, and why it does not prove `JC(2)`

The old sequence

```text
zero D2 plane -> zero D4 parabola -> zero D6 quartic -> nonzero D8 exit
```

was constructed from cokernel identities for pulled-back target two-forms.
Equation `(9)` identifies the same sequence with successive osculation
strata of an actual geometric object, the local labelled-triple
discriminant. This explains why the parabola and quartic reappeared in a
seemingly unrelated formal collision calculation: they are the same
equations, not numerical analogies.

The relation also sharpens the bridge to the Hamiltonian correction in
arXiv:2608.23777. In both settings a two-form residual is paid by a
branch-sensitive coordinate that does not belong to the smaller descended
core. Here that coordinate is the transverse `x` coefficient `chi_Q(s)`;
THM-4408 identifies the analogous role of Long's added Hamiltonian primitive.
At the exceptional quartic, the first unavoidable coefficient of this
sidecar is literally a rational multiple of the retained two-form
obstruction `kappa`.

This remains a local duality, not a planar counterexample or a proof of the
Jacobian conjecture. The labelled collision is a singular restriction curve,
not a polynomial Keller map. The formal compensator need not terminate or
descend from a target algebra. THM-4046 already proves that the associated
actual target lift stops at `J_8`. No passage from Long's four-dimensional
Poisson map to an injective or noninjective planar Keller map is supplied.

## 6. Verification

Reproduce with

```text
python3 -B 04-computation/jc2_exceptional_quartic_constant_normal_high_order_scout_s616.py
python3 -B -O 04-computation/jc2_exceptional_quartic_constant_normal_high_order_scout_s616.py
PYTHONHASHSEED=17 python3 -B 04-computation/jc2_exceptional_quartic_constant_normal_high_order_scout_s616.py

python3 -B 04-computation/jc2_constant_normal_debt_contact_independent_audit_s616.py
python3 -B -O 04-computation/jc2_constant_normal_debt_contact_independent_audit_s616.py
PYTHONHASHSEED=17 python3 -B 04-computation/jc2_constant_normal_debt_contact_independent_audit_s616.py
```

The primary script implements the quartic field and truncated formal series
using only `Fraction` tuples, derives each coefficient recursively, checks the
known quadratic second-order hostile, and verifies `(23)`. The independent
script imports neither it nor its arithmetic. It uses direct SymPy
substitution and follows the nested symbolic specializations

```text
(a,b,c) -> (p,r) -> r -> alpha mod F,
```

recovering `(11)`, `(14)`, `(17)`, and `(21)` independently. Both scripts
derive or check the determinant in `(5)`, keep all arithmetic exact, and
contain explicit hostile/scope gates.
